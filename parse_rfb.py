"""
Parser para o formato binário OmniPro V6 RFB (.rfb) da IBA Dosimetry.

Portado do RadPy (Stephen Terry et al., BSD-3-Clause) para Python 3 puro,
usando apenas stdlib (struct) + numpy.  Sem dependências externas.

Referência original: https://github.com/radpy/RadPy

Estrutura geral do arquivo
--------------------------
  PascalString         versão OmniPro
  [main_header]        parâmetros do feixe (energia, campo, SSD…)
    N × [measurement_data + 2 bytes padding]
  [3 bytes discriminador]
    → 0   : [additional_header]  (novo grupo com diferentes parâmetros)
    → N≠0 : [new_scan_type_header] (mesmo grupo, outro tipo de curva)
  [additional_headers repetidos até EOF]
  12 × 0x00            marcador de fim de arquivo

Todos os valores numéricos são little-endian (nativo Windows / IBA).
Posições no arquivo em mm; convertidas para cm na saída.

Uso:
    from parse_rfb import load_rfb_files
    scans = load_rfb_files(["arquivo.rfb"], index_start=0)
"""

import struct
import numpy as np
from pathlib import Path
from typing import List, Optional

from parse_snctxt import Scan


# ===========================================================================
# Mapeamentos
# ===========================================================================

_PARTICLE_MAP = {
    0: "Photon", 1: "Electron", 2: "Proton",
    3: "Neutron", 4: "Cobalt", 5: "Isotope",
}

_MEDIUM_MAP = {0: "Air", 1: "Water", 2: "Film"}

_DATA_TYPE_MAP = {
    1: "relative_optical_density",
    2: "relative_dose",
    3: "relative_ionization",
    4: "absolute_dose",
    5: "charge",
}

_WEDGE_MAP = {
    -1: "Open", 0: "Hard_Wedge", 1: "Dynamic_Wedge",
    2: "Enhanced_Wedge", 3: "Virtual_Wedge", 4: "Soft_Wedge",
}

_GANTRY_SCALE_MAP = {
    0: "CW_180_Down", 1: "CCW_180_Down",
    2: "CW_180_Up",   3: "CCW_180_Up",
}


# ===========================================================================
# Leitor binário
# ===========================================================================

class _R:
    """Leitor sequencial de buffer binário (little-endian)."""

    def __init__(self, data: bytes) -> None:
        self._d = data
        self._p = 0

    # ---- primitivos --------------------------------------------------------

    def remaining(self) -> int:
        return len(self._d) - self._p

    def at_end(self) -> bool:
        return self._p >= len(self._d)

    def read(self, n: int) -> bytes:
        b = self._d[self._p: self._p + n]
        self._p += n
        return b

    def peek(self, n: int) -> bytes:
        return self._d[self._p: self._p + n]

    def skip(self, n: int) -> None:
        self._p += n

    def u8(self)  -> int:   return struct.unpack_from('B', self._d, self._advance(1))[0]
    def s8(self)  -> int:   return struct.unpack_from('b', self._d, self._advance(1))[0]
    def u16(self) -> int:   return struct.unpack_from('<H', self._d, self._advance(2))[0]
    def s16(self) -> int:   return struct.unpack_from('<h', self._d, self._advance(2))[0]
    def s32(self) -> int:   return struct.unpack_from('<i', self._d, self._advance(4))[0]
    def f64(self) -> float: return struct.unpack_from('<d', self._d, self._advance(8))[0]

    def _advance(self, n: int) -> int:
        p = self._p; self._p += n; return p

    def pascal_string(self) -> str:
        n = self.u8()
        return self.read(n).decode('latin-1', errors='replace')

    def skip_to(self, marker: int) -> None:
        """Avança até encontrar o byte `marker` (inclusive)."""
        while not self.at_end():
            if self.u8() == marker:
                return

    def last_byte_as_s8(self, raw: bytes) -> int:
        """Retorna o último byte de `raw` como signed int8 (ScanTypeAdapter)."""
        return struct.unpack('b', bytes([raw[-1]]))[0]


# ===========================================================================
# Leitura de blocos estruturais
# ===========================================================================

def _read_shared_header_fields(r: _R) -> dict:
    """
    Lê os campos comuns entre main_header (header_data) e additional_header,
    de rad_device até jaw_crossplane_positive + gantry_scale + 1 byte padding.

    Retorna dict com os campos extraídos.
    """
    rad_device  = r.pascal_string()
    r.skip(2)                       # Padding(2)
    energy      = r.f64()
    particle    = r.u8()
    r.skip(1)                       # Padding(1) / UNInt8 raw
    wedge_type  = r.s16()
    r.skip(2)                       # Padding(2)
    wedge_angle = r.u8()
    r.skip(3)                       # Padding(3)
    gantry_angle     = r.u16()
    r.skip(2)
    collimator_angle = r.u16()
    r.skip(2)
    ssd = r.f64()
    r.skip(2)
    sad = r.f64()
    applicator  = r.pascal_string()
    medium      = r.s8()
    r.skip(1)                       # Padding(1)
    institution = r.pascal_string()
    address     = r.pascal_string()
    telephone   = r.pascal_string()
    email       = r.pascal_string()
    r.skip(2)
    inplane_jaw_neg  = r.f64()
    r.skip(2)
    inplane_jaw_pos  = r.f64()
    r.skip(2)
    crossplane_jaw_neg = r.f64()
    r.skip(2)
    crossplane_jaw_pos = r.f64()
    gantry_scale = r.u8()
    r.skip(1)                       # Padding(1)

    return {
        'rad_device':          rad_device,
        'energy':              energy,
        'particle':            particle,
        'wedge_type':          wedge_type,
        'wedge_angle':         wedge_angle,
        'gantry_angle':        gantry_angle,
        'collimator_angle':    collimator_angle,
        'SSD':                 ssd,
        'SAD':                 sad,
        'applicator':          applicator,
        'medium':              medium,
        'institution':         institution,
        'inplane_jaw_negative':    inplane_jaw_neg,
        'inplane_jaw_positive':    inplane_jaw_pos,
        'crossplane_jaw_negative': crossplane_jaw_neg,
        'crossplane_jaw_positive': crossplane_jaw_pos,
        'gantry_scale':        gantry_scale,
    }


def _read_main_header(r: _R) -> dict:
    """
    Lê main_header completo (header_data com padding inicial de 13 bytes
    e a lógica de num_scans/scan_type do primeiro bloco).
    """
    r.skip(13)                     # Padding(13) — "CBeam" prefix area
    h = _read_shared_header_fields(r)

    # ScanTypeAdapter: lê bytes até encontrar um não-zero;
    # o último byte (não-zero) é num_scans_with_this_header (signed int8)
    while True:
        b = r.u8()
        if b != 0:
            h['num_scans_with_this_header'] = struct.unpack('b', bytes([b]))[0]
            break

    r.skip(5)                      # Padding(5)

    # scan_type_field: Byte(length) + Padding(2) + MetaField(scan_type, length-1)
    st_len  = r.u8()
    r.skip(2)
    h['scan_type_str'] = r.read(st_len - 1).decode('latin-1', errors='replace')

    return h


def _read_additional_header(r: _R) -> dict:
    """
    Lê additional_header (novo grupo de scans com parâmetros diferentes).
    Deve ser chamado DEPOIS de _skip_to(0x80) do bloco anterior.
    """
    # Primeiro RepeatUntil(obj == '\x80', ...) — avança até 0x80
    r.skip_to(0x80)

    h = _read_shared_header_fields(r)

    # Byte("scan_type"): 0 = DepthDose, else = Profile
    scan_type_byte = r.u8()
    h['scan_type_byte'] = scan_type_byte
    h['scan_type_str']  = "DepthDose" if scan_type_byte == 0 else "Profile"
    h['num_scans_with_this_header'] = 0  # não fixo; controlado por delimitador

    # Segundo RepeatUntil(obj == '\x80', ...) — avança até 0x80
    r.skip_to(0x80)

    return h


def _read_measurement_data(r: _R) -> dict:
    """
    Lê um bloco measurement_data completo.
    Retorna dict com campos relevantes + arrays abscissa e ordinate.
    """
    measured_date  = r.s32()
    modified_date  = r.s32()
    data_type      = r.u8()
    r.skip(8)                  # chamber_radius
    r.skip(8)                  # calibration_factor
    r.skip(8)                  # unknown1
    r.skip(8)                  # unknown2
    calibration_date = r.pascal_string()
    r.skip(8)                  # peff_offset
    detector      = r.pascal_string()
    detector_type = r.u8()
    r.skip(1)                  # raw2

    operator            = r.pascal_string()
    measurement_comment = r.pascal_string()

    crossplane_servo = r.s16()
    inplane_servo    = r.s16()
    depth_servo      = r.s16()
    r.skip(2)                  # measurements_per_point
    r.skip(2)                  # raw3
    r.skip(8)                  # scan_speed
    r.skip(2)                  # servo_type
    r.skip(2)                  # measurement_mode_a
    r.skip(6)                  # raw4_b

    isocenter_crossplane = r.f64()
    isocenter_inplane    = r.f64()
    isocenter_depth      = r.f64()
    r.skip(8)                  # normalization_crossplane
    r.skip(8)                  # normalization_inplane
    r.skip(8)                  # normalization_depth

    # 6 float64 electrometer values
    r.skip(8 * 6)              # field_norm, ref_norm, field_dark, ref_dark, field_hv, ref_hv

    field_hv = 0.0             # extracted below via re-read of that area (optional)

    r.skip(2)                  # field_gain
    r.skip(2)                  # reference_gain
    r.pascal_string()          # field_range
    r.pascal_string()          # reference_range
    r.skip(8)                  # water_surface_correction
    r.skip(1)                  # measurement_mode_b
    r.skip(1)                  # raw5
    r.skip(1)                  # HV_connection
    r.skip(1)                  # raw6

    r.skip(8)                  # reference_maximum
    r.skip(8)                  # reference_minimum
    r.skip(8)                  # reference_average
    r.skip(8)                  # electrometer_sampling_time
    r.skip(2)                  # electrometer_type
    r.skip(8)                  # renormalization_factor
    r.skip(8)                  # curve_offset

    r.pascal_string()          # setup_comment
    r.skip(1)                  # ca24_calibration (Flag)
    r.skip(1)                  # Padding(1)

    # 4 positions A/B/C/D (3 float64 each = 12 × 8 = 96 bytes)
    r.skip(96)
    r.skip(10)                 # raw9

    scan_start_crossplane = r.f64()
    scan_start_inplane    = r.f64()
    scan_start_depth      = r.f64()
    scan_end_crossplane   = r.f64()
    scan_end_inplane      = r.f64()
    scan_end_depth        = r.f64()

    # data array: SNInt16(length) + length × 2 × NFloat64
    length    = r.s16()
    abscissa  = np.empty(length, dtype=np.float64)
    ordinate  = np.empty(length, dtype=np.float64)
    for i in range(length):
        abscissa[i] = r.f64()
        ordinate[i] = r.f64()

    return {
        'measured_date':          measured_date,
        'detector':               detector,
        'operator':               operator,
        'measurement_comment':    measurement_comment,
        'crossplane_servo':       crossplane_servo,
        'inplane_servo':          inplane_servo,
        'depth_servo':            depth_servo,
        'isocenter_crossplane':   isocenter_crossplane,
        'isocenter_inplane':      isocenter_inplane,
        'isocenter_depth':        isocenter_depth,
        'scan_start_crossplane':  scan_start_crossplane,
        'scan_start_inplane':     scan_start_inplane,
        'scan_start_depth':       scan_start_depth,
        'scan_end_crossplane':    scan_end_crossplane,
        'scan_end_inplane':       scan_end_inplane,
        'scan_end_depth':         scan_end_depth,
        'abscissa':               abscissa,
        'ordinate':               ordinate,
    }


# ===========================================================================
# Inferência de tipo de curva e construção de Scan
# ===========================================================================

def _infer_scan_type(m: dict) -> str:
    """
    Determina o tipo de curva comparando a variação nas 3 coordenadas entre
    os pontos de início e fim da varredura.

    Retorna: "Depth Scan" | "Crossline" | "Inline" | "Diagonal"
    """
    dx = abs(m['scan_end_crossplane'] - m['scan_start_crossplane'])
    dy = abs(m['scan_end_inplane']    - m['scan_start_inplane'])
    dz = abs(m['scan_end_depth']      - m['scan_start_depth'])

    max_d = max(dx, dy, dz)
    if max_d < 0.5:                       # scan praticamente em ponto único
        return 'Crossline'

    if dz == max_d:
        return 'Depth Scan'
    if dx > dy:
        return 'Crossline'
    if dy > dx:
        return 'Inline'
    return 'Diagonal'


def _build_scan(h: dict, m: dict, index: int) -> Optional[Scan]:
    """
    Constrói um Scan a partir dos dicts de main_header e measurement_data.
    Retorna None se os dados forem inválidos.
    """
    abscissa = m['abscissa']
    ordinate = m['ordinate']
    if len(abscissa) < 2:
        return None

    # posições mm → cm
    pos_cm = abscissa / 10.0
    # garantir que dose seja positiva
    if np.mean(ordinate) < 0:
        ordinate = -ordinate

    scan_type = _infer_scan_type(m)

    x = np.zeros_like(pos_cm)
    y = np.zeros_like(pos_cm)
    z = np.zeros_like(pos_cm)

    if scan_type == 'Depth Scan':
        z = pos_cm
    elif scan_type == 'Inline':
        y = pos_cm
    else:  # Crossline / Diagonal
        x = pos_cm

    # modalidade e energia
    particle  = h.get('particle', 0)
    beam_type = _PARTICLE_MAP.get(particle, 'Photon')
    energy    = str(h.get('energy', ''))
    # remover ".0" se for inteiro
    try:
        e_f = float(energy)
        energy = str(int(e_f)) if e_f == int(e_f) else energy
    except ValueError:
        pass

    # campo: jaws em mm → cm
    cp_neg = abs(h.get('crossplane_jaw_negative', 0))
    cp_pos = abs(h.get('crossplane_jaw_positive', 0))
    ip_neg = abs(h.get('inplane_jaw_negative',    0))
    ip_pos = abs(h.get('inplane_jaw_positive',    0))
    field_x = f"{(cp_neg + cp_pos) / 10:.4g}" if (cp_neg + cp_pos) > 0 else '?'
    field_y = f"{(ip_neg + ip_pos) / 10:.4g}" if (ip_neg + ip_pos) > 0 else '?'

    # SSD mm → cm
    ssd_val = h.get('SSD', 0)
    ssd_str = f"{ssd_val / 10:.4g}" if ssd_val > 0 else '?'

    medium   = _MEDIUM_MAP.get(h.get('medium', 1), 'Water')
    detector = m.get('detector', '')

    return Scan(
        index     = index,
        beam_type = beam_type,
        energy    = energy,
        field_x   = field_x,
        field_y   = field_y,
        scan_type = scan_type,
        ssd       = ssd_str,
        medium    = medium,
        detector  = detector,
        x = x, y = y, z = z,
        dose = ordinate,
    )


# ===========================================================================
# Parser principal
# ===========================================================================

def parse_rfb(filepath: str, index_start: int = 0) -> List[Scan]:
    """
    Lê um arquivo .rfb OmniPro V6 e retorna lista de Scan.

    Parâmetros
    ----------
    filepath    : caminho do arquivo .rfb
    index_start : índice base para os Scan retornados

    Lança IOError se o arquivo não puder ser lido ou não for um RFB válido.
    """
    with open(filepath, 'rb') as fh:
        data = fh.read()

    r = _R(data)
    scans: List[Scan] = []

    # ---- bloco 1: main_header + seus scans ---------------------------------

    # PascalString: versão OmniPro (primeira coisa no arquivo)
    omnipro_version = r.pascal_string()
    if not omnipro_version.startswith('OmniPro') and \
       not omnipro_version.startswith('Ver'):
        # tentativa de leitura mesmo assim; pode ser uma versão diferente
        pass

    try:
        main_hdr = _read_main_header(r)
    except Exception as exc:
        raise IOError(f"Falha ao ler main_header do RFB: {exc}") from exc

    n_main = abs(main_hdr.get('num_scans_with_this_header', 0))
    for _ in range(n_main):
        try:
            m = _read_measurement_data(r)
            r.skip(2)                  # main_measurement_data: Padding(2) ao final
            s = _build_scan(main_hdr, m, len(scans) + index_start)
            if s is not None:
                scans.append(s)
        except Exception:
            break

    # ---- blocos seguintes: discriminador → new_scan_type ou additional -----

    while r.remaining() > 12:
        try:
            # 3 bytes discriminadores (ScanTypeAdapter)
            disc_raw = r.read(3)
            disc = struct.unpack('b', bytes([disc_raw[-1]]))[0]
        except Exception:
            break

        if disc == 0:
            # --- additional_header: novo grupo com parâmetros diferentes ---
            try:
                add_hdr = _read_additional_header(r)
            except Exception:
                break

            # additional_measurement_data: repetição controlada por delimitador
            while r.remaining() > 12:
                try:
                    m = _read_measurement_data(r)
                    r.skip(1)          # Padding(1)
                    delim = r.u8()     # Byte("delimiter")

                    # Peek 3 bytes para delimiter2 (ScanTypeAdapter)
                    if r.remaining() >= 3:
                        peek = r.peek(3)
                        delim2 = struct.unpack('b', bytes([peek[-1]]))[0]
                    else:
                        delim2 = 0

                    s = _build_scan(add_hdr, m, len(scans) + index_start)
                    if s is not None:
                        scans.append(s)

                    # condição de parada do RepeatUntil
                    if delim != 0x80 and delim2 == 0:
                        break
                    # se delimiter == 0 e delimiter2 != 0: avançar até 0x80
                    if delim == 0 and delim2 != 0:
                        r.skip_to(0x80)

                except Exception:
                    break

        else:
            # --- new_scan_type_header: mesmo header, tipo de curva diferente ---
            try:
                r.skip(5)              # Padding(5)
                st_len = r.u8()
                r.skip(2)
                scan_type_str = r.read(st_len - 1).decode('latin-1', errors='replace')

                sub_hdr = dict(main_hdr)
                sub_hdr['scan_type_str'] = scan_type_str

                n_sub = abs(disc)
                for _ in range(n_sub):
                    try:
                        m = _read_measurement_data(r)
                        r.skip(2)
                        s = _build_scan(sub_hdr, m, len(scans) + index_start)
                        if s is not None:
                            scans.append(s)
                    except Exception:
                        break

            except Exception:
                break

    return scans


# ===========================================================================
# API pública
# ===========================================================================

def load_rfb_files(filepaths: list, index_start: int = 0) -> List[Scan]:
    """
    Carrega um ou mais arquivos .rfb e retorna lista unificada de Scan.

    Parâmetros
    ----------
    filepaths   : lista de caminhos de arquivos .rfb
    index_start : índice do primeiro Scan retornado
    """
    scans: List[Scan] = []
    idx = index_start
    for fp in filepaths:
        new = parse_rfb(fp, index_start=idx)
        scans.extend(new)
        idx += len(new)
    return scans

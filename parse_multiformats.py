"""
Parsers unificados para múltiplos formatos de perfis dosimétricos.
Auto-detecta o formato pelo conteúdo + extensão do arquivo.

Formatos suportados
-------------------
  Eclipse w2CAD / Line Profile Export    (.asc, .txt)
  OmniPro RFA300 ASCII BDS (IBA)         (.txt, .asc)
  RayStation Physics Export              (.csv)
  SNC Water Tank                         (.sncxml)
  Standard Imaging DV1D                  (.csv)

Todas as posições de saída estão em **cm**; dose relativa (u.a. ou % conforme
o arquivo).

Nota: arquivos .rfb (OmniPro V6 binário) são tratados pelo módulo
parse_rfb.py e não por este arquivo.

Uso:
    from parse_multiformats import load_auto
    scans = load_auto(["profile.asc", "export.csv"], index_start=0)
"""

import re
import csv
import xml.etree.ElementTree as ET
import numpy as np
from pathlib import Path
from typing import List, Optional

from parse_snctxt import Scan


# ===========================================================================
# Constantes de mapeamento compartilhadas
# ===========================================================================

_CURVETYPE_MAP = {
    "OPD":        "Depth Scan",   # w2CAD depth
    "OPP":        "Profile",      # w2CAD profile (eixo inferido depois)
    "DPR":        "Diagonal",
    "DPT":        "Depth Scan",   # IBA
    "PRO":        "Profile",      # IBA
    "DEPTHDOSE":  "Depth Scan",
    "DEPTH":      "Depth Scan",
    "CROSSLINE":  "Crossline",
    "INLINE":     "Inline",
    "PROFILE":    "Crossline",    # padrão ao usar "Profile" sem eixo
    "DIAGONAL":   "Diagonal",
}

_MODALITY_MAP = {
    "PHO": "Photon", "PHOTON": "Photon", "XRAY": "Photon", "X": "Photon",
    "ELE": "Electron", "ELECTRON": "Electron", "E": "Electron",
}

_DETECTOR_W2CAD = {"CHA": "Chamber", "DIO": "Diode", "DIA": "Diamond"}


# ===========================================================================
# Detecção automática de formato
# ===========================================================================

def _detect_format(path: Path, head: list[str]) -> str:
    """
    Retorna um identificador de formato baseado na extensão e nas primeiras
    linhas do arquivo.

    Retorna: 'w2cad' | 'ibatxt' | 'raystation' | 'sncxml' | 'dv1d' | 'unknown'
    """
    ext = path.suffix.lower()

    # --- SNC Water Tank XML ---
    if ext == '.sncxml':
        return 'sncxml'

    head_joined = '\n'.join(head[:40]).upper()

    # --- Standard Imaging DV1D ---
    if 'DOSEVIEW 1D' in head_joined or 'CHARGE TABLE' in head_joined:
        return 'dv1d'

    # --- RayStation CSV ---
    if '#EXPORTED' in head_joined and 'RAYSTATION' in head_joined:
        return 'raystation'
    if 'CURVETYPE' in head_joined and 'RADIATIONTYPE' in head_joined:
        return 'raystation'

    # --- Eclipse w2CAD / Line Profile ---
    if '%VERSION' in head_joined or '$STOM' in head_joined:
        return 'w2cad'
    if '%BMTY' in head_joined or '%FLSZ' in head_joined:
        return 'w2cad'

    # --- OmniPro RFA300 ASCII BDS ---
    if ':MSR' in head_joined or ':SYS BDS' in head_joined:
        return 'ibatxt'
    if '%SCN' in head_joined or '%BMT' in head_joined:
        return 'ibatxt'

    return 'unknown'


# ===========================================================================
# Construtor de Scan genérico
# ===========================================================================

def _make_scan(index: int, scan_type: str, pos_cm: np.ndarray,
               dose: np.ndarray, beam_type: str = "Photon",
               energy: str = "", field_x: str = "?", field_y: str = "?",
               ssd: str = "?", date: str = "", detector: str = "",
               medium: str = "Water") -> Scan:
    """Monta um Scan com as coordenadas corretas para o scan_type."""
    x = np.zeros_like(pos_cm)
    y = np.zeros_like(pos_cm)
    z = np.zeros_like(pos_cm)
    if scan_type == "Depth Scan":
        z = pos_cm
    elif scan_type == "Inline":
        y = pos_cm
    else:  # Crossline / Diagonal / genérico
        x = pos_cm
    return Scan(index=index, beam_type=beam_type, energy=energy,
                field_x=field_x, field_y=field_y, scan_type=scan_type,
                ssd=ssd, medium=medium, date=date, detector=detector,
                x=x, y=y, z=z, dose=dose)


def _infer_axis(x_arr: np.ndarray, y_arr: np.ndarray) -> str:
    """Infere qual eixo varia mais para classificar Crossline vs Inline."""
    return "Inline" if np.ptp(y_arr) > np.ptp(x_arr) else "Crossline"


# ===========================================================================
# 1. Eclipse w2CAD / Eclipse Line Profile Export
# ===========================================================================

def _parse_w2cad(lines: list[str], index_start: int) -> List[Scan]:
    """
    Parser para o formato Eclipse w2CAD (.asc) e Line Profile Export (.txt).

    Estrutura:
      %VERSION 2
      %BMTY PHO | ELE
      %FLSZ 100.0 * 100.0       (mm)
      %TYPE OPD | OPP | DPR
      %AXIS X | Y
      %SSD  1000.0               (mm)
      %DPTH 100.0                (mm)
      $STOM
      < X_mm  Y_mm  Z_mm  dose >
      ...
      $ENDOM
    """
    scans: List[Scan] = []
    idx = index_start

    # metadados globais (válidos até serem sobrescritos)
    g_meta: dict[str, str] = {}
    cur_meta: dict[str, str] = {}
    data_pts: list[tuple] = []      # (x_mm, y_mm, z_mm, dose)
    in_data = False

    def _flush(meta, pts):
        nonlocal idx
        if len(pts) < 2:
            return
        arr_x = np.array([p[0] for p in pts])
        arr_y = np.array([p[1] for p in pts])
        arr_z = np.array([p[2] for p in pts])
        arr_d = np.array([p[3] for p in pts])

        raw_type = meta.get('TYPE', meta.get('%TYPE', '')).upper()
        scan_type = _CURVETYPE_MAP.get(raw_type, 'Crossline')

        # se "Profile" genérico, inferir pelo eixo ou pela variação
        if scan_type == 'Profile':
            axis = meta.get('AXIS', '').upper()
            if axis == 'X':
                scan_type = 'Crossline'
            elif axis == 'Y':
                scan_type = 'Inline'
            else:
                scan_type = _infer_axis(arr_x, arr_y)

        # coordenada principal em cm
        if scan_type == 'Depth Scan':
            pos_cm = arr_z / 10.0
        elif scan_type == 'Inline':
            pos_cm = arr_y / 10.0
        else:
            pos_cm = arr_x / 10.0

        # modalidade
        raw_mod = meta.get('BMTY', '').upper().split()[0]
        beam_type = _MODALITY_MAP.get(raw_mod, 'Photon')

        # energy
        energy = meta.get('ENERGY', '')

        # campo: %FLSZ "100.0 * 100.0" mm → cm
        flsz = meta.get('FLSZ', '')
        parts_f = re.split(r'[\*x×\s]+', flsz.strip())
        parts_f = [p.strip() for p in parts_f if p.strip()]
        if len(parts_f) >= 2:
            try:
                field_x = f"{float(parts_f[0]) / 10:.4g}"
                field_y = f"{float(parts_f[1]) / 10:.4g}"
            except ValueError:
                field_x = field_y = '?'
        else:
            field_x = field_y = '?'

        # SSD mm → cm
        try:
            ssd_str = f"{float(meta.get('SSD', '1000')) / 10:.4g}"
        except ValueError:
            ssd_str = '?'

        detector = _DETECTOR_W2CAD.get(meta.get('DETY', '').upper(), '')
        date = meta.get('DATE', '')

        scans.append(_make_scan(idx, scan_type, pos_cm, arr_d,
                                beam_type=beam_type, energy=energy,
                                field_x=field_x, field_y=field_y,
                                ssd=ssd_str, date=date, detector=detector))
        idx += 1

    for line in lines:
        s = line.strip()

        if not s:
            continue

        if s.startswith('$STOM') or s.startswith('$FIELD'):
            in_data = True
            data_pts = []
            continue

        if s.startswith('$ENDOM') or s.startswith('$ENOM'):
            in_data = False
            _flush({**g_meta, **cur_meta}, data_pts)
            data_pts = []
            cur_meta = {}
            continue

        if in_data:
            # linhas de dado: < X Y Z dose >  (os <> são opcionais)
            cleaned = s.lstrip('<').rstrip('>')
            nums = cleaned.split()
            if len(nums) >= 4:
                try:
                    data_pts.append(tuple(float(n) for n in nums[:4]))
                except ValueError:
                    pass
            continue

        # linhas de metadados: %KEY value
        if s.startswith('%'):
            # remover o '%'
            rest = s[1:].strip()
            key, _, val = rest.partition(' ')
            key = key.strip().upper()
            val = val.strip()
            # metadados que se aplicam a todos os profiles
            if key in ('VERSION', 'DATE', 'DETY', 'BMTY'):
                g_meta[key] = val
            else:
                cur_meta[key] = val

    return scans


# ===========================================================================
# 2. OmniPro RFA300 ASCII BDS (IBA)
# ===========================================================================
#
# Estrutura:
#   :MSR <n>         — número total de perfis
#   :SYS BDS ...     — identificador do sistema
#
#   %VNR 2.0
#   %MOD PHO
#   %SCN DPT | PRO | DIA
#   %BMT PHO 6.0     — modality + energy
#   %FSZ 100.0 100.0 — field size mm
#   %SSD 1000.0      — mm
#   %PRD 100.0       — profundidade mm (para perfis laterais)
#   %PTS 200         — número de pontos
#   =                — marca início dos dados
#   X_mm Y_mm Z_mm   dose
#   ...
#   :EPF             — fim do perfil
#   %VNR ...         — início do próximo perfil
#   :EOM             — fim do arquivo
# ===========================================================================

def _parse_ibatxt(lines: list[str], index_start: int) -> List[Scan]:
    """Parser para OmniPro RFA300 ASCII BDS."""
    scans: List[Scan] = []
    idx = index_start

    meta: dict[str, str] = {}
    pts: list[tuple] = []
    in_data = False

    def _flush():
        nonlocal idx
        if len(pts) < 2:
            return

        arr_x = np.array([p[0] for p in pts])
        arr_y = np.array([p[1] for p in pts])
        arr_z = np.array([p[2] for p in pts])
        arr_d = np.array([p[3] for p in pts])

        raw_scn = meta.get('SCN', '').upper()
        scan_type = _CURVETYPE_MAP.get(raw_scn, 'Crossline')
        if scan_type == 'Profile':
            scan_type = _infer_axis(arr_x, arr_y)

        if scan_type == 'Depth Scan':
            pos_cm = arr_z / 10.0
        elif scan_type == 'Inline':
            pos_cm = arr_y / 10.0
        else:
            pos_cm = arr_x / 10.0

        # %BMT <MOD> <energy>
        bmt = meta.get('BMT', meta.get('MOD', '')).split()
        raw_mod = bmt[0].upper() if bmt else ''
        beam_type = _MODALITY_MAP.get(raw_mod, 'Photon')
        energy = bmt[1].rstrip('.0').rstrip('.') if len(bmt) >= 2 else meta.get('ENERGY', '')

        # campo %FSZ mm → cm
        fsz = meta.get('FSZ', meta.get('FLD', '')).split()
        if len(fsz) >= 2:
            try:
                field_x = f"{float(fsz[0]) / 10:.4g}"
                field_y = f"{float(fsz[1]) / 10:.4g}"
            except ValueError:
                field_x = field_y = '?'
        else:
            field_x = field_y = '?'

        try:
            ssd_str = f"{float(meta.get('SSD', '1000')) / 10:.4g}"
        except ValueError:
            ssd_str = '?'

        date = meta.get('DAT', '')
        scans.append(_make_scan(idx, scan_type, pos_cm, arr_d,
                                beam_type=beam_type, energy=energy,
                                field_x=field_x, field_y=field_y,
                                ssd=ssd_str, date=date))
        idx += 1

    for line in lines:
        s = line.strip()

        if not s:
            continue

        upper = s.upper()

        if upper.startswith(':MSR') or upper.startswith(':SYS'):
            continue

        if upper == ':EPF' or upper == ':EOM':
            if pts:
                _flush()
            meta = {}
            pts = []
            in_data = False
            continue

        if s == '=':
            in_data = True
            continue

        if in_data:
            nums = s.split()
            if len(nums) >= 4:
                try:
                    pts.append(tuple(float(n) for n in nums[:4]))
                except ValueError:
                    pass
            continue

        # metadados %KEY value1 [value2 ...]
        if s.startswith('%'):
            rest = s[1:].strip()
            key, _, val = rest.partition(' ')
            meta[key.strip().upper()] = val.strip()

    # flushar perfil sem :EPF final (arquivos mal terminados)
    if pts:
        _flush()

    return scans


# ===========================================================================
# 3. RayStation Physics Export (.csv)
# ===========================================================================
#
# Estrutura:
#   #Exported: RayStation 10B
#   #Time: 01/01/2024 12:00:00
#   #Machine Name: MACHINE
#   ...linhas # de comentário...
#
#   energy; 6
#   SSD; 1000           (mm)
#   Fieldsize; X1; X2; Y1; Y2   (mm)
#   CurveType; Crossline | Inline | Depth
#   RadiationType; Photon | Electron
#   Quantity; Relative Dose
#   StartPoint; X; Y; Z  (mm)
#   pos_mm; dose
#   ...
#   [linha em branco separa perfis]
# ===========================================================================

def _parse_raystation(lines: list[str], index_start: int) -> List[Scan]:
    """Parser para RayStation Physics Export CSV."""
    scans: List[Scan] = []
    idx = index_start

    meta: dict[str, str] = {}
    pts: list[tuple[float, float]] = []   # (pos_mm, dose)

    _RS_CURVETYPE = {
        'CROSSLINE': 'Crossline',
        'CROSSPLANE': 'Crossline',
        'INLINE':    'Inline',
        'INPLANE':   'Inline',
        'DEPTH':     'Depth Scan',
        'PDD':       'Depth Scan',
    }

    def _flush():
        nonlocal idx
        if len(pts) < 2:
            return
        raw_ct = meta.get('CURVETYPE', '').upper().replace(' ', '')
        scan_type = _RS_CURVETYPE.get(raw_ct, 'Crossline')

        pos_cm = np.array([p[0] for p in pts]) / 10.0
        dose   = np.array([p[1] for p in pts])

        raw_rad = meta.get('RADIATIONTYPE', 'PHOTON').upper()
        beam_type = _MODALITY_MAP.get(raw_rad[:3], 'Photon')

        energy = meta.get('ENERGY', '').rstrip('.0').rstrip('.')

        # Fieldsize: "X1; X2; Y1; Y2" mm → largura em cm
        fs_raw = meta.get('FIELDSIZE', '')
        fs_parts = [p.strip() for p in fs_raw.split(';') if p.strip()]
        if len(fs_parts) >= 4:
            try:
                fx = abs(float(fs_parts[0])) + abs(float(fs_parts[1]))
                fy = abs(float(fs_parts[2])) + abs(float(fs_parts[3]))
                field_x = f"{fx / 10:.4g}"
                field_y = f"{fy / 10:.4g}"
            except ValueError:
                field_x = field_y = '?'
        else:
            field_x = field_y = '?'

        try:
            ssd_str = f"{float(meta.get('SSD', '1000')) / 10:.4g}"
        except ValueError:
            ssd_str = '?'

        scans.append(_make_scan(idx, scan_type, pos_cm, dose,
                                beam_type=beam_type, energy=energy,
                                field_x=field_x, field_y=field_y,
                                ssd=ssd_str))
        idx += 1

    for line in lines:
        s = line.strip()

        if not s:
            # linha em branco separa perfis
            if pts:
                _flush()
                meta = {}
                pts = []
            continue

        if s.startswith('#'):
            continue

        # separar por ';'
        parts = [p.strip() for p in s.split(';')]

        if len(parts) >= 2:
            key = parts[0].strip()
            key_upper = key.upper().replace(' ', '').replace(':', '')

            # tentar parse como par posição; dose (dados numéricos)
            try:
                pos_val  = float(parts[0])
                dose_val = float(parts[1])
                pts.append((pos_val, dose_val))
                continue
            except ValueError:
                pass

            # linha de metadado
            if key_upper == 'ENERGY':
                meta['ENERGY'] = parts[1]
            elif key_upper == 'SSD':
                meta['SSD'] = parts[1]
            elif key_upper in ('FIELDSIZE', 'FIELDSIZEX', 'COLLIMATION'):
                meta['FIELDSIZE'] = '; '.join(parts[1:])
            elif key_upper == 'CURVETYPE':
                meta['CURVETYPE'] = parts[1]
            elif key_upper == 'RADIATIONTYPE':
                meta['RADIATIONTYPE'] = parts[1]
            elif key_upper == 'STARTPOINT':
                pass  # usado somente para orientação
            else:
                meta[key_upper] = '; '.join(parts[1:])

    if pts:
        _flush()

    return scans


# ===========================================================================
# 4. SNC Water Tank (.sncxml)
# ===========================================================================
#
# Formato XML.  As marcações podem vir com namespace (ex.: xmlns="...")
# portanto usamos localname para comparações.
#
# Estrutura esperada (simplificada):
#   <ScanDatabase>
#     <RadiationDevices>
#       <RadiationDevice Name="..." />
#     </RadiationDevices>
#     <BeamDataSets>
#       <BeamDataSet>
#         <BeamData Energy="..." BeamType="Photon|Electron" ...>
#           <FieldSize>
#             <JawLeft>-50</JawLeft>  <JawRight>50</JawRight>
#             <JawUpper>-50</JawUpper><JawLower>50</JawLower>
#           </FieldSize>
#           <Scans>
#             <Scan SSD="..." ScanType="..." ...>
#               <Datapoints>
#                 <Datapoint X="0" Y="0" Z="1.0" Dose="0.1"/>
#               </Datapoints>
#             </Scan>
#           </Scans>
#         </BeamData>
#       </BeamDataSet>
#     </BeamDataSets>
#   </ScanDatabase>
#
# Unidades: coordenadas em cm, dose relativa.
# ===========================================================================

def _localname(tag: str) -> str:
    """Extrai o nome local de um tag XML descartando o namespace '{uri}'."""
    return tag.split('}', 1)[-1] if '}' in tag else tag


def _xml_find_all(element, local: str) -> list:
    """Encontra todos os filhos com dado localname (ignora namespace)."""
    return [c for c in element.iter() if _localname(c.tag) == local]


def _xml_text(element, local: str, default: str = '') -> str:
    """Retorna o texto do primeiro filho com dado localname."""
    els = _xml_find_all(element, local)
    if els and els[0].text:
        return els[0].text.strip()
    return default


def _parse_sncxml(filepath: str, index_start: int) -> List[Scan]:
    """Parser para SNC Water Tank XML (.sncxml)."""
    scans: List[Scan] = []
    idx = index_start

    try:
        tree = ET.parse(filepath)
        root = tree.getroot()
    except ET.ParseError as e:
        raise ValueError(f"Erro ao ler XML: {e}")

    _SCANTYPE_XML = {
        'DEPTH': 'Depth Scan', 'PDD': 'Depth Scan', 'DEPTHDOSE': 'Depth Scan',
        'CROSSLINE': 'Crossline', 'CROSSPLANE': 'Crossline',
        'INLINE': 'Inline', 'INPLANE': 'Inline',
        'DIAGONAL': 'Diagonal',
    }

    # percorrer todos os BeamData
    for beam_el in _xml_find_all(root, 'BeamData'):
        # energia e tipo
        energy  = beam_el.get('Energy', _xml_text(beam_el, 'Energy', ''))
        raw_bt  = beam_el.get('BeamType', _xml_text(beam_el, 'BeamType', 'Photon'))
        beam_type = _MODALITY_MAP.get(raw_bt.upper()[:3], raw_bt.capitalize())

        # campo (jaw em mm → cm)
        def _jaw(name: str) -> Optional[float]:
            v = beam_el.get(name, _xml_text(beam_el, name, ''))
            try: return abs(float(v))
            except (ValueError, TypeError): return None

        jaw_l = _jaw('JawLeft')  or _jaw('MLCLeft')
        jaw_r = _jaw('JawRight') or _jaw('MLCRight')
        jaw_u = _jaw('JawUpper') or _jaw('MLCUpper')
        jaw_d = _jaw('JawLower') or _jaw('MLCLower')

        if jaw_l is not None and jaw_r is not None:
            fx = (jaw_l + jaw_r) / 10.0
        else:
            fx_raw = beam_el.get('FieldSizeX', _xml_text(beam_el, 'FieldSizeX', ''))
            try: fx = float(fx_raw) / 10.0
            except (ValueError, TypeError): fx = float('nan')
        if jaw_u is not None and jaw_d is not None:
            fy = (jaw_u + jaw_d) / 10.0
        else:
            fy_raw = beam_el.get('FieldSizeY', _xml_text(beam_el, 'FieldSizeY', ''))
            try: fy = float(fy_raw) / 10.0
            except (ValueError, TypeError): fy = fx

        field_x = f"{fx:.4g}" if not np.isnan(fx) else '?'
        field_y = f"{fy:.4g}" if not np.isnan(fy) else '?'

        # percorrer Scan nodes
        for scan_el in _xml_find_all(beam_el, 'Scan'):
            raw_type = scan_el.get('ScanType', scan_el.get('ScanAxis',
                                   _xml_text(scan_el, 'ScanType', ''))).upper()
            scan_type = _SCANTYPE_XML.get(raw_type.replace(' ', ''), 'Crossline')

            raw_ssd = scan_el.get('SSD', _xml_text(scan_el, 'SSD', '1000'))
            try:
                ssd_str = f"{float(raw_ssd) / 10:.4g}"
            except ValueError:
                ssd_str = '?'

            date = scan_el.get('MeasurementDate',
                               _xml_text(scan_el, 'MeasurementDate', ''))

            # pontos: <Datapoint X Y Z Dose /> em cm
            xs, ys, zs, ds = [], [], [], []
            for dp in _xml_find_all(scan_el, 'Datapoint'):
                try:
                    xs.append(float(dp.get('X', '0')))
                    ys.append(float(dp.get('Y', '0')))
                    zs.append(float(dp.get('Z', '0')))
                    ds.append(float(dp.get('Dose', dp.get('RelativeDose', '0'))))
                except (ValueError, TypeError):
                    pass

            if len(ds) < 2:
                continue

            arr_x = np.array(xs)
            arr_y = np.array(ys)
            arr_z = np.array(zs)
            arr_d = np.array(ds)

            # garantir dose positiva
            if np.mean(arr_d) < 0:
                arr_d = -arr_d

            if scan_type == 'Depth Scan':
                pos_cm = arr_z
            elif scan_type == 'Inline':
                pos_cm = arr_y
            else:
                pos_cm = arr_x

            scans.append(_make_scan(idx, scan_type, pos_cm, arr_d,
                                    beam_type=beam_type, energy=energy,
                                    field_x=field_x, field_y=field_y,
                                    ssd=ssd_str, date=date))
            idx += 1

    return scans


# ===========================================================================
# 5. Standard Imaging DoseView 1D (.csv)
# ===========================================================================
#
# Estrutura:
#   "DoseView 1D Software Export"
#   "timestamp","version","color","OS",...
#   ...
#   "Charge Table"
#   [header row]
#   row0, row1, ..., depth(col9), ch1(col10), ch2(col11), ..., duration(col21)
# ===========================================================================

def _parse_dv1d(lines: list[str], index_start: int) -> List[Scan]:
    """Parser simplificado para Standard Imaging DoseView 1D CSV."""
    scans: List[Scan] = []
    idx = index_start

    in_table = False
    skip_header = False
    data_rows: list[list[str]] = []

    for line in lines:
        s = line.strip().strip('"')
        if 'Charge Table' in s:
            in_table = True
            skip_header = True
            continue
        if not in_table:
            continue
        if skip_header:
            skip_header = False
            continue
        if not s:
            continue

        try:
            reader = csv.reader([line])
            parts = next(reader)
        except Exception:
            parts = line.split(',')

        parts = [p.strip().strip('"') for p in parts]
        if len(parts) > 11:
            data_rows.append(parts)

    if not data_rows:
        return scans

    depths, ch1, ch2, dur = [], [], [], []
    for row in data_rows:
        try:
            depths.append(float(row[9]))
            ch1.append(float(row[10]))
            ch2.append(float(row[11]))
            dur_val = float(row[21]) if len(row) > 21 and row[21] else 1.0
            dur.append(max(dur_val, 1e-6))
        except (ValueError, IndexError):
            continue

    if not depths:
        return scans

    arr_depth = np.array(depths)
    arr_ch1   = np.array(ch1)
    arr_dur   = np.array(dur)

    # normalizar por duração e garantir positivo
    signal = arr_ch1 / arr_dur
    if np.mean(signal) < 0:
        signal = -signal

    # se profundidades variam → PDD
    if np.ptp(arr_depth) > 0.5:
        # profundidades em mm → cm
        pos_cm = arr_depth / 10.0
        scan_type = 'Depth Scan'
    else:
        pos_cm = arr_depth / 10.0
        scan_type = 'Crossline'

    scans.append(_make_scan(idx, scan_type, pos_cm, signal))
    return scans


# ===========================================================================
# API pública
# ===========================================================================

def parse_file(filepath: str, index_start: int = 0) -> List[Scan]:
    """
    Lê um arquivo no formato suportado e retorna lista de Scan.

    Detecta automaticamente o formato pelo conteúdo + extensão.

    Lança ValueError se o formato não for reconhecido.
    """
    path = Path(filepath)

    # XML não precisa ler como texto primeiro
    if path.suffix.lower() == '.sncxml':
        return _parse_sncxml(filepath, index_start)

    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()

    head = [ln.rstrip() for ln in lines[:50]]
    fmt  = _detect_format(path, head)

    if fmt == 'w2cad':
        return _parse_w2cad(lines, index_start)
    elif fmt == 'ibatxt':
        return _parse_ibatxt(lines, index_start)
    elif fmt == 'raystation':
        return _parse_raystation(lines, index_start)
    elif fmt == 'dv1d':
        return _parse_dv1d(lines, index_start)
    else:
        raise ValueError(
            f"Formato não reconhecido: {path.name}\n"
            "Formatos suportados: w2CAD (.asc/.txt), OmniPro RFA300 (.txt/.asc), "
            "RayStation (.csv), SNC Water Tank (.sncxml), DoseView 1D (.csv)"
        )


def load_auto(filepaths: list[str], index_start: int = 0) -> List[Scan]:
    """
    Carrega uma lista de arquivos com auto-detecção de formato.

    Parâmetros
    ----------
    filepaths   : lista de caminhos
    index_start : índice do primeiro Scan retornado

    Retorna lista unificada de Scan, preservando a ordem dos arquivos.
    Lança ValueError se algum arquivo tiver formato não reconhecido.
    """
    scans: List[Scan] = []
    idx = index_start
    for fp in filepaths:
        new = parse_file(fp, index_start=idx)
        scans.extend(new)
        idx += len(new)
    return scans

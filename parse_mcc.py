"""
Parser para arquivos .mcc do PTW Mephisto (BeamScan / Mephisto MC2).

Formato: blocos BEGIN_SCAN / END_SCAN aninhados em BEGIN_SCAN_DATA /
END_SCAN_DATA.  Cada scan contém pares chave = valor e uma seção
BEGIN_DATA / END_DATA com linhas "pos(mm)  dose [ref]".

SCAN_CURVETYPE → scan_type do Scan:
  PDD                 → "Depth Scan"
  INPLANE_PROFILE     → "Inline"
  CROSSPLANE_PROFILE  → "Crossline"
  DIAGONAL_PROFILE    → "Diagonal"

Posições convertidas de mm para cm para compatibilidade com Scan.

Uso como módulo:
    from parse_mcc import load_mcc_files
    scans = load_mcc_files(["arquivo.mcc"], index_start=0)
"""

import numpy as np
from pathlib import Path
from typing import List

from parse_snctxt import Scan


# ---------------------------------------------------------------------------
# Mapeamentos
# ---------------------------------------------------------------------------

_CURVETYPE_MAP: dict[str, str] = {
    "PDD":                "Depth Scan",
    "INPLANE_PROFILE":    "Inline",
    "CROSSPLANE_PROFILE": "Crossline",
    "DIAGONAL_PROFILE":   "Diagonal",
    "CROSSPLANE":         "Crossline",
    "INPLANE":            "Inline",
}

_MODALITY_MAP: dict[str, str] = {
    "XRAY":     "Photon",
    "PHOTON":   "Photon",
    "ELECTRON": "Electron",
    "ELECTRONS":"Electron",
}


# ---------------------------------------------------------------------------
# Parser de baixo nível
# ---------------------------------------------------------------------------

def _parse_meta_value(lines: list[str], start: int, stop: int) -> dict[str, str]:
    """Lê pares CHAVE = VALOR entre as linhas start e stop (exclusive)."""
    meta: dict[str, str] = {}
    for ln in lines[start:stop]:
        s = ln.strip()
        if '=' in s and not s.startswith('#'):
            k, _, v = s.partition('=')
            meta[k.strip().upper()] = v.strip()
    return meta


def parse_mcc(filepath: str, index_start: int = 0) -> List[Scan]:
    """
    Lê um arquivo .mcc e retorna lista de Scan.

    Parâmetros
    ----------
    filepath    : caminho do arquivo .mcc
    index_start : índice base para os Scan retornados
    """
    path = Path(filepath)
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()

    scans: List[Scan] = []
    idx = index_start
    i = 0
    n = len(lines)

    while i < n:
        ln = lines[i].strip()

        if ln == 'BEGIN_SCAN':
            meta: dict[str, str] = {}
            raw_pts: list[tuple[float, float]] = []
            in_data = False
            i += 1

            while i < n:
                s = lines[i].strip()

                if s == 'END_SCAN':
                    break

                if s == 'BEGIN_DATA':
                    in_data = True
                elif s == 'END_DATA':
                    in_data = False
                elif in_data:
                    if s and not s.startswith('#'):
                        parts = s.split()
                        if len(parts) >= 2:
                            try:
                                raw_pts.append((float(parts[0]), float(parts[1])))
                            except ValueError:
                                pass
                elif '=' in s and not s.startswith('#'):
                    k, _, v = s.partition('=')
                    meta[k.strip().upper()] = v.strip()

                i += 1

            if raw_pts:
                scans.append(_build_scan(meta, raw_pts, idx, path.stem))
                idx += 1

        i += 1

    return scans


# ---------------------------------------------------------------------------
# Construção do Scan
# ---------------------------------------------------------------------------

def _build_scan(meta: dict[str, str], pts: list[tuple], index: int,
                filename: str) -> Scan:
    """Monta um Scan a partir dos metadados e pontos brutos de um bloco MCC."""

    # tipo de curva
    curvetype = meta.get('SCAN_CURVETYPE', '').upper()
    scan_type = _CURVETYPE_MAP.get(curvetype, curvetype or 'Unknown')

    # modalidade e energia
    modality  = meta.get('MODALITY', meta.get('BEAM_TYPE', '')).upper()
    beam_type = _MODALITY_MAP.get(modality, modality.capitalize() or 'Photon')
    energy    = meta.get('ENERGY', meta.get('ENERGY_MEV', ''))

    # campo: MCC armazena em mm → converter para cm
    def _mm_to_cm_str(key1: str, key2: str, fallback: str = '?') -> str:
        raw = meta.get(key1, meta.get(key2, ''))
        try:
            return f"{float(raw) / 10:.4g}"
        except (ValueError, TypeError):
            return fallback

    field_x = _mm_to_cm_str('FIELD_INPLANE',    'FIELD_X', '?')
    field_y = _mm_to_cm_str('FIELD_CROSSPLANE',  'FIELD_Y', '?')

    # SSD: mm → cm
    ssd_str = _mm_to_cm_str('SSD', 'ISOCENTER_DIST', '?')
    if ssd_str == '?':
        # alguns arquivos usam SAD em mm
        ssd_str = _mm_to_cm_str('SAD', 'SDD', '?')

    detector = meta.get('DETECTOR', '')

    # arrays (posição em mm → cm)
    pos_cm = np.array([p[0] for p in pts]) / 10.0
    dose   = np.array([p[1] for p in pts])

    x = np.zeros_like(pos_cm)
    y = np.zeros_like(pos_cm)
    z = np.zeros_like(pos_cm)

    if scan_type == 'Depth Scan':
        z = pos_cm
    elif scan_type == 'Inline':
        y = pos_cm
    elif scan_type in ('Crossline', 'Diagonal'):
        x = pos_cm

    return Scan(
        index    = index,
        beam_type= beam_type,
        energy   = energy,
        field_x  = field_x,
        field_y  = field_y,
        scan_type= scan_type,
        ssd      = ssd_str,
        medium   = 'Water',
        date     = meta.get('DATE', ''),
        detector = detector,
        x = x, y = y, z = z,
        dose = dose,
    )


# ---------------------------------------------------------------------------
# API pública
# ---------------------------------------------------------------------------

def load_mcc_files(filepaths: list[str], index_start: int = 0) -> List[Scan]:
    """
    Carrega um ou mais arquivos .mcc e retorna lista unificada de Scan.

    Parâmetros
    ----------
    filepaths   : lista de caminhos
    index_start : índice do primeiro Scan
    """
    scans: List[Scan] = []
    idx = index_start
    for fp in filepaths:
        new = parse_mcc(fp, index_start=idx)
        scans.extend(new)
        idx += len(new)
    return scans

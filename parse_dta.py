"""
Parser para arquivos de perfil 1D no formato DTA (IBA/Scanditronix e variantes).

Formatos suportados
-------------------
1. **IBA OmniPro Accept ASCII** (.dta / .asc): blocos delimitados por
   ":MSR <n> <type>" (início) e ":EOM" (fim), com chave-valor nos cabeçalhos
   e colunas posição(mm) + dose nas linhas de dados.

2. **Texto plano (fallback)**: um único perfil por arquivo, metadados em
   linhas "# CHAVE: VALOR" ou "% CHAVE = VALOR", seguido de duas colunas
   numéricas (posição mm, dose).

Scan_type inferido:
  - SCAN_TYPE explícito no cabeçalho  (PDD, INPLANE, CROSSLINE, PROFILE…)
  - Heurística: se todas as posições ≥ 0  → "Depth Scan"
                caso contrário            → "Crossline" (perfil genérico)

Posições convertidas de mm para cm para compatibilidade com Scan.

Uso como módulo:
    from parse_dta import load_dta_files
    scans = load_dta_files(["arquivo.dta"], index_start=0)
"""

import re
import numpy as np
from pathlib import Path
from typing import List

from parse_snctxt import Scan


# ---------------------------------------------------------------------------
# Mapeamentos
# ---------------------------------------------------------------------------

_SCANTYPE_MAP: dict[str, str] = {
    "PDD":        "Depth Scan",
    "DEPTH":      "Depth Scan",
    "DD":         "Depth Scan",
    "DOSEPROFILE":"Crossline",
    "PROFILE":    "Crossline",
    "CROSSLINE":  "Crossline",
    "CROSSPLANE": "Crossline",
    "INLINE":     "Inline",
    "INPLANE":    "Inline",
    "DIAGONAL":   "Diagonal",
}

_MODALITY_MAP: dict[str, str] = {
    "XRAY":    "Photon",
    "PHOTON":  "Photon",
    "X":       "Photon",
    "E":       "Electron",
    "ELECTRON":"Electron",
}

# Padrão de início de bloco IBA:  :MSR  3  DoseProfile
_IBA_BLOCK_RE = re.compile(r'^:MSR\s+(\d+)\s+(\S+)', re.IGNORECASE)
# Par chave-valor IBA:  \tEnergy\t=\t6.0
_IBA_KV_RE    = re.compile(r'^\t?(\w[\w\s]*?)\s*=\s*(.+)$')


# ---------------------------------------------------------------------------
# Detecção de formato
# ---------------------------------------------------------------------------

def _looks_like_iba(lines: list[str]) -> bool:
    """Retorna True se o arquivo parece ser no formato IBA OmniPro."""
    for ln in lines[:30]:
        if ln.strip().startswith(':') or ln.strip().startswith('%'):
            return True
    return False


# ---------------------------------------------------------------------------
# Parser IBA OmniPro ASCII
# ---------------------------------------------------------------------------

def _parse_iba(lines: list[str], index_start: int) -> List[Scan]:
    """Lê formato IBA OmniPro Accept ASCII com blocos :MSR … :EOM."""
    scans: List[Scan] = []
    idx = index_start
    i = 0
    n = len(lines)

    while i < n:
        ln = lines[i].strip()
        m = _IBA_BLOCK_RE.match(ln)
        if m:
            profile_type = m.group(2).upper()  # ex.: "DoseProfile", "PDD"
            meta: dict[str, str] = {}
            pts: list[tuple[float, float]] = []
            in_data = False
            i += 1

            while i < n:
                s = lines[i].rstrip('\n')
                stripped = s.strip()

                if stripped.upper() == ':EOM':
                    break

                # início da seção de dados (linha vazia ou "=")
                if stripped == '=':
                    in_data = True
                    i += 1
                    continue

                if in_data:
                    parts = stripped.split()
                    if len(parts) >= 2:
                        try:
                            pts.append((float(parts[0]), float(parts[1])))
                        except ValueError:
                            pass
                else:
                    # tentar extrair chave = valor (formato IBA usa tabulação)
                    kv = _IBA_KV_RE.match(s)
                    if kv:
                        meta[kv.group(1).strip().upper()] = kv.group(2).strip()
                    elif '\t' in s:
                        parts = [p.strip() for p in s.split('\t') if p.strip()]
                        if len(parts) == 2:
                            meta[parts[0].upper()] = parts[1]

                i += 1

            if pts:
                scans.append(_build_scan(meta, profile_type, pts, idx))
                idx += 1

        i += 1

    return scans


# ---------------------------------------------------------------------------
# Parser de texto plano (fallback)
# ---------------------------------------------------------------------------

_COMMENT_RE = re.compile(r'^[#%;!\*]')
_KV_RE      = re.compile(
    r'^[#%;!\*]?\s*([A-Za-z_][\w\s]*)[\s:=]+(.+)$'
)

def _parse_plain(lines: list[str], filepath: str, index_start: int) -> List[Scan]:
    """
    Fallback: um arquivo = (possivelmente vários) blocos de dados separados
    por linhas em branco ou cabeçalhos.  Cada bloco com ≥ 3 pontos numéricos
    vira um Scan.
    """
    scans: List[Scan] = []
    idx = index_start

    # coletar metadados globais (linhas de comentário no início)
    global_meta: dict[str, str] = {}
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        if _COMMENT_RE.match(s):
            kv = _KV_RE.match(s)
            if kv:
                global_meta[kv.group(1).strip().upper()] = kv.group(2).strip()
        else:
            break  # chegou nos dados

    # agrupar blocos numéricos separados por linhas em branco / comentário
    current_pts: list[tuple[float, float]] = []
    current_meta: dict[str, str] = dict(global_meta)

    def _flush():
        nonlocal current_pts, current_meta, idx
        if len(current_pts) >= 3:
            profile_type = ''  # sem informação explícita
            scans.append(_build_scan(current_meta, profile_type, current_pts, idx))
            idx += 1
        current_pts = []
        current_meta = dict(global_meta)

    for ln in lines:
        s = ln.strip()
        if not s:
            if current_pts:
                _flush()
            continue

        if _COMMENT_RE.match(s):
            kv = _KV_RE.match(s)
            if kv:
                current_meta[kv.group(1).strip().upper()] = kv.group(2).strip()
            continue

        parts = s.split()
        if len(parts) >= 2:
            try:
                current_pts.append((float(parts[0]), float(parts[1])))
                continue
            except ValueError:
                pass

        # linha não numérica dentro de bloco → novo cabeçalho
        if current_pts:
            _flush()

    if current_pts:
        _flush()

    return scans


# ---------------------------------------------------------------------------
# Construção do Scan
# ---------------------------------------------------------------------------

def _infer_scan_type(meta: dict[str, str], profile_type: str,
                     pos_mm: np.ndarray) -> str:
    """Determina scan_type a partir de metadados ou heurística de posição."""

    # 1. tipo explícito no cabeçalho
    for key in ('SCAN_TYPE', 'SCANTYPE', 'TYPE', 'MEASUREMENT_TYPE',
                'PROFILETYPE', 'CURVETYPE'):
        raw = meta.get(key, '').upper().replace(' ', '')
        if raw in _SCANTYPE_MAP:
            return _SCANTYPE_MAP[raw]

    # 2. tipo vindo do bloco IBA (:MSR n <type>)
    pt = profile_type.upper().replace(' ', '')
    if pt in _SCANTYPE_MAP:
        return _SCANTYPE_MAP[pt]

    # 3. heurística de posição
    if len(pos_mm) > 0 and np.all(pos_mm >= -0.5):
        return 'Depth Scan'

    return 'Crossline'


def _build_scan(meta: dict[str, str], profile_type: str,
                pts: list[tuple], index: int) -> Scan:
    """Constrói um Scan a partir dos metadados e pontos de um bloco DTA."""

    pos_mm = np.array([p[0] for p in pts])
    dose   = np.array([p[1] for p in pts])
    pos_cm = pos_mm / 10.0

    scan_type = _infer_scan_type(meta, profile_type, pos_mm)

    # modalidade / energia
    modality  = meta.get('MODALITY', meta.get('RADIATION', meta.get('BEAM', ''))).upper()
    beam_type = _MODALITY_MAP.get(modality.split()[0] if modality else '', 'Photon')
    energy    = meta.get('ENERGY', meta.get('ENERGY_MEV', meta.get('ENERGY_MEV_NOM', '')))
    # limpar "MeV" ou "MV" se vier junto
    energy = re.sub(r'[A-Za-z]', '', energy).strip()

    # campo: converter mm → cm
    def _try_cm(key1: str, key2: str = '') -> str:
        raw = meta.get(key1, meta.get(key2, '')) if key2 else meta.get(key1, '')
        try:
            return f"{float(raw) / 10:.4g}"
        except (ValueError, TypeError):
            return raw or '?'

    field_x = _try_cm('FIELD_INPLANE',   'FIELD_SIZE_X')
    field_y = _try_cm('FIELD_CROSSPLANE','FIELD_SIZE_Y')
    if field_x == '?' and field_y == '?':
        fs = _try_cm('FIELD_SIZE', 'FIELDSIZE')
        field_x = field_y = fs

    ssd_str = _try_cm('SSD', 'SID')
    if ssd_str == '?':
        ssd_str = _try_cm('ISOCENTER_DIST', 'ISODIST')

    detector = meta.get('DETECTOR', meta.get('DETECTOR_TYPE', ''))
    date     = meta.get('DATE', meta.get('MEAS_DATE', ''))

    # coordenadas
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
        date     = date,
        detector = detector,
        x = x, y = y, z = z,
        dose = dose,
    )


# ---------------------------------------------------------------------------
# API pública
# ---------------------------------------------------------------------------

def parse_dta(filepath: str, index_start: int = 0) -> List[Scan]:
    """
    Lê um arquivo .dta (ou .asc) e retorna lista de Scan.

    Detecta automaticamente o formato IBA OmniPro ou texto plano.

    Parâmetros
    ----------
    filepath    : caminho do arquivo
    index_start : índice base para os Scan retornados
    """
    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()

    if _looks_like_iba(lines):
        return _parse_iba(lines, index_start)
    else:
        return _parse_plain(lines, filepath, index_start)


def load_dta_files(filepaths: list[str], index_start: int = 0) -> List[Scan]:
    """
    Carrega um ou mais arquivos .dta / .asc e retorna lista unificada de Scan.

    Parâmetros
    ----------
    filepaths   : lista de caminhos
    index_start : índice do primeiro Scan
    """
    scans: List[Scan] = []
    idx = index_start
    for fp in filepaths:
        new = parse_dta(fp, index_start=idx)
        scans.extend(new)
        idx += len(new)
    return scans

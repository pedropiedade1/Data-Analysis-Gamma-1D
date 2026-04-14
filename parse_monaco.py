"""
Parser para arquivos de exportacao de dose do Monaco (Elekta TPS) no formato .ALL.

O arquivo .ALL contem uma grade 2D de dose (cGy) exportada num plano coronal ou
transversal do fantoma.  A partir da grade sao extraidos perfis 1D compatíveis
com o dataclass Scan do parse_snctxt.

Convencao de coordenadas:
  Arquivo Coronal  (PlaneDesc = "C: X.XX cm"):
      grade X × Y, onde X = crossline e Y = inline.
      Profundidade fixada no valor do arquivo (ex.: 50.00 mm = 5 cm).
      Extrai: Crossline em Y=0  e  Inline em X=0.

  Arquivo Transversal (PlaneDesc = "T: X.XX cm"):
      grade X × Y_monaco, onde X = crossline e Y_monaco = eixo de profundidade.
      Linhas superiores = Y alto (mais raso); linhas inferiores = Y baixo (mais fundo).
      Para que profundidade aumente com o índice, usamos z = -Y_monaco.
      Extrai: PDP (Depth Scan) em X=0.

Uso como modulo:
    from parse_monaco import load_monaco_files
    scans = load_monaco_files(["arquivo.ALL", ...], index_start=100)
"""

import re
import os
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Tuple


# ---------------------------------------------------------------------------
# Estrutura de dados Monaco
# ---------------------------------------------------------------------------

@dataclass
class MonacoFile:
    """Grade de dose 2D lida de um arquivo Monaco .ALL."""
    filepath: str
    patient_id: str
    plan_name: str        # ex.: "CCC1", "MontC05"
    algorithm: str        # ex.: "CCC", "MonteCarlo"
    grid_label: str       # ex.: "1mm", "0.5mm"
    plane_type: str       # "Coronal" ou "Transversal"
    position_mm: float    # profundidade (Coronal) ou posicao lateral (Transversal) em mm
    x_start: float        # mm — canto superior esquerdo, coordenada X
    y_start: float        # mm — canto superior esquerdo, coordenada Y (decresce por linha)
    nx: int
    ny: int
    res_mm: float
    dose_matrix: np.ndarray   # shape (ny, nx), cGy

    # --- coordenadas ---------------------------------------------------------

    @property
    def x_coords(self) -> np.ndarray:
        """Posicoes X (mm) correspondentes às colunas da grade."""
        return self.x_start + np.arange(self.nx) * self.res_mm

    @property
    def y_coords(self) -> np.ndarray:
        """Posicoes Y (mm) correspondentes às linhas da grade (decresce)."""
        return self.y_start - np.arange(self.ny) * self.res_mm

    # --- extracao de perfis --------------------------------------------------

    def profile_at_y(self, y_mm: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
        """Perfil ao longo de X na linha mais proxima de y_mm.
        Retorna (x_arr_mm, dose_arr_cGy)."""
        row = int(round((self.y_start - y_mm) / self.res_mm))
        row = max(0, min(self.ny - 1, row))
        return self.x_coords.copy(), self.dose_matrix[row, :].copy()

    def profile_at_x(self, x_mm: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
        """Perfil ao longo de Y na coluna mais proxima de x_mm.
        Retorna (y_arr_mm, dose_arr_cGy)."""
        col = int(round((x_mm - self.x_start) / self.res_mm))
        col = max(0, min(self.nx - 1, col))
        return self.y_coords.copy(), self.dose_matrix[:, col].copy()


# ---------------------------------------------------------------------------
# Parser do arquivo .ALL
# ---------------------------------------------------------------------------

def _parse_plan_name(plan_name: str) -> Tuple[str, str]:
    """Interpreta o nome do plano e retorna (algoritmo, label_grid).

    Exemplos:
        'CCC1'     -> ('CCC', '1mm')
        'CCC05'    -> ('CCC', '0.5mm')
        'CCC1G02'  -> ('CCC', '1mm')
        'MontC1'   -> ('MonteCarlo', '1mm')
        'MontC05'  -> ('MonteCarlo', '0.5mm')
        'MontC1G01'-> ('MonteCarlo', '1mm')
    """
    if plan_name.upper().startswith('MONTC'):
        alg = 'MonteCarlo'
        rest = plan_name[5:]
    elif plan_name.upper().startswith('CCC'):
        alg = 'CCC'
        rest = plan_name[3:]
    else:
        return plan_name, '?'

    # remove sufixo tipo G01, G02
    rest = re.sub(r'[Gg]\d+$', '', rest)

    try:
        if rest.startswith('0') and len(rest) > 1:
            # '05' -> 0.5 mm
            grid_mm = float('0.' + rest[1:])
        else:
            grid_mm = float(rest)
        label = f"{grid_mm:g}mm"
    except ValueError:
        label = rest or '?'

    return alg, label


def parse_monaco_all(filepath: str) -> MonacoFile:
    """Le um arquivo Monaco .ALL e retorna um MonacoFile.

    Formato esperado no cabecalho (linhas separadas por linha em branco):
        PatientID,<id>
        PlaneDesc,C:  5.00 cm      (C=Coronal, T=Transversal)
        Upperleft,<x0>,<y0>        (mm)
        DosePtsxy,<nx>,<ny>
        DoseResmm,<res>            (mm/pixel)

    Apos o cabecalho: ny linhas de dados com nx valores separados por virgula (cGy).
    Linhas em branco sao ignoradas durante a leitura dos dados.
    """
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        raw_lines = [ln.rstrip() for ln in f]

    # --- cabecalho -----------------------------------------------------------
    meta: dict[str, str] = {}
    data_start_line = None

    for i, line in enumerate(raw_lines):
        if not line:
            continue
        # linha de dado: comeca com digito ou '-' e contem virgulas com numeros
        if (line[0].isdigit() or line[0] == '-') and ',' in line:
            parts = line.split(',')
            try:
                float(parts[0])
                # verificacao extra: todos os itens sao numeros
                all_nums = all(
                    p.strip().lstrip('-').replace('.', '', 1).isdigit()
                    for p in parts if p.strip()
                )
                if all_nums and len(parts) > 10:
                    data_start_line = i
                    break
            except ValueError:
                pass
        # linha de cabecalho chave,valor
        if ',' in line:
            k, _, v = line.partition(',')
            meta[k.strip()] = v.strip()

    if data_start_line is None:
        raise ValueError(f"Nenhuma linha de dados encontrada em: {filepath}")

    # --- metadados -----------------------------------------------------------
    plane_desc = meta.get('PlaneDesc', 'C: 0.00 cm')
    plane_type = 'Coronal' if plane_desc.lstrip().startswith('C') else 'Transversal'

    ul = meta.get('Upperleft', '-155.0,150.0').split(',')
    x_start = float(ul[0].strip())
    y_start = float(ul[1].strip())

    dxy = meta.get('DosePtsxy', '311,301').split(',')
    nx = int(dxy[0].strip())
    ny = int(dxy[1].strip())

    res_mm = float(meta.get('DoseResmm', '1.0').strip())
    patient_id = meta.get('PatientID', 'Unknown').strip()

    # nome do plano e posicao extraidos do nome do arquivo
    basename = os.path.splitext(os.path.basename(filepath))[0]
    parts = basename.split('.')
    plan_name = parts[1] if len(parts) > 1 else 'Unknown'
    pos_str = '.'.join(parts[3:5]) if len(parts) >= 5 else '0.00'
    try:
        position_mm = float(pos_str)
    except ValueError:
        position_mm = 0.0

    algorithm, grid_label = _parse_plan_name(plan_name)

    # --- grade de dose -------------------------------------------------------
    data_rows: list = []
    i = data_start_line
    while i < len(raw_lines) and len(data_rows) < ny:
        line = raw_lines[i]
        if line:
            try:
                vals = [float(v) for v in line.split(',') if v.strip()]
                if len(vals) == nx:
                    data_rows.append(vals)
                elif len(vals) > nx:
                    data_rows.append(vals[:nx])
            except ValueError:
                pass
        i += 1

    if len(data_rows) == 0:
        raise ValueError(f"Grade de dose vazia em: {filepath}")

    dose_matrix = np.array(data_rows, dtype=float)  # (ny, nx)

    return MonacoFile(
        filepath=filepath,
        patient_id=patient_id,
        plan_name=plan_name,
        algorithm=algorithm,
        grid_label=grid_label,
        plane_type=plane_type,
        position_mm=position_mm,
        x_start=x_start,
        y_start=y_start,
        nx=nx,
        ny=ny,
        res_mm=res_mm,
        dose_matrix=dose_matrix,
    )


# ---------------------------------------------------------------------------
# Conversao para objetos Scan
# ---------------------------------------------------------------------------

def monaco_to_scans(mf: MonacoFile, index_start: int = 0,
                    pdp_surface_offset_cm: float = 10.0) -> list:
    """Converte um MonacoFile em objetos Scan compatíveis com o parse_snctxt.

    Coronal  -> Crossline (Y=0) + Inline (X=0)
    Transversal -> Depth Scan (X=0), com z = -Y_monaco (profundidade cresce)
    """
    from parse_snctxt import Scan

    scans = []
    tag = f"[TPS] {mf.algorithm} {mf.grid_label}"
    depth_cm = mf.position_mm / 10.0

    if mf.plane_type == 'Coronal':
        # ---- Crossline: perfil ao longo de X em Y=0 ------------------------
        x_arr, dose_cl = mf.profile_at_y(y_mm=0.0)
        # converte mm -> cm para posicao
        x_cm = x_arr / 10.0
        s_cl = Scan(
            index=index_start,
            scan_type='Crossline',
            beam_type='Photon',
            detector=tag,
            energy='',
            field_x='',
            field_y='',
            ssd='',
            medium='',
        )
        s_cl.x = x_cm
        s_cl.y = np.zeros_like(x_cm)
        s_cl.z = np.full_like(x_cm, depth_cm)
        s_cl.dose = dose_cl
        s_cl.name_override = f"{mf.algorithm} {mf.grid_label} CL d={depth_cm:.1f}cm"
        scans.append(s_cl)

        # ---- Inline: perfil ao longo de Y em X=0 ---------------------------
        y_arr, dose_il = mf.profile_at_x(x_mm=0.0)
        y_cm = y_arr / 10.0
        s_il = Scan(
            index=index_start + 1,
            scan_type='Inline',
            beam_type='Photon',
            detector=tag,
            energy='',
            field_x='',
            field_y='',
            ssd='',
            medium='',
        )
        s_il.x = np.zeros_like(y_cm)
        s_il.y = y_cm
        s_il.z = np.full_like(y_cm, depth_cm)
        s_il.dose = dose_il
        s_il.name_override = f"{mf.algorithm} {mf.grid_label} IL d={depth_cm:.1f}cm"
        scans.append(s_il)

    elif mf.plane_type == 'Transversal':
        # ---- Depth Scan: perfil ao longo de Y_monaco em X=0 ----------------
        # Y_monaco decresce da superficie para o fundo; z = -Y_monaco para que
        # profundidade fisica aumente com z.
        y_arr, dose_pdp = mf.profile_at_x(x_mm=0.0)
        z_cm = -y_arr / 10.0   # z aumenta com a profundidade

        s_pdp = Scan(
            index=index_start,
            scan_type='Depth Scan',
            beam_type='Photon',
            detector=tag,
            energy='',
            field_x='',
            field_y='',
            ssd='',
            medium='',
        )
        s_pdp.x = np.zeros_like(z_cm)
        s_pdp.y = np.zeros_like(z_cm)
        s_pdp.z = z_cm
        s_pdp.dose = dose_pdp
        lat_mm = mf.position_mm
        s_pdp.name_override = (
            f"{mf.algorithm} {mf.grid_label} PDP"
            + (f" x={lat_mm:.0f}mm" if lat_mm != 0.0 else "")
        )
        # Desloca origem para que z=0 coincida com a superfície da água.
        # O isocêntro Monaco está a pdp_surface_offset_cm da superfície (padrão 10cm).
        s_pdp._pos_offset = pdp_surface_offset_cm
        scans.append(s_pdp)

    return scans


# ---------------------------------------------------------------------------
# API principal
# ---------------------------------------------------------------------------

def load_monaco_files(filepaths, index_start: int = 0) -> list:
    """Carrega multiplos arquivos Monaco .ALL e retorna lista plana de Scan.

    Args:
        filepaths:    lista de caminhos de arquivo .ALL
        index_start:  indice inicial para os Scan gerados

    Returns:
        lista de Scan (crosslines, inlines, PDPs) prontos para o UI
    """
    all_scans = []
    idx = index_start
    errors = []
    for fp in filepaths:
        try:
            mf = parse_monaco_all(fp)
            new_scans = monaco_to_scans(mf, index_start=idx)
            all_scans.extend(new_scans)
            idx += len(new_scans)
        except Exception as e:
            errors.append(f"{os.path.basename(fp)}: {e}")
    if errors:
        import warnings
        warnings.warn("Erros ao carregar Monaco:\n" + "\n".join(errors))
    return all_scans


# ---------------------------------------------------------------------------
# Teste rapido de linha de comando
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import sys, glob

    patterns = sys.argv[1:] or [
        r'C:\RadCalc\Inst1\data\srs e2e\dados\TimeB\MonacoPhantom.CCC1.*.ALL'
    ]
    files = []
    for pat in patterns:
        files.extend(glob.glob(pat))

    if not files:
        print("Nenhum arquivo encontrado.")
        sys.exit(1)

    scans = load_monaco_files(sorted(files))
    print(f"\n{len(scans)} perfis extraidos:\n")
    for s in scans:
        pos = s.position
        det = s.detector_override or s.detector
        print(f"  [{s.index:>3}]  {s.scan_type:<12}  {s.name_override:<40}  "
              f"N={len(pos)}  pos=[{pos.min():.1f}, {pos.max():.1f}]cm  "
              f"dose_max={s.dose.max():.2f} cGy")

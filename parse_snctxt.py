"""
Parser e plotter para arquivos .snctxt do Sun Nuclear SNC
Adaptado de: https://github.com/mwgeurts/snc_extract (MATLAB)
Inclui cálculo de índice gama 1D

Uso:
    python parse_snctxt.py
    python parse_snctxt.py --file "Curso SRS.snctxt" --plot
    python parse_snctxt.py --gamma --ref 0 --eval 1  (scan 0 vs scan 1)
"""

import re
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Estrutura de dados
# ---------------------------------------------------------------------------

@dataclass
class Scan:
    index: int
    beam_type: str = ""
    energy: str = ""
    field_x: str = ""
    field_y: str = ""
    scan_type: str = ""
    ssd: str = ""
    medium: str = ""
    date: str = ""
    action: str = ""
    detector: str = ""
    # arrays numpy (dados originais, nunca alterados)
    x: np.ndarray = field(default_factory=lambda: np.array([]))
    y: np.ndarray = field(default_factory=lambda: np.array([]))
    z: np.ndarray = field(default_factory=lambda: np.array([]))
    dose: np.ndarray = field(default_factory=lambda: np.array([]))
    # campos editáveis
    name_override: str = ""
    energy_override: str = ""
    detector_override: str = ""
    field_x_override: str = ""
    field_y_override: str = ""
    # edições aplicadas à dose (None = usar original)
    _dose_edit: Optional[np.ndarray] = field(default=None, repr=False)
    # deslocamento de posição para centralização (cm)
    _pos_offset: float = 0.0

    # ------------------------------------------------------------------
    # Labels e display
    # ------------------------------------------------------------------

    @property
    def display_name(self) -> str:
        if self.name_override:
            return self.name_override
        fs = f"{self.display_field_x}x{self.display_field_y}"
        e  = self.energy_override or self.energy
        b  = self.beam_type
        return f"[{self.index}] {b} {e}MeV/MV  {fs}cm  {self.scan_type}"

    @property
    def label(self) -> str:
        return self.display_name

    @property
    def display_energy(self) -> str:
        return self.energy_override or self.energy

    @property
    def display_detector(self) -> str:
        return self.detector_override or self.detector

    @property
    def display_field_x(self) -> str:
        return self.field_x_override or self.field_x

    @property
    def display_field_y(self) -> str:
        return self.field_y_override or self.field_y

    # ------------------------------------------------------------------
    # Posição e dose (com edições aplicadas)
    # ------------------------------------------------------------------

    @property
    def _sort_order(self) -> np.ndarray:
        if self.scan_type == "Depth Scan":
            return np.argsort(self.z)
        elif self.scan_type == "Crossline":
            return np.argsort(self.x)
        elif self.scan_type == "Inline":
            return np.argsort(self.y)
        return np.arange(len(self.z))

    @property
    def position(self) -> np.ndarray:
        """Posição ordenada + offset de centralização."""
        if self.scan_type == "Depth Scan":
            return self.z[self._sort_order] + self._pos_offset
        elif self.scan_type == "Crossline":
            return self.x[self._sort_order] + self._pos_offset
        elif self.scan_type == "Inline":
            return self.y[self._sort_order] + self._pos_offset
        return self.z + self._pos_offset

    @property
    def dose_sorted(self) -> np.ndarray:
        """Dose original ordenada (sem edições)."""
        return self.dose[self._sort_order]

    @property
    def dose_display(self) -> np.ndarray:
        """Dose com edições (suavização/renormalização) aplicadas."""
        if self._dose_edit is not None:
            return self._dose_edit
        return self.dose_sorted

    # ------------------------------------------------------------------
    # Edição de dose
    # ------------------------------------------------------------------

    def reset_edits(self):
        """Remove todas as edições; restaura dose original e offset de posição."""
        self._dose_edit = None
        self._pos_offset = 0.0
        self.field_x_override = ""
        self.field_y_override = ""


# ---------------------------------------------------------------------------
# Métricas dosimétricase
# ---------------------------------------------------------------------------

def _interp_crossing(pos: np.ndarray, dose: np.ndarray,
                     level: float, side: str = "left") -> Optional[float]:
    """Encontra a posição onde dose == level por interpolação linear.
    side='left'  : primeira travessia da esquerda (subida)
    side='right' : última travessia da direita (descida)
    side='falling': primeira travessia na parte descendente (após pico)
    """
    crossings = []
    for i in range(len(dose) - 1):
        d0, d1 = dose[i], dose[i + 1]
        if (d0 - level) * (d1 - level) < 0:
            t = (level - d0) / (d1 - d0)
            crossings.append(pos[i] + t * (pos[i + 1] - pos[i]))
    if not crossings:
        return None
    if side == "left":
        return crossings[0]
    if side == "right":
        return crossings[-1]
    if side == "falling":
        # após o pico
        idx_peak = np.argmax(dose)
        pos_peak = pos[idx_peak]
        after = [c for c in crossings if c > pos_peak]
        return after[0] if after else None
    return crossings[0]


def compute_profile_metrics(pos: np.ndarray, dose: np.ndarray) -> dict:
    """Métricas de perfil transversal (Crossline / Inline).

    Retorna dicionário com:
      fwhm        : largura a 50% (cm)
      center      : centro geométrico do campo (média dos lados 50%) (cm)
      l50, r50    : posições esquerda/direita a 50%
      l80, r80    : posições esquerda/direita a 80%
      l20, r20    : posições esquerda/direita a 20%
      penumbra_l  : penumbra esquerda 80%→20% (cm)
      penumbra_r  : penumbra direita  20%→80% (cm)
      flatness    : (Dmax-Dmin)/Dmax × 100 na região central 80% do campo
      symmetry    : (Dl - Dr) / Dcenter × 100 (simetria)
    """
    d_norm = dose / np.max(dose) * 100.0

    l50 = _interp_crossing(pos, d_norm, 50.0, "left")
    r50 = _interp_crossing(pos, d_norm, 50.0, "right")
    l80 = _interp_crossing(pos, d_norm, 80.0, "left")
    r80 = _interp_crossing(pos, d_norm, 80.0, "right")
    l20 = _interp_crossing(pos, d_norm, 20.0, "left")
    r20 = _interp_crossing(pos, d_norm, 20.0, "right")

    fwhm   = (r50 - l50) if (l50 is not None and r50 is not None) else None
    center = ((l50 + r50) / 2) if (l50 is not None and r50 is not None) else None

    pen_l = (l80 - l20) if (l80 is not None and l20 is not None) else None
    pen_r = (r20 - r80) if (r20 is not None and r80 is not None) else None

    # Planura: região central de 80% do FWHM
    flatness = symmetry = None
    if fwhm is not None and center is not None:
        half = fwhm * 0.4   # ±40% do centro = 80% do campo
        mask = (pos >= center - half) & (pos <= center + half)
        if mask.sum() >= 2:
            d_inner = d_norm[mask]
            flatness = (d_inner.max() - d_inner.min()) / d_inner.max() * 100

        # simetria: espelha em torno do centro
        sym_pts = []
        for xi, di in zip(pos, d_norm):
            mirror = center - (xi - center)
            di_mirror = np.interp(mirror, pos, d_norm)
            if not np.isnan(di_mirror):
                sym_pts.append(abs(di - di_mirror))
        if sym_pts:
            d_cax = np.interp(center, pos, d_norm) if center is not None else 1
            symmetry = max(sym_pts) / d_cax * 100 if d_cax else None

    return dict(fwhm=fwhm, center=center,
                l50=l50, r50=r50,
                l80=l80, r80=r80,
                l20=l20, r20=r20,
                penumbra_l=pen_l, penumbra_r=pen_r,
                flatness=flatness, symmetry=symmetry)


def compute_pdp_metrics(pos: np.ndarray, dose: np.ndarray) -> dict:
    """Métricas de perfil de profundidade (Depth Scan).

    pos deve estar ordenado ascendente (superfície → fundo).

    Retorna dicionário com:
      zmax   : profundidade do máximo (cm)
      dmax   : valor máximo de dose (% relativa, antes de normalizar)
      r90    : profundidade onde dose cai a 90% do dmax (borda distal)
      r80    : profundidade onde dose cai a 80% do dmax
      r50    : profundidade onde dose cai a 50% do dmax
      r10    : profundidade onde dose cai a 10% do dmax
      r5     : profundidade onde dose cai a 5% do dmax
      rp     : alcance prático (intersecção da tangente máxima com a linha de base)
    """
    d_norm = dose / np.max(dose) * 100.0
    idx_peak = int(np.argmax(d_norm))
    zmax  = float(pos[idx_peak])
    dmax  = float(np.max(dose))

    def _after_peak(level):
        """Crossings somente na parte descendente (pos > zmax)."""
        return _interp_crossing(pos, d_norm, level, "falling")

    r90 = _after_peak(90.0)
    r80 = _after_peak(80.0)
    r50 = _after_peak(50.0)
    r10 = _after_peak(10.0)
    r5  = _after_peak(5.0)

    # Alcance prático Rp: tangente no ponto de inflexão da queda
    rp = None
    try:
        # região descendente após o pico
        mask_fall = pos > zmax
        if mask_fall.sum() > 5:
            p_fall = pos[mask_fall]
            d_fall = d_norm[mask_fall]
            # gradiente
            grad = np.gradient(d_fall, p_fall)
            idx_infl = int(np.argmin(grad))   # gradiente mais negativo = ponto de inflexão
            slope = grad[idx_infl]
            if slope != 0:
                intercept = d_fall[idx_infl] - slope * p_fall[idx_infl]
                rp = -intercept / slope   # cruzamento com dose=0
    except Exception:
        pass

    return dict(zmax=zmax, dmax=dmax,
                r90=r90, r80=r80, r50=r50,
                r10=r10, r5=r5, rp=rp)


def format_metrics(scan: Scan) -> str:
    """Retorna string formatada com as métricas do scan (usa dose_display)."""
    pos  = scan.position
    dose = scan.dose_display
    lines = []

    if scan.scan_type == "Depth Scan":
        m = compute_pdp_metrics(pos, dose)
        lines.append(f"  Zmax        = {m['zmax']:.3f} cm")
        lines.append(f"  Dmax (rel.) = {m['dmax']:.2f} %")
        for key, label in [("r90","R90"), ("r80","R80"), ("r50","R50"),
                            ("r10","R10"), ("r5","R5"),  ("rp","Rp")]:
            val = m[key]
            lines.append(f"  {label:<12}= {val:.3f} cm" if val is not None
                         else f"  {label:<12}= n.d.")
    else:
        m = compute_profile_metrics(pos, dose)
        def _fmt(v): return f"{v:.3f} cm" if v is not None else "n.d."
        lines.append(f"  FWHM (50%)  = {_fmt(m['fwhm'])}")
        lines.append(f"  Centro      = {_fmt(m['center'])}")
        lines.append(f"  L50 / R50   = {_fmt(m['l50'])} / {_fmt(m['r50'])}")
        lines.append(f"  L80 / R80   = {_fmt(m['l80'])} / {_fmt(m['r80'])}")
        lines.append(f"  L20 / R20   = {_fmt(m['l20'])} / {_fmt(m['r20'])}")
        lines.append(f"  Penumbra E  = {_fmt(m['penumbra_l'])}")
        lines.append(f"  Penumbra D  = {_fmt(m['penumbra_r'])}")
        if m['flatness'] is not None:
            lines.append(f"  Planura     = {m['flatness']:.2f} %")
        if m['symmetry'] is not None:
            lines.append(f"  Simetria    = {m['symmetry']:.2f} %")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Funções de edição de curva
# ---------------------------------------------------------------------------

def apply_smooth(scan: Scan,
                 method: str = "savgol",
                 window: int = 11,
                 polyorder: int = 3,
                 sigma: float = 1.0,
                 region: Optional[tuple[float, float]] = None) -> None:
    """Suaviza a dose do scan in-place.

    Parâmetros
    ----------
    method  : "savgol" | "moving_avg" | "gaussian"
    window  : tamanho da janela (savgol / moving_avg); deve ser ímpar
    polyorder: grau do polinômio (savgol)
    sigma   : desvio-padrão em pontos (gaussian)
    region  : (pos_start, pos_end) em cm; None = curva toda
    """
    from scipy.signal import savgol_filter
    from scipy.ndimage import gaussian_filter1d

    pos  = scan.position
    base = scan._dose_edit if scan._dose_edit is not None else scan.dose_sorted
    dose = base.copy()

    # índices da região
    if region is not None:
        mask = (pos >= region[0]) & (pos <= region[1])
    else:
        mask = np.ones(len(pos), dtype=bool)

    idx = np.where(mask)[0]
    if len(idx) < max(window, 5):
        raise ValueError(f"Região muito pequena ({len(idx)} pts) para suavizar.")

    segment = dose[idx]

    if method == "savgol":
        w = window if window % 2 == 1 else window + 1
        w = max(w, polyorder + 1)
        segment = savgol_filter(segment, w, polyorder)
    elif method == "moving_avg":
        kernel = np.ones(window) / window
        pad = window // 2
        padded = np.pad(segment, pad, mode="edge")
        segment = np.convolve(padded, kernel, mode="valid")[:len(idx)]
    elif method == "gaussian":
        segment = gaussian_filter1d(segment.astype(float), sigma)
    else:
        raise ValueError(f"Método desconhecido: {method}")

    dose[idx] = segment
    scan._dose_edit = dose


def apply_centering(scan: Scan, method: str = "fwhm") -> float:
    """Centraliza o eixo de posição in-place.

    Para perfis (Cross/Inline): calcula o centro geométrico do campo e
    aplica um offset para que o centro fique em posição=0.

    Para PDDs: opcionalmente desloca Z para que dmax fique em Z=0
    (desativado por padrão; use method="dmax" para ativar).

    Retorna o offset aplicado (cm).
    """
    pos  = scan.position - scan._pos_offset   # posição sem offset atual
    dose = scan.dose_display

    if scan.scan_type == "Depth Scan":
        if method == "dmax":
            offset = -pos[np.argmax(dose)]
        else:
            offset = 0.0  # sem sentido centralizar PDP em 0

    else:  # Cross ou Inline
        dose_norm = dose / np.max(dose) * 100.0

        if method == "fwhm":
            # média das duas posições onde dose = 50%
            half = 50.0
            crossings = []
            for i in range(len(dose_norm) - 1):
                if (dose_norm[i] - half) * (dose_norm[i+1] - half) < 0:
                    t = (half - dose_norm[i]) / (dose_norm[i+1] - dose_norm[i])
                    crossings.append(pos[i] + t * (pos[i+1] - pos[i]))
            if len(crossings) >= 2:
                center = (crossings[0] + crossings[-1]) / 2
            else:
                center = pos[np.argmax(dose_norm)]
        elif method == "cax":
            center = pos[np.argmax(dose_norm)]
        else:
            center = 0.0

        offset = -center

    scan._pos_offset = offset
    return offset


def apply_renormalize(scan: Scan,
                      method: str = "dmax",
                      value: float = 100.0,
                      ref_pos: Optional[float] = None,
                      region: Optional[tuple[float, float]] = None) -> None:
    """Renormaliza a dose in-place.

    method:
      "dmax"   – normaliza pelo pico (dose_max → value)
      "cax"    – normaliza pelo ponto mais próximo de pos=0
      "point"  – normaliza pelo ponto mais próximo de ref_pos
      "region" – normaliza pela média da região [region[0], region[1]]
    value : valor percentual a atribuir ao ponto de referência (padrão 100)
    """
    pos  = scan.position
    base = scan._dose_edit if scan._dose_edit is not None else scan.dose_sorted
    dose = base.copy()

    if method == "dmax":
        ref_val = np.max(dose)
    elif method == "cax":
        idx = np.argmin(np.abs(pos))
        ref_val = dose[idx]
    elif method == "point":
        if ref_pos is None:
            raise ValueError("ref_pos obrigatório para method='point'")
        idx = np.argmin(np.abs(pos - ref_pos))
        ref_val = dose[idx]
    elif method == "region":
        if region is None:
            raise ValueError("region obrigatório para method='region'")
        mask = (pos >= region[0]) & (pos <= region[1])
        ref_val = np.mean(dose[mask]) if mask.any() else np.max(dose)
    else:
        raise ValueError(f"Método desconhecido: {method}")

    if ref_val == 0:
        raise ValueError("Valor de referência = 0; não é possível renormalizar.")

    scan._dose_edit = dose / ref_val * value


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_snctxt(filepath: str) -> list[Scan]:
    """Lê um arquivo .snctxt e retorna lista de Scan."""
    path = Path(filepath)
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        lines = f.readlines()

    scans = []
    scan: Optional[Scan] = None
    in_dose_table = False
    dose_rows = []
    scan_idx = 0

    # mapeamento de campos do header para atributos do Scan
    header_map = {
        "Summary Beam Type":           "beam_type",
        "Energy (MV / MeV)":           "energy",
        "Summary FieldSize X (cm)":    "field_x",
        "Summary FieldSize Y (cm)":    "field_y",
        "Summary Scan Type":           "scan_type",
        "Source to Surface Distance (cm)": "ssd",
        "Scan Medium":                 "medium",
        "Scan Date":                   "date",
        "Action":                      "action",
        "Field Detector Model #":      "detector",
    }

    for line in lines:
        line = line.rstrip("\n").rstrip("\r")

        if line.strip() == "BEGIN SCAN":
            if scan is not None:
                _finish_scan(scan, dose_rows, scans)
            scan = Scan(index=scan_idx)
            scan_idx += 1
            dose_rows = []
            in_dose_table = False
            continue

        if line.strip() == "BEGIN DOSE TABLE":
            in_dose_table = True
            continue

        if line.strip() == "END DOSE TABLE":
            in_dose_table = False
            continue

        if scan is None:
            continue

        # dentro da tabela de dose
        if in_dose_table:
            parts = line.split("\t")
            # linha de cabeçalho das colunas
            if "X (cm)" in line:
                continue
            # linha de dados: começa com tab vazio + 4 valores
            if len(parts) >= 5:
                try:
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    d = float(parts[4])
                    dose_rows.append((x, y, z, d))
                except ValueError:
                    pass
            continue

        # fora da tabela: parse de campos de metadados
        parts = line.split("\t")
        if len(parts) >= 2:
            key = parts[0].strip()
            val = parts[1].strip() if len(parts) > 1 else ""
            attr = header_map.get(key)
            if attr:
                setattr(scan, attr, val)

    # último scan
    if scan is not None:
        _finish_scan(scan, dose_rows, scans)

    return scans


def _finish_scan(scan: Scan, dose_rows: list, scans: list):
    if dose_rows:
        arr = np.array(dose_rows)
        scan.x = arr[:, 0]
        scan.y = arr[:, 1]
        scan.z = arr[:, 2]
        scan.dose = arr[:, 3]
    scans.append(scan)


# ---------------------------------------------------------------------------
# Índice Gama 1D
# ---------------------------------------------------------------------------

def gamma_1d(
    pos_ref: np.ndarray,
    dose_ref: np.ndarray,
    pos_eval: np.ndarray,
    dose_eval: np.ndarray,
    dd: float = 3.0,     # critério dose-diferença em %
    dta: float = 0.3,    # critério DTA em cm (3 mm)
    norm: str = "max",   # "max" ou "local"
    threshold: float = 10.0,  # % do máximo abaixo do qual não calcula
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calcula o índice gama 1D entre curva de referência e avaliação.

    Algoritmo:
        Para cada ponto de avaliação (interpolado na grade de referência),
        γ = min sqrt( (ΔD/DD)² + (Δd/DTA)² )
        onde Δd percorre todos os pontos de referência.

    Retorna:
        pos_out : posições onde γ foi calculado (acima do threshold)
        gamma   : valores de γ
    """
    # Normalização
    if norm == "max":
        norm_ref  = np.max(dose_ref)
        norm_eval = np.max(dose_eval)
    else:
        norm_ref = norm_eval = 1.0  # local: usa dose local

    dose_ref_n  = dose_ref  / norm_ref  * 100
    dose_eval_n = dose_eval / norm_eval * 100

    # Threshold baseado no máximo da referência
    thresh_val = threshold  # em % normalizado

    # Interpolação da referência numa grade fina
    interp_pos = np.linspace(pos_ref.min(), pos_ref.max(), len(pos_ref) * 10)
    interp_dose = np.interp(interp_pos, pos_ref, dose_ref_n)

    # Calcula gama para cada ponto de avaliação dentro do range de referência
    pos_out = []
    gamma_out = []

    for i, p in enumerate(pos_eval):
        d_eval = dose_eval_n[i]
        if d_eval < thresh_val:
            continue
        if p < pos_ref.min() or p > pos_ref.max():
            continue

        if norm == "local":
            d_ref_local = np.interp(p, interp_pos, interp_dose)
            dd_val = dd * d_ref_local / 100
        else:
            dd_val = dd

        # Vetor de γ² para todos os pontos de referência
        delta_d = (interp_dose - d_eval) / dd_val
        delta_pos = (interp_pos - p) / dta
        gamma_sq = delta_d**2 + delta_pos**2
        gamma_val = np.sqrt(np.min(gamma_sq))

        pos_out.append(p)
        gamma_out.append(gamma_val)

    return np.array(pos_out), np.array(gamma_out)


def gamma_pass_rate(gamma: np.ndarray) -> float:
    """Retorna a taxa de aprovação (γ ≤ 1) em %."""
    if len(gamma) == 0:
        return float("nan")
    return 100.0 * np.sum(gamma <= 1.0) / len(gamma)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def _normalize(dose: np.ndarray, pos: np.ndarray, scan_type: str) -> np.ndarray:
    """
    Normaliza a curva:
      - Depth Scan  -> 100% no dmax (pico)
      - Cross/Inline -> 100% no centro (posição mais próxima de 0)
    dose e pos já devem estar ordenados (use scan.dose_sorted / scan.position).
    """
    if scan_type == "Depth Scan":
        peak = np.max(dose)
        return dose / peak * 100.0 if peak != 0 else dose
    else:
        idx_center = np.argmin(np.abs(pos))
        cax = dose[idx_center]
        return dose / cax * 100.0 if cax != 0 else dose


def plot_scans(scans: list[Scan], title: str = "", normalize: bool = True):
    """Plota todos os scans agrupados por beam_type e scan_type.

    Melhorias:
      - PDD: eixo X começa na superfície (z crescente -> profundidade)
      - Depth Scan normalizado ao dmax (100%)
      - Cross/Inline normalizados ao CAX (100% no centro)
      - Cores consistentes por tamanho de campo entre subplots
    """
    from collections import defaultdict

    groups = defaultdict(list)
    for s in scans:
        key = (s.beam_type, s.scan_type)
        groups[key].append(s)

    n_groups = len(groups)
    if n_groups == 0:
        print("Nenhum scan para plotar.")
        return

    ncols = min(3, n_groups)
    nrows = (n_groups + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows))
    axes = np.array(axes).flatten()

    color_map = plt.colormaps["tab10"]

    # Mapeamento global: tamanho de campo -> cor
    all_fs = sorted({f"{s.field_x}x{s.field_y}" for s in scans})
    fs_color = {fs: color_map(i % 10) for i, fs in enumerate(all_fs)}

    for ax_idx, ((beam, scan_t), scan_list) in enumerate(sorted(groups.items())):
        ax = axes[ax_idx]

        for s in scan_list:
            fs = f"{s.field_x}x{s.field_y}"
            color = fs_color[fs]
            pos  = s.position          # já ordenado + offset de centralização
            dose_s = s.dose_display    # com edições aplicadas
            dose = _normalize(dose_s, pos, scan_t) if normalize else dose_s

            if scan_t == "Depth Scan":
                # Z já ordenado: Z_min (superfície) na esquerda, Z_max (fundo) na direita
                ax.plot(pos, dose, label=f"{fs} cm  [{s.index}]",
                        color=color, linewidth=1.5)
                xlabel = "Profundidade Z (cm)"
            else:
                ax.plot(pos, dose, label=f"{fs} cm  [{s.index}]",
                        color=color, linewidth=1.5)
                ax.axvline(0, color="gray", linewidth=0.8, linestyle=":")
                xlabel = "Posição X (cm)" if scan_t == "Crossline" else "Posição Y (cm)"

        ylabel = "Dose Relativa (%)" + (" – norm. dmax/CAX" if normalize else "")
        ax.set_title(f"{beam}  –  {scan_t}", fontsize=10, fontweight="bold")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=7, loc="best")
        ax.grid(True, alpha=0.3)
        if normalize:
            ax.axhline(100, color="gray", linewidth=0.7, linestyle="--")
            ax.axhline(50,  color="gray", linewidth=0.5, linestyle=":")

    for ax in axes[n_groups:]:
        ax.set_visible(False)

    fig.suptitle(title or "Scans SNC .snctxt", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.show()


def plot_gamma(
    pos_ref, dose_ref,
    pos_eval, dose_eval,
    pos_g, gamma,
    label_ref="Referencia",
    label_eval="Avaliacao",
    dd=3.0, dta=0.3,
    scan_type="",
):
    """Plota curvas normalizadas + diferenca de dose + indice gama 1D."""
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True,
                              gridspec_kw={"height_ratios": [2.5, 1, 1.5]})
    ax1, ax2, ax3 = axes

    # --- Painel 1: curvas normalizadas ---
    dose_ref_n  = dose_ref  / np.max(dose_ref)  * 100
    dose_eval_n = dose_eval / np.max(dose_eval) * 100

    # pos já vem ordenado de scan.position (Z_min na esquerda para Depth)
    x_ref  = pos_ref
    x_eval = pos_eval
    xlabel = "Profundidade Z (cm)" if "Depth" in scan_type else "Posicao (cm)"

    ax1.plot(x_ref,  dose_ref_n,  "b-",  linewidth=2,   label=label_ref)
    ax1.plot(x_eval, dose_eval_n, "r--", linewidth=1.5, label=label_eval)
    ax1.axhline(100, color="gray", linewidth=0.7, linestyle="--")
    ax1.axhline(50,  color="gray", linewidth=0.5, linestyle=":")
    ax1.set_ylabel("Dose Relativa (%)")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_title(f"Indice Gama 1D  -  criterio {dd}% / {dta*10:.0f}mm", fontsize=11)

    # --- Painel 2: diferença de dose (%) ---
    dose_eval_interp = np.interp(pos_ref, pos_eval, dose_eval_n)
    diff = dose_eval_interp - dose_ref_n
    ax2.bar(x_ref, diff, width=np.diff(x_ref).mean() * 0.9,
            color=np.where(np.abs(diff) <= dd, "steelblue", "tomato"),
            alpha=0.7)
    ax2.axhline(0,   color="black", linewidth=0.8)
    ax2.axhline( dd, color="orange", linewidth=1, linestyle="--", label=f"+{dd}%")
    ax2.axhline(-dd, color="orange", linewidth=1, linestyle="--", label=f"-{dd}%")
    ax2.set_ylabel("Diff dose (%)")
    ax2.legend(fontsize=7, loc="upper right")
    ax2.grid(True, alpha=0.3)

    # --- Painel 3: gama ---
    x_g = pos_g

    colors_g = ["limegreen" if v <= 1.0 else "crimson" for v in gamma]
    ax3.scatter(x_g, gamma, c=colors_g, s=10, zorder=3)
    ax3.axhline(1.0, color="black", linewidth=1.2, linestyle="--")
    ax3.set_ylabel("gamma")
    ax3.set_xlabel(xlabel)
    ax3.set_ylim(0, max(2.0, np.nanmax(gamma) * 1.15) if len(gamma) > 0 else 2)
    ax3.grid(True, alpha=0.3)

    pr = gamma_pass_rate(gamma)
    n_pass = int(np.sum(gamma <= 1.0))
    ax3.set_title(
        f"Pass rate (gamma<=1): {pr:.1f}%   "
        f"({n_pass}/{len(gamma)} pts)   "
        f"gamma_max={np.max(gamma):.3f}   gamma_mean={np.mean(gamma):.3f}",
        fontsize=9
    )

    plt.tight_layout()
    plt.show()
    return pr


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Parser e análise de arquivos .snctxt do Sun Nuclear SNC"
    )
    parser.add_argument(
        "--file", default="Curso SRS.snctxt",
        help="Caminho do arquivo .snctxt"
    )
    parser.add_argument("--list", action="store_true",
                        help="Lista todos os scans e sai")
    parser.add_argument("--plot", action="store_true",
                        help="Plota todos os scans")
    parser.add_argument("--no-normalize", action="store_true",
                        help="Desliga normalizacao dmax/CAX nos plots")
    parser.add_argument("--gamma", action="store_true",
                        help="Calcula índice gama entre dois scans")
    parser.add_argument("--ref",  type=int, default=0,
                        help="Índice do scan de referência (padrão: 0)")
    parser.add_argument("--eval", type=int, default=1,
                        help="Índice do scan de avaliação (padrão: 1)")
    parser.add_argument("--dd",  type=float, default=3.0,
                        help="Critério dose-diferença em %% (padrão: 3)")
    parser.add_argument("--dta", type=float, default=0.3,
                        help="Critério DTA em cm (padrão: 0.3 = 3mm)")
    parser.add_argument("--threshold", type=float, default=10.0,
                        help="Threshold %% do máximo (padrão: 10)")
    parser.add_argument("--norm", choices=["max", "local"], default="max",
                        help="Tipo de normalização (padrão: max)")
    args = parser.parse_args()

    # resolução do caminho
    fpath = Path(args.file)
    if not fpath.is_absolute():
        fpath = Path(__file__).parent / fpath

    print(f"Lendo: {fpath}")
    scans = parse_snctxt(str(fpath))
    print(f"  -> {len(scans)} scans encontrados\n")

    # lista de scans
    print(f"{'#':>3}  {'Feixe':<10} {'Energia':>8}  {'Campo':>8}  {'Tipo':<15}  {'Pts':>5}  {'SSD':>6}  Data")
    print("-" * 80)
    for s in scans:
        npts = len(s.dose)
        print(f"{s.index:>3}  {s.beam_type:<10} {s.energy:>8}  "
              f"{s.field_x+'x'+s.field_y:>8}  {s.scan_type:<15}  {npts:>5}  "
              f"{s.ssd:>6}  {s.date}")

    if args.list:
        return

    if args.plot:
        fname = fpath.stem
        plot_scans(scans, title=fname, normalize=not args.no_normalize)

    if args.gamma:
        if args.ref >= len(scans) or args.eval >= len(scans):
            print(f"Erro: índices inválidos (máx {len(scans)-1})")
            sys.exit(1)

        s_ref  = scans[args.ref]
        s_eval = scans[args.eval]

        print(f"\nGama 1D:")
        print(f"  Referencia : {s_ref.label}")
        print(f"  Avaliacao  : {s_eval.label}")
        print(f"  Criterio   : {args.dd}% / {args.dta*10:.0f}mm  threshold={args.threshold}%")

        pos_g, g = gamma_1d(
            s_ref.position, s_ref.dose_sorted,
            s_eval.position, s_eval.dose_sorted,
            dd=args.dd,
            dta=args.dta,
            norm=args.norm,
            threshold=args.threshold,
        )
        pr = gamma_pass_rate(g)
        print(f"  Pass rate (gamma<=1): {pr:.2f}%  ({int(np.sum(g<=1))}/{len(g)} pts)")
        print(f"  gamma_max = {np.max(g):.3f}  gamma_mean = {np.mean(g):.3f}")

        plot_gamma(
            s_ref.position, s_ref.dose_sorted,
            s_eval.position, s_eval.dose_sorted,
            pos_g, g,
            label_ref=f"Ref [{args.ref}] {s_ref.scan_type} {s_ref.field_x}x{s_ref.field_y}cm",
            label_eval=f"Eval [{args.eval}] {s_eval.scan_type} {s_eval.field_x}x{s_eval.field_y}cm",
            dd=args.dd, dta=args.dta,
            scan_type=s_ref.scan_type,
        )


if __name__ == "__main__":
    main()

"""
Interface gráfica para análise e edição de arquivos .snctxt do Sun Nuclear SNC.

Uso:
    py -3 srs_ui.py
    py -3 srs_ui.py --file "Curso SRS.snctxt"
"""

import argparse
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from parse_snctxt import (
    parse_snctxt, Scan, gamma_1d, gamma_pass_rate, _normalize,
    apply_smooth, apply_centering, apply_renormalize,
    format_metrics, compute_profile_metrics, compute_pdp_metrics,
)
from parse_monaco import load_monaco_files
from parse_mcc          import load_mcc_files
from parse_dta          import load_dta_files
from parse_multiformats import load_auto as load_multi
from parse_rfb          import load_rfb_files

# ---------------------------------------------------------------------------
# Paleta: cor por detector
# ---------------------------------------------------------------------------
_CMAP = plt.colormaps["tab10"]
_DET_COLORS: dict[str, tuple] = {}

def det_color(det: str) -> tuple:
    key = det or "?"
    if key not in _DET_COLORS:
        _DET_COLORS[key] = _CMAP(len(_DET_COLORS) % 10)
    return _DET_COLORS[key]


# ---------------------------------------------------------------------------
# Helpers de widget
# ---------------------------------------------------------------------------

def labeled_entry(parent, row, text, default="", width=14, col=0):
    ttk.Label(parent, text=text).grid(row=row, column=col,   sticky=tk.W, pady=1, padx=2)
    var = tk.StringVar(value=default)
    ttk.Entry(parent, textvariable=var, width=width).grid(
        row=row, column=col+1, sticky=tk.EW, pady=1, padx=2)
    return var

def labeled_combo(parent, row, text, values, default=None, width=14, col=0):
    ttk.Label(parent, text=text).grid(row=row, column=col,   sticky=tk.W, pady=1, padx=2)
    cb = ttk.Combobox(parent, values=values, state="readonly", width=width)
    cb.grid(row=row, column=col+1, sticky=tk.EW, pady=1, padx=2)
    if default and default in values:
        cb.set(default)
    elif values:
        cb.current(0)
    return cb

def scrolled_listbox(parent, height=8, font=("Courier", 9), mode=tk.EXTENDED):
    frame = ttk.Frame(parent)
    lb = tk.Listbox(frame, selectmode=mode, font=font,
                    activestyle="dotbox", height=height,
                    exportselection=False)
    sy = ttk.Scrollbar(frame, orient=tk.VERTICAL,   command=lb.yview)
    sx = ttk.Scrollbar(frame, orient=tk.HORIZONTAL, command=lb.xview)
    lb.config(yscrollcommand=sy.set, xscrollcommand=sx.set)
    lb.grid(row=0, column=0, sticky=tk.NSEW)
    sy.grid(row=0, column=1, sticky=tk.NS)
    sx.grid(row=1, column=0, sticky=tk.EW)
    frame.rowconfigure(0, weight=1)
    frame.columnconfigure(0, weight=1)
    return frame, lb


# ---------------------------------------------------------------------------
# Aplicação
# ---------------------------------------------------------------------------

class SRSApp(tk.Tk):
    def __init__(self, filepath=""):
        super().__init__()
        self.title("SNC .snctxt — Análise de Curvas")
        self.geometry("1420x860")
        self.resizable(True, True)

        self.scans: list[Scan] = []
        self.filepath = ""
        self._scans_filtered: list[Scan] = []  # lista atual filtrada (aba Scans)

        self._build_ui()
        if filepath:
            self._load_file(filepath)

    # ======================================================================
    # Layout geral
    # ======================================================================

    def _build_ui(self):
        # Barra de arquivo
        top = ttk.Frame(self, padding=(4, 3))
        top.pack(side=tk.TOP, fill=tk.X)
        ttk.Button(top, text="Abrir arquivo...", command=self._open_file   ).pack(side=tk.LEFT, padx=2)
        ttk.Button(top, text="Importar Monaco...", command=self._import_monaco).pack(side=tk.LEFT, padx=2)
        ttk.Button(top, text="Importar MCC...",   command=self._import_mcc   ).pack(side=tk.LEFT, padx=2)
        ttk.Button(top, text="Importar DTA...",   command=self._import_dta   ).pack(side=tk.LEFT, padx=2)
        ttk.Button(top, text="Importar RFB...",   command=self._import_rfb  ).pack(side=tk.LEFT, padx=2)
        ttk.Button(top, text="Importar outros...",command=self._import_multi).pack(side=tk.LEFT, padx=2)
        ttk.Button(top, text="Reset tudo",       command=self._reset_all  ).pack(side=tk.LEFT, padx=2)
        self.lbl_file = ttk.Label(top, text="Nenhum arquivo", foreground="gray")
        self.lbl_file.pack(side=tk.LEFT, padx=8)

        # Painel H: controles | plot
        pane = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        pane.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        left = ttk.Frame(pane, width=360)
        pane.add(left, weight=0)
        self._build_left(left)

        right = ttk.Frame(pane)
        pane.add(right, weight=1)
        self._build_plot(right)

        self.status = ttk.Label(self, text="Pronto.", relief=tk.SUNKEN, anchor=tk.W)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

    def _build_left(self, parent):
        nb = ttk.Notebook(parent)
        nb.pack(fill=tk.BOTH, expand=True)

        t_scans = ttk.Frame(nb, padding=4)
        t_edit  = ttk.Frame(nb, padding=4)
        t_gamma = ttk.Frame(nb, padding=4)

        nb.add(t_scans, text="  Scans  ")
        nb.add(t_edit,  text="  Editar  ")
        nb.add(t_gamma, text="  Gama  ")

        self._build_tab_scans(t_scans)
        self._build_tab_edit(t_edit)
        self._build_tab_gamma(t_gamma)

    # ======================================================================
    # Aba: Scans
    # ======================================================================

    def _build_tab_scans(self, p):
        # --- Filtros ---
        ff = ttk.LabelFrame(p, text="Filtros", padding=4)
        ff.pack(fill=tk.X, pady=(0, 4))
        ff.columnconfigure(1, weight=1)

        ttk.Label(ff, text="Feixe:").grid(row=0, column=0, sticky=tk.W)
        self.cb_filt_beam = ttk.Combobox(ff, state="readonly", width=13)
        self.cb_filt_beam.grid(row=0, column=1, sticky=tk.EW, padx=2, pady=1)
        self.cb_filt_beam.bind("<<ComboboxSelected>>", lambda _: self._refresh_scans_list())

        ttk.Label(ff, text="Tipo:").grid(row=1, column=0, sticky=tk.W)
        self.cb_filt_type = ttk.Combobox(ff, state="readonly", width=13)
        self.cb_filt_type.grid(row=1, column=1, sticky=tk.EW, padx=2, pady=1)
        self.cb_filt_type.bind("<<ComboboxSelected>>", lambda _: self._refresh_scans_list())

        ttk.Label(ff, text="Detector:").grid(row=2, column=0, sticky=tk.W)
        self.cb_filt_det = ttk.Combobox(ff, state="readonly", width=13)
        self.cb_filt_det.grid(row=2, column=1, sticky=tk.EW, padx=2, pady=1)
        self.cb_filt_det.bind("<<ComboboxSelected>>", lambda _: self._refresh_scans_list())

        ttk.Button(ff, text="Limpar filtros",
                   command=self._clear_filters).grid(row=3, column=0, columnspan=2, pady=2)

        # --- Lista de scans ---
        lf = ttk.LabelFrame(p, text="Scans (Ctrl+clique = múltipla seleção)", padding=3)
        lf.pack(fill=tk.BOTH, expand=True)
        lf.rowconfigure(0, weight=1)
        lf.columnconfigure(0, weight=1)

        frm, self.lb_scans = scrolled_listbox(lf, height=14)
        frm.grid(row=0, column=0, sticky=tk.NSEW)
        self.lb_scans.bind("<<ListboxSelect>>", self._on_scans_select)

        btn_f = ttk.Frame(p)
        btn_f.pack(fill=tk.X, pady=2)
        ttk.Button(btn_f, text="Selec. todos", command=lambda: self.lb_scans.select_set(0, tk.END)).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=1)
        ttk.Button(btn_f, text="Limpar seleção", command=lambda: self.lb_scans.selection_clear(0, tk.END)).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=1)

        opt = ttk.LabelFrame(p, text="Plot", padding=4)
        opt.pack(fill=tk.X, pady=2)
        self.var_norm = tk.BooleanVar(value=True)
        ttk.Checkbutton(opt, text="Normalizar (dmax / CAX = 100%)", variable=self.var_norm).pack(anchor=tk.W)
        self.var_color_by = tk.StringVar(value="detector")
        ttk.Label(opt, text="Colorir por:").pack(side=tk.LEFT)
        ttk.Radiobutton(opt, text="Detector", variable=self.var_color_by, value="detector").pack(side=tk.LEFT)
        ttk.Radiobutton(opt, text="Campo",    variable=self.var_color_by, value="field"   ).pack(side=tk.LEFT)

        ttk.Button(p, text=">> Plotar selecionados", command=self._plot_selected).pack(fill=tk.X, pady=3)

    # ======================================================================
    # Aba: Editar
    # ======================================================================

    def _build_tab_edit(self, p):
        # Envolve tudo num canvas rolável para não cortar conteúdo
        canvas = tk.Canvas(p, borderwidth=0, highlightthickness=0)
        vsb = ttk.Scrollbar(p, orient=tk.VERTICAL, command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        # frame interno que recebe todos os widgets
        p = ttk.Frame(canvas)
        win_id = canvas.create_window((0, 0), window=p, anchor="nw")
        def _on_resize(e):
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas.itemconfig(win_id, width=e.width)
        p.bind("<Configure>", _on_resize)
        canvas.bind("<Configure>", lambda e: canvas.itemconfig(win_id, width=e.width))
        # scroll com roda do mouse
        def _on_wheel(e):
            canvas.yview_scroll(int(-1 * (e.delta / 120)), "units")
        canvas.bind_all("<MouseWheel>", _on_wheel)

        # --- Seletor de scan (próprio da aba) ---
        sf = ttk.LabelFrame(p, text="Scan para editar", padding=3)
        sf.pack(fill=tk.X, pady=(0, 4))
        sf.columnconfigure(0, weight=1)

        frm, self.lb_edit = scrolled_listbox(sf, height=6)
        frm.pack(fill=tk.BOTH, expand=True)
        self.lb_edit.bind("<<ListboxSelect>>", self._on_edit_select)

        ttk.Label(sf, text="(Ctrl+clique = múltiplos; edições aplicam a todos selecionados)",
                  font=("", 7), foreground="gray").pack(anchor=tk.W)

        # Info do scan corrente
        self.lbl_edit_info = ttk.Label(p, text="", foreground="navy",
                                        font=("", 8), wraplength=340)
        self.lbl_edit_info.pack(anchor=tk.W, pady=(0, 4))

        # ---- Renomear / Metadados ----
        mf = ttk.LabelFrame(p, text="Renomear / Corrigir informações", padding=4)
        mf.pack(fill=tk.X, pady=(0, 4))
        mf.columnconfigure(1, weight=1)

        self.var_name     = labeled_entry(mf, 0, "Nome:",          width=22)
        self.var_beam     = labeled_entry(mf, 1, "Feixe:",         width=12)
        self.var_energy   = labeled_entry(mf, 2, "Energia:",       width=10)
        self.var_detector = labeled_entry(mf, 3, "Detector:",      width=18)

        # campo nominal: X e Y na mesma linha
        ttk.Label(mf, text="Campo (cm):").grid(row=4, column=0, sticky=tk.W, pady=1, padx=2)
        field_row = ttk.Frame(mf)
        field_row.grid(row=4, column=1, sticky=tk.EW, pady=1, padx=2)
        self.var_field_x = tk.StringVar()
        self.var_field_y = tk.StringVar()
        ttk.Entry(field_row, textvariable=self.var_field_x, width=6).pack(side=tk.LEFT)
        ttk.Label(field_row, text=" x ").pack(side=tk.LEFT)
        ttk.Entry(field_row, textvariable=self.var_field_y, width=6).pack(side=tk.LEFT)

        ttk.Button(mf, text="Aplicar informações",
                   command=self._apply_meta).grid(row=5, column=0, columnspan=2,
                                                   sticky=tk.EW, pady=4)

        # ---- Suavização ----
        sm = ttk.LabelFrame(p, text="Suavização", padding=4)
        sm.pack(fill=tk.X, pady=(0, 4))
        sm.columnconfigure(1, weight=1)

        self.cb_sm_method = labeled_combo(sm, 0, "Método:",
                                           ["savgol", "moving_avg", "gaussian"], "savgol")
        self.var_sm_win    = labeled_entry(sm, 1, "Janela (pts):",   "11", width=6)
        self.var_sm_poly   = labeled_entry(sm, 2, "Grau (savgol):",  "3",  width=6)
        self.var_sm_sigma  = labeled_entry(sm, 3, "Sigma (gauss):",  "1.0",width=6)

        rng = ttk.Frame(sm)
        rng.grid(row=4, column=0, columnspan=2, sticky=tk.EW, pady=1)
        ttk.Label(rng, text="Região (cm):").pack(side=tk.LEFT)
        self.var_sm_start = tk.StringVar()
        self.var_sm_end   = tk.StringVar()
        ttk.Entry(rng, textvariable=self.var_sm_start, width=7).pack(side=tk.LEFT, padx=2)
        ttk.Label(rng, text="até").pack(side=tk.LEFT)
        ttk.Entry(rng, textvariable=self.var_sm_end,   width=7).pack(side=tk.LEFT, padx=2)
        ttk.Label(rng, text="(vazio=tudo)", font=("", 7), foreground="gray").pack(side=tk.LEFT)

        ttk.Button(sm, text="Suavizar",
                   command=self._apply_smooth).grid(row=5, column=0, columnspan=2,
                                                     sticky=tk.EW, pady=3)

        # ---- Centralização ----
        cf = ttk.LabelFrame(p, text="Centralizar perfil / PDP", padding=4)
        cf.pack(fill=tk.X, pady=(0, 4))
        cf.columnconfigure(1, weight=1)

        self.cb_center = labeled_combo(cf, 0, "Método:",
                                        ["fwhm — ponto médio 50%",
                                         "cax — posição do pico",
                                         "dmax — deslocar dmax para Z=0 (PDP)"],
                                        "fwhm — ponto médio 50%")
        ttk.Button(cf, text="Centralizar",
                   command=self._apply_center).grid(row=1, column=0, columnspan=2,
                                                     sticky=tk.EW, pady=3)
        self.lbl_center_result = ttk.Label(cf, text="", foreground="gray", font=("", 8))
        self.lbl_center_result.grid(row=2, column=0, columnspan=2, sticky=tk.W)

        # ---- Deslocar posição manualmente ----
        df = ttk.LabelFrame(p, text="Deslocar posição (eixo X)", padding=4)
        df.pack(fill=tk.X, pady=(0, 4))
        df.columnconfigure(1, weight=1)

        self.var_pos_delta = labeled_entry(df, 0, "Delta (cm):", "0.0", width=8)
        ttk.Button(df, text="Aplicar deslocamento",
                   command=self._apply_pos_delta).grid(row=1, column=0, columnspan=2,
                                                        sticky=tk.EW, pady=3)
        self.lbl_pos_delta_result = ttk.Label(df, text="", foreground="gray", font=("", 8))
        self.lbl_pos_delta_result.grid(row=2, column=0, columnspan=2, sticky=tk.W)

        # ---- Renormalização ----
        rf = ttk.LabelFrame(p, text="Renormalizar", padding=4)
        rf.pack(fill=tk.X, pady=(0, 4))
        rf.columnconfigure(1, weight=1)

        self.cb_renorm = labeled_combo(rf, 0, "Método:",
                                        ["dmax", "cax", "point", "region"], "dmax")
        self.var_rn_val   = labeled_entry(rf, 1, "Valor ref. (%):", "100", width=7)
        self.var_rn_pos   = labeled_entry(rf, 2, "Posição (cm):",   "",    width=7)

        rng2 = ttk.Frame(rf)
        rng2.grid(row=3, column=0, columnspan=2, sticky=tk.EW, pady=1)
        ttk.Label(rng2, text="Região (cm):").pack(side=tk.LEFT)
        self.var_rn_start = tk.StringVar()
        self.var_rn_end   = tk.StringVar()
        ttk.Entry(rng2, textvariable=self.var_rn_start, width=7).pack(side=tk.LEFT, padx=2)
        ttk.Label(rng2, text="até").pack(side=tk.LEFT)
        ttk.Entry(rng2, textvariable=self.var_rn_end,   width=7).pack(side=tk.LEFT, padx=2)

        ttk.Button(rf, text="Renormalizar",
                   command=self._apply_renorm).grid(row=4, column=0, columnspan=2,
                                                     sticky=tk.EW, pady=3)

        # ---- Métricas ----
        mf = ttk.LabelFrame(p, text="Métricas (curva editada)", padding=4)
        mf.pack(fill=tk.X, pady=(0, 4))

        self.txt_metrics = tk.Text(mf, height=10, font=("Courier", 8),
                                    state=tk.DISABLED, wrap=tk.NONE,
                                    background="#f5f5f5", relief=tk.FLAT)
        sb = ttk.Scrollbar(mf, orient=tk.VERTICAL, command=self.txt_metrics.yview)
        self.txt_metrics.config(yscrollcommand=sb.set)
        self.txt_metrics.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        # ---- Aceitar / Reset ----
        btns = ttk.LabelFrame(p, text="Confirmar edição", padding=4)
        btns.pack(fill=tk.X, pady=4)

        ttk.Button(btns, text="✔  Aceitar — consolida edição como novo original",
                   command=self._accept_edit).pack(fill=tk.X, pady=1)
        ttk.Separator(btns).pack(fill=tk.X, pady=2)
        ttk.Button(btns, text="✖  Rejeitar / Reset selecionados",
                   command=self._reset_selected).pack(fill=tk.X, pady=1)
        ttk.Button(btns, text="✖  Reset TODOS",
                   command=self._reset_all).pack(fill=tk.X, pady=1)

    # ======================================================================
    # Aba: Gama
    # ======================================================================

    def _build_tab_gamma(self, p):
        g = ttk.LabelFrame(p, text="Índice Gama 1D", padding=6)
        g.pack(fill=tk.BOTH, expand=True)
        g.columnconfigure(1, weight=1)

        ttk.Label(g, text="Referência:").grid(row=0, column=0, sticky=tk.W)
        self.cb_gref = ttk.Combobox(g, state="readonly", width=24)
        self.cb_gref.grid(row=0, column=1, sticky=tk.EW, padx=2, pady=2)

        ttk.Label(g, text="Avaliação:").grid(row=1, column=0, sticky=tk.W)
        self.cb_geval = ttk.Combobox(g, state="readonly", width=24)
        self.cb_geval.grid(row=1, column=1, sticky=tk.EW, padx=2, pady=2)

        self.var_gdd    = labeled_entry(g, 2, "DD (%):",        "3.0", width=7)
        self.var_gdta   = labeled_entry(g, 3, "DTA (cm):",      "0.3", width=7)
        self.var_gthresh= labeled_entry(g, 4, "Threshold (%):", "10.0",width=7)

        self.cb_gnorm = labeled_combo(g, 5, "Normaliz.:", ["max", "local"], "max", width=8)

        ttk.Button(g, text=">> Calcular Gama",
                   command=self._calc_gamma).grid(row=6, column=0, columnspan=2,
                                                   sticky=tk.EW, pady=6)

        self.lbl_gresult = ttk.Label(g, text="", foreground="navy",
                                      font=("", 9, "bold"), wraplength=320)
        self.lbl_gresult.grid(row=7, column=0, columnspan=2, sticky=tk.W)

    # ======================================================================
    # Área de plot
    # ======================================================================

    def _build_plot(self, parent):
        self.fig = Figure(figsize=(9, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, parent).update()

    # ======================================================================
    # Carga de arquivo
    # ======================================================================

    def _open_file(self):
        p = filedialog.askopenfilename(
            title="Abrir .snctxt",
            filetypes=[("SNC Text", "*.snctxt"), ("Todos", "*.*")])
        if p:
            self._load_file(p)

    def _load_file(self, path):
        try:
            self.scans = parse_snctxt(path)
            self.filepath = path
            self.lbl_file.config(text=Path(path).name, foreground="black")
            self._populate_filters()
            self._refresh_scans_list()
            self._refresh_edit_list()
            self._refresh_gamma_combos()
            self._st(f"Carregado: {len(self.scans)} scans | {Path(path).name}")
        except Exception as e:
            messagebox.showerror("Erro", str(e))

    def _import_monaco(self):
        """Importa um ou mais arquivos Monaco .ALL e adiciona os perfis à lista."""
        paths = filedialog.askopenfilenames(
            title="Importar Monaco .ALL",
            filetypes=[("Monaco ALL", "*.ALL"), ("Todos", "*.*")])
        if not paths:
            return
        idx_start = max((s.index for s in self.scans), default=-1) + 1
        try:
            new_scans = load_monaco_files(list(paths), index_start=idx_start)
        except Exception as e:
            messagebox.showerror("Erro ao importar Monaco", str(e))
            return
        if not new_scans:
            messagebox.showinfo("Monaco", "Nenhum perfil extraido dos arquivos selecionados.")
            return
        self.scans.extend(new_scans)
        self._populate_filters()
        self._refresh_scans_list()
        self._refresh_edit_list()
        self._refresh_gamma_combos()
        nf = len(paths)
        ns = len(new_scans)
        self._st(f"Monaco: {nf} arquivo(s) importado(s) → {ns} perfil(is) adicionado(s).")

    def _import_mcc(self):
        """Importa um ou mais arquivos PTW Mephisto .mcc."""
        paths = filedialog.askopenfilenames(
            title="Importar PTW Mephisto .mcc",
            filetypes=[("PTW Mephisto", "*.mcc"), ("Todos", "*.*")])
        if not paths:
            return
        idx_start = max((s.index for s in self.scans), default=-1) + 1
        try:
            new_scans = load_mcc_files(list(paths), index_start=idx_start)
        except Exception as e:
            messagebox.showerror("Erro ao importar MCC", str(e))
            return
        if not new_scans:
            messagebox.showinfo("MCC", "Nenhum perfil extraído dos arquivos selecionados.")
            return
        self.scans.extend(new_scans)
        self._populate_filters()
        self._refresh_scans_list()
        self._refresh_edit_list()
        self._refresh_gamma_combos()
        nf = len(paths)
        ns = len(new_scans)
        self._st(f"MCC: {nf} arquivo(s) importado(s) → {ns} perfil(is) adicionado(s).")

    def _import_dta(self):
        """Importa um ou mais arquivos DTA (IBA/Scanditronix) .dta / .asc."""
        paths = filedialog.askopenfilenames(
            title="Importar DTA / ASC",
            filetypes=[("DTA / ASC", "*.dta *.asc"), ("DTA", "*.dta"),
                       ("ASC", "*.asc"), ("Todos", "*.*")])
        if not paths:
            return
        idx_start = max((s.index for s in self.scans), default=-1) + 1
        try:
            new_scans = load_dta_files(list(paths), index_start=idx_start)
        except Exception as e:
            messagebox.showerror("Erro ao importar DTA", str(e))
            return
        if not new_scans:
            messagebox.showinfo("DTA", "Nenhum perfil extraído dos arquivos selecionados.")
            return
        self.scans.extend(new_scans)
        self._populate_filters()
        self._refresh_scans_list()
        self._refresh_edit_list()
        self._refresh_gamma_combos()
        nf = len(paths)
        ns = len(new_scans)
        self._st(f"DTA: {nf} arquivo(s) importado(s) → {ns} perfil(is) adicionado(s).")

    def _import_rfb(self):
        """Importa arquivos OmniPro V6 binários .rfb (IBA Dosimetry)."""
        paths = filedialog.askopenfilenames(
            title="Importar OmniPro V6 RFB",
            filetypes=[("OmniPro RFB", "*.rfb"), ("Todos", "*.*")])
        if not paths:
            return
        idx_start = max((s.index for s in self.scans), default=-1) + 1
        try:
            new_scans = load_rfb_files(list(paths), index_start=idx_start)
        except Exception as e:
            messagebox.showerror("Erro ao importar RFB", str(e))
            return
        if not new_scans:
            messagebox.showinfo("RFB", "Nenhum perfil extraído dos arquivos selecionados.")
            return
        self.scans.extend(new_scans)
        self._populate_filters()
        self._refresh_scans_list()
        self._refresh_edit_list()
        self._refresh_gamma_combos()
        nf, ns = len(paths), len(new_scans)
        self._st(f"RFB: {nf} arquivo(s) importado(s) → {ns} perfil(is) adicionado(s).")

    def _import_multi(self):
        """Importa arquivos em múltiplos formatos com auto-detecção:
        Eclipse w2CAD (.asc/.txt), OmniPro RFA300 (.txt/.asc),
        RayStation (.csv), SNC Water Tank (.sncxml), DoseView 1D (.csv)."""
        paths = filedialog.askopenfilenames(
            title="Importar — w2CAD / OmniPro / RayStation / SNC XML / DV1D",
            filetypes=[
                ("Todos os suportados", "*.asc *.txt *.csv *.sncxml"),
                ("Eclipse w2CAD / OmniPro (.asc)", "*.asc"),
                ("Eclipse Line Profile / OmniPro (.txt)", "*.txt"),
                ("RayStation / DoseView 1D (.csv)", "*.csv"),
                ("SNC Water Tank (.sncxml)", "*.sncxml"),
                ("Todos", "*.*"),
            ])
        if not paths:
            return
        idx_start = max((s.index for s in self.scans), default=-1) + 1
        errors = []
        new_scans = []
        for p in paths:
            try:
                new_scans.extend(load_multi([p], index_start=idx_start + len(new_scans)))
            except Exception as e:
                errors.append(f"{Path(p).name}: {e}")
        if errors:
            messagebox.showwarning("Aviso — erros de importação",
                                   "\n".join(errors))
        if not new_scans:
            messagebox.showinfo("Importar outros", "Nenhum perfil extraído.")
            return
        self.scans.extend(new_scans)
        self._populate_filters()
        self._refresh_scans_list()
        self._refresh_edit_list()
        self._refresh_gamma_combos()
        nf = len(paths)
        ns = len(new_scans)
        self._st(f"Outros: {nf} arquivo(s) → {ns} perfil(is) adicionado(s)."
                 + (f"  {len(errors)} erro(s)." if errors else ""))

    # ======================================================================
    # Filtros e listas
    # ======================================================================

    def _populate_filters(self):
        beams = sorted({s.beam_type for s in self.scans})
        types = sorted({s.scan_type for s in self.scans})
        dets  = sorted({(s.detector_override or s.detector) for s in self.scans})

        self.cb_filt_beam["values"] = ["(todos)"] + beams
        self.cb_filt_type["values"] = ["(todos)"] + types
        self.cb_filt_det["values"]  = ["(todos)"] + [d for d in dets if d]

        self.cb_filt_beam.set("(todos)")
        self.cb_filt_type.set("(todos)")
        self.cb_filt_det.set("(todos)")

    def _clear_filters(self):
        self.cb_filt_beam.set("(todos)")
        self.cb_filt_type.set("(todos)")
        self.cb_filt_det.set("(todos)")
        self._refresh_scans_list()

    def _apply_filters(self) -> list[Scan]:
        beam = self.cb_filt_beam.get()
        stype= self.cb_filt_type.get()
        det  = self.cb_filt_det.get()
        out  = self.scans
        if beam  and beam  != "(todos)": out = [s for s in out if s.beam_type == beam]
        if stype and stype != "(todos)": out = [s for s in out if s.scan_type == stype]
        if det   and det   != "(todos)": out = [s for s in out if (s.detector_override or s.detector) == det]
        return out

    def _scan_row(self, s: Scan) -> str:
        det   = s.detector_override or s.detector or "—"
        name  = s.name_override if s.name_override else f"{s.beam_type} {s.display_energy}MeV/MV"
        beam  = s.beam_override if hasattr(s, "beam_override") and s.beam_override else s.beam_type
        edited= "*" if (s._dose_edit is not None or s._pos_offset != 0.0
                        or s.name_override or s.energy_override
                        or s.detector_override) else " "
        fs = f"{s.display_field_x}x{s.display_field_y}"
        return (f"[{s.index:>2}]{edited}{name:<17} "
                f"{fs:>7}cm  "
                f"{s.scan_type:<12}  "
                f"{det}")

    def _refresh_scans_list(self):
        selected_indices = {self._scans_filtered[i].index
                            for i in self.lb_scans.curselection()
                            if i < len(self._scans_filtered)}
        self._scans_filtered = self._apply_filters()
        self.lb_scans.delete(0, tk.END)
        for pos, s in enumerate(self._scans_filtered):
            self.lb_scans.insert(tk.END, self._scan_row(s))
            if s.index in selected_indices:
                self.lb_scans.selection_set(pos)

    def _refresh_edit_list(self):
        # preserva seleção pelos scan.index antes de reconstruir a lista
        selected_indices = {self.scans[i].index
                            for i in self.lb_edit.curselection()
                            if i < len(self.scans)}
        self.lb_edit.delete(0, tk.END)
        for pos, s in enumerate(self.scans):
            self.lb_edit.insert(tk.END, self._scan_row(s))
            if s.index in selected_indices:
                self.lb_edit.selection_set(pos)

    def _refresh_gamma_combos(self):
        rows = [self._scan_row(s) for s in self.scans]
        self.cb_gref["values"]  = rows
        self.cb_geval["values"] = rows
        if len(rows) > 0: self.cb_gref.current(0)
        if len(rows) > 1: self.cb_geval.current(1)

    def _selected_in_scans(self) -> list[Scan]:
        return [self._scans_filtered[i] for i in self.lb_scans.curselection()]

    def _selected_in_edit(self) -> list[Scan]:
        return [self.scans[i] for i in self.lb_edit.curselection()]

    def _on_scans_select(self, _=None):
        sel = self._selected_in_scans()
        if sel:
            self._st(f"{len(sel)} scan(s) selecionado(s)")

    def _on_edit_select(self, _=None):
        sel = self._selected_in_edit()
        if not sel:
            self.lbl_edit_info.config(text="")
            return
        s = sel[-1]
        det = s.detector_override or s.detector or "—"
        self.lbl_edit_info.config(
            text=(f"[{s.index}]  {s.beam_type} {s.display_energy}MeV/MV  "
                  f"{s.field_x}x{s.field_y}cm  {s.scan_type}  "
                  f"SSD={s.ssd}cm  N={len(s.dose)}  Det: {det}  "
                  f"Pos: {s.position.min():.2f}…{s.position.max():.2f}cm")
        )
        self.var_name.set(s.name_override)
        self.var_beam.set(s.beam_type)
        self.var_energy.set(s.energy_override or s.energy)
        self.var_detector.set(s.detector_override or s.detector)
        self.var_field_x.set(s.field_x_override or s.field_x)
        self.var_field_y.set(s.field_y_override or s.field_y)
        # mostra curva e métricas imediatamente
        self._plot_edit_preview()
        self._update_metrics()

    # ======================================================================
    # Plot
    # ======================================================================

    def _plot_selected(self):
        sel = self._selected_in_scans()
        if not sel:
            messagebox.showinfo("Aviso", "Selecione ao menos uma curva.")
            return

        norm     = self.var_norm.get()
        color_by = self.var_color_by.get()

        self.fig.clear()
        from collections import defaultdict
        groups: dict[str, list[Scan]] = defaultdict(list)
        for s in sel:
            groups[s.scan_type].append(s)

        n     = len(groups)
        ncols = min(3, n)
        nrows = (n + ncols - 1) // ncols
        axes  = self.fig.subplots(nrows, ncols)
        if n == 1: axes = np.array([[axes]])
        elif nrows == 1: axes = np.array([axes])
        axes_flat = axes.flatten()

        for ax_idx, (scan_t, slist) in enumerate(sorted(groups.items())):
            ax = axes_flat[ax_idx]
            for s in slist:
                pos  = s.position
                dose = s.dose_display
                if norm:
                    dose = _normalize(dose, pos, scan_t)

                if color_by == "detector":
                    color = det_color(s.detector_override or s.detector)
                else:
                    fs = f"{s.field_x}x{s.field_y}"
                    color = det_color(fs)

                det   = s.detector_override or s.detector or "—"
                lbl   = s.name_override or f"[{s.index}] {s.field_x}x{s.field_y}cm"
                if s._dose_edit is not None or s._pos_offset != 0.0:
                    lbl += "*"
                lbl += f"  {det}"

                ax.plot(pos, dose, label=lbl, color=color, linewidth=1.6)

            xlabel = {"Depth Scan": "Profundidade Z (cm)",
                      "Crossline":  "Posição X (cm)",
                      "Inline":     "Posição Y (cm)"}.get(scan_t, "Posição (cm)")
            ax.set_title(f"{scan_t}", fontsize=9, fontweight="bold")
            ax.set_xlabel(xlabel, fontsize=8)
            ax.set_ylabel("Dose (%)" + (" norm." if norm else ""), fontsize=8)
            ax.legend(fontsize=7, loc="best")
            ax.grid(True, alpha=0.3)
            if norm:
                ax.axhline(100, color="gray", lw=0.7, ls="--")
                ax.axhline(50,  color="gray", lw=0.5, ls=":")
            if scan_t != "Depth Scan":
                ax.axvline(0, color="gray", lw=0.7, ls=":")

        for ax in axes_flat[n:]: ax.set_visible(False)
        self.fig.suptitle(Path(self.filepath).stem, fontsize=10)
        self.fig.tight_layout()
        self.canvas.draw()
        self._st(f"Plotadas {len(sel)} curva(s).")

    # ======================================================================
    # Preview de edição (original vs editado)
    # ======================================================================

    def _annotate_metrics_on_ax(self, ax, pos, dose_n, scan_type, color, label=""):
        """Adiciona marcadores gráficos e caixa de métricas ao eixo do preview."""
        try:
            if scan_type == "Depth Scan":
                m = compute_pdp_metrics(pos, dose_n)

                # linha vertical no Zmax
                ax.axvline(m['zmax'], color=color, lw=1.0, ls=':', alpha=0.8)
                ax.text(m['zmax'] + 0.05, 103,
                        f"Zmax={m['zmax']:.2f}cm",
                        ha='left', va='bottom', fontsize=5.5, color=color)

                # linhas nos pontos de queda + marcadores verticais
                for level, key, lbl in [(90,'r90','R90'), (50,'r50','R50'), (10,'r10','R10')]:
                    if m[key] is not None:
                        ax.plot(m[key], level, 'o', color=color, ms=4, zorder=5)
                        ax.text(m[key] + 0.1, level, f"{lbl}={m[key]:.2f}",
                                ha='left', va='center', fontsize=5.5, color=color)

                # caixa de texto consolidada
                def _f(v): return "n.d." if v is None else f"{v:.2f} cm"
                rows = [
                    f"  Zmax={_f(m['zmax'])}",
                    f"  R90 ={_f(m['r90'])}",
                    f"  R80 ={_f(m['r80'])}",
                    f"  R50 ={_f(m['r50'])}",
                    f"  R10 ={_f(m['r10'])}",
                    f"  Rp  ={_f(m['rp'])}",
                ]
                hdr = label if label else "Métricas"
                txt = hdr + "\n" + "\n".join(rows)

            else:  # perfil
                m = compute_profile_metrics(pos, dose_n)

                # linhas verticais nos bordos 50% (FWHM)
                if m['l50'] is not None and m['r50'] is not None:
                    ax.axvline(m['l50'], color=color, lw=1.0, ls='--', alpha=0.7)
                    ax.axvline(m['r50'], color=color, lw=1.0, ls='--', alpha=0.7)
                    mid = (m['l50'] + m['r50']) / 2.0
                    # seta dupla no nível 50%
                    ax.annotate('', xy=(m['r50'], 50), xytext=(m['l50'], 50),
                                arrowprops=dict(arrowstyle='<->', color=color,
                                                lw=1.0, mutation_scale=8))
                    ax.text(mid, 52, f"FWHM={m['fwhm']:.2f}cm",
                            ha='center', va='bottom', fontsize=5.5, color=color)

                # linhas 80%/20% para penumbra
                for lvl in (80, 20):
                    ax.axhline(lvl, color='gray', lw=0.4, ls=':', alpha=0.6)

                def _fc(v): return "n.d." if v is None else f"{v:.3f} cm"
                def _fm(v): return "n.d." if v is None else f"{v*10:.2f} mm"
                def _fp(v): return "n.d." if v is None else f"{v:.2f} %"
                rows = [
                    f"  FWHM  ={_fc(m['fwhm'])}",
                    f"  Centro={_fc(m['center'])}",
                    f"  PenL  ={_fm(m['penumbra_l'])}",
                    f"  PenR  ={_fm(m['penumbra_r'])}",
                    f"  Flat  ={_fp(m['flatness'])}",
                    f"  Sim   ={_fp(m['symmetry'])}",
                ]
                hdr = label if label else "Métricas"
                txt = hdr + "\n" + "\n".join(rows)

            ax.text(0.99, 0.99, txt, transform=ax.transAxes,
                    fontsize=5.5, va='top', ha='right', family='monospace',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.82,
                              edgecolor=color, lw=0.7))
        except Exception:
            pass  # nunca interrompe o plot por erro em métricas

    def _plot_edit_preview(self):
        """Plota original (cinza tracejado) vs editado (colorido sólido) para os
        scans selecionados na aba Editar. Painel de diferença aparece se há edição."""
        sel = self._selected_in_edit()
        if not sel:
            return

        self.fig.clear()

        # Agrupa por scan_type para separar subplots
        from collections import defaultdict
        groups: dict[str, list[Scan]] = defaultdict(list)
        for s in sel:
            groups[s.scan_type].append(s)

        n_types = len(groups)
        # Layout: 2 linhas (curvas | diferença) × n_types colunas
        ncols = min(3, n_types)
        nrows_g = (n_types + ncols - 1) // ncols
        # gridspec: (nrows_g curvas + nrows_g diffs) × ncols
        total_rows = nrows_g * 2
        gs = self.fig.add_gridspec(total_rows, ncols, hspace=0.45, wspace=0.35,
                                   height_ratios=[2.5, 1] * nrows_g)

        xlabel_map = {"Depth Scan": "Profundidade Z (cm)",
                      "Crossline":  "Posição X (cm)",
                      "Inline":     "Posição Y (cm)"}

        for col_idx, (scan_t, slist) in enumerate(sorted(groups.items())):
            row_curve = (col_idx // ncols) * 2
            row_diff  = row_curve + 1
            c         = col_idx % ncols
            ax_curve  = self.fig.add_subplot(gs[row_curve, c])
            ax_diff   = self.fig.add_subplot(gs[row_diff,  c], sharex=ax_curve)

            has_any_edit = any(s._dose_edit is not None or s._pos_offset != 0.0
                               for s in slist)

            for s in slist:
                pos_orig = s.position     # já ordenado, com offset
                d_orig   = s.dose_sorted  # dose original (sem edições)
                d_edit   = s.dose_display # dose com edições (= dose_sorted se sem edição)

                d_orig_n = _normalize(d_orig, pos_orig, scan_t)
                d_edit_n = _normalize(d_edit, pos_orig, scan_t)

                det = s.detector_override or s.detector or "—"
                lbl = s.name_override or f"[{s.index}] {s.field_x}x{s.field_y}cm"
                color = det_color(det)

                # original: sempre cinza tracejado
                ax_curve.plot(pos_orig, d_orig_n,
                              color="gray", lw=1.2, ls="--", alpha=0.6,
                              label=f"{lbl} original")
                # editado: colorido sólido (sobreposto; idêntico se sem edição)
                ax_curve.plot(pos_orig, d_edit_n,
                              color=color, lw=2.0, ls="-",
                              label=f"{lbl} editado  [{det}]")

                # anotações de métricas sobre a curva editada
                short = s.name_override or f"[{s.index}]"
                self._annotate_metrics_on_ax(ax_curve, pos_orig, d_edit_n,
                                             scan_t, color, label=short)

                # diferença (editado - original)
                diff = d_edit_n - d_orig_n
                ax_diff.plot(pos_orig, diff, color=color, lw=1.2)

            ax_curve.axhline(100, color="gray", lw=0.6, ls="--")
            ax_curve.axhline(50,  color="gray", lw=0.5, ls=":")
            if scan_t != "Depth Scan":
                ax_curve.axvline(0, color="gray", lw=0.6, ls=":")
            ax_curve.set_ylabel("Dose (%) norm.")
            ax_curve.set_title(
                scan_t + (" — com edição*" if has_any_edit else " — sem edição"),
                fontsize=9, fontweight="bold")
            ax_curve.legend(fontsize=6, loc="best")
            ax_curve.grid(True, alpha=0.3)

            ax_diff.axhline(0, color="black", lw=0.8)
            ax_diff.set_ylabel("Δ (%)")
            ax_diff.set_xlabel(xlabel_map.get(scan_t, "Posição (cm)"), fontsize=8)
            ax_diff.grid(True, alpha=0.3)
            if not has_any_edit:
                ax_diff.set_title("(sem diferença)", fontsize=7, color="gray")

        self.fig.suptitle("Preview: original (---) vs editado (—)", fontsize=9)
        self.canvas.draw()

    # ======================================================================
    # Edição
    # ======================================================================

    def _edit_targets(self) -> list[Scan]:
        sel = self._selected_in_edit()
        if not sel:
            messagebox.showinfo("Aviso",
                "Selecione ao menos um scan na lista 'Scan para editar' (aba Editar).")
        return sel

    # ---- Metadados ----

    def _apply_meta(self):
        targets = self._edit_targets()
        if not targets: return
        name     = self.var_name.get().strip()
        beam     = self.var_beam.get().strip()
        energy   = self.var_energy.get().strip()
        detector = self.var_detector.get().strip()
        field_x  = self.var_field_x.get().strip()
        field_y  = self.var_field_y.get().strip()
        for s in targets:
            s.name_override     = name
            s.beam_type         = beam or s.beam_type
            s.energy_override   = energy
            s.detector_override = detector
            s.field_x_override  = field_x
            s.field_y_override  = field_y
        self._refresh_scans_list()
        self._refresh_edit_list()
        self._refresh_gamma_combos()
        self._populate_filters()
        self._plot_edit_preview()
        self._st(f"Informações atualizadas em {len(targets)} scan(s).")

    # ---- Suavização ----

    def _apply_smooth(self):
        targets = self._edit_targets()
        if not targets: return
        method = self.cb_sm_method.get()
        try:
            window  = int(self.var_sm_win.get())
            polyord = int(self.var_sm_poly.get())
            sigma   = float(self.var_sm_sigma.get())
        except ValueError:
            messagebox.showerror("Erro", "Janela/grau/sigma inválidos."); return

        region = None
        s0, s1 = self.var_sm_start.get().strip(), self.var_sm_end.get().strip()
        if s0 and s1:
            try: region = (float(s0), float(s1))
            except ValueError:
                messagebox.showerror("Erro", "Região inválida."); return

        errs = []
        for s in targets:
            try:
                apply_smooth(s, method=method, window=window,
                             polyorder=polyord, sigma=sigma, region=region)
            except Exception as e:
                errs.append(f"[{s.index}]: {e}")

        self._refresh_scans_list()
        self._refresh_edit_list()
        if errs: messagebox.showwarning("Aviso", "\n".join(errs))
        else: self._st(f"Suavização ({method}) aplicada em {len(targets)} scan(s).")
        self._plot_edit_preview()
        self._update_metrics()

    # ---- Centralização ----

    def _apply_center(self):
        targets = self._edit_targets()
        if not targets: return
        raw = self.cb_center.get()
        method = "dmax" if "dmax" in raw else ("cax" if "cax" in raw else "fwhm")

        offsets, errs = [], []
        for s in targets:
            try:
                off = apply_centering(s, method=method)
                offsets.append(f"[{s.index}] {off:+.4f}cm")
            except Exception as e:
                errs.append(f"[{s.index}]: {e}")

        self._refresh_scans_list()
        self._refresh_edit_list()
        result = "  ".join(offsets)
        self.lbl_center_result.config(text=result)
        if errs: messagebox.showwarning("Aviso", "\n".join(errs))
        else: self._st("Centralização: " + result)
        self._plot_edit_preview()
        self._update_metrics()

    # ---- Deslocamento manual de posição ----

    def _apply_pos_delta(self):
        targets = self._edit_targets()
        if not targets: return
        try:
            delta = float(self.var_pos_delta.get())
        except ValueError:
            messagebox.showerror("Erro", "Valor de delta inválido."); return

        for s in targets:
            s._pos_offset += delta

        self._refresh_scans_list()
        self._refresh_edit_list()
        results = "  ".join(f"[{s.index}] offset={s._pos_offset:+.3f}cm" for s in targets)
        self.lbl_pos_delta_result.config(text=results)
        self._plot_edit_preview()
        self._update_metrics()
        self._st(f"Delta {delta:+.3f}cm aplicado em {len(targets)} scan(s).")

    # ---- Renormalização ----

    def _apply_renorm(self):
        targets = self._edit_targets()
        if not targets: return
        method = self.cb_renorm.get()
        try:
            value = float(self.var_rn_val.get())
        except ValueError:
            messagebox.showerror("Erro", "Valor inválido."); return

        ref_pos = None
        txt = self.var_rn_pos.get().strip()
        if txt:
            try: ref_pos = float(txt)
            except ValueError:
                messagebox.showerror("Erro", "Posição inválida."); return

        region = None
        s0, s1 = self.var_rn_start.get().strip(), self.var_rn_end.get().strip()
        if s0 and s1:
            try: region = (float(s0), float(s1))
            except ValueError:
                messagebox.showerror("Erro", "Região inválida."); return

        errs = []
        for s in targets:
            try:
                apply_renormalize(s, method=method, value=value,
                                  ref_pos=ref_pos, region=region)
            except Exception as e:
                errs.append(f"[{s.index}]: {e}")

        self._refresh_scans_list()
        self._refresh_edit_list()
        if errs: messagebox.showwarning("Aviso", "\n".join(errs))
        else: self._st(f"Renormalização ({method}) aplicada em {len(targets)} scan(s).")
        self._plot_edit_preview()
        self._update_metrics()

    # ---- Reset ----

    def _accept_edit(self):
        """Consolida a dose editada como novo 'original' do scan.
        A partir daqui, reset não volta ao dado do arquivo — apenas à última
        versão aceita. Útil para edições iterativas (suavizar → aceitar → renormalizar)."""
        targets = self._edit_targets()
        if not targets: return
        count = 0
        for s in targets:
            if s._dose_edit is not None:
                # substitui dose_sorted pelo estado editado atual
                order = s._sort_order
                # reordena o array de dose original para alinhamento
                new_dose = np.empty_like(s.dose)
                new_dose[order] = s._dose_edit   # mapeia de volta à ordem original
                s.dose = new_dose
                s._dose_edit = None              # limpa edição pendente
                count += 1
        self._refresh_scans_list()
        self._refresh_edit_list()
        self._plot_edit_preview()
        self._update_metrics()
        self._st(f"Edição aceita e consolidada em {count} scan(s). Reset agora parte deste estado.")

    def _reset_selected(self):
        targets = self._edit_targets()
        for s in targets:
            s.reset_edits()
            s.name_override = s.energy_override = s.detector_override = ""
            s.field_x_override = s.field_y_override = ""
        self._refresh_scans_list()
        self._refresh_edit_list()
        self._refresh_gamma_combos()
        self._plot_edit_preview()
        self._update_metrics()
        self._st(f"Reset em {len(targets)} scan(s).")

    def _reset_all(self):
        for s in self.scans:
            s.reset_edits()
            s.name_override = s.energy_override = s.detector_override = ""
            s.field_x_override = s.field_y_override = ""
        self._refresh_scans_list()
        self._refresh_edit_list()
        self._refresh_gamma_combos()
        self._plot_edit_preview()
        self._st("Reset completo.")

    # ======================================================================
    # Gama
    # ======================================================================

    def _calc_gamma(self):
        if not self.scans:
            messagebox.showinfo("Aviso", "Carregue um arquivo primeiro."); return
        ri = self.cb_gref.current()
        ei = self.cb_geval.current()
        if ri < 0 or ei < 0:
            messagebox.showinfo("Aviso", "Selecione referência e avaliação."); return
        if ri == ei:
            messagebox.showinfo("Aviso", "Ref e eval devem ser diferentes."); return
        try:
            dd    = float(self.var_gdd.get())
            dta   = float(self.var_gdta.get())
            thr   = float(self.var_gthresh.get())
        except ValueError:
            messagebox.showerror("Erro", "DD/DTA/Threshold inválidos."); return

        sr = self.scans[ri]
        se = self.scans[ei]
        if sr.scan_type != se.scan_type:
            if not messagebox.askyesno("Tipos diferentes",
                f"Ref: {sr.scan_type}\nEval: {se.scan_type}\nContinuar?"):
                return

        pos_g, g = gamma_1d(sr.position, sr.dose_display,
                             se.position, se.dose_display,
                             dd=dd, dta=dta,
                             norm=self.cb_gnorm.get(), threshold=thr)

        pr   = gamma_pass_rate(g)
        gmax = np.max(g)  if len(g) else float("nan")
        gmu  = np.mean(g) if len(g) else float("nan")
        res  = f"Pass rate: {pr:.1f}%   γ_max={gmax:.3f}   γ_mean={gmu:.3f}   N={len(g)}"
        self.lbl_gresult.config(text=res)
        self._st(f"Gama: {res}")
        self._plot_gamma(sr, se, pos_g, g, dd, dta, pr, gmax, gmu)

    def _plot_gamma(self, sr, se, pos_g, g, dd, dta, pr, gmax, gmu):
        scan_t = sr.scan_type
        self.fig.clear()
        gs  = self.fig.add_gridspec(3, 1, height_ratios=[2.5, 1, 1.5], hspace=0.38)
        ax1 = self.fig.add_subplot(gs[0])
        ax2 = self.fig.add_subplot(gs[1], sharex=ax1)
        ax3 = self.fig.add_subplot(gs[2], sharex=ax1)

        pr_  = sr.position;  dr_  = _normalize(sr.dose_display,  pr_,  scan_t)
        pe_  = se.position;  de_  = _normalize(se.dose_display,  pe_,  scan_t)

        xlabel = {"Depth Scan": "Profundidade Z (cm)",
                  "Crossline":  "Posição X (cm)",
                  "Inline":     "Posição Y (cm)"}.get(scan_t, "Posição (cm)")

        det_r = sr.detector_override or sr.detector or "—"
        det_e = se.detector_override or se.detector or "—"
        lbl_r = (sr.name_override or f"[{sr.index}] {sr.field_x}x{sr.field_y}cm") + f"  {det_r}"
        lbl_e = (se.name_override or f"[{se.index}] {se.field_x}x{se.field_y}cm") + f"  {det_e}"

        ax1.plot(pr_, dr_, "b-",  lw=2,   label=lbl_r)
        ax1.plot(pe_, de_, "r--", lw=1.5, label=lbl_e)
        ax1.axhline(100, color="gray", lw=0.7, ls="--")
        ax1.axhline(50,  color="gray", lw=0.5, ls=":")
        ax1.set_ylabel("Dose (%) norm.")
        ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3)
        ax1.set_title(
            f"Gama 1D  |  {sr.beam_type} {sr.display_energy}MeV/MV  "
            f"|  {dd}% / {dta*10:.0f}mm  |  {scan_t}",
            fontsize=9, fontweight="bold")

        de_i = np.interp(pr_, pe_, de_)
        diff = de_i - dr_
        bw   = np.diff(pr_).mean() * 0.9 if len(pr_) > 1 else 0.1
        ax2.bar(pr_, diff, width=bw,
                color=["steelblue" if abs(d) <= dd else "tomato" for d in diff],
                alpha=0.75)
        ax2.axhline(0,   color="black",  lw=0.8)
        ax2.axhline( dd, color="orange", lw=1, ls="--")
        ax2.axhline(-dd, color="orange", lw=1, ls="--")
        ax2.set_ylabel("Diff (%)"); ax2.grid(True, alpha=0.3)

        if len(g):
            ax3.scatter(pos_g, g,
                        c=["limegreen" if v <= 1 else "crimson" for v in g],
                        s=10, zorder=3)
        ax3.axhline(1.0, color="black", lw=1.2, ls="--")
        ax3.set_ylabel("γ"); ax3.set_xlabel(xlabel)
        ax3.set_ylim(0, max(2.0, gmax * 1.15) if len(g) else 2)
        ax3.grid(True, alpha=0.3)
        n_pass = int(np.sum(g <= 1.0)) if len(g) else 0
        ax3.set_title(f"Pass rate: {pr:.1f}%  ({n_pass}/{len(g)})  "
                      f"γ_max={gmax:.3f}  γ_mean={gmu:.3f}", fontsize=8)
        self.canvas.draw()

    # ======================================================================
    # Utilidades
    # ======================================================================

    def _update_metrics(self):
        """Atualiza o painel de métricas para os scans selecionados na aba Editar."""
        sel = self._selected_in_edit()
        self.txt_metrics.config(state=tk.NORMAL)
        self.txt_metrics.delete("1.0", tk.END)
        if not sel:
            self.txt_metrics.config(state=tk.DISABLED)
            return
        for s in sel:
            name = s.name_override or f"[{s.index}] {s.field_x}x{s.field_y}cm {s.scan_type}"
            det  = s.detector_override or s.detector or "—"
            self.txt_metrics.insert(tk.END, f"{name}  |  {det}\n", "header")
            try:
                self.txt_metrics.insert(tk.END, format_metrics(s) + "\n\n")
            except Exception as e:
                self.txt_metrics.insert(tk.END, f"  Erro: {e}\n\n")
        self.txt_metrics.tag_config("header", font=("Courier", 8, "bold"),
                                     foreground="navy")
        self.txt_metrics.config(state=tk.DISABLED)

    def _st(self, msg):
        self.status.config(text=msg)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="")
    args = ap.parse_args()
    fp = args.file
    if fp and not Path(fp).is_absolute():
        fp = str(Path(__file__).parent / fp)
    SRSApp(filepath=fp).mainloop()

if __name__ == "__main__":
    main()

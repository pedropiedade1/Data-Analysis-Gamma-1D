# SRS Dosimetry Analysis Tool

Ferramenta gráfica em Python para análise, edição e comparação de curvas
dosimétricas de fantoma de água.  Desenvolvida para uso clínico em física
médica: comissionamento de aceleradores lineares, QA de SRS/SBRT e
comparação medida × cálculo.

---

## Índice

1. [Requisitos](#1-requisitos)
2. [Execução](#2-execução)
3. [Formatos de arquivo suportados](#3-formatos-de-arquivo-suportados)
4. [Interface gráfica — visão geral](#4-interface-gráfica--visão-geral)
   - 4.1 [Barra de ferramentas](#41-barra-de-ferramentas)
   - 4.2 [Aba Scans](#42-aba-scans)
   - 4.3 [Aba Editar](#43-aba-editar)
   - 4.4 [Aba Gama](#44-aba-gama)
   - 4.5 [Área de gráfico](#45-área-de-gráfico)
5. [Módulos — referência técnica](#5-módulos--referência-técnica)
   - [parse_snctxt.py](#parse_snctxtpy)
   - [parse_monaco.py](#parse_monacopy)
   - [parse_mcc.py](#parse_mccpy)
   - [parse_dta.py](#parse_dtapy)
   - [parse_rfb.py](#parse_rfbpy)
   - [parse_multiformats.py](#parse_multiformatspy)
   - [srs_ui.py](#srs_uipy)
6. [Estrutura de dados — classe `Scan`](#6-estrutura-de-dados--classe-scan)
7. [Métricas calculadas](#7-métricas-calculadas)
8. [Índice Gama 1D](#8-índice-gama-1d)
9. [Fluxo de trabalho típico](#9-fluxo-de-trabalho-típico)

---

## 1. Requisitos

| Pacote | Versão mínima |
|--------|--------------|
| Python | 3.10 |
| numpy  | 1.23 |
| matplotlib | 3.5 |
| scipy  | 1.9 (para suavização Savitzky-Golay) |
| tkinter | incluído no Python padrão |

Instalação das dependências:

```bash
pip install numpy matplotlib scipy
```

---

## 2. Execução

### Executável Windows (sem instalar Python)

Um executável pré-compilado para Windows está disponível na pasta `dist/`
do repositório e também na página de
[Releases](https://github.com/pedropiedade1/Data-Analysis-Gamma-1D/releases).

Basta baixar e executar `SRS_Analysis.exe` — sem instalação, sem Python.

### Via Python (desenvolvimento)

```bash
# sem arquivo pré-carregado
py -3 srs_ui.py

# com arquivo snctxt inicial
py -3 srs_ui.py --file "Curso SRS.snctxt"
```

### Gerar o executável localmente

```bash
pip install pyinstaller
pyinstaller --onefile --windowed --name "SRS_Analysis" srs_ui.py
# o executável fica em dist/SRS_Analysis.exe
```

---

## 3. Formatos de arquivo suportados

| Botão | Formato | Extensões | Módulo parser |
|-------|---------|-----------|---------------|
| **Abrir arquivo...** | Sun Nuclear SNC Water Tank | `.snctxt` | `parse_snctxt.py` |
| **Importar Monaco...** | Monaco TPS — grade de dose 2D | `.ALL` | `parse_monaco.py` |
| **Importar MCC...** | PTW Mephisto / BeamScan ASCII | `.mcc` | `parse_mcc.py` |
| **Importar DTA...** | IBA OmniPro Accept / texto plano | `.dta` `.asc` | `parse_dta.py` |
| **Importar RFB...** | OmniPro V6 binário (IBA Dosimetry) | `.rfb` | `parse_rfb.py` |
| **Importar outros...** | Múltiplos formatos (auto-detecção) | `.asc` `.txt` `.csv` `.sncxml` | `parse_multiformats.py` |

### Formatos disponíveis em "Importar outros..."

| Formato | Extensão | Notas |
|---------|----------|-------|
| Eclipse w2CAD | `.asc` | Exportação direta do TPS Eclipse |
| Eclipse Line Profile Export | `.txt` | Mesmo parser do w2CAD |
| OmniPro RFA300 ASCII BDS | `.txt` `.asc` | Detectado por `:MSR` / `%SCN` |
| RayStation Physics Export | `.csv` | Detectado por `#Exported: RayStation` |
| SNC Water Tank XML | `.sncxml` | Formato XML nativo do SNC |
| Standard Imaging DoseView 1D | `.csv` | Detectado por `DoseView 1D Software Export` |

> **Formato não suportado no momento:**
> SNC IC Profiler (`.prm`) — formato não documentado publicamente.

---

## 4. Interface gráfica — visão geral

```
┌──────────────────────────────────────────────────────────────────────────────────┐
│ [Abrir] [Monaco] [MCC] [DTA] [RFB] [Outros] [Reset tudo]  arquivo.snctxt        │
├──────────────────────┬───────────────────────────────────────────────────────────┤
│ ┌────────────────┐   │                                                        │
│ │  Scans  Editar │   │                                                        │
│ │  Gama          │   │          Área de gráfico (matplotlib)                  │
│ │                │   │                                                        │
│ │                │   │                                                        │
│ └────────────────┘   │                                                        │
├──────────────────────┴───────────────────────────────────────────────────────┤
│ Barra de status                                                               │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 4.1 Barra de ferramentas

| Botão | Ação |
|-------|------|
| **Abrir arquivo...** | Carrega um `.snctxt` como conjunto principal (substitui scans existentes) |
| **Importar Monaco...** | Adiciona perfis extraídos de arquivos Monaco `.ALL` |
| **Importar MCC...** | Adiciona scans de arquivos PTW Mephisto `.mcc` |
| **Importar DTA...** | Adiciona scans de arquivos DTA/ASC (IBA ou texto plano) |
| **Importar RFB...** | Adiciona scans de arquivos OmniPro V6 binários `.rfb` |
| **Importar outros...** | Adiciona scans de qualquer formato suportado com auto-detecção |
| **Reset tudo** | Remove todas as edições (suavização, renorm., offset) de todos os scans |

> Todos os botões de **importação adicionam** scans à lista atual (não substituem).
> Use "Reset tudo" + recarregamento manual para recomeçar do zero.

---

### 4.2 Aba Scans

**Filtros**

Permite filtrar a lista por:
- **Feixe** — Photon, Electron, etc.
- **Tipo** — Depth Scan, Crossline, Inline, Diagonal
- **Detector** — nome do detector lido do arquivo (ou editado)

O botão **Limpar filtros** restaura a lista completa.

**Lista de scans**

Exibe todos os scans que passam pelos filtros.  Cada linha mostra:
```
[idx]* BeamType Energy  FieldXxFieldY cm   ScanType       Detector
```
O `*` indica que o scan possui alguma edição pendente.

Suporta seleção múltipla com `Ctrl+clique` ou `Shift+clique`.

**Opções de plot**

- **Normalizar** — normaliza cada curva a 100% no dmax (PDP) ou no CAX (perfil)
- **Colorir por Detector / Campo** — controla a cor das curvas no gráfico

Botão **">> Plotar selecionados"** gera o gráfico agrupado por tipo de curva.

---

### 4.3 Aba Editar

O painel de edição possui scroll vertical.  Todas as edições são **não
destrutivas** — o dado original fica intacto até que o botão
**"✔ Aceitar"** seja pressionado.

#### Seletor de scan

Lista todos os scans (sem filtro).  Múltiplos scans podem ser selecionados;
edições são aplicadas a **todos** os selecionados simultaneamente.

#### Renomear / Corrigir informações

Campos editáveis:
- **Nome** — rótulo exibido nos gráficos
- **Feixe** — modalidade (Photon, Electron…)
- **Energia** — valor numérico (ex.: `6`, `15`, `9`)
- **Detector** — nome do detector
- **Campo (cm)** — largura X e Y do campo colimado

Botão **"Aplicar informações"** confirma as alterações de metadados.

#### Suavização

| Parâmetro | Descrição |
|-----------|-----------|
| Método | `savgol` (Savitzky-Golay), `moving_avg` (média móvel), `gaussian` (filtro gaussiano) |
| Janela (pts) | Número de pontos da janela; deve ser ímpar para savgol |
| Grau (savgol) | Ordem do polinômio (padrão 3) |
| Sigma (gauss) | Desvio padrão do kernel gaussiano (em pontos) |
| Região (cm) | Aplica a suavização somente entre `start` e `end` (vazio = tudo) |

#### Centralizar perfil / PDP

| Método | Ação |
|--------|------|
| `fwhm` | Desloca a posição de modo que o ponto médio do FWHM fique em 0 |
| `cax` | Desloca de modo que o pico (CAX) fique em 0 |
| `dmax` | Desloca de modo que o dmax fique em Z = 0 (para PDPs) |

Exibe o offset calculado para cada scan.

#### Deslocar posição (eixo X)

Aplica um deslocamento manual em cm ao eixo de posição. Útil para alinhar
curvas medidas com eixos de referência.

#### Renormalizar

| Método | Ponto de referência |
|--------|---------------------|
| `dmax` | Pico da curva (dose máxima) |
| `cax` | Ponto central (posição 0) |
| `point` | Posição específica definida em "Posição (cm)" |
| `region` | Média da dose dentro da região [start, end] |

#### Métricas

Caixa de texto com as métricas calculadas para o(s) scan(s) selecionado(s).
Ver [Seção 7](#7-métricas-calculadas) para a lista completa.

#### Aceitar / Reset

| Botão | Efeito |
|-------|--------|
| **✔ Aceitar** | Consolida a dose editada como novo "original"; reset futuro parte deste estado |
| **✖ Rejeitar / Reset selecionados** | Descarta todas as edições dos scans selecionados |
| **✖ Reset TODOS** | Descarta edições de todos os scans na sessão |

---

### 4.4 Aba Gama

Calcula o **índice gama 1D** entre dois scans.

| Campo | Descrição |
|-------|-----------|
| **Referência** | Scan que serve como referência (medição ou TPS) |
| **Avaliação** | Scan a ser comparado contra a referência |
| **DD (%)** | Critério de diferença de dose (ex.: 3%) |
| **DTA (cm)** | Critério de distância de concordância (ex.: 0.3 = 3 mm) |
| **Threshold (%)** | Dose mínima (% do máximo) para incluir ponto no cálculo |
| **Normaliz.** | `max` = normaliza pelo máximo global; `local` = ponto a ponto |

O botão **">> Calcular Gama"** gera um gráfico com 3 painéis:
1. Curvas sobrepostas (referência azul, avaliação vermelho)
2. Diferença ponto-a-ponto (barras coloridas: azul ≤ DD, vermelho > DD)
3. Índice γ por posição (verde ≤ 1, vermelho > 1)

Exibe: **pass rate (%)**, **γ_max**, **γ_mean** e contagem de pontos.

---

### 4.5 Área de gráfico

Gráfico interativo Matplotlib com barra de navegação padrão:
- Pan / Zoom com mouse
- Salvar figura (PNG, PDF, SVG…)
- Reset de zoom

No modo **"Plotar selecionados"**, as curvas são agrupadas por tipo
(Depth Scan, Crossline, Inline) em subplots separados.

No modo **"Preview de edição"** (acionado ao selecionar na aba Editar), cada
subplot mostra original (---cinza) vs editado (—colorido) mais um painel de
diferença abaixo.

---

## 5. Módulos — referência técnica

### `parse_snctxt.py`

Parser principal para o formato Sun Nuclear SNC Water Tank (`.snctxt`).

**Funções públicas:**

| Função | Descrição |
|--------|-----------|
| `parse_snctxt(filepath)` | Lê o arquivo e retorna `list[Scan]` |
| `gamma_1d(pos_r, dose_r, pos_e, dose_e, dd, dta, norm, threshold)` | Calcula índice γ 1D; retorna `(positions, gamma_array)` |
| `gamma_pass_rate(gamma_array)` | Retorna taxa de aprovação em % (γ ≤ 1) |
| `compute_profile_metrics(pos, dose)` | Calcula métricas de perfil lateral |
| `compute_pdp_metrics(pos, dose)` | Calcula métricas de PDP/PDD |
| `format_metrics(scan)` | Retorna string formatada com todas as métricas do scan |
| `apply_smooth(scan, method, window, polyorder, sigma, region)` | Suavização in-place |
| `apply_centering(scan, method)` | Centralização in-place; retorna offset em cm |
| `apply_renormalize(scan, method, value, ref_pos, region)` | Renormalização in-place |

**Classe `Scan`** — ver [Seção 6](#6-estrutura-de-dados--classe-scan).

---

### `parse_monaco.py`

Parser para exportação de dose 2D do Monaco (Elekta TPS) no formato `.ALL`.

**Funções públicas:**

| Função | Descrição |
|--------|-----------|
| `parse_monaco_all(filepath)` | Lê um arquivo `.ALL` e retorna `MonacoFile` |
| `load_monaco_files(filepaths, index_start)` | Carrega vários arquivos e extrai `list[Scan]` |

**Classe `MonacoFile`:**

| Campo | Tipo | Descrição |
|-------|------|-----------|
| `filepath` | str | Caminho do arquivo |
| `plan_name` | str | Nome do plano (ex.: "CCC1", "MontC05") |
| `algorithm` | str | Algoritmo: "CCC" ou "MonteCarlo" |
| `plane_type` | str | "Coronal" ou "Transversal" |
| `position_mm` | float | Profundidade (Coronal) ou posição lateral (Transversal) em mm |
| `dose_matrix` | ndarray | Grade de dose (ny × nx) em cGy |

**Convenção de coordenadas:**
- **Coronal** (PlaneDesc = `C: X.XX cm`): extrai Crossline em Y=0 e Inline em X=0
- **Transversal** (PlaneDesc = `T: X.XX cm`): extrai PDP em X=0

---

### `parse_mcc.py`

Parser para o formato PTW Mephisto ASCII (`.mcc`).

**Funções públicas:**

| Função | Descrição |
|--------|-----------|
| `parse_mcc(filepath, index_start)` | Lê um arquivo `.mcc` e retorna `list[Scan]` |
| `load_mcc_files(filepaths, index_start)` | Carrega vários arquivos; retorna `list[Scan]` |

**Estrutura do arquivo `.mcc`:**

```
BEGIN_SCAN_DATA
  BEGIN_SCAN
    SCAN_CURVETYPE = PDD | INPLANE_PROFILE | CROSSPLANE_PROFILE | DIAGONAL_PROFILE
    MODALITY       = XRAY | ELECTRON
    ENERGY         = 6.0
    FIELD_INPLANE  = 100.0     (mm)
    FIELD_CROSSPLANE = 100.0   (mm)
    SSD            = 1000.0    (mm)
    DETECTOR       = SEMIFLEX
    BEGIN_DATA
      pos_mm  dose
      ...
    END_DATA
  END_SCAN
END_SCAN_DATA
```

**Mapeamento `SCAN_CURVETYPE` → `scan_type`:**

| MCC | Scan.scan_type |
|-----|----------------|
| PDD | Depth Scan |
| INPLANE_PROFILE | Inline |
| CROSSPLANE_PROFILE | Crossline |
| DIAGONAL_PROFILE | Diagonal |

Posições convertidas de **mm → cm**.

---

### `parse_dta.py`

Parser para arquivos de perfil 1D `.dta` / `.asc`.  Detecta automaticamente
dois sub-formatos:

**Sub-formato IBA OmniPro Accept ASCII:**
- Blocos delimitados por `:MSR` / `:EOM`
- Cabeçalho com pares `chave = valor` (separados por tabulação)
- Dados após a linha `=`

**Sub-formato texto plano (fallback):**
- Comentários com `#`, `%`, `;`
- Dados em duas colunas (posição mm, dose)
- Suporta múltiplos blocos por arquivo separados por linhas em branco

**Funções públicas:**

| Função | Descrição |
|--------|-----------|
| `parse_dta(filepath, index_start)` | Lê um arquivo DTA/ASC; auto-detecta formato |
| `load_dta_files(filepaths, index_start)` | Carrega vários arquivos |

---

### `parse_rfb.py`

Parser para o formato binário OmniPro V6 RFB (`.rfb`) da IBA Dosimetry.
Portado do projeto RadPy (Stephen Terry et al., BSD-3-Clause) para Python 3 puro
— sem dependências externas além de `numpy`.

**Funções públicas:**

| Função | Descrição |
|--------|-----------|
| `parse_rfb(filepath, index_start)` | Lê um arquivo `.rfb` e retorna `list[Scan]` |
| `load_rfb_files(filepaths, index_start)` | Carrega vários arquivos; retorna `list[Scan]` |

**Estrutura binária do arquivo (little-endian / Windows nativo):**

```
PascalString          versão OmniPro (ex.: "OmniPro 6.4")
[main_header]
  Padding(13)
  PascalString          rad_device (linac)
  float64               energy (MeV/MV)
  uint8                 particle  (0=Photon, 1=Electron, 2=Proton…)
  int16                 wedge_type
  uint8                 wedge_angle
  uint16                gantry_angle, collimator_angle
  float64               SSD, SAD  (mm)
  PascalString          applicator
  int8                  medium (0=Air, 1=Water, 2=Film)
  PascalStrings         institution, address, telephone, email
  float64 ×4            jaw positions (inplane±, crossplane±)  (mm)
  uint8                 gantry_scale
  [bytes até != 0x00]   num_scans_with_this_header
  Padding(5)
  [scan_type_field]     string com tipo da varredura
N × [measurement_data + 2 bytes padding]
[3 bytes discriminador]
  → 0:   [additional_header]  — novo grupo (parâmetros diferentes)
  → N≠0: [new_scan_type_header] — mesmo header, outro tipo de curva
[...repete até EOF]
12 × 0x00              marcador de fim de arquivo
```

**measurement_data (bloco por varredura):**

```
int32           measured_date (Unix timestamp)
uint8           data_type (1=rel_opt_density, 2=rel_dose, 3=rel_ionization…)
float64 ×2      chamber_radius, calibration_factor
float64 ×2      unknown1, unknown2
PascalString    calibration_date
PascalString    detector
...             (campos de eletrômetro, servo, correntes, HV)
float64 ×3      isocenter (crossplane, inplane, depth)  (mm)
float64 ×3      scan_start (crossplane, inplane, depth) (mm)
float64 ×3      scan_end   (crossplane, inplane, depth) (mm)
int16           N  (número de pontos)
N × float64 ×2  (abscissa_mm, ordinate)
```

**Determinação do tipo de curva** — comparação da variação nas 3 coordenadas:
- `|Δz| > |Δx|, |Δy|` → Depth Scan
- `|Δx| > |Δy|` → Crossline
- `|Δy| > |Δx|` → Inline

**Conversão de unidades:** posições mm → cm (`÷ 10`), jaws mm → cm.

**Créditos:** RadPy © 2011 Stephen Terry et al. — BSD-3-Clause.

---

### `parse_multiformats.py`

Parser unificado com auto-detecção de formato.  Cobre os formatos do
`ParseProfile.m` (water_tank, mwgeurts).

**Formatos e critérios de detecção:**

| Formato | Critério de detecção |
|---------|---------------------|
| Eclipse w2CAD | Linhas com `%VERSION`, `%BMTY`, `%FLSZ` ou `$STOM` |
| OmniPro RFA300 ASCII BDS | Linha `:MSR` ou `:SYS BDS` ou `%SCN` |
| RayStation Physics Export | `#Exported:` contendo "RayStation" ou campos `CurveType` + `RadiationType` |
| SNC Water Tank XML | Extensão `.sncxml` |
| Standard Imaging DV1D | Linha contendo "DoseView 1D" ou "Charge Table" |

**Funções públicas:**

| Função | Descrição |
|--------|-----------|
| `parse_file(filepath, index_start)` | Lê um arquivo com auto-detecção; lança `ValueError` se não reconhecido |
| `load_auto(filepaths, index_start)` | Carrega lista de arquivos; retorna `list[Scan]` unificada |

**Sub-parsers internos:**

| Função | Formato |
|--------|---------|
| `_parse_w2cad(lines, index_start)` | Eclipse w2CAD / Line Profile |
| `_parse_ibatxt(lines, index_start)` | OmniPro RFA300 ASCII BDS |
| `_parse_raystation(lines, index_start)` | RayStation CSV |
| `_parse_sncxml(filepath, index_start)` | SNC Water Tank XML (usa `xml.etree.ElementTree`) |
| `_parse_dv1d(lines, index_start)` | Standard Imaging DoseView 1D |

---

### `srs_ui.py`

Aplicação principal — janela Tkinter `SRSApp(tk.Tk)`.

**Classe `SRSApp`:**

| Atributo | Tipo | Descrição |
|----------|------|-----------|
| `scans` | `list[Scan]` | Lista mestra de todos os scans carregados |
| `filepath` | str | Caminho do arquivo snctxt principal |
| `_scans_filtered` | `list[Scan]` | Subconjunto filtrado exibido na aba Scans |

**Métodos principais:**

| Método | Descrição |
|--------|-----------|
| `_build_ui()` | Constrói toda a interface (barra, painéis, abas) |
| `_open_file()` | Diálogo + carregamento de `.snctxt` |
| `_load_file(path)` | Parseia snctxt e atualiza todos os controles |
| `_import_monaco()` | Importa arquivos Monaco `.ALL` |
| `_import_mcc()` | Importa arquivos PTW Mephisto `.mcc` |
| `_import_dta()` | Importa arquivos DTA/ASC |
| `_import_rfb()` | Importa arquivos OmniPro V6 `.rfb` |
| `_import_multi()` | Importa múltiplos formatos com auto-detecção |
| `_populate_filters()` | Preenche os comboboxes de filtro com valores únicos |
| `_apply_filters()` | Retorna lista filtrada conforme seleções |
| `_scan_row(scan)` | Formata a string de exibição de um scan na listbox |
| `_refresh_scans_list()` | Reconstrói a listbox da aba Scans |
| `_refresh_edit_list()` | Reconstrói a listbox da aba Editar |
| `_refresh_gamma_combos()` | Atualiza os comboboxes da aba Gama |
| `_plot_selected()` | Plota os scans selecionados na aba Scans |
| `_plot_edit_preview()` | Plota original vs editado na aba Editar |
| `_annotate_metrics_on_ax(ax, ...)` | Adiciona caixa de métricas ao subplot |
| `_apply_meta()` | Aplica edições de metadados aos scans selecionados |
| `_apply_smooth()` | Aplica suavização |
| `_apply_center()` | Aplica centralização |
| `_apply_pos_delta()` | Aplica deslocamento manual de posição |
| `_apply_renorm()` | Aplica renormalização |
| `_accept_edit()` | Consolida edição como novo original |
| `_reset_selected()` | Reverte edições dos scans selecionados |
| `_reset_all()` | Reverte todas as edições |
| `_calc_gamma()` | Calcula índice γ 1D e plota resultado |
| `_plot_gamma(...)` | Renderiza o gráfico de gama (3 painéis) |
| `_update_metrics()` | Atualiza painel de métricas na aba Editar |

---

## 6. Estrutura de dados — classe `Scan`

Definida em `parse_snctxt.py`.  Representa um único perfil 1D.

**Campos de dados (imutáveis após leitura):**

| Campo | Tipo | Descrição |
|-------|------|-----------|
| `index` | int | Índice único (atribuído na ordem de carregamento) |
| `beam_type` | str | "Photon" ou "Electron" |
| `energy` | str | Energia nominal (ex.: "6", "15", "9") |
| `field_x` | str | Largura do campo em X (cm) |
| `field_y` | str | Largura do campo em Y (cm) |
| `scan_type` | str | "Depth Scan", "Crossline", "Inline", "Diagonal" |
| `ssd` | str | SSD em cm |
| `medium` | str | Meio de medição (ex.: "Water") |
| `date` | str | Data da medição (conforme o arquivo) |
| `detector` | str | Identificador do detector |
| `x` | ndarray | Coordenada X (cm) — array de zeros exceto para Crossline |
| `y` | ndarray | Coordenada Y (cm) — array de zeros exceto para Inline |
| `z` | ndarray | Coordenada Z (cm) — array de zeros exceto para Depth Scan |
| `dose` | ndarray | Dose bruta (unidades do arquivo) |

**Campos de edição (sobreposições não destrutivas):**

| Campo | Tipo | Descrição |
|-------|------|-----------|
| `name_override` | str | Rótulo personalizado para gráficos |
| `energy_override` | str | Energia editada |
| `detector_override` | str | Detector editado |
| `field_x_override` | str | Campo X editado |
| `field_y_override` | str | Campo Y editado |
| `_dose_edit` | ndarray\|None | Dose editada (None = usar original) |
| `_pos_offset` | float | Deslocamento de posição em cm |

**Propriedades computadas:**

| Propriedade | Retorna |
|-------------|---------|
| `position` | Array de posição ordenado + offset aplicado |
| `dose_sorted` | Dose original ordenada (sem edições) |
| `dose_display` | Dose com edições aplicadas |
| `display_name` | Nome formatado para exibição |
| `label` | Alias de `display_name` |

---

## 7. Métricas calculadas

### Perfis laterais (`compute_profile_metrics`)

| Métrica | Descrição |
|---------|-----------|
| `fwhm` | Largura a meia altura (cm) |
| `center` | Centro geométrico do FWHM (cm) |
| `l50`, `r50` | Posições 50% esquerda e direita (cm) |
| `l80`, `r80` | Posições 80% esquerda e direita (cm) |
| `l20`, `r20` | Posições 20% esquerda e direita (cm) |
| `penumbra_l` | Penumbra esquerda: distância 20%→80% (cm) |
| `penumbra_r` | Penumbra direita: distância 80%→20% (cm) |
| `flatness` | Planura dentro do campo: `(Dmax - Dmin)/(Dmax + Dmin) × 100` (%) |
| `symmetry` | Simetria: máxima assimetria de pares de pontos equidistantes do centro (%) |

### PDPs e PDDs (`compute_pdp_metrics`)

| Métrica | Descrição |
|---------|-----------|
| `zmax` | Profundidade do dose máxima (cm) |
| `r90`, `r80`, `r50`, `r10` | Profundidade em que a dose cai para 90/80/50/10% (cm) |
| `rp` | Range prático (elétrons): projeção da tangente em R50 até zero (cm) |

---

## 8. Índice Gama 1D

Implementado em `parse_snctxt.gamma_1d`.

**Fórmula:**

```
γ(xᵣ) = min_xₑ √[ (D_e(xₑ) - D_r(xᵣ))² / DD²  +  (xₑ - xᵣ)² / DTA² ]
```

**Parâmetros:**

| Parâmetro | Padrão | Descrição |
|-----------|--------|-----------|
| `dd` | 3.0 | Critério DD em % |
| `dta` | 0.3 | Critério DTA em cm (0.3 = 3 mm) |
| `norm` | "max" | "max" = normaliza pelo máximo global; "local" = ponto a ponto |
| `threshold` | 10.0 | Threshold em % do máximo; pontos abaixo são ignorados |

A função interpola a curva de avaliação nos pontos da referência antes do
cálculo.  A taxa de aprovação (`gamma_pass_rate`) conta pontos com γ ≤ 1.

---

## 8.1 Melhoria planejada do Índice Gama — inspirada no PyMedPhys

> **Status:** documentada, aguardando implementação.
>
> Referência: [PyMedPhys Gamma Analysis](https://github.com/pymedphys/pymedphys)
> — Wendling et al. 2007 (doi:10.1118/1.2721657)

### Motivação

A implementação atual realiza uma busca de força bruta: para cada ponto de
avaliação, o γ mínimo é calculado contra **todos** os pontos da grade de
referência interpolada.  Para perfis 1D com poucas centenas de pontos isso é
suficiente, mas a abordagem tem limitações que o PyMedPhys resolve com um
algoritmo mais robusto.

### Diferenças principais em relação à implementação atual

| Aspecto | Implementação atual | Melhoria planejada |
|---------|--------------------|--------------------|
| **Interpolação** | Grade de referência × 10 fixa | Grade de avaliação com passo = `DTA / interp_fraction` |
| **Threshold padrão** | 10 % do máximo | 20 % do máximo (padrão clínico AAPM TG-218) |
| **Normalização local** | Parcialmente implementada | Corretamente calculada sobre a dose de referência no ponto avaliado |
| **Terminação antecipada** | Não | Parâmetro `max_gamma`; interrompe busca ao ultrapassar o limite |
| **Algoritmo de busca** | Força bruta sobre toda a grade | Expansão por cascas (shells) a partir do ponto de referência |
| **Fração de interpolação** | Fixa em 10× | Configurável via `interp_fraction` (padrão 10) |

### Algoritmo de cascas (shell-based search)

Em vez de varrer todos os pontos, o algoritmo expande progressivamente a busca:

```
distância = 0
enquanto distância ≤ distância_máxima_de_teste:
    interpolar dose de avaliação nos pontos à distância atual
    calcular dose_diff mínima para cada ponto de referência ainda em busca
    atualizar γ = √[ (dose_diff/DD)² + (distância/DTA)² ]
    remover pontos onde γ < 1 (já aprovados)
    distância += passo_adaptativo
```

Para 1D, as "cascas" são simplesmente dois pontos simétricos (xᵣ ± d).
Isso elimina varreduras desnecessárias e permite terminação antecipada.

### Normalização global vs. local (detalhada)

**Global (atual padrão — "max"):**
```
DD_abs = (dd/100) × max(dose_referência)
γ = √[ ((D_eval - D_ref) / DD_abs)² + (Δx / DTA)² ]
```

**Local ("local") — melhoria:**
```
DD_abs = (dd/100) × D_ref(xᵣ)   ← dose de referência no ponto xᵣ
γ = √[ ((D_eval - D_ref) / DD_abs)² + (Δx / DTA)² ]
```

A normalização local é mais sensível em regiões de dose baixa e deve ser
usada com cuidado em penumbras — o padrão clínico recomendado é global.

### Novos parâmetros propostos

| Parâmetro | Padrão proposto | Descrição |
|-----------|----------------|-----------|
| `dd` | 3.0 | Critério DD em % (sem alteração) |
| `dta` | 0.3 | Critério DTA em cm (sem alteração) |
| `norm` | `"global"` | `"global"` ou `"local"` |
| `threshold` | 20.0 | Alterado de 10 % para 20 % (TG-218) |
| `interp_fraction` | 10 | Resolução = `dta / interp_fraction` |
| `max_gamma` | 2.0 | Interrompe busca acima deste valor |

### Impacto esperado

- **Pass rate:** o threshold de 20 % exclui mais pontos de baixa dose,
  podendo aumentar ligeiramente a taxa de aprovação em comparação com 10 %.
- **Velocidade:** irrelevante para 1D (< 1 ms), mas a arquitetura modular
  facilita extensão futura para 2D.
- **Normalização local:** pode reduzir a pass rate em penumbras onde há
  gradientes de dose altos — resultado esperado e clinicamente correto.

### Compatibilidade retroativa

A assinatura proposta é compatível com a atual:
```python
# atual
gamma_1d(pos_ref, dose_ref, pos_eval, dose_eval,
          dd=3.0, dta=0.3, norm="max", threshold=10.0)

# proposta
gamma_1d(pos_ref, dose_ref, pos_eval, dose_eval,
          dd=3.0, dta=0.3, norm="global", threshold=20.0,
          interp_fraction=10, max_gamma=2.0)
```

A interface gráfica (aba Gama) receberá um novo campo **"Fração interp."**
e o padrão do threshold será atualizado de 10 % para 20 %.

---

## 9. Fluxo de trabalho típico

### Comparação medida × TPS

1. **Abrir arquivo...** → selecionar `.snctxt` com dados medidos
2. **Importar Monaco...** → selecionar arquivos `.ALL` com dados do TPS
3. Aba **Scans** → selecionar curvas do mesmo tipo e plotar
4. Aba **Editar** → selecionar um scan medido → **Centralizar** → **Renormalizar**
5. Aba **Gama** → selecionar referência (TPS) e avaliação (medido) → calcular

### Importação de múltiplos formatos

1. **Importar outros...** → selecionar arquivos (`.asc`, `.csv`, `.sncxml`, etc.)
2. O programa detecta automaticamente o formato de cada arquivo
3. Erros por arquivo são listados em um aviso; arquivos válidos são importados normalmente
4. Usar a aba **Editar** para corrigir metadados se necessário

### Edição iterativa

```
Suavizar → visualizar preview → Aceitar → Renormalizar → Aceitar → ...
```

O botão **"✔ Aceitar"** consolida o estado atual como novo ponto de partida
para edições subsequentes sem perda de informação de iterações anteriores.

---

*Desenvolvido para uso clínico em física médica.*
*Baseado em: [mwgeurts/snc_extract](https://github.com/mwgeurts/snc_extract)
e [mwgeurts/water_tank](https://github.com/mwgeurts/water_tank) (referências MATLAB).*

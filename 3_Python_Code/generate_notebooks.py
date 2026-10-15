#!/usr/bin/env python3
# ==============================================================================
# Script generador de los 5 Notebooks interactivos de Jupyter por cada etapa
# ==============================================================================

import os
import nbformat as nbf

NOTEBOOKS_DIR = "/Users/efrainlav/Attenuation_Q/4_Notebooks"
os.makedirs(NOTEBOOKS_DIR, exist_ok=True)

def create_nb_01():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Etapa 1: Carga, Visualización y Control de Calidad de Datos Sísmicos\n\n"
            "Este notebook implementa el **Paso 1** del flujo de estimación de atenuación sísmica. "
            "Su objetivo es importar la información sísmica (formato SEG-Y o matrices de datos), "
            "estructurar el muestreo temporal y realizar controles de calidad geofísicos esenciales:\n\n"
            "* Análisis de traza individual en amplitud y tiempo.\n"
            "* Despliegue en *Wiggle* con relleno positivo para gathers y secciones sísmicas.\n"
            "* Espectro de amplitud y potencia en decibelios (dB) para identificar la banda pasante."
        ),
        nbf.v4.new_code_cell(
            "import os\n"
            "import sys\n"
            "import numpy as np\n"
            "import matplotlib.pyplot as plt\n\n"
            "# Añadir la carpeta de módulos en Python al PATH\n"
            "MODULE_DIR = os.path.abspath('../3_Python_Code')\n"
            "if MODULE_DIR not in sys.path:\n"
            "    sys.path.insert(0, MODULE_DIR)\n\n"
            "from io_seismic import generate_synthetic_seismic_trace, save_mat_seismic\n"
            "from config import SeismicSurveyConfig\n\n"
            "print('Módulos importados correctamente.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. Definición de Parámetros y Carga de la Señal Sísmica\n\n"
            "Definimos una configuración de levantamiento con una tasa de muestreo estándar de "
            "$\\Delta t = 2\\text{ ms}$ ($0.002\\text{ s}$) y duración de $2.0\\text{ s}$. "
            "Generamos una señal sísmica multicomponente con reflectores y atenuación."
        ),
        nbf.v4.new_code_cell(
            "survey_cfg = SeismicSurveyConfig(\n"
            "    survey_name='Seismic_Survey_Area',\n"
            "    line_name='Seismic_Line_01',\n"
            "    well_name='Well_01',\n"
            "    sample_rate_s=0.002,\n"
            "    t_max_s=2.0\n"
            ")\n\n"
            "# Generamos una traza sísmica con reflectores y factor Q = 60\n"
            "t_vec, main_trace = generate_synthetic_seismic_trace(\n"
            "    sample_rate_s=survey_cfg.sample_rate_s,\n"
            "    t_max_s=survey_cfg.t_max_s,\n"
            "    peak_frequencies=(25.0, 38.0, 50.0),\n"
            "    arrival_times_s=(0.35, 0.85, 1.35),\n"
            "    amplitudes=(1.0, 0.72, 0.48),\n"
            "    true_q=60.0,\n"
            "    noise_level=0.015,\n"
            "    random_seed=42\n"
            ")\n\n"
            "# Generamos un gather sísmico sintético de 15 trazas adyacentes\n"
            "n_traces = 15\n"
            "gather = np.zeros((len(t_vec), n_traces))\n"
            "for i in range(n_traces):\n"
            "    offset_factor = 1.0 - 0.015 * abs(i - n_traces // 2)\n"
            "    jitter = np.random.normal(0, 0.01, size=len(t_vec))\n"
            "    gather[:, i] = main_trace * offset_factor + jitter\n\n"
            "print(f'Traza principal: {len(main_trace)} muestras.')\n"
            "print(f'Gather sísmico: {gather.shape[1]} trazas x {gather.shape[0]} muestras.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. Visualización y Control de Calidad (QC)\n\n"
            "Generamos los gráficos clásicos del procesamiento sísmico:"
        ),
        nbf.v4.new_code_cell(
            "fig, axes = plt.subplots(1, 3, figsize=(16, 5))\n\n"
            "# Panel 1: Traza individual\n"
            "axes[0].plot(t_vec, main_trace, color='#1f77b4', linewidth=1.2)\n"
            "axes[0].set_title('A) Traza Sísmica Individual', fontweight='bold')\n"
            "axes[0].set_xlabel('Tiempo (s)')\n"
            "axes[0].set_ylabel('Amplitud Normalizada')\n"
            "axes[0].set_xlim(0, 2.0)\n"
            "axes[0].grid(True, linestyle=':', alpha=0.6)\n\n"
            "# Panel 2: Despliegue en Wiggle\n"
            "scale = 0.8\n"
            "for i in range(n_traces):\n"
            "    tr = gather[:, i]\n"
            "    tr_norm = tr / np.max(np.abs(gather)) * scale\n"
            "    axes[1].plot(i + tr_norm, t_vec, color='black', linewidth=0.7)\n"
            "    axes[1].fill_betweenx(t_vec, i, i + tr_norm, where=(tr_norm > 0), color='#2ca02c', alpha=0.6)\n"
            "axes[1].set_title('B) Sección Sísmica / Wiggle Display', fontweight='bold')\n"
            "axes[1].set_xlabel('Número de Traza / CDP')\n"
            "axes[1].set_ylabel('Tiempo Doble (s)')\n"
            "axes[1].set_ylim(2.0, 0)  # Invertir eje de tiempo\n"
            "axes[1].set_xlim(-1, n_traces)\n"
            "axes[1].grid(True, linestyle=':', alpha=0.5)\n\n"
            "# Panel 3: Espectro de amplitud\n"
            "n_fft = 2048\n"
            "fft_vals = np.abs(np.fft.rfft(main_trace, n=n_fft))\n"
            "freqs = np.fft.rfftfreq(n_fft, d=survey_cfg.sample_rate_s)\n"
            "fft_db = 20 * np.log10(fft_vals / np.max(fft_vals) + 1e-6)\n"
            "axes[2].plot(freqs, fft_db, color='#d62728', linewidth=1.5)\n"
            "axes[2].set_title('C) Espectro de Amplitud (dB)', fontweight='bold')\n"
            "axes[2].set_xlabel('Frecuencia (Hz)')\n"
            "axes[2].set_ylabel('Amplitud (dB)')\n"
            "axes[2].set_xlim(0, 100)\n"
            "axes[2].set_ylim(-40, 5)\n"
            "axes[2].grid(True, linestyle=':', alpha=0.6)\n\n"
            "plt.tight_layout()\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 3. Exportación de Datos para la Etapa de Descomposición\n\n"
            "Guardamos la traza procesada y sus metadatos en un archivo `.mat` y `.npz` "
            "para que puedan ser cargados por los siguientes notebooks."
        ),
        nbf.v4.new_code_cell(
            "os.makedirs('data_processed', exist_ok=True)\n"
            "save_mat_seismic(\n"
            "    'data_processed/traza_preprocesada.mat',\n"
            "    {'trace': main_trace, 'time': t_vec, 'dt': survey_cfg.sample_rate_s}\n"
            ")\n"
            "np.savez('data_processed/traza_preprocesada.npz', trace=main_trace, time=t_vec, dt=survey_cfg.sample_rate_s)\n"
            "print('Datos exportados exitosamente a data_processed/traza_preprocesada.mat y .npz')"
        )
    ]
    return nb

def create_nb_02():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Etapa 2: Generación del Diccionario de Ondículas de Morlet\n\n"
            "Este notebook implementa el **Paso 2** del flujo de trabajo: la síntesis de un diccionario redundante "
            "tridimensional de ondículas de Morlet en el espacio $[Frecuencia \\times Fase \\times Tiempo]$.\n\n"
            "### Formulación Matemática\n"
            "Cada ondícula compleja de Morlet $m(t; f, \\phi)$ se evalúa como:\n\n"
            "$$m(t; f, \\phi) = \\exp\\left( -t^2 \\cdot f^2 \\cdot \\ln(2/k) \\right) \\cdot \\exp(i 2\\pi f t) \\cdot \\exp(i \\phi)$$\n\n"
            "donde:\n"
            "* $f$: Frecuencia central (Hz).\n"
            "* $\\phi$: Ángulo de fase (grados o radianes).\n"
            "* $t$: Coordenada temporal relativa a la envolvente (s).\n"
            "* $k$: Parámetro de ciclos de la ondícula (típicamente $k = 0.2$)."
        ),
        nbf.v4.new_code_cell(
            "import os\n"
            "import sys\n"
            "import numpy as np\n"
            "import matplotlib.pyplot as plt\n\n"
            "MODULE_DIR = os.path.abspath('../3_Python_Code')\n"
            "if MODULE_DIR not in sys.path:\n"
            "    sys.path.insert(0, MODULE_DIR)\n\n"
            "from config import MorletDictionaryConfig\n"
            "from morlet_dictionary import compute_morlet_atom, build_morlet_dictionary\n\n"
            "print('Módulo de diccionario cargado.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. Configuración de Parámetros del Diccionario\n\n"
            "Definimos el rango de frecuencias exploradas (de $5$ a $120\\text{ Hz}$), "
            "el paso angular de fase ($-180^\\circ$ a $+180^\\circ$) y la ventana temporal."
        ),
        nbf.v4.new_code_cell(
            "dict_cfg = MorletDictionaryConfig(\n"
            "    f_min=5.0,\n"
            "    f_max=100.0,\n"
            "    df=2.0,\n"
            "    phase_min_deg=-180.0,\n"
            "    phase_max_deg=180.0,\n"
            "    dphase_deg=10.0,\n"
            "    wavelet_half_duration_s=0.25,\n"
            "    k_param=0.2\n"
            ")\n"
            "sample_rate_s = 0.002\n\n"
            "morlet_dict, f_vec, phi_vec, t_atom = build_morlet_dictionary(dict_cfg, sample_rate_s)\n"
            "print(f'Dimensiones del diccionario tridimensional: {morlet_dict.shape}')\n"
            "print(f'Frecuencias evaluadas: {len(f_vec)} valores ({f_vec[0]} a {f_vec[-1]} Hz).')\n"
            "print(f'Fases evaluadas: {len(phi_vec)} valores ({phi_vec[0]}° a {phi_vec[-1]}°).')\n"
            "print(f'Muestras temporales: {len(t_atom)} muestras ({t_atom[0]*1000:.0f} a {t_atom[-1]*1000:.0f} ms).')"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. Inspección Visual de Cortes del Diccionario\n\n"
            "Visualizamos la dependencia con la frecuencia, el comportamiento ante rotaciones de fase "
            "y la huella energética espectral en el plano tiempo-frecuencia."
        ),
        nbf.v4.new_code_cell(
            "fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))\n\n"
            "# Corte A: Frecuencias crecientes a fase 0°\n"
            "freqs_show = [15, 30, 50, 75]\n"
            "colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']\n"
            "for f, c in zip(freqs_show, colors):\n"
            "    atom = compute_morlet_atom(t_atom, frequency_hz=f, phase_deg=0.0, k_param=dict_cfg.k_param)\n"
            "    axes[0].plot(t_atom * 1000, np.real(atom), label=f'{f} Hz', color=c, linewidth=1.5)\n"
            "axes[0].set_title('A) Ondículas vs Frecuencia (Fase = 0°)', fontweight='bold')\n"
            "axes[0].set_xlabel('Tiempo (ms)')\n"
            "axes[0].set_ylabel('Amplitud Real')\n"
            "axes[0].legend(loc='upper right')\n"
            "axes[0].grid(True, linestyle=':', alpha=0.6)\n\n"
            "# Corte B: Rotación de fase a 35 Hz\n"
            "phases_show = [-90, -45, 0, 45, 90]\n"
            "phase_colors = ['#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']\n"
            "for ph, c in zip(phases_show, phase_colors):\n"
            "    atom = compute_morlet_atom(t_atom, frequency_hz=35.0, phase_deg=ph, k_param=dict_cfg.k_param)\n"
            "    axes[1].plot(t_atom * 1000, np.real(atom), label=f'$\\phi = {ph}^\\circ$', color=c, linewidth=1.4)\n"
            "axes[1].set_title('B) Rotación de Fase (f = 35 Hz)', fontweight='bold')\n"
            "axes[1].set_xlabel('Tiempo (ms)')\n"
            "axes[1].legend(loc='upper right')\n"
            "axes[1].grid(True, linestyle=':', alpha=0.6)\n\n"
            "# Corte C: Huella espectral del átomo de Morlet en (t, f)\n"
            "f_axis = np.linspace(5, 100, 200)\n"
            "t_axis = np.linspace(-150, 150, 200)\n"
            "T, F = np.meshgrid(t_axis / 1000, f_axis)\n"
            "f0 = 40.0\n"
            "sigma_t = np.sqrt(np.log(2.0 / dict_cfg.k_param)) / (2.0 * np.pi * f0)\n"
            "sigma_f = 1.0 / (2.0 * np.pi * sigma_t)\n"
            "envelope_2d = np.exp(-0.5 * (T / sigma_t)**2 - 0.5 * ((F - f0) / sigma_f)**2)\n"
            "im = axes[2].contourf(t_axis, f_axis, envelope_2d, levels=25, cmap='magma')\n"
            "axes[2].set_title('C) Huella Tiempo-Frecuencia del Átomo', fontweight='bold')\n"
            "axes[2].set_xlabel('Tiempo (ms)')\n"
            "axes[2].set_ylabel('Frecuencia (Hz)')\n"
            "plt.colorbar(im, ax=axes[2], label='Energía Relativa')\n\n"
            "plt.tight_layout()\n"
            "plt.show()"
        )
    ]
    return nb

def create_nb_03():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Etapa 3: Descomposición Espectral Matching Pursuit (MPD)\n\n"
            "Este notebook implementa el **Paso 3**: la descomposición voraz (*greedy*) de la traza sísmica "
            "en átomos elementales de Morlet y el cálculo de la **Distribución Tiempo-Frecuencia (TFD)**.\n\n"
            "### Algoritmo MPD\n"
            "1. Se calcula la energía inicial de la señal residual: $E_0 = \\|R_0\\|^2$.\n"
            "2. En cada iteración $n$, mediante la señal analítica de Hilbert se estiman los atributos instantáneos:\n"
            "   * Tiempo central de llegada $t_0 = \\arg\\max |z(t)|$.\n"
            "   * Amplitud envolvente $A_n = |z(t_0)|$.\n"
            "   * Fase instantánea $\\phi_n = \\arg(z(t_0))$.\n"
            "   * Frecuencia instantánea $f_n = \\frac{1}{2\\pi} \\frac{d\\theta}{dt}(t_0)$.\n"
            "3. Se selecciona y sustrae el átomo de Morlet óptimo: $R_n(t) = R_{n-1}(t) - A_n m_n(t)$.\n"
            "4. Se acumula en la traza modelada $\\hat{s}(t) = \\sum A_n m_n(t)$ hasta alcanzar el umbral de energía (`cutoff`)."
        ),
        nbf.v4.new_code_cell(
            "import os\n"
            "import sys\n"
            "import numpy as np\n"
            "import matplotlib.pyplot as plt\n"
            "from matplotlib.gridspec import GridSpec\n\n"
            "MODULE_DIR = os.path.abspath('../3_Python_Code')\n"
            "if MODULE_DIR not in sys.path:\n"
            "    sys.path.insert(0, MODULE_DIR)\n\n"
            "from config import MPDConfig, MorletDictionaryConfig\n"
            "from matching_pursuit import decompose_trace_mpd, compute_time_frequency_distribution\n"
            "from io_seismic import load_mat_seismic\n\n"
            "# Cargamos los datos preprocesados de la Etapa 1\n"
            "mat_data, _ = load_mat_seismic('data_processed/traza_preprocesada.mat')\n"
            "trace = mat_data.squeeze()\n"
            "sample_rate_s = 0.002\n"
            "t_vec = np.arange(len(trace)) * sample_rate_s\n"
            "print(f'Traza cargada con {len(trace)} muestras.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. Ejecución de la Descomposición Matching Pursuit"
        ),
        nbf.v4.new_code_cell(
            "mpd_cfg = MPDConfig(energy_cutoff=0.003, max_iterations=150)\n"
            "dict_cfg = MorletDictionaryConfig(f_min=5.0, f_max=100.0, df=1.0, wavelet_half_duration_s=0.25)\n\n"
            "mpd_results = decompose_trace_mpd(\n"
            "    trace,\n"
            "    sample_rate_s=sample_rate_s,\n"
            "    mpd_config=mpd_cfg,\n"
            "    dict_config=dict_cfg\n"
            ")\n\n"
            "n_atoms = len(mpd_results['atoms'])\n"
            "final_energy = mpd_results['energy_history'][-1] * 100\n"
            "print(f'Iteraciones completadas: {mpd_results[\"iterations\"]}')\n"
            "print(f'Átomos de Morlet extraídos: {n_atoms}')\n"
            "print(f'Energía residual final: {final_energy:.3f}% (explicó {100-final_energy:.3f}% de la señal).')"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. Reconstrucción y Señal Residual\n\n"
            "Comparamos la traza original frente a la reconstruida por la suma de átomos y evaluamos el decaimiento de energía."
        ),
        nbf.v4.new_code_cell(
            "fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=False)\n\n"
            "# A: Traza original vs modelada\n"
            "axes[0].plot(t_vec, trace, 'k-', label='Traza Original', linewidth=1.3, alpha=0.8)\n"
            "axes[0].plot(t_vec, mpd_results['modeled_trace'], 'r--', label='Traza Reconstruida (MPD)', linewidth=1.2)\n"
            "axes[0].set_title('A) Comparación: Traza Original vs Reconstrucción MPD', fontweight='bold')\n"
            "axes[0].set_ylabel('Amplitud')\n"
            "axes[0].legend(loc='upper right')\n"
            "axes[0].grid(True, linestyle=':', alpha=0.6)\n\n"
            "# B: Traza residual\n"
            "axes[1].plot(t_vec, mpd_results['residual_trace'], color='#7f7f7f', linewidth=1.0)\n"
            "axes[1].set_title('B) Señal Residual $R_N(t)$ (Ruido no Explicado)', fontweight='bold')\n"
            "axes[1].set_ylabel('Amplitud')\n"
            "axes[1].set_xlabel('Tiempo (s)')\n"
            "axes[1].grid(True, linestyle=':', alpha=0.6)\n\n"
            "# C: Convergencia de energía\n"
            "iters = np.arange(len(mpd_results['energy_history']))\n"
            "axes[2].plot(iters, np.array(mpd_results['energy_history']) * 100, color='#2ca02c', linewidth=2.0)\n"
            "axes[2].axhline(mpd_cfg.energy_cutoff * 100, color='red', linestyle=':', label='Umbral de corte')\n"
            "axes[2].set_title('C) Convergencia de Energía Residual (%)', fontweight='bold')\n"
            "axes[2].set_xlabel('Número de Iteración / Átomo')\n"
            "axes[2].set_ylabel('Energía (%)')\n"
            "axes[2].set_yscale('log')\n"
            "axes[2].legend()\n"
            "axes[2].grid(True, linestyle=':', which='both', alpha=0.6)\n\n"
            "plt.tight_layout()\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 3. Construcción del Espectrograma Tiempo-Frecuencia (TFD)\n\n"
            "Sintetizamos el mapa continuo de densidad de energía tiempo-frecuencia proyectando los átomos de Morlet."
        ),
        nbf.v4.new_code_cell(
            "f_eval = np.arange(2.0, 90.0, 1.0)\n"
            "tfdist = compute_time_frequency_distribution(\n"
            "    mpd_results['atoms'],\n"
            "    time_vector=t_vec,\n"
            "    freq_vector=f_eval,\n"
            "    k_param=0.2\n"
            ")\n\n"
            "plt.figure(figsize=(12, 6))\n"
            "im = plt.imshow(\n"
            "    tfdist,\n"
            "    extent=[f_eval[0], f_eval[-1], t_vec[-1], t_vec[0]],\n"
            "    aspect='auto',\n"
            "    cmap='turbo'\n"
            ")\n"
            "plt.title('Distribución Tiempo-Frecuencia de Alta Resolución (MPD)', fontsize=13, fontweight='bold')\n"
            "plt.xlabel('Frecuencia (Hz)', fontsize=11)\n"
            "plt.ylabel('Tiempo (s)', fontsize=11)\n"
            "plt.colorbar(im, label='Densidad de Energía')\n"
            "plt.tight_layout()\n"
            "plt.show()\n\n"
            "# Guardar distribución tiempo-frecuencia para etapas 4 y 5\n"
            "np.savez('data_processed/tfdist_results.npz', tfdist=tfdist, t_vec=t_vec, f_eval=f_eval)\n"
            "print('Distribución TF exportada a data_processed/tfdist_results.npz')"
        )
    ]
    return nb

def create_nb_04():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Etapa 4: Modelado Teórico de Atenuación y Catálogo de Ganancias (Wang)\n\n"
            "Este notebook implementa el **Paso 4**: el cálculo analítico de las curvas de atenuación anelástica "
            "y la generación del catálogo teórico de compensación de amplitudes estabilizadas de Wang.\n\n"
            "### Formulación Física y Matemática\n"
            "El decaimiento teórico de la amplitud con la distancia y el tiempo está regido por:\n\n"
            "$$\\beta(X, Q) = \\exp\\left( -\\frac{\\pi X}{Q} \\right) = \\exp\\left( -\\frac{\\omega t}{2Q} \\right)$$\n\n"
            "donde $X = f \\cdot t$ es la variable adimensional del producto tiempo por frecuencia.\n\n"
            "Para compensar la amplitud sin amplificar de forma descontrolada el ruido en altas frecuencias, "
            "se define el parámetro de regularización de Wang con límite de ganancia $G_{\\lim} = 20\\text{ dB}$:\n\n"
            "$$\\sigma = \\exp\\left( -(0.23 \\cdot G_{\\lim} + 1.63) \\right)$$\n\n"
            "$$\\text{Ganancia}(X, Q) = \\frac{\\beta(X, Q) + \\sigma}{\\beta^2(X, Q) + \\sigma}$$"
        ),
        nbf.v4.new_code_cell(
            "import os\n"
            "import sys\n"
            "import numpy as np\n"
            "import matplotlib.pyplot as plt\n\n"
            "MODULE_DIR = os.path.abspath('../3_Python_Code')\n"
            "if MODULE_DIR not in sys.path:\n"
            "    sys.path.insert(0, MODULE_DIR)\n\n"
            "from config import QInversionConfig\n"
            "from q_theoretical import compute_theoretical_attenuation\n\n"
            "print('Módulo teórico cargado.')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. Generación de la Base Teórica en función de Q\n\n"
            "Exploramos el catálogo continuo para factores de calidad desde $Q = 10$ hasta $Q = 140$ "
            "y valores del producto tiempo-frecuencia $X$ de $1$ a $250$."
        ),
        nbf.v4.new_code_cell(
            "q_cfg = QInversionConfig(\n"
            "    q_min=10.0,\n"
            "    q_max=140.0,\n"
            "    dq=0.5,\n"
            "    x_min=1.0,\n"
            "    x_max=220.0,\n"
            "    dx=1.0,\n"
            "    gain_limit_db=20.0\n"
            ")\n\n"
            "beta, gan_teorica, x_vec, q_vec = compute_theoretical_attenuation(q_cfg)\n"
            "print(f'Matriz Beta(X, Q): {beta.shape}')\n"
            "print(f'Matriz de Ganancia Teórica: {gan_teorica.shape}')\n"
            "print(f'Factor de estabilización desvi (Glim=20 dB): {q_cfg.desvi:.6f}')"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. Inspección del Catálogo Teórico\n\n"
            "Comparamos cómo varía la atenuación espectral $\\beta(X, Q)$ y la ganancia estabilizada resultante "
            "para medios de alta disipación ($Q=15$) frente a medios poco atenuantes ($Q=100$)."
        ),
        nbf.v4.new_code_cell(
            "fig, axes = plt.subplots(1, 2, figsize=(15, 5))\n\n"
            "q_targets = [15, 30, 50, 75, 100]\n"
            "for q_t in q_targets:\n"
            "    idx_q = int(np.argmin(np.abs(q_vec - q_t)))\n"
            "    axes[0].plot(x_vec, beta[:, idx_q], label=f'Q = {q_vec[idx_q]:.0f}', linewidth=1.8)\n"
            "axes[0].set_title('A) Decaimiento Espectral Teórico $\\\\beta(X, Q)$', fontweight='bold')\n"
            "axes[0].set_xlabel('X (frecuencia × tiempo)')\n"
            "axes[0].set_ylabel('Amplitud $\\\\beta$')\n"
            "axes[0].legend(loc='upper right')\n"
            "axes[0].grid(True, linestyle=':', alpha=0.6)\n\n"
            "for q_t in q_targets:\n"
            "    idx_q = int(np.argmin(np.abs(q_vec - q_t)))\n"
            "    axes[1].plot(x_vec, gan_teorica[:, idx_q], label=f'Q = {q_vec[idx_q]:.0f}', linewidth=1.8)\n"
            "axes[1].set_title('B) Ganancia Estabilizada de Wang ($G_{lim}=20$ dB)', fontweight='bold')\n"
            "axes[1].set_xlabel('X (frecuencia × tiempo)')\n"
            "axes[1].set_ylabel('Ganancia Teórica')\n"
            "axes[1].legend(loc='upper left')\n"
            "axes[1].grid(True, linestyle=':', alpha=0.6)\n\n"
            "plt.tight_layout()\n"
            "plt.show()"
        )
    ]
    return nb

def create_nb_05():
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell(
            "# Etapa 5: Cálculo e Inversión del Factor de Calidad (Q)\n\n"
            "Este notebook implementa el **Paso 5**: el ajuste final y la estimación del **Factor de Calidad ($Q$)** "
            "a partir de la distribución tiempo-frecuencia calculada por MPD.\n\n"
            "### Metodología de Inversión\n"
            "1. **Isolíneas Chi ($\\chi = f \\cdot t$):** Se proyecta la distribución TF sobre isolíneas de $X = t \\cdot f$.\n"
            "2. **Integración de Energía:** Se integran las amplitudes al cuadrado a lo largo de cada trayectoria de $X$.\n"
            "3. **Ganancia Observada:** Se normaliza la curva de amplitud y se aplica la fórmula de compensación de Wang.\n"
            "4. **Minimización de Residual:** Se compara la ganancia observada contra el catálogo teórico, evaluando el error residual:\n"
            "   $$Q_{\\text{opt}} = \\arg\\min_Q \\sum_X \\left| \\text{Gan}_{\\text{datos}}(X) - \\text{Gan}_{\\text{teórica}}(X, Q) \\right|$$"
        ),
        nbf.v4.new_code_cell(
            "import os\n"
            "import sys\n"
            "import numpy as np\n"
            "import matplotlib.pyplot as plt\n\n"
            "MODULE_DIR = os.path.abspath('../3_Python_Code')\n"
            "if MODULE_DIR not in sys.path:\n"
            "    sys.path.insert(0, MODULE_DIR)\n\n"
            "from config import QInversionConfig\n"
            "from q_estimation import estimate_q_from_tfdist, extract_chi_isoline_curve\n\n"
            "# Cargar la distribución tiempo-frecuencia generada en la Etapa 3\n"
            "data = np.load('data_processed/tfdist_results.npz')\n"
            "tfdist = data['tfdist']\n"
            "t_vec = data['t_vec']\n"
            "f_eval = data['f_eval']\n"
            "print(f'Distribución TF cargada: {tfdist.shape}')"
        ),
        nbf.v4.new_markdown_cell(
            "## 1. Inversión del Factor Q"
        ),
        nbf.v4.new_code_cell(
            "q_cfg = QInversionConfig(\n"
            "    q_min=10.0,\n"
            "    q_max=120.0,\n"
            "    dq=0.5,\n"
            "    x_min=1.0,\n"
            "    x_max=220.0,\n"
            "    dx=1.0,\n"
            "    gain_limit_db=20.0\n"
            ")\n\n"
            "results = estimate_q_from_tfdist(\n"
            "    tfdist=tfdist,\n"
            "    time_vector=t_vec,\n"
            "    freq_vector=f_eval,\n"
            "    q_config=q_cfg\n"
            ")\n\n"
            "q_est = results['q_estimated']\n"
            "print('=' * 60)\n"
            "print(f'  FACTOR DE CALIDAD ESTIMADO: Q = {q_est:.1f}')\n"
            "print('=' * 60)"
        ),
        nbf.v4.new_markdown_cell(
            "## 2. Panel Diagnóstico de Control de Calidad de la Inversión\n\n"
            "Analizamos de forma integrada las 4 componentes del proceso de inversión:\n"
            "1. Espectrograma TF con energía de átomos.\n"
            "2. Curva observada de amplitudes sobre isolíneas Chi.\n"
            "3. Comparación entre la ganancia de los datos y la ganancia teórica óptima.\n"
            "4. Curva de error residual mostrando la convexidad y la ubicación del mínimo global."
        ),
        nbf.v4.new_code_cell(
            "fig, axes = plt.subplots(2, 2, figsize=(14, 9))\n\n"
            "# Panel 1: Espectrograma MPD\n"
            "im = axes[0, 0].imshow(\n"
            "    tfdist,\n"
            "    extent=[f_eval[0], f_eval[-1], t_vec[-1], t_vec[0]],\n"
            "    aspect='auto',\n"
            "    cmap='turbo'\n"
            ")\n"
            "axes[0, 0].set_title('A) Distribución Tiempo - Frecuencia (MPD)', fontweight='bold')\n"
            "axes[0, 0].set_xlabel('Frecuencia (Hz)')\n"
            "axes[0, 0].set_ylabel('Tiempo (s)')\n"
            "plt.colorbar(im, ax=axes[0, 0], label='Energía')\n\n"
            "# Panel 2: Amplitud integrada en isolíneas Chi\n"
            "x_vec = results['x_vector']\n"
            "axes[0, 1].plot(x_vec, results['amp_x_norm'], color='navy', linewidth=1.8)\n"
            "axes[0, 1].set_title('B) Amplitud Integrada sobre Isolíneas $\\\\chi = f \\\\cdot t$', fontweight='bold')\n"
            "axes[0, 1].set_xlabel('X (frecuencia × tiempo)')\n"
            "axes[0, 1].set_ylabel('Amplitud Normalizada')\n"
            "axes[0, 1].grid(True, linestyle=':', alpha=0.6)\n\n"
            "# Panel 3: Ajuste de ganancias\n"
            "axes[1, 0].plot(x_vec, results['gan_datos'], 'r-', label='Ganancia Observada (Datos)', linewidth=2.0)\n"
            "axes[1, 0].plot(x_vec, results['gan_teorica_opt'], 'b--', label=f'Teórica Óptima (Q={q_est:.1f})', linewidth=1.8)\n"
            "axes[1, 0].set_title('C) Ajuste de Ganancia: Observada vs Teórica', fontweight='bold')\n"
            "axes[1, 0].set_xlabel('X (frecuencia × tiempo)')\n"
            "axes[1, 0].set_ylabel('Ganancia Estabilizada')\n"
            "axes[1, 0].legend()\n"
            "axes[1, 0].grid(True, linestyle=':', alpha=0.6)\n\n"
            "# Panel 4: Función objetivo residual vs Q\n"
            "q_vec = results['q_vector']\n"
            "axes[1, 1].plot(q_vec, results['residuals'], color='darkgreen', linewidth=2.0)\n"
            "axes[1, 1].axvline(q_est, color='red', linestyle=':', label=f'Mínimo: Q = {q_est:.1f}')\n"
            "axes[1, 1].set_title('D) Función Objetivo Residual vs Factor Q', fontweight='bold')\n"
            "axes[1, 1].set_xlabel('Factor de Calidad (Q)')\n"
            "axes[1, 1].set_ylabel('Error Residual Acumulado')\n"
            "axes[1, 1].legend()\n"
            "axes[1, 1].grid(True, linestyle=':', alpha=0.6)\n\n"
            "plt.tight_layout()\n"
            "plt.show()"
        ),
        nbf.v4.new_markdown_cell(
            "## 3. Conclusiones e Interpretación Geofísica\n\n"
            "* El algoritmo Matching Pursuit permite aislar con precisión la dispersión y decaimiento de energía en el plano $(t, f)$.\n"
            "* La regularización de Wang evita que las altas frecuencias dominen el cálculo de ganancia, manteniendo el proceso estable ante ruido.\n"
            "* La curva residual presenta un mínimo global bien definido, proporcionando una estimación robusta y reproducible del factor $Q$."
        )
    ]
    return nb

if __name__ == "__main__":
    notebooks = {
        "01_Carga_y_Control_Calidad_Datos.ipynb": create_nb_01(),
        "02_Diccionario_Ondiculas_Morlet.ipynb": create_nb_02(),
        "03_Descomposicion_Matching_Pursuit_MPD.ipynb": create_nb_03(),
        "04_Modelado_Teorico_Atenuacion_Q.ipynb": create_nb_04(),
        "05_Estimacion_e_Inversion_Q.ipynb": create_nb_05()
    }
    
    for filename, nb in notebooks.items():
        filepath = os.path.join(NOTEBOOKS_DIR, filename)
        with open(filepath, "w", encoding="utf-8") as f:
            nbf.write(nb, f)
        print(f"Notebook creado: {filepath}")

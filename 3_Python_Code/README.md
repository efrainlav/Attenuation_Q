# Pipeline de Estimación del Factor de Calidad Sísmico ($Q$) y MPD en Python

Implementación modular y autónoma en **Python** del flujo de estimación de atenuación anelástica ($Q$) mediante descomposición espectral *Matching Pursuit* (MPD) y compensación estabilizada de amplitudes (Wang).

---

## 📁 Estructura del Módulo Python

```text
3_Python_Code/
├── config.py             # Dataclasses de configuración (parámetros geofísicos, rangos y muestreo)
├── io_seismic.py         # Carga/guardado de matrices .mat, SEG-Y y generador de trazas sintéticas
├── morlet_dictionary.py  # Síntesis vectorizada del diccionario 3D de ondículas de Morlet
├── matching_pursuit.py   # Algoritmo MPD con atributos instantáneos de Hilbert y distribución TF
├── q_theoretical.py      # Catálogo teórico analítico Beta(X, Q) y ganancia estabilizada de Wang
├── q_estimation.py       # Integración sobre isolíneas Chi (f·t), cálculo de ganancia e inversión de Q
├── pipeline_demo.py      # Script ejecutable end-to-end con datos sintéticos y generación de figuras
├── test_pipeline.py      # Suite de pruebas unitarias matemáticas y funcionales
├── requirements.txt      # Dependencias del entorno
├── output/               # Gráficos de control de calidad generados y matrices exportadas
└── README.md             # Esta documentación
```

---

## ⚙️ Instalación y Dependencias

El entorno requiere Python 3.8+ y las siguientes librerías científicas estándar:

```bash
pip install -r 3_Python_Code/requirements.txt
```

*Contenido de dependencias:*
* `numpy`: Operaciones matriciales y cálculo vectorizado.
* `scipy`: Transformada de Hilbert (`scipy.signal.hilbert`) e I/O de archivos MATLAB (`scipy.io`).
* `matplotlib`: Visualización de espectrogramas y paneles de control de calidad.
* *(Opcional)* `segyio`: Requerido únicamente si se leen directamente archivos binarios SEG-Y desde Python (`pip install segyio`).

---

## 🚀 Ejecución Rápida

### 1. Ejecutar el Pipeline Completo (Demo Autónoma)
El script `pipeline_demo.py` ejecuta los 5 pasos del flujo, generando una traza sísmica con atenuación conocida, descomponiéndola con MPD y recuperando el valor de $Q$:

```bash
python3 3_Python_Code/pipeline_demo.py
```

**Salidas generadas en `3_Python_Code/output/`:**
* `02_morlet_dictionary_sample.png`: Cortes representativos del diccionario de Morlet (frecuencia, fase y plano tiempo-frecuencia).
* `04_theoretical_gain_catalog.png`: Curvas de ganancia estabilizada teórica en función de $Q$.
* `05_q_inversion_results_panel.png`: Panel integral con la distribución TF, curva de amplitudes sobre isolíneas $X$, ajuste de ganancias y función objetivo residual vs $Q$.
* `resultados_inversion_q_demo.mat`: Archivo compatible con MATLAB con las variables resultantes.

### 2. Ejecutar Pruebas Unitarias
Para validar la coherencia matemática y el correcto funcionamiento de los algoritmos:

```bash
python3 -m unittest 3_Python_Code/test_pipeline.py
```

---

## 🔬 Equivalencia y Mapeo: MATLAB vs Python

| Fase del Flujo | Script Original en MATLAB | Módulo Equivalente en Python |
| :--- | :--- | :--- |
| **Configuración** | Constantes dispersas en scripts | `config.py` (`dataclasses` centralizadas y estructuradas) |
| **Paso 1: I/O de Datos** | `1_Carga_Datos.m` | `io_seismic.py` (`load_mat_seismic`, `load_segy_traces`, `generate_synthetic_seismic_trace`) |
| **Paso 2: Diccionario** | `2_Diccionario de Ondiculas.m` | `morlet_dictionary.py` (`build_morlet_dictionary`, `compute_morlet_atom`) |
| **Paso 3: MPD y TF** | `3_Implementacion_MPD_Quifa18.m` | `matching_pursuit.py` (`decompose_trace_mpd`, `compute_time_frequency_distribution`) |
| **Paso 4: $Q$ Teórico** | `4_Calculo_Teorico_Q.m` | `q_theoretical.py` (`compute_theoretical_attenuation`, `plot_theoretical_gain_catalog`) |
| **Paso 5: Inversión $Q$** | `5_Calculo_Datos_Q.m` | `q_estimation.py` (`estimate_q_from_tfdist`, `extract_chi_isoline_curve`) |

---

## 💻 Ejemplo de Uso Modular en Código Propio

```python
import numpy as np
from config import QInversionConfig, MorletDictionaryConfig, MPDConfig
from io_seismic import load_mat_seismic
from matching_pursuit import decompose_trace_mpd, compute_time_frequency_distribution
from q_estimation import estimate_q_from_tfdist

# 1. Cargar una traza o sección sísmica (previamente exportada a .mat)
# data, meta = load_mat_seismic("data/seismic_line_trace.mat")
# trace = data[:, 0]  # primera traza

# 2. Descomposición con MPD
mpd_res = decompose_trace_mpd(trace, sample_rate_s=0.002)

# 3. Construir distribución tiempo-frecuencia
t_vec = np.arange(len(trace)) * 0.002
f_eval = np.arange(2.0, 120.0, 1.0)
tfdist = compute_time_frequency_distribution(mpd_res["atoms"], t_vec, f_eval)

# 4. Estimar el Factor de Calidad Q
q_results = estimate_q_from_tfdist(tfdist, t_vec, f_eval)
print(f"Factor Q estimado: {q_results['q_estimated']}")
```

# Notebooks Interactivos de Jupyter por Etapas (Python)

Colección de **Jupyter Notebooks (`.ipynb`)** interactivos y auto-contenidos que desarrollan paso a paso cada fase del flujo de estimación del **Factor de Calidad Sísmico ($Q$)** mediante **Matching Pursuit Decomposition (MPD)**.

Los notebooks incluyen explicaciones teóricas, ecuaciones en LaTeX, bloques de código ejecutables y gráficos interactivos pre-renderizados para exploración inmediata.

---

## 📚 Índice de Notebooks

| Notebook | Etapa del Flujo | Descripción y Contenido |
| :--- | :--- | :--- |
| [`01_Carga_y_Control_Calidad_Datos.ipynb`](01_Carga_y_Control_Calidad_Datos.ipynb) | **Paso 1: Datos y QC** | Lectura/generación de trazas, estructuración de tiempo ($dt=2\text{ ms}$), despliegue en *Wiggle* con relleno positivo y análisis espectral (dB). |
| [`02_Diccionario_Ondiculas_Morlet.ipynb`](02_Diccionario_Ondiculas_Morlet.ipynb) | **Paso 2: Diccionario** | Síntesis del diccionario 3D $[Frecuencia \times Fase \times Tiempo]$, análisis de cortes de frecuencia, rotación de fase y huella espectral en el plano $(t, f)$. |
| [`03_Descomposicion_Matching_Pursuit_MPD.ipynb`](03_Descomposicion_Matching_Pursuit_MPD.ipynb) | **Paso 3: MPD y TF** | Algoritmo voraz (*greedy*) con atributos de Hilbert, comparación traza original vs modelada, señal residual $R_N(t)$, convergencia de energía y distribución TF de alta resolución. |
| [`04_Modelado_Teorico_Atenuacion_Q.ipynb`](04_Modelado_Teorico_Atenuacion_Q.ipynb) | **Paso 4: Base Teórica** | Modelado analítico de atenuación espectral $\beta(X, Q)$ y catálogo de curvas de ganancia estabilizada de Wang ($G_{\lim}=20\text{ dB}$) para $Q \in [10, 140]$. |
| [`05_Estimacion_e_Inversion_Q.ipynb`](05_Estimacion_e_Inversion_Q.ipynb) | **Paso 5: Inversión de $Q$** | Proyección e integración en isolíneas $\chi = f \cdot t$, ganancia observada, minimización de residuales $L_1$ contra catálogo teórico y diagnóstico del mínimo de $Q$. |

---

## 🚀 Cómo Ejecutar los Notebooks

### Opción A: Desde Jupyter Lab o Jupyter Notebook
Desde la terminal, sitúese en la raíz del repositorio o dentro de la carpeta `4_Notebooks/` y ejecute:

```bash
cd 4_Notebooks
jupyter notebook
```
o bien:
```bash
jupyter lab
```

### Opción B: Desde Visual Studio Code / IDE
Abra la carpeta del repositorio en su editor. Abra cualquier archivo `.ipynb` y seleccione el kernel de Python instalado en su entorno virtual o sistema.

---

## 💾 Persistencia de Datos Intermedios
Los notebooks comparten información entre etapas mediante la carpeta local `4_Notebooks/data_processed/`:
* `traza_preprocesada.mat` / `.npz`: Traza sísmica y vector de tiempos generados en el Notebook 01.
* `tfdist_results.npz`: Matriz continua de densidad de energía tiempo-frecuencia calculada por MPD en el Notebook 03, consumida por el Notebook 05 para la inversión de $Q$.

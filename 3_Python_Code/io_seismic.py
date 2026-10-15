# ==============================================================================
# Módulo de Entrada/Salida y Generación de Datos Sísmicos
# Provee lectura de SEG-Y, matrices .mat de MATLAB, archivos NumPy .npz
# y un generador de datos sintéticos con atenuación Q para pruebas independientes.
# ==============================================================================

import os
from typing import Optional, Tuple, Dict, Any
import numpy as np
import scipy.io as sio

def load_mat_seismic(filepath: str, var_name: Optional[str] = None) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Carga trazas o secciones sísmicas desde un archivo .mat de MATLAB.
    
    Args:
        filepath: Ruta al archivo .mat
        var_name: Nombre específico de la variable a extraer (opcional).
        
    Returns:
        data: Matriz NumPy con los datos sísmicos (muestras x trazas o vector 1D).
        meta: Diccionario con variables adicionales presentes en el archivo.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Archivo no encontrado: {filepath}")
    
    mat_data = sio.loadmat(filepath)
    clean_dict = {k: v for k, v in mat_data.items() if not k.startswith('__')}
    
    if var_name and var_name in clean_dict:
        data = np.asarray(clean_dict[var_name], dtype=np.float64)
        return data, clean_dict
    
    # Si no se especifica, toma la primera matriz numérica encontrada
    for k, v in clean_dict.items():
        if isinstance(v, np.ndarray) and v.ndim in (1, 2) and np.issubdtype(v.dtype, np.number):
            return np.asarray(v, dtype=np.float64), clean_dict
            
    raise ValueError(f"No se encontró una matriz sísmica válida en {filepath}")

def save_mat_seismic(filepath: str, data_dict: Dict[str, Any]) -> None:
    """Guarda variables y matrices procesadas en formato .mat compatible con MATLAB."""
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    sio.savemat(filepath, data_dict)

def load_segy_traces(filepath: str, trace_indices: Optional[slice] = None) -> Tuple[np.ndarray, float]:
    """
    Lee trazas de un archivo SEG-Y estándar utilizando segyio (si está instalado).
    
    Args:
        filepath: Ruta al archivo .segy o .sgy
        trace_indices: Rango de trazas a leer.
        
    Returns:
        traces: Matriz de datos sísmicos (muestras x trazas).
        sample_rate_s: Tasa de muestreo en segundos.
    """
    try:
        import segyio
    except ImportError:
        raise ImportError(
            "La librería 'segyio' no está instalada. Puede instalarla con 'pip install segyio' "
            "o utilizar archivos en formato .mat / .npz con las funciones de este módulo."
        )
        
    with segyio.open(filepath, "r", ignore_geometry=True) as f:
        dt_us = segyio.tools.dt(f)
        sample_rate_s = dt_us / 1e6
        
        if trace_indices is None:
            raw_traces = [f.trace[i] for i in range(f.tracecount)]
        else:
            raw_traces = [f.trace[i] for i in range(*trace_indices.indices(f.tracecount))]
            
        traces = np.column_stack(raw_traces).astype(np.float64)
        return traces, sample_rate_s

def generate_synthetic_seismic_trace(
    sample_rate_s: float = 0.002,
    t_max_s: float = 2.0,
    peak_frequencies: Tuple[float, ...] = (30.0, 45.0, 25.0),
    arrival_times_s: Tuple[float, ...] = (0.35, 0.75, 1.25),
    amplitudes: Tuple[float, ...] = (1.0, 0.65, 0.40),
    true_q: Optional[float] = 60.0,
    noise_level: float = 0.02,
    random_seed: int = 42
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Genera una traza sísmica sintética con reflectores discretos y atenuación dependiente de Q.
    Permite validar de forma autónoma el pipeline de inversión de Q y MPD.
    
    Args:
        sample_rate_s: Intervalo de muestreo (s).
        t_max_s: Duración total (s).
        peak_frequencies: Frecuencias centrales de las ondículas de entrada (Hz).
        arrival_times_s: Tiempos de arribo de los reflectores (s).
        amplitudes: Amplitudes relativas de cada reflector.
        true_q: Factor de calidad Q aplicado a la propagación (None para señal elástica).
        noise_level: Desviación estándar de ruido gaussiano añadido.
        random_seed: Semilla para reproducibilidad.
        
    Returns:
        t: Vector de tiempo (s).
        trace: Amplitud de la traza sísmica (1D array).
    """
    np.random.seed(random_seed)
    n_samples = int(np.round(t_max_s / sample_rate_s)) + 1
    t = np.linspace(0.0, t_max_s, n_samples)
    trace = np.zeros(n_samples, dtype=np.float64)
    
    for f0, t0, amp in zip(peak_frequencies, arrival_times_s, amplitudes):
        # Atenuación anelástica de la ondícula en función de Q y tiempo t0:
        # La frecuencia dominante decae: f_eff = f0 / (1 + pi*f0*t0 / Q)
        # La amplitud experimenta decaimiento exponencial exp(-pi*f0*t0 / Q)
        if true_q is not None and true_q > 0:
            attenuation_factor = np.exp(-np.pi * f0 * t0 / true_q)
            f_eff = f0 / (1.0 + (np.pi * f0 * t0) / (2.0 * true_q))
            eff_amp = amp * attenuation_factor
        else:
            eff_amp = amp
            f_eff = f0
            
        # Ondícula de Ricker analítica
        tau = t - t0
        pift = (np.pi * f_eff * tau) ** 2
        ricker = eff_amp * (1.0 - 2.0 * pift) * np.exp(-pift)
        trace += ricker
        
    if noise_level > 0:
        trace += np.random.normal(0, noise_level * np.max(np.abs(trace)), size=n_samples)
        
    return t, trace

# ==============================================================================
# Módulo de Descomposición Matching Pursuit (MPD) y Análisis Tiempo-Frecuencia
# Implementa el algoritmo iterativo voraz mediante atributos complejos de Hilbert,
# sustracción de ondículas y generación de la distribución tiempo-frecuencia.
# ==============================================================================

from typing import List, Tuple, Dict, Any, Optional
import numpy as np
from scipy.signal import hilbert

try:
    from .config import MPDConfig, MorletDictionaryConfig
    from .morlet_dictionary import compute_morlet_atom
except ImportError:
    from config import MPDConfig, MorletDictionaryConfig
    from morlet_dictionary import compute_morlet_atom

def energy_calc(signal: np.ndarray) -> float:
    """Calcula la energía cuadrática total de una señal discreta."""
    return float(np.sum(signal ** 2))

def compute_complex_attributes(
    trace: np.ndarray,
    sample_rate_s: float,
    f_min: float = 2.0,
    f_max: float = 200.0
) -> Tuple[int, float, float, float]:
    """
    Estima los atributos complejos instantáneos del átomo dominante mediante la señal analítica:
    - Posición temporal del pico de envolvente (muestra t0_idx).
    - Amplitud envolvente en t0 (amp).
    - Fase instantánea en t0 (phi_deg).
    - Frecuencia instantánea en t0 (freq_hz).
    
    Args:
        trace: Traza sísmica residual (1D array).
        sample_rate_s: Tasa de muestreo (s).
        f_min: Frecuencia mínima admisible (Hz).
        f_max: Frecuencia máxima admisible (Hz).
        
    Returns:
        (t0_idx, amp, phi_deg, freq_hz)
    """
    analytic_signal = hilbert(trace)
    envelope = np.abs(analytic_signal)
    
    # 1. Posición del máximo de la envolvente
    t0_idx = int(np.argmax(envelope))
    amp = float(envelope[t0_idx])
    
    # 2. Fase instantánea en el punto central (en grados, -180 a 180)
    inst_phase = np.angle(analytic_signal)
    phi_deg = float(np.rad2deg(inst_phase[t0_idx]))
    
    # 3. Frecuencia instantánea derivada de la fase desenrollada
    unwrapped_phase = np.unwrap(inst_phase)
    # Diferenciación centrada o de extremos
    if 0 < t0_idx < len(trace) - 1:
        dphi_dt = (unwrapped_phase[t0_idx + 1] - unwrapped_phase[t0_idx - 1]) / (2.0 * sample_rate_s)
    elif t0_idx == 0 and len(trace) > 1:
        dphi_dt = (unwrapped_phase[1] - unwrapped_phase[0]) / sample_rate_s
    elif t0_idx == len(trace) - 1 and len(trace) > 1:
        dphi_dt = (unwrapped_phase[-1] - unwrapped_phase[-2]) / sample_rate_s
    else:
        dphi_dt = 2.0 * np.pi * 30.0  # Fallback
        
    freq_hz = float(dphi_dt / (2.0 * np.pi))
    
    # Acotamiento dentro del rango geofísico válido
    freq_hz = np.clip(freq_hz, f_min, f_max)
    
    return t0_idx, amp, phi_deg, freq_hz

def wave_sub_add(
    trace: np.ndarray,
    modeled_trace: np.ndarray,
    t0_idx: int,
    atom: np.ndarray,
    t_atom_len: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Resta el átomo seleccionado de la traza residual y lo acumula en la traza modelada.
    """
    half_len = (t_atom_len - 1) // 2
    n_samples = len(trace)
    
    start_trace = max(0, t0_idx - half_len)
    end_trace = min(n_samples, t0_idx + half_len + 1)
    
    start_atom = start_trace - (t0_idx - half_len)
    end_atom = start_atom + (end_trace - start_trace)
    
    real_atom = np.real(atom[start_atom:end_atom])
    
    trace[start_trace:end_trace] -= real_atom
    modeled_trace[start_trace:end_trace] += real_atom
    
    return trace, modeled_trace

def decompose_trace_mpd(
    trace: np.ndarray,
    sample_rate_s: float = 0.002,
    mpd_config: Optional[MPDConfig] = None,
    dict_config: Optional[MorletDictionaryConfig] = None
) -> Dict[str, Any]:
    """
    Descompone una traza sísmica individual mediante Matching Pursuit Decomposition (MPD).
    
    Args:
        trace: Traza sísmica 1D.
        sample_rate_s: Intervalo de muestreo (s).
        mpd_config: Parámetros del algoritmo MPD.
        dict_config: Parámetros de la ondícula de Morlet.
        
    Returns:
        Diccionario con:
        - 'modeled_trace': Traza reconstruida sumando los átomos de Morlet.
        - 'residual_trace': Residuo no explicado.
        - 'atoms': Matriz de parámetros [t0_s, amp, phi_deg, freq_hz].
        - 'iterations': Número de iteraciones ejecutadas.
        - 'energy_history': Lista con evolución del porcentaje de energía residual.
    """
    if mpd_config is None:
        mpd_config = MPDConfig()
    if dict_config is None:
        dict_config = MorletDictionaryConfig()
        
    residual = trace.copy().astype(np.float64)
    modeled = np.zeros_like(residual)
    initial_energy = energy_calc(residual)
    
    if initial_energy == 0.0:
        return {
            "modeled_trace": modeled,
            "residual_trace": residual,
            "atoms": np.empty((0, 4)),
            "iterations": 0,
            "energy_history": [0.0]
        }
        
    ts = dict_config.wavelet_half_duration_s
    t_atom = np.arange(-ts, ts + sample_rate_s / 2.0, sample_rate_s)
    t_atom_len = len(t_atom)
    
    atoms = []
    energy_history = [1.0]
    curr_rel_energy = 1.0
    iteration = 0
    
    while curr_rel_energy > mpd_config.energy_cutoff and iteration < mpd_config.max_iterations:
        t0_idx, amp, phi_deg, freq_hz = compute_complex_attributes(
            residual,
            sample_rate_s,
            f_min=dict_config.f_min,
            f_max=dict_config.f_max
        )
        
        # Sintetizar el átomo de Morlet óptimo
        atom = amp * compute_morlet_atom(
            t_atom,
            frequency_hz=freq_hz,
            phase_deg=phi_deg,
            k_param=dict_config.k_param
        )
        
        # Sustraer de la señal residual y sumar a la traza modelada
        residual, modeled = wave_sub_add(residual, modeled, t0_idx, atom, t_atom_len)
        
        t0_s = t0_idx * sample_rate_s
        atoms.append([t0_s, amp, phi_deg, freq_hz])
        
        new_energy = energy_calc(residual)
        curr_rel_energy = new_energy / initial_energy
        energy_history.append(curr_rel_energy)
        iteration += 1
        
    return {
        "modeled_trace": modeled,
        "residual_trace": residual,
        "atoms": np.array(atoms, dtype=np.float64),
        "iterations": iteration,
        "energy_history": energy_history
    }

def compute_time_frequency_distribution(
    atoms: np.ndarray,
    time_vector: np.ndarray,
    freq_vector: np.ndarray,
    k_param: float = 0.2
) -> np.ndarray:
    """
    Construye la distribución espectral continua de tiempo-frecuencia (tfdist)
    a partir de los átomos de ondículas descompuestos por MPD.
    
    Args:
        atoms: Matriz de dimensiones (N_atomos, 4) con [t0_s, amp, phi_deg, freq_hz].
        time_vector: Vector de tiempo de la traza (s).
        freq_vector: Vector de frecuencias de evaluación (Hz).
        k_param: Parámetro de ciclos de Morlet.
        
    Returns:
        tfdist: Matriz 2D de distribución de energía tiempo-frecuencia (n_times x n_freqs).
    """
    n_times = len(time_vector)
    n_freqs = len(freq_vector)
    tfdist = np.zeros((n_times, n_freqs), dtype=np.float64)
    
    if len(atoms) == 0:
        return tfdist
        
    T = time_vector[:, np.newaxis]  # (Nt, 1)
    F = freq_vector[np.newaxis, :]  # (1, Nf)
    
    for atom in atoms:
        t0, amp, _, f0 = atom
        if f0 <= 0:
            continue
            
        # Varianza en tiempo y frecuencia de la ondícula de Morlet
        sigma_t = np.sqrt(np.log(2.0 / k_param)) / (2.0 * np.pi * f0)
        sigma_f = 1.0 / (2.0 * np.pi * sigma_t)
        
        # Envolvente gaussiana 2D en el plano tiempo-frecuencia
        t_gauss = np.exp(-0.5 * ((T - t0) / sigma_t) ** 2)
        f_gauss = np.exp(-0.5 * ((F - f0) / sigma_f) ** 2)
        
        # Acumular la energía del átomo
        atom_energy_plane = (amp ** 2) * (t_gauss @ f_gauss)
        tfdist += atom_energy_plane
        
    return tfdist

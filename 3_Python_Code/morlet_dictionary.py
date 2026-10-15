# ==============================================================================
# Módulo de Construcción del Diccionario de Ondículas de Morlet
# Genera el diccionario redundante 3D (Frecuencia x Fase x Tiempo)
# optimizado mediante operaciones vectorizadas en NumPy.
# ==============================================================================

from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt

try:
    from .config import MorletDictionaryConfig
except ImportError:
    from config import MorletDictionaryConfig

def compute_morlet_atom(
    t: np.ndarray,
    frequency_hz: float,
    phase_deg: float,
    k_param: float = 0.2
) -> np.ndarray:
    """
    Calcula una ondícula de Morlet compleja individual según la formulación del proyecto.
    
    Ecuación:
        w(t) = exp(-t^2 * f^2 * ln(2/k)) * exp(i * 2pi * f * t) * exp(i * phase_rad)
        
    Args:
        t: Vector temporal centrado en 0 (s).
        frequency_hz: Frecuencia central de la ondícula (Hz).
        phase_deg: Desfase en grados (-180 a 180).
        k_param: Parámetro de ciclos de la ondícula (default = 0.2).
        
    Returns:
        Ondícula compleja evaluada en los puntos t.
    """
    phase_rad = np.deg2rad(phase_deg)
    gaussian_envelope = np.exp(- (t ** 2) * (frequency_hz ** 2) * np.log(2.0 / k_param))
    oscillation = np.exp(1j * 2.0 * np.pi * frequency_hz * t)
    phase_shift = np.exp(1j * phase_rad)
    return gaussian_envelope * oscillation * phase_shift

def build_morlet_dictionary(
    config: MorletDictionaryConfig,
    sample_rate_s: float = 0.002
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Construye el diccionario tridimensional de ondículas de Morlet.
    
    Args:
        config: Objeto con parámetros de frecuencia, fase y duración.
        sample_rate_s: Intervalo de muestreo temporal (s).
        
    Returns:
        morlet_dict: Tensor 3D complejo de dimensiones (n_freqs, n_phases, n_time_samples).
        f_vec: Vector 1D de frecuencias (Hz).
        phi_vec: Vector 1D de fases (grados).
        t_vec: Vector 1D de tiempos de la ondícula (s).
    """
    ts = config.wavelet_half_duration_s
    t_vec = np.arange(-ts, ts + sample_rate_s / 2.0, sample_rate_s)
    f_vec = np.arange(config.f_min, config.f_max + config.df / 2.0, config.df)
    phi_vec = np.arange(config.phase_min_deg, config.phase_max_deg + config.dphase_deg / 2.0, config.dphase_deg)
    
    # Construcción vectorizada usando broadcasting tridimensional:
    # F: (Nf, 1, 1), PHI: (1, Nphi, 1), T: (1, 1, Nt)
    F = f_vec[:, np.newaxis, np.newaxis]
    PHI = np.deg2rad(phi_vec)[np.newaxis, :, np.newaxis]
    T = t_vec[np.newaxis, np.newaxis, :]
    
    k = config.k_param
    gaussian = np.exp(- (T ** 2) * (F ** 2) * np.log(2.0 / k))
    carrier = np.exp(1j * 2.0 * np.pi * F * T)
    phase_term = np.exp(1j * PHI)
    
    morlet_dict = gaussian * carrier * phase_term
    return morlet_dict, f_vec, phi_vec, t_vec

def plot_morlet_dictionary_sample(
    morlet_dict: np.ndarray,
    f_vec: np.ndarray,
    phi_vec: np.ndarray,
    t_vec: np.ndarray,
    save_path: str = None
) -> None:
    """Visualiza cortes representativos del diccionario de ondículas."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Corte 1: Variación de frecuencia a fase constante (phi=0)
    idx_phi_0 = int(np.argmin(np.abs(phi_vec)))
    freqs_to_show = [15, 30, 60, 100]
    for f in freqs_to_show:
        idx_f = int(np.argmin(np.abs(f_vec - f)))
        atom = np.real(morlet_dict[idx_f, idx_phi_0, :])
        axes[0].plot(t_vec * 1000, atom, label=f"{f_vec[idx_f]:.0f} Hz")
    axes[0].set_title("Ondículas de Morlet vs Frecuencia (Fase = 0°)")
    axes[0].set_xlabel("Tiempo (ms)")
    axes[0].set_ylabel("Amplitud Real")
    axes[0].legend(loc="upper right")
    axes[0].grid(True, alpha=0.3)
    
    # Corte 2: Variación de fase a frecuencia constante (f=40 Hz)
    idx_f_40 = int(np.argmin(np.abs(f_vec - 40.0)))
    phases_to_show = [-90, -45, 0, 45, 90]
    for ph in phases_to_show:
        idx_ph = int(np.argmin(np.abs(phi_vec - ph)))
        atom = np.real(morlet_dict[idx_f_40, idx_ph, :])
        axes[1].plot(t_vec * 1000, atom, label=f"{phi_vec[idx_ph]:.0f}°")
    axes[1].set_title("Rotación de Fase (f = 40 Hz)")
    axes[1].set_xlabel("Tiempo (ms)")
    axes[1].legend(loc="upper right")
    axes[1].grid(True, alpha=0.3)
    
    # Corte 3: Espectrograma de amplitudes del plano (Frecuencia vs Tiempo) en fase 0
    slice_ft = np.abs(morlet_dict[:, idx_phi_0, :])
    im = axes[2].imshow(
        slice_ft,
        extent=[t_vec[0] * 1000, t_vec[-1] * 1000, f_vec[0], f_vec[-1]],
        aspect="auto",
        origin="lower",
        cmap="viridis"
    )
    axes[2].set_title("Magnitud Envolvente en Plano (t, f)")
    axes[2].set_xlabel("Tiempo (ms)")
    axes[2].set_ylabel("Frecuencia (Hz)")
    plt.colorbar(im, ax=axes[2], label="|Amplitud|")
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200)
    plt.close()

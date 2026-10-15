# ==============================================================================
# Módulo de Configuración - Parámetros de Atenuación Sísmica (Q) y MPD
# ==============================================================================

from dataclasses import dataclass, field
from typing import Tuple, List
import numpy as np

@dataclass
class SeismicSurveyConfig:
    """Configuración general del levantamiento sísmico y parámetros temporales."""
    survey_name: str = "Seismic_Survey_Area"
    line_name: str = "Seismic_Line_A"
    well_name: str = "Well_01"
    
    # Parámetros temporales y de muestreo
    sample_rate_s: float = 0.002  # Tasa de muestreo (2 ms = 0.002 s)
    t_max_s: float = 4.0          # Tiempo máximo de registro (4.0 s)
    
    # Rango de CDPs o trazas de interés para análisis
    selected_cdp_range: Tuple[int, int] = (100, 300)

@dataclass
class MorletDictionaryConfig:
    """Parámetros para la construcción del diccionario de ondículas de Morlet."""
    f_min: float = 2.0       # Frecuencia mínima (Hz)
    f_max: float = 200.0     # Frecuencia máxima (Hz)
    df: float = 1.0          # Paso de frecuencia (Hz)
    
    phase_min_deg: float = -180.0  # Fase inicial (grados)
    phase_max_deg: float = 180.0   # Fase final (grados)
    dphase_deg: float = 1.0        # Paso de fase (grados)
    
    wavelet_half_duration_s: float = 0.5  # Semiduración temporal de la ondícula (ts = 0.5 s)
    k_param: float = 0.2                  # Parámetro de ciclos de la ondícula de Morlet

@dataclass
class MPDConfig:
    """Configuración del algoritmo Matching Pursuit Decomposition."""
    energy_cutoff: float = 0.001  # Umbral de energía residual (0.001 = 0.1% energía restante)
    max_iterations: int = 500     # Número máximo de átomos por traza
    verbose: bool = True

@dataclass
class QInversionConfig:
    """Configuración para el modelado teórico e inversión del Factor de Calidad Q."""
    # Rango de búsqueda para Q
    q_min: float = 1.0
    q_max: float = 150.0
    dq: float = 0.1
    
    # Espacio X = w * t o f * t
    x_min: float = 1.0
    x_max: float = 450.0
    dx: float = 1.0
    
    # Límite de ganancia estabilizada de Wang (en dB)
    gain_limit_db: float = 20.0
    
    @property
    def q_vector(self) -> np.ndarray:
        return np.arange(self.q_min, self.q_max + self.dq, self.dq)
    
    @property
    def x_vector(self) -> np.ndarray:
        return np.arange(self.x_min, self.x_max + self.dx, self.dx)
    
    @property
    def desvi(self) -> float:
        """Factor de regularización/estabilización de Wang."""
        return float(np.exp(-(0.23 * self.gain_limit_db + 1.63)))

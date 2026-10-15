# ==============================================================================
# Módulo de Modelado Teórico del Factor de Calidad (Q)
# Calcula las matrices analíticas Beta(X, Q) y las curvas de Ganancia Estabilizada
# de compensación de amplitudes según la formulación de Wang.
# ==============================================================================

from typing import Tuple, Optional
import numpy as np
import matplotlib.pyplot as plt

try:
    from .config import QInversionConfig
except ImportError:
    from config import QInversionConfig

def compute_theoretical_attenuation(
    config: Optional[QInversionConfig] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Genera las matrices teóricas de decaimiento Beta(X, Q) y ganancia estabilizada.
    
    Ecuaciones:
        Beta(X, Q) = exp( - (2 * pi * X) / (2 * Q) ) = exp( - pi * X / Q )
        desvi = exp( - (0.23 * Glim + 1.63) )
        B1 = Beta + desvi
        B2 = Beta^2 + desvi
        Ganancia_Teorica = B1 / B2
        
    Args:
        config: Configuración de los rangos de X, Q y límite de ganancia.
        
    Returns:
        beta_q: Matriz de decaimiento teórico (len(x) x len(q)).
        gan_teorica: Matriz de ganancia estabilizada teórica (len(x) x len(q)).
        x_vec: Vector 1D del parámetro X = w*t o f*t.
        q_vec: Vector 1D del factor de calidad Q.
    """
    if config is None:
        config = QInversionConfig()
        
    x_vec = config.x_vector  # (Nx,)
    q_vec = config.q_vector  # (Nq,)
    
    # Broadcasting 2D: X como columnas (Nx, 1), Q como filas (1, Nq)
    X = x_vec[:, np.newaxis]
    Q = q_vec[np.newaxis, :]
    
    # 1. Matriz de atenuación teórica Beta(X, Q)
    beta_q = np.exp(- (2.0 * np.pi * X) / (2.0 * Q))
    
    # 2. Factor de estabilización / desviación de ganancia de Wang
    desvi = config.desvi
    
    # 3. Ganancia teórica estabilizada
    b1 = beta_q + desvi
    b2 = (beta_q ** 2) + desvi
    gan_teorica = b1 / b2
    
    return beta_q, gan_teorica, x_vec, q_vec

def plot_theoretical_gain_catalog(
    gan_teorica: np.ndarray,
    x_vec: np.ndarray,
    q_vec: np.ndarray,
    q_values_to_plot: Tuple[float, ...] = (10, 25, 50, 75, 100, 150),
    save_path: Optional[str] = None
) -> None:
    """Grafica el catálogo de curvas teóricas de ganancia para distintos valores de Q."""
    plt.figure(figsize=(9, 5))
    
    for q_target in q_values_to_plot:
        idx_q = int(np.argmin(np.abs(q_vec - q_target)))
        plt.plot(x_vec, gan_teorica[:, idx_q], label=f"Q = {q_vec[idx_q]:.1f}", linewidth=1.8)
        
    plt.title("Curvas Teóricas de Ganancia Estabilizada de Wang vs X (ω·t)", fontsize=12)
    plt.xlabel("X (frecuencia × tiempo acumulado)", fontsize=10)
    plt.ylabel("Ganancia Teórica", fontsize=10)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(title="Factor Q", loc="upper left")
    plt.xlim([x_vec[0], x_vec[-1]])
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=200)
    plt.close()

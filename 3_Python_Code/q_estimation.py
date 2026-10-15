# ==============================================================================
# Módulo de Estimación e Inversión del Factor de Calidad (Q) sobre Datos
# Realiza la proyección sobre isolíneas Chi (X = t * f), cálculo de ganancia
# observada y optimización por ajuste residual contra el catálogo teórico.
# ==============================================================================

from typing import Tuple, Dict, Any, Optional
import numpy as np
import matplotlib.pyplot as plt

try:
    from .config import QInversionConfig
    from .q_theoretical import compute_theoretical_attenuation
except ImportError:
    from config import QInversionConfig
    from q_theoretical import compute_theoretical_attenuation

def extract_chi_isoline_curve(
    tfdist: np.ndarray,
    time_vector: np.ndarray,
    freq_vector: np.ndarray,
    x_bins: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integra la energía de la distribución tiempo-frecuencia a lo largo de isolíneas
    del producto tiempo-frecuencia: X = t * f (Chi).
    
    Args:
        tfdist: Matriz tiempo-frecuencia de energía (len(t) x len(f)).
        time_vector: Vector de tiempos (s).
        freq_vector: Vector de frecuencias (Hz).
        x_bins: Valores discretos de X para agrupar e integrar.
        
    Returns:
        x_valid: Vector de valores de X evaluados.
        amp_x: Amplitud acumulada a lo largo de cada isolínea de X.
    """
    # Matriz 2D de coordenadas X = tiempo x frecuencia
    T = time_vector[:, np.newaxis]
    F = freq_vector[np.newaxis, :]
    X_grid = T * F  # Dimensiones (Nt, Nf)
    
    dx = x_bins[1] - x_bins[0] if len(x_bins) > 1 else 1.0
    half_dx = dx / 2.0
    
    amp_x = np.zeros(len(x_bins), dtype=np.float64)
    
    for i, x_val in enumerate(x_bins):
        # Máscara booleana para seleccionar celdas en el rango [x_val - half_dx, x_val + half_dx)
        mask = (X_grid >= (x_val - half_dx)) & (X_grid < (x_val + half_dx))
        if np.any(mask):
            amp_x[i] = np.sum(tfdist[mask] ** 2)
            
    return x_bins, amp_x

def estimate_q_from_tfdist(
    tfdist: np.ndarray,
    time_vector: np.ndarray,
    freq_vector: np.ndarray,
    q_config: Optional[QInversionConfig] = None
) -> Dict[str, Any]:
    """
    Estima el Factor de Calidad (Q) óptimo a partir de la distribución tiempo-frecuencia.
    
    Procedimiento:
    1. Proyecta las amplitudes sobre isolíneas X = t * f.
    2. Normaliza la curva observada y aplica la función de ganancia estabilizada de Wang.
    3. Compara la ganancia de los datos frente al catálogo teórico.
    4. Encuentra el Q óptimo que minimiza la norma L1 del error residual.
    
    Args:
        tfdist: Matriz 2D de distribución tiempo-frecuencia.
        time_vector: Vector temporal (s).
        freq_vector: Vector de frecuencias (Hz).
        q_config: Configuración de los rangos de Q, X y límite de ganancia.
        
    Returns:
        Diccionario con resultados:
        - 'q_estimated': Valor numérico de Q estimado.
        - 'residuals': Vector de residuales para cada Q evaluado.
        - 'q_vector': Vector de búsqueda de Q.
        - 'x_vector': Vector de X evaluado.
        - 'gan_datos': Curva de ganancia calculada sobre los datos reales.
        - 'gan_teorica_opt': Curva teórica correspondiente al Q óptimo.
        - 'amp_x_norm': Curva de amplitudes observada normalizada.
    """
    if q_config is None:
        q_config = QInversionConfig()
        
    # 1. Catálogo teórico
    _, gan_teorica, x_vec, q_vec = compute_theoretical_attenuation(q_config)
    
    # 2. Curva observada de amplitudes sobre isolíneas X
    _, amp_x = extract_chi_isoline_curve(tfdist, time_vector, freq_vector, x_vec)
    
    max_amp = np.max(amp_x)
    if max_amp == 0:
        max_amp = 1.0
        
    loc_max = int(np.argmax(amp_x))
    amp_x_norm = amp_x / max_amp
    
    # 3. Cálculo de la curva de ganancia de los datos (formulación de Wang)
    desvi = q_config.desvi
    a_term = amp_x_norm + desvi
    b_term = (amp_x_norm ** 2) + desvi
    gan_datos = a_term / b_term
    
    # Fijar ganancia unitaria antes del punto máximo de amplitud
    gan_datos[:loc_max + 1] = 1.0
    
    # 4. Minimización del error residual entre datos observados y catálogo teórico:
    # residual(Q) = sum | Gan_Teorica(X, Q) - Gan_Datos(X) |
    residuals = np.sum(np.abs(gan_teorica - gan_datos[:, np.newaxis]), axis=0)
    
    opt_idx = int(np.argmin(residuals))
    q_estimated = float(q_vec[opt_idx])
    gan_teorica_opt = gan_teorica[:, opt_idx]
    
    return {
        "q_estimated": q_estimated,
        "residuals": residuals,
        "q_vector": q_vec,
        "x_vector": x_vec,
        "gan_datos": gan_datos,
        "gan_teorica_opt": gan_teorica_opt,
        "amp_x_norm": amp_x_norm
    }

def plot_q_inversion_results(
    results: Dict[str, Any],
    tfdist: np.ndarray,
    time_vector: np.ndarray,
    freq_vector: np.ndarray,
    save_path: Optional[str] = None
) -> None:
    """Genera un panel multipanel con el análisis integral de inversión de Q."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: Distribución Tiempo-Frecuencia (Espectrograma)
    im = axes[0, 0].imshow(
        tfdist,
        extent=[freq_vector[0], freq_vector[-1], time_vector[-1], time_vector[0]],
        aspect="auto",
        cmap="jet"
    )
    axes[0, 0].set_title("Distribución Tiempo - Frecuencia (MPD)")
    axes[0, 0].set_xlabel("Frecuencia (Hz)")
    axes[0, 0].set_ylabel("Tiempo (s)")
    plt.colorbar(im, ax=axes[0, 0], label="Energía")
    
    # Panel 2: Curva de Amplitud Normalizada sobre isolíneas Chi
    x_vec = results["x_vector"]
    axes[0, 1].plot(x_vec, results["amp_x_norm"], color="navy", linewidth=2.0)
    axes[0, 1].set_title("Amplitud Integrada sobre Isolíneas X (ω·t)")
    axes[0, 1].set_xlabel("X (frecuencia × tiempo)")
    axes[0, 1].set_ylabel("Amplitud Normalizada")
    axes[0, 1].grid(True, linestyle="--", alpha=0.5)
    
    # Panel 3: Comparación Ganancia Observada vs Ganancia Teórica Óptima
    q_opt = results["q_estimated"]
    axes[1, 0].plot(x_vec, results["gan_datos"], "r-", label="Ganancia Observada (Datos)", linewidth=2.0)
    axes[1, 0].plot(x_vec, results["gan_teorica_opt"], "b--", label=f"Ganancia Teórica Óptima (Q = {q_opt:.1f})", linewidth=2.0)
    axes[1, 0].set_title("Ajuste de Ganancia: Observada vs Teórica")
    axes[1, 0].set_xlabel("X (frecuencia × tiempo)")
    axes[1, 0].set_ylabel("Ganancia Estabilizada")
    axes[1, 0].legend()
    axes[1, 0].grid(True, linestyle="--", alpha=0.5)
    
    # Panel 4: Curva de Error Residual vs Factor Q
    q_vec = results["q_vector"]
    axes[1, 1].plot(q_vec, results["residuals"], color="darkgreen", linewidth=2.0)
    axes[1, 1].axvline(q_opt, color="red", linestyle=":", label=f"Q Mínimo = {q_opt:.1f}")
    axes[1, 1].set_title("Función Objetivo Residual vs Factor Q")
    axes[1, 1].set_xlabel("Factor de Calidad (Q)")
    axes[1, 1].set_ylabel("Error Residual Acumulado")
    axes[1, 1].legend()
    axes[1, 1].grid(True, linestyle="--", alpha=0.5)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200)
    plt.close()

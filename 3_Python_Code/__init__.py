"""
Paquete Python para Estimación del Factor de Calidad Sísmico (Q) mediante Matching Pursuit Decomposition (MPD).
Todos los nombres de proyectos, áreas y líneas han sido anonimizados.
"""

from .config import (
    SeismicSurveyConfig,
    MorletDictionaryConfig,
    MPDConfig,
    QInversionConfig
)
from .io_seismic import (
    load_mat_seismic,
    save_mat_seismic,
    generate_synthetic_seismic_trace
)
from .morlet_dictionary import (
    compute_morlet_atom,
    build_morlet_dictionary
)
from .matching_pursuit import (
    decompose_trace_mpd,
    compute_time_frequency_distribution
)
from .q_theoretical import (
    compute_theoretical_attenuation
)
from .q_estimation import (
    extract_chi_isoline_curve,
    estimate_q_from_tfdist
)

__all__ = [
    "SeismicSurveyConfig",
    "MorletDictionaryConfig",
    "MPDConfig",
    "QInversionConfig",
    "load_mat_seismic",
    "save_mat_seismic",
    "generate_synthetic_seismic_trace",
    "compute_morlet_atom",
    "build_morlet_dictionary",
    "decompose_trace_mpd",
    "compute_time_frequency_distribution",
    "compute_theoretical_attenuation",
    "extract_chi_isoline_curve",
    "estimate_q_from_tfdist",
]

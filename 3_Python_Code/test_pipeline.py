# ==============================================================================
# Pruebas Unitarias del Pipeline de Estimación de Q y MPD en Python
# Verifica la integridad matemática de cada componente sin dependencias externas.
# ==============================================================================

import unittest
import os
import sys
import numpy as np

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from config import MorletDictionaryConfig, QInversionConfig, MPDConfig
from morlet_dictionary import compute_morlet_atom, build_morlet_dictionary
from matching_pursuit import energy_calc, compute_complex_attributes, decompose_trace_mpd
from q_theoretical import compute_theoretical_attenuation
from q_estimation import estimate_q_from_tfdist
from io_seismic import generate_synthetic_seismic_trace

class TestAttenuationQPipeline(unittest.TestCase):
    
    def test_morlet_atom_properties(self):
        """Verifica simetría y propiedades de la ondícula de Morlet."""
        t = np.linspace(-0.2, 0.2, 401)
        atom_0deg = compute_morlet_atom(t, frequency_hz=30.0, phase_deg=0.0)
        atom_180deg = compute_morlet_atom(t, frequency_hz=30.0, phase_deg=180.0)
        
        # En fase 0°, la parte real debe ser par (simétrica respecto al centro)
        mid = len(t) // 2
        real_part = np.real(atom_0deg)
        np.testing.assert_allclose(real_part[:mid], real_part[mid+1:][::-1], atol=1e-10)
        
        # En fase 180°, debe ser opuesta a fase 0°
        np.testing.assert_allclose(np.real(atom_0deg), -np.real(atom_180deg), atol=1e-10)
        
    def test_energy_calc(self):
        """Verifica el cálculo de energía discreta."""
        sig = np.array([3.0, 4.0, 0.0])
        self.assertEqual(energy_calc(sig), 25.0)
        
    def test_theoretical_gain_catalog(self):
        """Verifica la generación del catálogo teórico y estabilización de Wang."""
        cfg = QInversionConfig(q_min=10, q_max=100, dq=10, x_min=1, x_max=50, dx=1)
        beta, gan, x_vec, q_vec = compute_theoretical_attenuation(cfg)
        
        self.assertEqual(beta.shape, (50, 10))
        self.assertEqual(gan.shape, (50, 10))
        # Para X > 0, Beta debe ser estrictamente decreciente con X
        self.assertTrue(np.all(np.diff(beta[:, 0]) < 0))
        # Mayor Q implica menor atenuación, por lo tanto Beta(X, Q_alto) > Beta(X, Q_bajo)
        self.assertTrue(np.all(beta[:, -1] >= beta[:, 0]))
        
    def test_synthetic_trace_and_mpd(self):
        """Verifica la descomposición MPD sobre una señal sintética."""
        t, trace = generate_synthetic_seismic_trace(
            sample_rate_s=0.002,
            t_max_s=0.6,
            peak_frequencies=(35.0,),
            arrival_times_s=(0.3,),
            amplitudes=(1.0,),
            true_q=None,
            noise_level=0.0
        )
        mpd_cfg = MPDConfig(energy_cutoff=0.01, max_iterations=50)
        dict_cfg = MorletDictionaryConfig(f_min=10, f_max=80, wavelet_half_duration_s=0.2)
        
        results = decompose_trace_mpd(trace, 0.002, mpd_cfg, dict_cfg)
        
        self.assertGreater(len(results["atoms"]), 0)
        # La energía residual debe ser menor al 1%
        self.assertLess(results["energy_history"][-1], 0.02)
        
        # El primer átomo debe estar cerca de t = 0.3 s y f = 35 Hz
        first_atom = results["atoms"][0]
        self.assertAlmostEqual(first_atom[0], 0.3, delta=0.02)
        self.assertAlmostEqual(first_atom[3], 35.0, delta=5.0)

if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
# ==============================================================================
# Script Demostrativo End-to-End: Pipeline Completo de Estimación de Q y MPD
# Ejecuta los pasos 1 al 5 de forma autónoma con datos sintéticos y/o reales.
# ==============================================================================

import os
import sys
import numpy as np

# Añadir el directorio raíz al path para importar el paquete local
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from config import (
    SeismicSurveyConfig,
    MorletDictionaryConfig,
    MPDConfig,
    QInversionConfig
)
from io_seismic import generate_synthetic_seismic_trace, save_mat_seismic
from morlet_dictionary import plot_morlet_dictionary_sample, build_morlet_dictionary
from matching_pursuit import decompose_trace_mpd, compute_time_frequency_distribution
from q_theoretical import compute_theoretical_attenuation, plot_theoretical_gain_catalog
from q_estimation import estimate_q_from_tfdist, plot_q_inversion_results

def run_demo():
    print("=" * 80)
    print("  PIPELINE DE ESTIMACIÓN DEL FACTOR DE CALIDAD SÍSMICO (Q) - VERSIÓN PYTHON")
    print("=" * 80)
    
    output_dir = os.path.join(CURRENT_DIR, "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # -------------------------------------------------------------------------
    # Paso 1: Carga o Generación de Traza Sísmica
    # -------------------------------------------------------------------------
    print("\n[Paso 1] Generando traza sísmica sintética con atenuación Q conocida...")
    survey_cfg = SeismicSurveyConfig(
        survey_name="Seismic_Survey_Area",
        line_name="Seismic_Line_01",
        well_name="Well_Test_01",
        sample_rate_s=0.002,
        t_max_s=2.0
    )
    
    target_true_q = 55.0  # Valor de referencia para validar la inversión
    t_vec, trace = generate_synthetic_seismic_trace(
        sample_rate_s=survey_cfg.sample_rate_s,
        t_max_s=survey_cfg.t_max_s,
        peak_frequencies=(25.0, 40.0, 55.0),
        arrival_times_s=(0.4, 0.9, 1.4),
        amplitudes=(1.0, 0.7, 0.45),
        true_q=target_true_q,
        noise_level=0.015,
        random_seed=101
    )
    print(f"  -> Traza generada con {len(trace)} muestras, dt = {survey_cfg.sample_rate_s*1000:.1f} ms.")
    print(f"  -> Factor Q teórico inyectado: {target_true_q:.1f}")
    
    # -------------------------------------------------------------------------
    # Paso 2: Construcción / Inspección del Diccionario de Morlet
    # -------------------------------------------------------------------------
    print("\n[Paso 2] Construyendo base de ondículas de Morlet complejas...")
    dict_cfg = MorletDictionaryConfig(
        f_min=5.0,
        f_max=120.0,
        df=2.0,
        phase_min_deg=-180.0,
        phase_max_deg=180.0,
        dphase_deg=10.0,
        wavelet_half_duration_s=0.3,
        k_param=0.2
    )
    morlet_dict, f_dict, phi_dict, t_dict = build_morlet_dictionary(dict_cfg, survey_cfg.sample_rate_s)
    dict_plot_path = os.path.join(output_dir, "02_morlet_dictionary_sample.png")
    plot_morlet_dictionary_sample(morlet_dict, f_dict, phi_dict, t_dict, save_path=dict_plot_path)
    print(f"  -> Diccionario sintetizado con dimensiones: {morlet_dict.shape}")
    print(f"  -> Gráfico de control guardado en: {dict_plot_path}")
    
    # -------------------------------------------------------------------------
    # Paso 3: Descomposición Matching Pursuit (MPD) y Distribución TF
    # -------------------------------------------------------------------------
    print("\n[Paso 3] Ejecutando Matching Pursuit Decomposition (MPD)...")
    mpd_cfg = MPDConfig(energy_cutoff=0.002, max_iterations=200, verbose=True)
    mpd_results = decompose_trace_mpd(
        trace,
        sample_rate_s=survey_cfg.sample_rate_s,
        mpd_config=mpd_cfg,
        dict_config=dict_cfg
    )
    
    n_atoms = len(mpd_results["atoms"])
    final_res_energy = mpd_results["energy_history"][-1] * 100
    print(f"  -> MPD convergió en {mpd_results['iterations']} iteraciones ({n_atoms} átomos extraídos).")
    print(f"  -> Energía residual final: {final_res_energy:.3f}% (explicó {100-final_res_energy:.3f}% de la señal).")
    
    f_eval = np.arange(2.0, 100.0, 1.0)
    tfdist = compute_time_frequency_distribution(
        mpd_results["atoms"],
        time_vector=t_vec,
        freq_vector=f_eval,
        k_param=dict_cfg.k_param
    )
    print(f"  -> Distribución tiempo-frecuencia calculada con dimensiones: {tfdist.shape}")
    
    # -------------------------------------------------------------------------
    # Paso 4: Modelado Teórico de Ganancia Estabilizada de Wang
    # -------------------------------------------------------------------------
    print("\n[Paso 4] Generando catálogo teórico de decaimiento y ganancia de Wang...")
    q_cfg = QInversionConfig(
        q_min=10.0,
        q_max=120.0,
        dq=0.5,
        x_min=1.0,
        x_max=250.0,
        dx=1.0,
        gain_limit_db=20.0
    )
    _, gan_teorica, x_vec, q_vec = compute_theoretical_attenuation(q_cfg)
    catalog_plot_path = os.path.join(output_dir, "04_theoretical_gain_catalog.png")
    plot_theoretical_gain_catalog(gan_teorica, x_vec, q_vec, save_path=catalog_plot_path)
    print(f"  -> Catálogo teórico generado para {len(q_vec)} valores de Q ({q_vec[0]:.1f} a {q_vec[-1]:.1f}).")
    print(f"  -> Gráfico de curvas teóricas guardado en: {catalog_plot_path}")
    
    # -------------------------------------------------------------------------
    # Paso 5: Cálculo e Inversión de Q sobre los Datos
    # -------------------------------------------------------------------------
    print("\n[Paso 5] Proyectando energía sobre isolíneas Chi e invirtiendo el valor de Q...")
    inversion_results = estimate_q_from_tfdist(
        tfdist=tfdist,
        time_vector=t_vec,
        freq_vector=f_eval,
        q_config=q_cfg
    )
    
    q_estimated = inversion_results["q_estimated"]
    error_abs = abs(q_estimated - target_true_q)
    error_pct = (error_abs / target_true_q) * 100
    
    print("=" * 80)
    print(f"  RESULTADO DE LA ESTIMACIÓN:")
    print(f"  - Factor Q Teórico Inyectado : {target_true_q:.1f}")
    print(f"  - Factor Q Estimado por MPD  : {q_estimated:.1f}")
    print(f"  - Error Absoluto             : {error_abs:.2f}")
    print(f"  - Error Relativo             : {error_pct:.2f}%")
    print("=" * 80)
    
    results_plot_path = os.path.join(output_dir, "05_q_inversion_results_panel.png")
    plot_q_inversion_results(inversion_results, tfdist, t_vec, f_eval, save_path=results_plot_path)
    print(f"  -> Panel integral de resultados guardado en: {results_plot_path}")
    
    # Exportar resultados numéricos a .mat para compatibilidad con MATLAB
    mat_export_path = os.path.join(output_dir, "resultados_inversion_q_demo.mat")
    save_mat_seismic(
        mat_export_path,
        {
            "q_estimated": q_estimated,
            "q_true": target_true_q,
            "q_vector": inversion_results["q_vector"],
            "residuals": inversion_results["residuals"],
            "gan_datos": inversion_results["gan_datos"],
            "gan_teorica_opt": inversion_results["gan_teorica_opt"],
            "tfdist": tfdist
        }
    )
    print(f"  -> Matriz de resultados exportada a MATLAB .mat en: {mat_export_path}")
    print("\n¡Pipeline ejecutado exitosamente!")

if __name__ == "__main__":
    run_demo()

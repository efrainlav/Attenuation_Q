#!/usr/bin/env python3
# ==============================================================================
# Script de Generación de Figuras de Calidad para la Documentación
# Genera imágenes científicas en alta resolución de cada paso del flujo.
# ==============================================================================

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(CURRENT_DIR)
ASSETS_DIR = os.path.join(REPO_ROOT, "assets", "images")
os.makedirs(ASSETS_DIR, exist_ok=True)

if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from config import MorletDictionaryConfig, MPDConfig, QInversionConfig
from io_seismic import generate_synthetic_seismic_trace
from morlet_dictionary import compute_morlet_atom, build_morlet_dictionary
from matching_pursuit import decompose_trace_mpd, compute_time_frequency_distribution
from q_theoretical import compute_theoretical_attenuation
from q_estimation import estimate_q_from_tfdist

plt.rcParams.update({
    'font.sans-serif': 'DejaVu Sans',
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 13
})

def generate_figure_1_seismic_qc():
    """Figura Paso 1: Carga y control de calidad de trazas sísmicas."""
    dt = 0.002
    t, main_trace = generate_synthetic_seismic_trace(
        sample_rate_s=dt,
        t_max_s=2.0,
        peak_frequencies=(25.0, 38.0, 50.0),
        arrival_times_s=(0.35, 0.85, 1.35),
        amplitudes=(1.0, 0.72, 0.48),
        true_q=60.0,
        noise_level=0.015,
        random_seed=42
    )
    
    # Generar un gather sintético adyacente (15 trazas)
    n_traces = 15
    gather = np.zeros((len(t), n_traces))
    for i in range(n_traces):
        offset_factor = 1.0 - 0.015 * abs(i - n_traces // 2)
        jitter = np.random.normal(0, 0.01, size=len(t))
        gather[:, i] = main_trace * offset_factor + jitter
    
    fig = plt.figure(figsize=(15, 5))
    gs = GridSpec(1, 3, width_ratios=[1.2, 1.4, 1.0], wspace=0.28)
    
    # Panel A: Traza individual
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(t, main_trace, color='#1f77b4', linewidth=1.2)
    ax1.set_title("A) Traza Sísmica Individual", fontweight='bold')
    ax1.set_xlabel("Tiempo (s)")
    ax1.set_ylabel("Amplitud Normalizada")
    ax1.set_xlim(0, 2.0)
    ax1.grid(True, linestyle=":", alpha=0.6)
    
    # Panel B: Despliegue Wiggle / Gather sísmico
    ax2 = fig.add_subplot(gs[1])
    scale = 0.8
    for i in range(n_traces):
        tr = gather[:, i]
        tr_norm = tr / np.max(np.abs(gather)) * scale
        ax2.plot(i + tr_norm, t, color='black', linewidth=0.7)
        ax2.fill_betweenx(t, i, i + tr_norm, where=(tr_norm > 0), color='#2ca02c', alpha=0.6)
    ax2.set_title("B) Sección Sísmica / Wiggle Display", fontweight='bold')
    ax2.set_xlabel("Número de Traza / CDP")
    ax2.set_ylabel("Tiempo Doble (s)")
    ax2.set_ylim(2.0, 0)  # Invertir eje de tiempo (estándar sísmico)
    ax2.set_xlim(-1, n_traces)
    ax2.grid(True, linestyle=":", alpha=0.5)
    
    # Panel C: Espectro de amplitud
    ax3 = fig.add_subplot(gs[2])
    n_fft = 2048
    fft_vals = np.abs(np.fft.rfft(main_trace, n=n_fft))
    freqs = np.fft.rfftfreq(n_fft, d=dt)
    fft_db = 20 * np.log10(fft_vals / np.max(fft_vals) + 1e-6)
    ax3.plot(freqs, fft_db, color='#d62728', linewidth=1.5)
    ax3.set_title("C) Espectro de Amplitud (dB)", fontweight='bold')
    ax3.set_xlabel("Frecuencia (Hz)")
    ax3.set_ylabel("Amplitud (dB)")
    ax3.set_xlim(0, 100)
    ax3.set_ylim(-40, 5)
    ax3.grid(True, linestyle=":", alpha=0.6)
    
    fig.suptitle("Paso 1: Carga, Visualización y Control de Calidad de Datos Sísmicos", fontsize=14, y=1.02)
    out_path = os.path.join(ASSETS_DIR, "01_seismic_trace_qc.png")
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Figura 1 guardada en: {out_path}")
    return t, main_trace

def generate_figure_2_morlet_dictionary():
    """Figura Paso 2: Base y diccionario de ondículas de Morlet."""
    dt = 0.002
    ts = 0.25
    t_atom = np.arange(-ts, ts + dt/2, dt)
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    
    # Panel A: Frecuencias variables a fase 0
    freqs = [15, 30, 50, 75]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    for f, c in zip(freqs, colors):
        atom = compute_morlet_atom(t_atom, frequency_hz=f, phase_deg=0.0, k_param=0.2)
        axes[0].plot(t_atom * 1000, np.real(atom), label=f"{f} Hz", color=c, linewidth=1.5)
    axes[0].set_title("A) Ondículas vs Frecuencia (Fase = 0°)", fontweight='bold')
    axes[0].set_xlabel("Tiempo (ms)")
    axes[0].set_ylabel("Amplitud Real")
    axes[0].legend(loc="upper right", framealpha=0.9)
    axes[0].grid(True, linestyle=":", alpha=0.6)
    
    # Panel B: Rotación de fase a 35 Hz
    phases = [-90, -45, 0, 45, 90]
    phase_colors = ['#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']
    for ph, c in zip(phases, phase_colors):
        atom = compute_morlet_atom(t_atom, frequency_hz=35.0, phase_deg=ph, k_param=0.2)
        axes[1].plot(t_atom * 1000, np.real(atom), label=f"$\\phi = {ph}^\\circ$", color=c, linewidth=1.4)
    axes[1].set_title("B) Rotación de Fase (f = 35 Hz)", fontweight='bold')
    axes[1].set_xlabel("Tiempo (ms)")
    axes[1].set_ylabel("Amplitud Real")
    axes[1].legend(loc="upper right", framealpha=0.9)
    axes[1].grid(True, linestyle=":", alpha=0.6)
    
    # Panel C: Mapa de energía del átomo en plano (tiempo, frecuencia)
    f_axis = np.linspace(5, 100, 200)
    t_axis = np.linspace(-150, 150, 200)
    T, F = np.meshgrid(t_axis / 1000, f_axis)
    f0 = 40.0
    sigma_t = np.sqrt(np.log(2.0 / 0.2)) / (2.0 * np.pi * f0)
    sigma_f = 1.0 / (2.0 * np.pi * sigma_t)
    envelope_2d = np.exp(-0.5 * (T / sigma_t)**2 - 0.5 * ((F - f0) / sigma_f)**2)
    
    im = axes[2].contourf(t_axis, f_axis, envelope_2d, levels=25, cmap="magma")
    axes[2].set_title("C) Huella Tiempo-Frecuencia del Átomo", fontweight='bold')
    axes[2].set_xlabel("Tiempo (ms)")
    axes[2].set_ylabel("Frecuencia (Hz)")
    plt.colorbar(im, ax=axes[2], label="Energía Relativa")
    
    fig.suptitle("Paso 2: Base Tridimensional de Ondículas Complejas de Morlet", fontsize=14, y=1.02)
    out_path = os.path.join(ASSETS_DIR, "02_morlet_wavelet_dictionary.png")
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Figura 2 guardada en: {out_path}")

def generate_figure_3_mpd_and_tfd(t_vec, trace):
    """Figura Paso 3: Descomposición MPD y Distribución Tiempo-Frecuencia."""
    dt = t_vec[1] - t_vec[0]
    mpd_cfg = MPDConfig(energy_cutoff=0.003, max_iterations=120)
    dict_cfg = MorletDictionaryConfig(f_min=5.0, f_max=100.0, df=1.0, wavelet_half_duration_s=0.25)
    
    mpd_res = decompose_trace_mpd(trace, sample_rate_s=dt, mpd_config=mpd_cfg, dict_config=dict_cfg)
    
    f_eval = np.arange(2.0, 90.0, 1.0)
    tfdist = compute_time_frequency_distribution(
        mpd_res["atoms"],
        time_vector=t_vec,
        freq_vector=f_eval,
        k_param=0.2
    )
    
    fig = plt.figure(figsize=(15, 9))
    gs = GridSpec(2, 2, height_ratios=[1.0, 1.2], hspace=0.32, wspace=0.25)
    
    # Subplot 1: Traza Original vs Reconstruida
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(t_vec, trace, 'k-', label="Traza Original", linewidth=1.3, alpha=0.8)
    ax1.plot(t_vec, mpd_res["modeled_trace"], 'r--', label="Traza Reconstruida (MPD)", linewidth=1.2)
    ax1.set_title("A) Comparación: Traza Original vs Modelo MPD", fontweight='bold')
    ax1.set_xlabel("Tiempo (s)")
    ax1.set_ylabel("Amplitud")
    ax1.legend(loc="upper right")
    ax1.grid(True, linestyle=":", alpha=0.6)
    
    # Subplot 2: Residuo y Convergencia de Energía
    ax2 = fig.add_subplot(gs[0, 1])
    iterations = np.arange(len(mpd_res["energy_history"]))
    energy_pct = np.array(mpd_res["energy_history"]) * 100
    ax2.plot(iterations, energy_pct, color='#2ca02c', linewidth=2.0)
    ax2.axhline(mpd_cfg.energy_cutoff * 100, color='red', linestyle=':', label=f"Umbral Cutoff ({mpd_cfg.energy_cutoff*100:.1f}%)")
    ax2.set_title("B) Curva de Convergencia de Energía Residual", fontweight='bold')
    ax2.set_xlabel("Iteración MPD (Átomo extraído)")
    ax2.set_ylabel("Energía Residual (%)")
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, linestyle=":", which="both", alpha=0.6)
    
    # Subplot 3: Traza Residual
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(t_vec, mpd_res["residual_trace"], color='#7f7f7f', linewidth=1.0)
    ax3.set_title("C) Señal Residual $R_N(t)$ no Explicada", fontweight='bold')
    ax3.set_xlabel("Tiempo (s)")
    ax3.set_ylabel("Amplitud Residual")
    ax3.grid(True, linestyle=":", alpha=0.6)
    
    # Subplot 4: Espectrograma Tiempo-Frecuencia de Alta Resolución
    ax4 = fig.add_subplot(gs[1, 1])
    im = ax4.imshow(
        tfdist,
        extent=[f_eval[0], f_eval[-1], t_vec[-1], t_vec[0]],
        aspect="auto",
        cmap="turbo"
    )
    ax4.set_title("D) Distribución Tiempo-Frecuencia de Alta Resolución", fontweight='bold')
    ax4.set_xlabel("Frecuencia (Hz)")
    ax4.set_ylabel("Tiempo (s)")
    plt.colorbar(im, ax=ax4, label="Densidad de Energía")
    
    fig.suptitle("Paso 3: Descomposición Espectral Matching Pursuit y Distribución TF", fontsize=14, y=0.98)
    out_path = os.path.join(ASSETS_DIR, "03_mpd_decomposition_and_tfd.png")
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Figura 3 guardada en: {out_path}")
    return tfdist, f_eval

def generate_figure_4_and_5_q_workflow(t_vec, tfdist, f_eval):
    """Figuras Pasos 4 y 5: Modelado teórico e inversión del factor Q."""
    q_cfg = QInversionConfig(
        q_min=10.0,
        q_max=120.0,
        dq=0.5,
        x_min=1.0,
        x_max=220.0,
        dx=1.0,
        gain_limit_db=20.0
    )
    beta, gan_teorica, x_vec, q_vec = compute_theoretical_attenuation(q_cfg)
    
    # Figura 4: Catálogo teórico
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    q_targets = [15, 30, 50, 75, 100]
    for q_t in q_targets:
        idx_q = int(np.argmin(np.abs(q_vec - q_t)))
        axes[0].plot(x_vec, beta[:, idx_q], label=f"Q = {q_vec[idx_q]:.0f}", linewidth=1.8)
    axes[0].set_title("A) Decaimiento Espectral Teórico $\\beta(X, Q)$", fontweight='bold')
    axes[0].set_xlabel("X (frecuencia × tiempo)")
    axes[0].set_ylabel("Amplitud $\\beta$")
    axes[0].legend(loc="upper right")
    axes[0].grid(True, linestyle=":", alpha=0.6)
    
    for q_t in q_targets:
        idx_q = int(np.argmin(np.abs(q_vec - q_t)))
        axes[1].plot(x_vec, gan_teorica[:, idx_q], label=f"Q = {q_vec[idx_q]:.0f}", linewidth=1.8)
    axes[1].set_title("B) Curvas de Ganancia Estabilizada de Wang ($G_{lim}=20$ dB)", fontweight='bold')
    axes[1].set_xlabel("X (frecuencia × tiempo)")
    axes[1].set_ylabel("Ganancia Teórica")
    axes[1].legend(loc="upper left")
    axes[1].grid(True, linestyle=":", alpha=0.6)
    
    fig.suptitle("Paso 4: Modelado Teórico de Atenuación y Catálogo de Ganancias de Wang", fontsize=14, y=1.02)
    out_path_4 = os.path.join(ASSETS_DIR, "04_q_theoretical_curves.png")
    plt.savefig(out_path_4, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Figura 4 guardada en: {out_path_4}")
    
    # Inversión
    inv_res = estimate_q_from_tfdist(tfdist, t_vec, f_eval, q_cfg)
    q_est = inv_res["q_estimated"]
    
    # Figura 5: Inversión sobre datos
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    
    # A: Curva Chi
    axes[0].plot(x_vec, inv_res["amp_x_norm"], color='navy', linewidth=1.8)
    axes[0].set_title("A) Amplitud Integrada en Isolíneas $\\chi = f \\cdot t$", fontweight='bold')
    axes[0].set_xlabel("X (frecuencia × tiempo)")
    axes[0].set_ylabel("Amplitud Normalizada")
    axes[0].grid(True, linestyle=":", alpha=0.6)
    
    # B: Ajuste de ganancias
    axes[1].plot(x_vec, inv_res["gan_datos"], 'r-', label="Ganancia Datos", linewidth=2.0)
    axes[1].plot(x_vec, inv_res["gan_teorica_opt"], 'b--', label=f"Teórica Óptima (Q={q_est:.1f})", linewidth=1.8)
    axes[1].set_title("B) Ajuste de Ganancia Observada vs Teórica", fontweight='bold')
    axes[1].set_xlabel("X (frecuencia × tiempo)")
    axes[1].set_ylabel("Ganancia Estabilizada")
    axes[1].legend()
    axes[1].grid(True, linestyle=":", alpha=0.6)
    
    # C: Minimización residual
    axes[2].plot(q_vec, inv_res["residuals"], color='darkgreen', linewidth=2.0)
    axes[2].axvline(q_est, color='red', linestyle=':', label=f"Mínimo Global: Q = {q_est:.1f}")
    axes[2].set_title("C) Función Objetivo Residual vs Factor Q", fontweight='bold')
    axes[2].set_xlabel("Factor de Calidad (Q)")
    axes[2].set_ylabel("Error Residual Acumulado")
    axes[2].legend()
    axes[2].grid(True, linestyle=":", alpha=0.6)
    
    fig.suptitle("Paso 5: Inversión y Estimación del Factor de Calidad Sísmico (Q)", fontsize=14, y=1.02)
    out_path_5 = os.path.join(ASSETS_DIR, "05_q_inversion_workflow.png")
    plt.savefig(out_path_5, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Figura 5 guardada en: {out_path_5}")

if __name__ == "__main__":
    t_vec, trace = generate_figure_1_seismic_qc()
    generate_figure_2_morlet_dictionary()
    tfdist, f_eval = generate_figure_3_mpd_and_tfd(t_vec, trace)
    generate_figure_4_and_5_q_workflow(t_vec, tfdist, f_eval)
    print("Todas las figuras generadas exitosamente en assets/images/")

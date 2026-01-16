#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced Parameter Space Analysis for PNAS-style Publication
重点：多参数空间的系统分析和微生物生态学意义

Author: Jian Wang
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import seaborn as sns

# PNAS figure settings
plt.rcParams['font.size'] = 8
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

COLORS = {
    'S': '#3498db',
    'M': '#e74c3c',
    'G': '#2ecc71',
    'coexist': '#f39c12'
}


def compute_invasion_landscape(omega_range, sigma_MS_range, params):
    """
    Compute G invasion rate across (ω, σ_MS) parameter space

    Returns invasion rate matrix and critical boundary
    """
    n_omega = len(omega_range)
    n_sigma = len(sigma_MS_range)

    invasion_rate = np.zeros((n_sigma, n_omega))
    s_SM_matrix = np.zeros((n_sigma, n_omega))
    m_SM_matrix = np.zeros((n_sigma, n_omega))

    sigma_SM = params['sigma_SM']
    sigma_GS = params['sigma_GS']
    sigma_GM = params['sigma_GM']
    alpha_GS = params['alpha_GS']
    alpha_GM = params['alpha_GM']

    for i, sigma_MS in enumerate(sigma_MS_range):
        for j, omega in enumerate(omega_range):
            # S-M equilibrium
            denom = 1 - sigma_MS * sigma_SM

            if abs(denom) > 1e-10 and sigma_MS > 1:
                s_SM = (1 - sigma_SM) / denom
                m_SM = (sigma_MS - 1) / denom

                # Net interaction parameters
                c = (1 - omega) * sigma_GS - omega * alpha_GS
                e = omega * sigma_GM - (1 - omega) * alpha_GM
                d = 2 * omega - 1

                # Invasion rate
                lambda_G = d + c * s_SM + e * m_SM

                invasion_rate[i, j] = lambda_G
                s_SM_matrix[i, j] = s_SM
                m_SM_matrix[i, j] = m_SM
            else:
                invasion_rate[i, j] = np.nan
                s_SM_matrix[i, j] = np.nan
                m_SM_matrix[i, j] = np.nan

    return invasion_rate, s_SM_matrix, m_SM_matrix


def compute_coexistence_region(omega_range, sigma_SM_range, sigma_MS_fixed, params):
    """
    Compute coexistence region in (ω, σ_SM) space with fixed σ_MS
    """
    n_omega = len(omega_range)
    n_sigma_SM = len(sigma_SM_range)

    invasion_rate = np.zeros((n_sigma_SM, n_omega))
    stability_SM = np.zeros((n_sigma_SM, n_omega), dtype=bool)

    for i, sigma_SM in enumerate(sigma_SM_range):
        for j, omega in enumerate(omega_range):
            # S-M equilibrium
            denom = 1 - sigma_MS_fixed * sigma_SM

            # Check S-M stability condition
            if abs(denom) > 1e-10 and sigma_MS_fixed > 1 and sigma_MS_fixed * sigma_SM < 1:
                stability_SM[i, j] = True

                s_SM = (1 - sigma_SM) / denom
                m_SM = (sigma_MS_fixed - 1) / denom

                # G invasion rate
                c = (1 - omega) * params['sigma_GS'] - omega * params['alpha_GS']
                e = omega * params['sigma_GM'] - (1 - omega) * params['alpha_GM']
                d = 2 * omega - 1

                invasion_rate[i, j] = d + c * s_SM + e * m_SM
            else:
                stability_SM[i, j] = False
                invasion_rate[i, j] = np.nan

    return invasion_rate, stability_SM


def plot_enhanced_parameter_analysis():
    """
    Generate comprehensive parameter space analysis (PNAS Figure 2)
    """
    fig = plt.figure(figsize=(7.5, 9))
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.3,
                  left=0.1, right=0.95, top=0.95, bottom=0.06)

    # Fixed parameters
    params = {
        'sigma_SM': 0.5,
        'sigma_GS': 0.4,
        'sigma_GM': 0.4,
        'alpha_GS': 0.3,
        'alpha_GM': 0.3
    }

    # ========================================================================
    # Panel A: Invasion rate in (ω, σ_MS) space - HIGH RESOLUTION
    # ========================================================================
    ax1 = fig.add_subplot(gs[0, :])

    omega_range = np.linspace(0, 1, 200)
    sigma_MS_range = np.linspace(1.0, 2.5, 200)

    invasion_rate, s_SM_matrix, m_SM_matrix = compute_invasion_landscape(
        omega_range, sigma_MS_range, params
    )

    # Create beautiful contour plot
    Omega, Sigma_MS = np.meshgrid(omega_range, sigma_MS_range)

    # Filled contours with better colormap
    levels_fill = np.linspace(-1.0, 1.5, 30)
    cf = ax1.contourf(Omega, Sigma_MS, invasion_rate,
                     levels=levels_fill, cmap='RdBu_r',
                     extend='both', alpha=0.9)

    # Critical boundary (invasion rate = 0)
    ct_zero = ax1.contour(Omega, Sigma_MS, invasion_rate,
                         levels=[0], colors='black',
                         linewidths=2.5, linestyles='-')

    # Additional contour lines
    levels_lines = [-0.5, -0.25, 0.25, 0.5, 0.75, 1.0]
    ct_lines = ax1.contour(Omega, Sigma_MS, invasion_rate,
                          levels=levels_lines, colors='gray',
                          linewidths=0.8, alpha=0.6)
    ax1.clabel(ct_lines, inline=True, fontsize=6, fmt='%0.2f')

    # Colorbar
    cbar = plt.colorbar(cf, ax=ax1, label='G invasion rate ($\\lambda_G$)')
    cbar.ax.tick_params(labelsize=7)

    # Add labels to regions
    ax1.text(0.25, 2.2, 'G excluded', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7),
            ha='center', weight='bold')
    ax1.text(0.75, 1.3, 'G invades', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7),
            ha='center', weight='bold')

    # Mark critical line
    ax1.plot([], [], 'k-', linewidth=2.5, label='Invasion threshold ($\\lambda_G = 0$)')

    ax1.set_xlabel('Generalist pathway parameter ($\\omega$)', fontweight='bold')
    ax1.set_ylabel('S$\\to$M benefit ($\\sigma_{MS}$)', fontweight='bold')
    ax1.set_title('A', loc='left', fontweight='bold', fontsize=12)
    ax1.legend(loc='upper left', framealpha=0.9, fontsize=7)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(1.0, 2.5)

    # ========================================================================
    # Panel B: S-M equilibrium densities across parameter space
    # ========================================================================
    ax2 = fig.add_subplot(gs[1, 0])

    # Show how S-M equilibrium changes
    cf2 = ax2.contourf(Omega, Sigma_MS, s_SM_matrix,
                      levels=20, cmap='Blues', alpha=0.8)
    ct2 = ax2.contour(Omega, Sigma_MS, invasion_rate, levels=[0],
                     colors='red', linewidths=2, linestyles='--')

    cbar2 = plt.colorbar(cf2, ax=ax2, label='$s^*_{SM}$')
    cbar2.ax.tick_params(labelsize=7)

    ax2.set_xlabel('$\\omega$', fontweight='bold')
    ax2.set_ylabel('$\\sigma_{MS}$', fontweight='bold')
    ax2.set_title('B', loc='left', fontweight='bold', fontsize=12)
    ax2.set_xlim(0, 1)
    ax2.set_ylim(1.0, 2.5)

    # ========================================================================
    # Panel C: M equilibrium density
    # ========================================================================
    ax3 = fig.add_subplot(gs[1, 1])

    cf3 = ax3.contourf(Omega, Sigma_MS, m_SM_matrix,
                      levels=20, cmap='Oranges', alpha=0.8)
    ct3 = ax3.contour(Omega, Sigma_MS, invasion_rate, levels=[0],
                     colors='red', linewidths=2, linestyles='--')

    cbar3 = plt.colorbar(cf3, ax=ax3, label='$m^*_{SM}$')
    cbar3.ax.tick_params(labelsize=7)

    ax3.set_xlabel('$\\omega$', fontweight='bold')
    ax3.set_ylabel('$\\sigma_{MS}$', fontweight='bold')
    ax3.set_title('C', loc='left', fontweight='bold', fontsize=12)
    ax3.set_xlim(0, 1)
    ax3.set_ylim(1.0, 2.5)

    # ========================================================================
    # Panel D: Alternative parameter space (ω, σ_SM)
    # ========================================================================
    ax4 = fig.add_subplot(gs[2, 0])

    sigma_SM_range = np.linspace(0.1, 1.2, 150)
    sigma_MS_fixed = 1.5

    invasion_rate_2, stability_SM_2 = compute_coexistence_region(
        omega_range, sigma_SM_range, sigma_MS_fixed, params
    )

    Omega2, Sigma_SM = np.meshgrid(omega_range, sigma_SM_range)

    cf4 = ax4.contourf(Omega2, Sigma_SM, invasion_rate_2,
                      levels=levels_fill, cmap='RdBu_r',
                      extend='both', alpha=0.9)
    ct4_zero = ax4.contour(Omega2, Sigma_SM, invasion_rate_2,
                          levels=[0], colors='black',
                          linewidths=2.5)

    # S-M stability boundary
    sigma_SM_crit = 1 / sigma_MS_fixed
    ax4.axhline(sigma_SM_crit, color='orange', linestyle='--',
               linewidth=2, label=f'$\\sigma_{{SM}} = 1/\\sigma_{{MS}}$ = {sigma_SM_crit:.2f}')

    # Shade unstable S-M region
    ax4.fill_between([0, 1], [sigma_SM_crit, sigma_SM_crit], [1.2, 1.2],
                    alpha=0.2, color='gray', label='S-M unstable')

    cbar4 = plt.colorbar(cf4, ax=ax4, label='$\\lambda_G$')
    cbar4.ax.tick_params(labelsize=7)

    ax4.set_xlabel('$\\omega$', fontweight='bold')
    ax4.set_ylabel('M$\\to$S benefit ($\\sigma_{SM}$)', fontweight='bold')
    ax4.set_title('D', loc='left', fontweight='bold', fontsize=12)
    ax4.legend(loc='upper right', framealpha=0.9, fontsize=7)
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0.1, 1.2)

    # ========================================================================
    # Panel E: Cooperation vs Competition balance
    # ========================================================================
    ax5 = fig.add_subplot(gs[2, 1])

    # Vary σ_GS (cooperation) vs α_GS (competition)
    sigma_GS_range = np.linspace(0.1, 0.8, 100)
    alpha_GS_range = np.linspace(0.1, 0.8, 100)

    # Fixed omega
    omega_fixed = 0.5
    sigma_MS_fixed2 = 1.5
    sigma_SM_fixed = 0.5

    # S-M equilibrium
    denom = 1 - sigma_MS_fixed2 * sigma_SM_fixed
    s_SM_val = (1 - sigma_SM_fixed) / denom
    m_SM_val = (sigma_MS_fixed2 - 1) / denom

    InvasionRate5 = np.zeros((len(alpha_GS_range), len(sigma_GS_range)))

    for i, alpha_GS in enumerate(alpha_GS_range):
        for j, sigma_GS in enumerate(sigma_GS_range):
            c = (1 - omega_fixed) * sigma_GS - omega_fixed * alpha_GS
            e = omega_fixed * params['sigma_GM'] - (1 - omega_fixed) * params['alpha_GM']
            d = 2 * omega_fixed - 1

            InvasionRate5[i, j] = d + c * s_SM_val + e * m_SM_val

    Sigma_GS, Alpha_GS = np.meshgrid(sigma_GS_range, alpha_GS_range)

    cf5 = ax5.contourf(Sigma_GS, Alpha_GS, InvasionRate5,
                      levels=20, cmap='RdBu_r', alpha=0.9)
    ct5 = ax5.contour(Sigma_GS, Alpha_GS, InvasionRate5,
                     levels=[0], colors='black', linewidths=2.5)

    # Diagonal line (cooperation = competition)
    ax5.plot([0.1, 0.8], [0.1, 0.8], 'k--', linewidth=1.5, alpha=0.5,
            label='$\\sigma_{GS} = \\alpha_{GS}$')

    cbar5 = plt.colorbar(cf5, ax=ax5, label='$\\lambda_G$')
    cbar5.ax.tick_params(labelsize=7)

    ax5.set_xlabel('S$\\to$G cooperation ($\\sigma_{GS}$)', fontweight='bold')
    ax5.set_ylabel('S$\\to$G competition ($\\alpha_{GS}$)', fontweight='bold')
    ax5.set_title('E', loc='left', fontweight='bold', fontsize=12)
    ax5.legend(loc='upper left', framealpha=0.9, fontsize=7)

    plt.savefig('figures/Figure2_parameter_space_analysis.png',
                dpi=300, bbox_inches='tight')
    print("✓ Enhanced Figure 2 saved")

    return fig


def plot_pairwise_publication():
    """Generate publication-quality pairwise systems figure"""

    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.5))

    # Parameters
    sigma_SM, sigma_MS = 0.5, 1.5

    # Panel A: S-M
    ax = axes[0]

    s_range = np.linspace(0, 2.5, 300)
    m_range = np.linspace(0, 2.5, 300)

    s_null = 1 + sigma_SM * m_range
    m_null = sigma_MS * s_range - 1

    ax.plot(s_null, m_range, color=COLORS['S'], linewidth=2,
            label='$ds/dt=0$', alpha=0.8)
    ax.plot(s_range, m_null, color=COLORS['M'], linewidth=2,
            label='$dm/dt=0$', alpha=0.8)

    # Equilibrium
    denom = 1 - sigma_MS * sigma_SM
    s_eq = (1 - sigma_SM) / denom
    m_eq = (sigma_MS - 1) / denom

    ax.plot(s_eq, m_eq, 'o', color='black', markersize=8,
           markeredgewidth=2, markerfacecolor=COLORS['coexist'],
           zorder=10, label='Stable equilibrium')

    ax.set_xlabel('$s$ (S density)', fontweight='bold')
    ax.set_ylabel('$m$ (M density)', fontweight='bold')
    ax.set_title('A  S-M mutualism', loc='left', fontweight='bold')
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 2.5)
    ax.legend(fontsize=6, framealpha=0.9)
    ax.grid(True, alpha=0.2)

    # Panel B: S-G
    ax = axes[1]

    a, c, d = 0.2, 0.3, 0.0

    s_null_sg = 1 + a * m_range
    g_null_sg = d + c * s_range

    ax.plot(s_null_sg, m_range, color=COLORS['S'], linewidth=2,
            label='$ds/dt=0$', alpha=0.8)
    ax.plot(s_range, g_null_sg, color=COLORS['G'], linewidth=2,
            label='$dg/dt=0$', alpha=0.8)

    # Equilibrium
    g_eq = (d + c) / (1 - a*c)
    s_eq = (1 + a*d) / (1 - a*c)

    if g_eq > 0 and s_eq > 0:
        ax.plot(s_eq, g_eq, 'o', color='black', markersize=8,
               markeredgewidth=2, markerfacecolor=COLORS['coexist'],
               zorder=10, label='Stable equilibrium')

    ax.set_xlabel('$s$ (S density)', fontweight='bold')
    ax.set_ylabel('$g$ (G density)', fontweight='bold')
    ax.set_title('B  S-G interaction', loc='left', fontweight='bold')
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 2.5)
    ax.legend(fontsize=6, framealpha=0.9)
    ax.grid(True, alpha=0.2)

    # Panel C: M-G
    ax = axes[2]

    b, e, d = 0.3, 0.2, 0.0

    m_null_mg = b * m_range - 1
    g_null_mg = d + e * s_range

    ax.plot(m_null_mg, m_range, color=COLORS['M'], linewidth=2,
            label='$dm/dt=0$', alpha=0.8)
    ax.plot(s_range, g_null_mg, color=COLORS['G'], linewidth=2,
            label='$dg/dt=0$', alpha=0.8)

    ax.set_xlabel('$m$ (M density)', fontweight='bold')
    ax.set_ylabel('$g$ (G density)', fontweight='bold')
    ax.set_title('C  M-G interaction', loc='left', fontweight='bold')
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 2.5)
    ax.legend(fontsize=6, framealpha=0.9)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig('figures/Figure1_pairwise_PNAS.png', dpi=300, bbox_inches='tight')
    print("✓ Figure 1 (PNAS style) saved")

    return fig


def plot_bifurcation_PNAS():
    """PNAS-style bifurcation diagram"""

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 6))

    # Parameters
    sigma_SM, sigma_MS = 0.5, 1.5
    sigma_GS, sigma_GM = 0.4, 0.4
    alpha_GS, alpha_GM = 0.3, 0.3

    omega_range = np.linspace(0, 1, 150)

    # ... (bifurcation computation code similar to before)

    plt.tight_layout()
    plt.savefig('figures/Figure3_bifurcation_PNAS.png', dpi=300, bbox_inches='tight')
    print("✓ Figure 3 (PNAS style) saved")

    return fig


if __name__ == "__main__":
    print("\nGenerating PNAS-style figures...")
    print("="*60)

    plot_pairwise_publication()
    plot_enhanced_parameter_analysis()

    print("="*60)
    print("✅ All figures generated\n")

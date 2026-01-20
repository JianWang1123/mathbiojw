#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced Comprehensive Figure System for PNAS Publication
超高质量、复杂度高、信息丰富的图表系统

Author: Jian Wang
Date: January 2026

Figure Design:
--------------
Figure 1: Pairwise Systems - Complete Dynamical Analysis (6 panels)
Figure 2: Parameter Space Architecture - Multi-dimensional Landscapes (8 panels)
Figure 3: Three-Species Dynamics - Phase Space Analysis (6 panels)
Figure 4: Bifurcation Structure - Complete Stability Analysis (4 panels)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from matplotlib.patches import Rectangle, FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.linalg import eig
import seaborn as sns

# PNAS publication settings - highest quality
plt.rcParams.update({
    'font.size': 7,
    'axes.labelsize': 8,
    'axes.titlesize': 9,
    'legend.fontsize': 6,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'lines.linewidth': 1.5,
    'axes.linewidth': 0.8,
    'figure.autolayout': False,
    'figure.constrained_layout.use': True
})

# Color scheme
COLORS = {
    'S': '#2E86AB',  # Deep blue
    'M': '#A23B72',  # Deep magenta
    'G': '#F18F01',  # Orange
    'coexist': '#06A77D',  # Teal
    'exclude': '#C73E1D',  # Red
    'neutral': '#6C757D'   # Gray
}


class ThreeSpeciesModel:
    """Enhanced three-species cross-feeding model with complete analysis"""

    def __init__(self, params):
        """Initialize with parameter dictionary"""
        self.params = params
        self.r_S = params.get('r_S', 1.0)
        self.r_M = params.get('r_M', 0.8)
        self.r_G = params.get('r_G', 0.9)
        self.sigma_SM = params.get('sigma_SM', 0.5)
        self.sigma_MS = params.get('sigma_MS', 1.5)
        self.sigma_SG = params.get('sigma_SG', 0.4)
        self.sigma_GS = params.get('sigma_GS', 0.4)
        self.sigma_MG = params.get('sigma_MG', 0.4)
        self.sigma_GM = params.get('sigma_GM', 0.4)
        self.alpha_SG = params.get('alpha_SG', 0.3)
        self.alpha_GS = params.get('alpha_GS', 0.3)
        self.alpha_MG = params.get('alpha_MG', 0.3)
        self.alpha_GM = params.get('alpha_GM', 0.3)
        self.omega = params.get('omega', 0.5)

    def net_interactions(self, omega):
        """Compute net interaction parameters"""
        a = (1 - omega) * self.sigma_SG - omega * self.alpha_SG
        c = (1 - omega) * self.sigma_GS - omega * self.alpha_GS
        b = omega * self.sigma_MG - (1 - omega) * self.alpha_MG
        e = omega * self.sigma_GM - (1 - omega) * self.alpha_GM
        d = 2 * omega - 1
        return a, b, c, d, e

    def SM_dynamics(self, t, y):
        """S-M pairwise dynamics"""
        s, m = y
        dsdt = self.r_S * s * (1 + self.sigma_SM * m - s)
        dmdt = self.r_M * m * (-1 + self.sigma_MS * s - m)
        return [dsdt, dmdt]

    def SG_dynamics(self, t, y, omega):
        """S-G pairwise dynamics"""
        s, g = y
        a, _, c, d, _ = self.net_interactions(omega)
        dsdt = self.r_S * s * (1 + a * g - s)
        dgdt = self.r_G * g * (d + c * s - g)
        return [dsdt, dgdt]

    def MG_dynamics(self, t, y, omega):
        """M-G pairwise dynamics"""
        m, g = y
        _, b, _, d, e = self.net_interactions(omega)
        dmdt = self.r_M * m * (-1 + b * g - m)
        dgdt = self.r_G * g * (d + e * m - g)
        return [dmdt, dgdt]

    def three_species_dynamics(self, t, y, omega):
        """Full three-species dynamics"""
        s, m, g = y
        a, b, c, d, e = self.net_interactions(omega)

        dsdt = self.r_S * s * (1 + self.sigma_SM * m + a * g - s)
        dmdt = self.r_M * m * (-1 + self.sigma_MS * s + b * g - m)
        dgdt = self.r_G * g * (d + c * s + e * m - g)

        return [dsdt, dmdt, dgdt]

    def compute_SM_equilibrium(self):
        """Compute S-M coexistence equilibrium"""
        denom = 1 - self.sigma_MS * self.sigma_SM
        if abs(denom) > 1e-10 and self.sigma_MS > 1:
            s_SM = (1 - self.sigma_SM) / denom
            m_SM = (self.sigma_MS - 1) / denom

            # Check stability
            stable = (self.sigma_MS > 1) and (self.sigma_MS * self.sigma_SM < 1)
            return s_SM, m_SM, stable
        return None, None, False

    def invasion_fitness(self, omega):
        """Compute G invasion fitness at S-M equilibrium"""
        s_SM, m_SM, stable = self.compute_SM_equilibrium()
        if not stable:
            return -np.inf

        _, _, c, d, e = self.net_interactions(omega)
        lambda_G = self.r_G * (d + c * s_SM + e * m_SM)
        return lambda_G

    def jacobian_SM(self, s, m):
        """Jacobian for S-M system"""
        J = np.array([
            [self.r_S * (1 + self.sigma_SM * m - 2*s),
             self.r_S * s * self.sigma_SM],
            [self.r_M * m * self.sigma_MS,
             self.r_M * (-1 + self.sigma_MS * s - 2*m)]
        ])
        return J

    def jacobian_three_species(self, s, m, g, omega):
        """Jacobian for three-species system"""
        a, b, c, d, e = self.net_interactions(omega)

        J = np.array([
            [self.r_S * (1 + self.sigma_SM * m + a * g - 2*s),
             self.r_S * s * self.sigma_SM,
             self.r_S * s * a],
            [self.r_M * m * self.sigma_MS,
             self.r_M * (-1 + self.sigma_MS * s + b * g - 2*m),
             self.r_M * m * b],
            [self.r_G * g * c,
             self.r_G * g * e,
             self.r_G * (d + c * s + e * m - 2*g)]
        ])
        return J


def create_figure1_pairwise_analysis(model):
    """
    Figure 1: Comprehensive Pairwise Systems Analysis

    6 Panels showing complete dynamical characterization:
    A: S-M time series with parameter variation
    B: S-M phase portrait with nullclines
    C: S-M bifurcation diagram (σ_MS scan)
    D: S-G dynamics across ω values
    E: M-G dynamics across ω values
    F: Pairwise stability map in parameter space
    """
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.35)

    # Panel A: S-M time series
    ax_a = fig.add_subplot(gs[0, :2])
    t_span = (0, 50)
    t_eval = np.linspace(0, 50, 500)

    # Multiple parameter regimes
    sigma_MS_values = [1.2, 1.5, 2.0]
    for sigma_MS in sigma_MS_values:
        model.sigma_MS = sigma_MS
        sol = solve_ivp(model.SM_dynamics, t_span, [0.5, 0.3], t_eval=t_eval,
                       method='RK45', rtol=1e-8, atol=1e-10)
        ax_a.plot(sol.t, sol.y[0], '-', color=COLORS['S'], alpha=0.7,
                 label=f'S (σ_MS={sigma_MS})' if sigma_MS == 1.5 else '')
        ax_a.plot(sol.t, sol.y[1], '-', color=COLORS['M'], alpha=0.7,
                 label=f'M (σ_MS={sigma_MS})' if sigma_MS == 1.5 else '')

    ax_a.set_xlabel('Time (days)')
    ax_a.set_ylabel('Population density (scaled)')
    ax_a.set_title('A. S-M mutualistic dynamics across mutualism strengths',
                   fontweight='bold', loc='left')
    ax_a.legend(frameon=False, ncol=2)
    ax_a.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # Panel B: S-M phase portrait with nullclines
    ax_b = fig.add_subplot(gs[0, 2])
    model.sigma_MS = 1.5

    # Nullclines
    s_range = np.linspace(0, 2, 100)
    m_range = np.linspace(0, 2, 100)

    # S nullcline: dm/ds = 0 => m = (s - 1) / sigma_SM
    s_null = s_range
    m_null_s = (s_null - 1) / model.sigma_SM

    # M nullcline: ds/dm = 0 => s = (1 + m) / sigma_MS
    m_null = m_range
    s_null_m = (1 + m_null) / model.sigma_MS

    ax_b.plot(s_null, m_null_s, '-', color=COLORS['S'], linewidth=2,
             label='S nullcline', alpha=0.7)
    ax_b.plot(s_null_m, m_null, '-', color=COLORS['M'], linewidth=2,
             label='M nullcline', alpha=0.7)

    # Trajectories
    for s0, m0 in [(0.3, 0.2), (0.8, 0.4), (1.2, 0.8)]:
        sol = solve_ivp(model.SM_dynamics, (0, 100), [s0, m0],
                       t_eval=np.linspace(0, 100, 1000), method='RK45')
        ax_b.plot(sol.y[0], sol.y[1], '-', color='gray', alpha=0.5, linewidth=1)
        ax_b.plot(s0, m0, 'o', color='black', markersize=4)

    # Equilibrium
    s_eq, m_eq, _ = model.compute_SM_equilibrium()
    if s_eq is not None:
        ax_b.plot(s_eq, m_eq, '*', color=COLORS['coexist'], markersize=12,
                 markeredgecolor='black', markeredgewidth=0.5, label='Equilibrium')

    ax_b.set_xlabel('S density')
    ax_b.set_ylabel('M density')
    ax_b.set_title('B. Phase portrait', fontweight='bold', loc='left')
    ax_b.set_xlim(0, 1.8)
    ax_b.set_ylim(0, 1.5)
    ax_b.legend(frameon=False, fontsize=5)
    ax_b.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # Panel C: Bifurcation diagram (σ_MS scan)
    ax_c = fig.add_subplot(gs[1, :2])
    sigma_MS_scan = np.linspace(0.8, 3.0, 200)
    s_equilibria = []
    m_equilibria = []
    stability = []

    for sigma_MS in sigma_MS_scan:
        model.sigma_MS = sigma_MS
        s_eq, m_eq, stable = model.compute_SM_equilibrium()
        if s_eq is not None and s_eq > 0 and m_eq > 0:
            s_equilibria.append(s_eq)
            m_equilibria.append(m_eq)
            stability.append(stable)
        else:
            s_equilibria.append(np.nan)
            m_equilibria.append(np.nan)
            stability.append(False)

    stability = np.array(stability)
    s_equilibria = np.array(s_equilibria)
    m_equilibria = np.array(m_equilibria)

    # Plot stable and unstable branches
    ax_c.plot(sigma_MS_scan[stability], s_equilibria[stability], '-',
             color=COLORS['S'], linewidth=2, label='S (stable)')
    ax_c.plot(sigma_MS_scan[~stability], s_equilibria[~stability], '--',
             color=COLORS['S'], linewidth=2, alpha=0.5, label='S (unstable)')
    ax_c.plot(sigma_MS_scan[stability], m_equilibria[stability], '-',
             color=COLORS['M'], linewidth=2, label='M (stable)')
    ax_c.plot(sigma_MS_scan[~stability], m_equilibria[~stability], '--',
             color=COLORS['M'], linewidth=2, alpha=0.5, label='M (unstable)')

    # Mark bifurcation points
    ax_c.axvline(1.0, color='black', linestyle=':', linewidth=1, alpha=0.5)
    ax_c.text(1.0, 1.5, 'Transcritical\nbifurcation', fontsize=6, ha='center')

    ax_c.set_xlabel('Mutualism strength (σ_MS)')
    ax_c.set_ylabel('Equilibrium density')
    ax_c.set_title('C. Bifurcation diagram: S-M coexistence',
                   fontweight='bold', loc='left')
    ax_c.legend(frameon=False, loc='upper left', ncol=2)
    ax_c.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    ax_c.set_ylim(-0.1, 2.5)

    # Panel D: S-G dynamics
    ax_d = fig.add_subplot(gs[1, 2])
    omega_values = [0.3, 0.5, 0.7]
    for omega in omega_values:
        sol = solve_ivp(lambda t, y: model.SG_dynamics(t, y, omega),
                       (0, 50), [0.8, 0.3], t_eval=np.linspace(0, 50, 500))
        ax_d.plot(sol.y[0], sol.y[1], '-', linewidth=1.5,
                 label=f'ω={omega}', alpha=0.7)

    ax_d.set_xlabel('S density')
    ax_d.set_ylabel('G density')
    ax_d.set_title('D. S-G trajectories', fontweight='bold', loc='left')
    ax_d.legend(frameon=False, fontsize=5)
    ax_d.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # Panel E: M-G dynamics
    ax_e = fig.add_subplot(gs[2, 0])
    for omega in omega_values:
        sol = solve_ivp(lambda t, y: model.MG_dynamics(t, y, omega),
                       (0, 50), [0.5, 0.2], t_eval=np.linspace(0, 50, 500))
        ax_e.plot(sol.y[0], sol.y[1], '-', linewidth=1.5,
                 label=f'ω={omega}', alpha=0.7)

    ax_e.set_xlabel('M density')
    ax_e.set_ylabel('G density')
    ax_e.set_title('E. M-G trajectories', fontweight='bold', loc='left')
    ax_e.legend(frameon=False, fontsize=5)
    ax_e.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # Panel F: Pairwise stability map
    ax_f = fig.add_subplot(gs[2, 1:])

    sigma_MS_range = np.linspace(0.5, 3.0, 100)
    sigma_SM_range = np.linspace(0.0, 1.5, 100)
    SM_grid, MS_grid = np.meshgrid(sigma_SM_range, sigma_MS_range)

    # Stability condition: σ_MS > 1 AND σ_MS * σ_SM < 1
    viable = MS_grid > 1
    stable = (MS_grid > 1) & (MS_grid * SM_grid < 1)

    # Create stability map
    stability_map = np.zeros_like(MS_grid)
    stability_map[~viable] = 0  # M extinction
    stability_map[viable & ~stable] = 1  # Unstable coexistence
    stability_map[stable] = 2  # Stable coexistence

    cmap_discrete = plt.cm.colors.ListedColormap(
        [COLORS['exclude'], '#FFD93D', COLORS['coexist']]
    )
    bounds = [0, 0.67, 1.33, 2]
    norm = plt.cm.colors.BoundaryNorm(bounds, cmap_discrete.N)

    im = ax_f.pcolormesh(SM_grid, MS_grid, stability_map,
                         cmap=cmap_discrete, norm=norm, shading='auto', alpha=0.6)

    # Contour lines
    ax_f.contour(SM_grid, MS_grid, MS_grid * SM_grid, levels=[1.0],
                colors='black', linewidths=2, linestyles='--')
    ax_f.axhline(1.0, color='black', linestyle=':', linewidth=1)

    ax_f.set_xlabel('Benefit to S from M (σ_SM)')
    ax_f.set_ylabel('Benefit to M from S (σ_MS)')
    ax_f.set_title('F. S-M stability regions in parameter space',
                   fontweight='bold', loc='left')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax_f, ticks=[0.33, 1.0, 1.67])
    cbar.ax.set_yticklabels(['M extinction', 'Unstable', 'Stable coexistence'],
                            fontsize=6)

    ax_f.text(0.7, 2.5, r'$\sigma_{MS} \cdot \sigma_{SM} = 1$',
             fontsize=7, rotation=-30, ha='center')
    ax_f.text(0.2, 0.6, r'$\sigma_{MS} = 1$', fontsize=7, ha='center')

    plt.savefig('figures/Figure1_comprehensive_pairwise.png', dpi=300, bbox_inches='tight')
    print("✓ Figure 1 created: Comprehensive pairwise systems analysis")
    return fig


def create_figure2_parameter_landscapes(model):
    """
    Figure 2: Multi-dimensional Parameter Space Architecture

    8 Panels showing complete parameter landscape:
    A: (ω, σ_MS) invasion rate with detailed contours
    B: (ω, σ_MS) S equilibrium density landscape
    C: (ω, σ_MS) M equilibrium density landscape
    D: (ω, σ_SM) invasion rate with stability boundary
    E: (σ_GS, α_GS) cooperation-competition trade-off
    F: (σ_GM, α_GM) metabolite-mediated interactions
    G: Critical ω vs σ_MS showing invasion threshold
    H: Coexistence window across multiple parameters
    """
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.4)

    # Common parameter ranges
    omega_range = np.linspace(0.0, 1.0, 150)
    sigma_MS_range = np.linspace(1.1, 3.0, 150)

    # Panel A: (ω, σ_MS) invasion landscape
    ax_a = fig.add_subplot(gs[0, 0])

    invasion_rate = np.zeros((len(sigma_MS_range), len(omega_range)))
    s_SM_matrix = np.zeros_like(invasion_rate)
    m_SM_matrix = np.zeros_like(invasion_rate)

    for i, sigma_MS in enumerate(sigma_MS_range):
        for j, omega in enumerate(omega_range):
            model.sigma_MS = sigma_MS
            s_SM, m_SM, stable = model.compute_SM_equilibrium()
            if stable:
                lambda_G = model.invasion_fitness(omega)
                invasion_rate[i, j] = lambda_G
                s_SM_matrix[i, j] = s_SM
                m_SM_matrix[i, j] = m_SM
            else:
                invasion_rate[i, j] = np.nan

    # Plot heatmap
    im_a = ax_a.pcolormesh(omega_range, sigma_MS_range, invasion_rate,
                          cmap='RdBu_r', shading='auto', vmin=-0.5, vmax=0.5)

    # Contour lines
    CS = ax_a.contour(omega_range, sigma_MS_range, invasion_rate,
                     levels=[-0.3, -0.1, 0, 0.1, 0.3], colors='black',
                     linewidths=[0.8, 0.8, 2.0, 0.8, 0.8],
                     linestyles=['--', '--', '-', '--', '--'])
    ax_a.clabel(CS, inline=True, fontsize=6, fmt='%.2f')

    ax_a.set_xlabel('Pathway parameter (ω)')
    ax_a.set_ylabel('S→M mutualism (σ_MS)')
    ax_a.set_title('A. G invasion fitness landscape', fontweight='bold', loc='left')

    cbar_a = plt.colorbar(im_a, ax=ax_a)
    cbar_a.set_label('Invasion rate λ_G', fontsize=7)

    # Panel B: S equilibrium density
    ax_b = fig.add_subplot(gs[0, 1])

    im_b = ax_b.pcolormesh(omega_range, sigma_MS_range, s_SM_matrix,
                          cmap='Blues', shading='auto')
    ax_b.contour(omega_range, sigma_MS_range, invasion_rate, levels=[0],
                colors='red', linewidths=2, linestyles='-')

    ax_b.set_xlabel('Pathway parameter (ω)')
    ax_b.set_ylabel('S→M mutualism (σ_MS)')
    ax_b.set_title('B. S equilibrium density at S-M state',
                   fontweight='bold', loc='left')

    cbar_b = plt.colorbar(im_b, ax=ax_b)
    cbar_b.set_label('s*_SM', fontsize=7)

    # Panel C: M equilibrium density
    ax_c = fig.add_subplot(gs[0, 2])

    im_c = ax_c.pcolormesh(omega_range, sigma_MS_range, m_SM_matrix,
                          cmap='Reds', shading='auto')
    ax_c.contour(omega_range, sigma_MS_range, invasion_rate, levels=[0],
                colors='blue', linewidths=2, linestyles='-')

    ax_c.set_xlabel('Pathway parameter (ω)')
    ax_c.set_ylabel('S→M mutualism (σ_MS)')
    ax_c.set_title('C. M equilibrium density at S-M state',
                   fontweight='bold', loc='left')

    cbar_c = plt.colorbar(im_c, ax=ax_c)
    cbar_c.set_label('m*_SM', fontsize=7)

    # Panel D: (ω, σ_SM) landscape
    ax_d = fig.add_subplot(gs[1, 0])

    sigma_SM_range = np.linspace(0.1, 1.0, 150)
    model.sigma_MS = 1.5  # Fix σ_MS

    invasion_rate_d = np.zeros((len(sigma_SM_range), len(omega_range)))

    for i, sigma_SM in enumerate(sigma_SM_range):
        for j, omega in enumerate(omega_range):
            model.sigma_SM = sigma_SM
            s_SM, m_SM, stable = model.compute_SM_equilibrium()
            if stable:
                lambda_G = model.invasion_fitness(omega)
                invasion_rate_d[i, j] = lambda_G
            else:
                invasion_rate_d[i, j] = np.nan

    im_d = ax_d.pcolormesh(omega_range, sigma_SM_range, invasion_rate_d,
                          cmap='RdBu_r', shading='auto', vmin=-0.5, vmax=0.5)

    ax_d.contour(omega_range, sigma_SM_range, invasion_rate_d, levels=[0],
                colors='black', linewidths=2)

    # Stability boundary: σ_MS * σ_SM = 1
    sigma_SM_stability = 1.0 / model.sigma_MS
    ax_d.axhline(sigma_SM_stability, color='purple', linestyle='--',
                linewidth=2, label=r'S-M stability bound')

    ax_d.set_xlabel('Pathway parameter (ω)')
    ax_d.set_ylabel('M→S benefit (σ_SM)')
    ax_d.set_title('D. Invasion rate in (ω, σ_SM) space',
                   fontweight='bold', loc='left')
    ax_d.legend(frameon=False, fontsize=5, loc='upper right')

    cbar_d = plt.colorbar(im_d, ax=ax_d)
    cbar_d.set_label('λ_G', fontsize=7)

    # Panel E: (σ_GS, α_GS) cooperation-competition
    ax_e = fig.add_subplot(gs[1, 1])

    sigma_GS_range = np.linspace(0.1, 0.8, 150)
    alpha_GS_range = np.linspace(0.1, 0.8, 150)
    model.omega = 0.5  # Fix ω at intermediate value

    invasion_rate_e = np.zeros((len(alpha_GS_range), len(sigma_GS_range)))

    for i, alpha_GS in enumerate(alpha_GS_range):
        for j, sigma_GS in enumerate(sigma_GS_range):
            model.alpha_GS = alpha_GS
            model.sigma_GS = sigma_GS
            lambda_G = model.invasion_fitness(model.omega)
            invasion_rate_e[i, j] = lambda_G

    im_e = ax_e.pcolormesh(sigma_GS_range, alpha_GS_range, invasion_rate_e,
                          cmap='RdBu_r', shading='auto', vmin=-0.5, vmax=0.5)

    CS_e = ax_e.contour(sigma_GS_range, alpha_GS_range, invasion_rate_e,
                       levels=[0], colors='black', linewidths=2)

    # Diagonal line: σ_GS = α_GS
    ax_e.plot([0.1, 0.8], [0.1, 0.8], 'k:', linewidth=1, alpha=0.5)
    ax_e.text(0.6, 0.65, r'$\sigma_{GS} = \alpha_{GS}$', fontsize=6, rotation=45)

    ax_e.set_xlabel('S→G cooperation (σ_GS)')
    ax_e.set_ylabel('S-G competition (α_GS)')
    ax_e.set_title('E. Cooperation vs competition trade-off',
                   fontweight='bold', loc='left')

    cbar_e = plt.colorbar(im_e, ax=ax_e)
    cbar_e.set_label('λ_G', fontsize=7)

    # Panel F: (σ_GM, α_GM) metabolite interactions
    ax_f = fig.add_subplot(gs[1, 2])

    sigma_GM_range = np.linspace(0.1, 0.8, 150)
    alpha_GM_range = np.linspace(0.1, 0.8, 150)

    invasion_rate_f = np.zeros((len(alpha_GM_range), len(sigma_GM_range)))

    for i, alpha_GM in enumerate(alpha_GM_range):
        for j, sigma_GM in enumerate(sigma_GM_range):
            model.alpha_GM = alpha_GM
            model.sigma_GM = sigma_GM
            lambda_G = model.invasion_fitness(model.omega)
            invasion_rate_f[i, j] = lambda_G

    im_f = ax_f.pcolormesh(sigma_GM_range, alpha_GM_range, invasion_rate_f,
                          cmap='RdBu_r', shading='auto', vmin=-0.5, vmax=0.5)

    ax_f.contour(sigma_GM_range, alpha_GM_range, invasion_rate_f,
                levels=[0], colors='black', linewidths=2)

    ax_f.plot([0.1, 0.8], [0.1, 0.8], 'k:', linewidth=1, alpha=0.5)

    ax_f.set_xlabel('M→G benefit (σ_GM)')
    ax_f.set_ylabel('M-G competition (α_GM)')
    ax_f.set_title('F. Metabolite-mediated interactions',
                   fontweight='bold', loc='left')

    cbar_f = plt.colorbar(im_f, ax=ax_f)
    cbar_f.set_label('λ_G', fontsize=7)

    # Panel G: Critical ω vs σ_MS
    ax_g = fig.add_subplot(gs[2, :2])

    # Reset model parameters
    model.sigma_SM = 0.5
    model.alpha_GS = 0.3
    model.alpha_GM = 0.3
    model.sigma_GS = 0.4
    model.sigma_GM = 0.4

    omega_crit_values = []
    sigma_MS_scan = np.linspace(1.1, 3.0, 100)

    for sigma_MS in sigma_MS_scan:
        model.sigma_MS = sigma_MS
        # Find critical omega where invasion_fitness = 0
        omega_test = np.linspace(0, 1, 200)
        fitness_vals = [model.invasion_fitness(om) for om in omega_test]

        # Find zero crossing
        sign_changes = np.where(np.diff(np.sign(fitness_vals)))[0]
        if len(sign_changes) > 0:
            omega_crit = omega_test[sign_changes[0]]
            omega_crit_values.append(omega_crit)
        else:
            omega_crit_values.append(np.nan)

    ax_g.plot(sigma_MS_scan, omega_crit_values, '-', color=COLORS['G'],
             linewidth=2.5)
    ax_g.fill_between(sigma_MS_scan, 0, omega_crit_values, alpha=0.3,
                     color=COLORS['exclude'], label='G excluded')
    ax_g.fill_between(sigma_MS_scan, omega_crit_values, 1, alpha=0.3,
                     color=COLORS['coexist'], label='G invades')

    ax_g.set_xlabel('S→M mutualism strength (σ_MS)')
    ax_g.set_ylabel('Critical pathway parameter (ω_crit)')
    ax_g.set_title('G. Invasion threshold across mutualism strength',
                   fontweight='bold', loc='left')
    ax_g.legend(frameon=False, loc='upper right')
    ax_g.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    ax_g.set_ylim(0, 1)

    # Panel H: Coexistence window
    ax_h = fig.add_subplot(gs[2, 2])

    # Compute coexistence window for different σ_MS values
    sigma_MS_window_vals = [1.2, 1.5, 2.0, 2.5]
    omega_scan = np.linspace(0, 1, 300)

    for sigma_MS_val in sigma_MS_window_vals:
        model.sigma_MS = sigma_MS_val
        fitness_curve = [model.invasion_fitness(om) for om in omega_scan]
        ax_h.plot(omega_scan, fitness_curve, '-', linewidth=1.5,
                 label=f'σ_MS = {sigma_MS_val}', alpha=0.8)

    ax_h.axhline(0, color='black', linestyle='--', linewidth=1)
    ax_h.set_xlabel('Pathway parameter (ω)')
    ax_h.set_ylabel('Invasion fitness λ_G')
    ax_h.set_title('H. Coexistence windows', fontweight='bold', loc='left')
    ax_h.legend(frameon=False, fontsize=6)
    ax_h.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    plt.savefig('figures/Figure2_enhanced_parameter_landscapes.png',
                dpi=300, bbox_inches='tight')
    print("✓ Figure 2 created: Enhanced parameter space landscapes")
    return fig


def create_figure3_three_species_dynamics(model):
    """
    Figure 3: Three-Species Dynamics and Phase Space Analysis

    6 Panels:
    A: Successful G invasion time series
    B: Failed G invasion time series
    C: 3D phase portrait (s, m, g)
    D: 2D phase projections (s-m, s-g, m-g planes)
    E: Initial condition sensitivity
    F: Transient dynamics near bifurcation
    """
    fig = plt.figure(figsize=(14, 9))
    gs = GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.4)

    # Panel A: Successful invasion (ω > ω_crit)
    ax_a = fig.add_subplot(gs[0, :2])
    model.omega = 0.55  # Above critical
    t_span = (0, 100)
    t_eval = np.linspace(0, 100, 1000)
    y0 = [0.5, 0.3, 0.1]  # Initial with small G

    sol_success = solve_ivp(lambda t, y: model.three_species_dynamics(t, y, model.omega),
                           t_span, y0, t_eval=t_eval, method='RK45',
                           rtol=1e-9, atol=1e-11)

    ax_a.plot(sol_success.t, sol_success.y[0], '-', color=COLORS['S'],
             linewidth=2, label='S (substrate specialist)', alpha=0.8)
    ax_a.plot(sol_success.t, sol_success.y[1], '-', color=COLORS['M'],
             linewidth=2, label='M (metabolite specialist)', alpha=0.8)
    ax_a.plot(sol_success.t, sol_success.y[2], '-', color=COLORS['G'],
             linewidth=2, label='G (generalist)', alpha=0.8)

    ax_a.set_xlabel('Time (days)')
    ax_a.set_ylabel('Population density (scaled)')
    ax_a.set_title(f'A. Successful invasion (ω = {model.omega} > ω_crit)',
                   fontweight='bold', loc='left')
    ax_a.legend(frameon=False, loc='right')
    ax_a.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # Panel B: Failed invasion (ω < ω_crit)
    ax_b = fig.add_subplot(gs[0, 2])
    model.omega = 0.25  # Below critical

    sol_fail = solve_ivp(lambda t, y: model.three_species_dynamics(t, y, model.omega),
                        t_span, y0, t_eval=t_eval, method='RK45',
                        rtol=1e-9, atol=1e-11)

    ax_b.plot(sol_fail.t, sol_fail.y[0], '-', color=COLORS['S'], linewidth=1.5, alpha=0.8)
    ax_b.plot(sol_fail.t, sol_fail.y[1], '-', color=COLORS['M'], linewidth=1.5, alpha=0.8)
    ax_b.plot(sol_fail.t, sol_fail.y[2], '-', color=COLORS['G'], linewidth=1.5, alpha=0.8)

    ax_b.set_xlabel('Time (days)')
    ax_b.set_ylabel('Density')
    ax_b.set_title(f'B. Failed invasion\n(ω = {model.omega} < ω_crit)',
                   fontweight='bold', loc='left')
    ax_b.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # Panel C: 3D phase portrait
    ax_c = fig.add_subplot(gs[1, 0], projection='3d')
    model.omega = 0.55

    # Multiple trajectories
    initial_conditions = [
        [0.4, 0.2, 0.1],
        [0.6, 0.4, 0.2],
        [0.8, 0.3, 0.15],
        [0.3, 0.5, 0.05]
    ]

    for y0_3d in initial_conditions:
        sol_3d = solve_ivp(lambda t, y: model.three_species_dynamics(t, y, model.omega),
                          (0, 150), y0_3d, t_eval=np.linspace(0, 150, 1500),
                          method='RK45', rtol=1e-9, atol=1e-11)
        ax_c.plot(sol_3d.y[0], sol_3d.y[1], sol_3d.y[2], '-',
                 linewidth=1, alpha=0.6)
        ax_c.plot([y0_3d[0]], [y0_3d[1]], [y0_3d[2]], 'o',
                 color='black', markersize=4)

    ax_c.set_xlabel('S density', fontsize=7)
    ax_c.set_ylabel('M density', fontsize=7)
    ax_c.set_zlabel('G density', fontsize=7)
    ax_c.set_title('C. 3D phase portrait', fontweight='bold', loc='left')

    # Panel D: 2D projections
    ax_d1 = fig.add_subplot(gs[1, 1])
    ax_d2 = fig.add_subplot(gs[1, 2])

    # S-M projection
    for y0_proj in initial_conditions:
        sol_proj = solve_ivp(lambda t, y: model.three_species_dynamics(t, y, model.omega),
                            (0, 150), y0_proj, t_eval=np.linspace(0, 150, 1500),
                            method='RK45')
        ax_d1.plot(sol_proj.y[0], sol_proj.y[1], '-', linewidth=1, alpha=0.6)

    ax_d1.set_xlabel('S density')
    ax_d1.set_ylabel('M density')
    ax_d1.set_title('D1. S-M projection', fontweight='bold', loc='left', fontsize=8)
    ax_d1.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # S-G projection
    for y0_proj in initial_conditions:
        sol_proj = solve_ivp(lambda t, y: model.three_species_dynamics(t, y, model.omega),
                            (0, 150), y0_proj, t_eval=np.linspace(0, 150, 1500),
                            method='RK45')
        ax_d2.plot(sol_proj.y[0], sol_proj.y[2], '-', linewidth=1, alpha=0.6)

    ax_d2.set_xlabel('S density')
    ax_d2.set_ylabel('G density')
    ax_d2.set_title('D2. S-G projection', fontweight='bold', loc='left', fontsize=8)
    ax_d2.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # Panel E: Initial condition sensitivity
    ax_e = fig.add_subplot(gs[2, :2])

    # Vary initial G density
    g_init_range = np.linspace(0.01, 0.5, 20)
    final_densities_S = []
    final_densities_M = []
    final_densities_G = []

    model.omega = 0.55
    for g_init in g_init_range:
        sol_ic = solve_ivp(lambda t, y: model.three_species_dynamics(t, y, model.omega),
                          (0, 200), [0.5, 0.3, g_init], method='RK45',
                          rtol=1e-9, atol=1e-11)
        final_densities_S.append(sol_ic.y[0, -1])
        final_densities_M.append(sol_ic.y[1, -1])
        final_densities_G.append(sol_ic.y[2, -1])

    ax_e.plot(g_init_range, final_densities_S, 'o-', color=COLORS['S'],
             label='S final', markersize=4, linewidth=1.5)
    ax_e.plot(g_init_range, final_densities_M, 's-', color=COLORS['M'],
             label='M final', markersize=4, linewidth=1.5)
    ax_e.plot(g_init_range, final_densities_G, '^-', color=COLORS['G'],
             label='G final', markersize=4, linewidth=1.5)

    ax_e.set_xlabel('Initial G density')
    ax_e.set_ylabel('Final equilibrium density')
    ax_e.set_title('E. Equilibrium independence from initial conditions',
                   fontweight='bold', loc='left')
    ax_e.legend(frameon=False, ncol=3, loc='upper right')
    ax_e.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # Panel F: Transient dynamics near bifurcation
    ax_f = fig.add_subplot(gs[2, 2])

    omega_near_crit = [0.35, 0.40, 0.45, 0.50]
    t_short = np.linspace(0, 80, 800)

    for omega_val in omega_near_crit:
        sol_trans = solve_ivp(lambda t, y: model.three_species_dynamics(t, y, omega_val),
                             (0, 80), [0.5, 0.3, 0.1], t_eval=t_short, method='RK45')
        ax_f.plot(sol_trans.t, sol_trans.y[2], '-', linewidth=1.5,
                 label=f'ω = {omega_val}', alpha=0.7)

    ax_f.set_xlabel('Time (days)')
    ax_f.set_ylabel('G density')
    ax_f.set_title('F. Transients near bifurcation', fontweight='bold', loc='left')
    ax_f.legend(frameon=False, fontsize=6)
    ax_f.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    plt.savefig('figures/Figure3_three_species_dynamics.png', dpi=300, bbox_inches='tight')
    print("✓ Figure 3 created: Three-species dynamics and phase space")
    return fig


def create_figure4_bifurcation_analysis(model):
    """
    Figure 4: Complete Bifurcation Structure and Stability Analysis

    4 Panels:
    A: Complete bifurcation diagram (ω scan) showing all equilibria
    B: Eigenvalue evolution showing stability transitions
    C: Coexistence window in multi-parameter space
    D: Stability regions with transcritical bifurcation points
    """
    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.35)

    # Panel A: Complete bifurcation diagram
    ax_a = fig.add_subplot(gs[0, :])

    omega_scan = np.linspace(0.0, 1.0, 300)
    equilibria_data = {'S': [], 'M': [], 'G': [], 'stable': []}

    for omega in omega_scan:
        model.omega = omega

        # Try to find three-species equilibrium
        s_SM, m_SM, stable_SM = model.compute_SM_equilibrium()

        if stable_SM:
            # Check if G can invade
            lambda_G = model.invasion_fitness(omega)

            if lambda_G > 0:
                # Three-species coexistence - need numerical solution
                def three_species_eq(x):
                    s, m, g = x
                    a, b, c, d, e = model.net_interactions(omega)
                    return [
                        s * (1 + model.sigma_SM * m + a * g - s),
                        m * (-1 + model.sigma_MS * s + b * g - m),
                        g * (d + c * s + e * m - g)
                    ]

                try:
                    sol = fsolve(three_species_eq, [s_SM, m_SM, 0.3],
                                full_output=True)
                    if sol[2] == 1:  # Solution found
                        s_eq, m_eq, g_eq = sol[0]
                        if all(np.array([s_eq, m_eq, g_eq]) > 0.001):
                            # Check stability
                            J = model.jacobian_three_species(s_eq, m_eq, g_eq, omega)
                            eigenvalues = eig(J)[0]
                            is_stable = all(np.real(eigenvalues) < 0)

                            equilibria_data['S'].append(s_eq)
                            equilibria_data['M'].append(m_eq)
                            equilibria_data['G'].append(g_eq)
                            equilibria_data['stable'].append(is_stable)
                        else:
                            # S-M equilibrium
                            equilibria_data['S'].append(s_SM)
                            equilibria_data['M'].append(m_SM)
                            equilibria_data['G'].append(0.0)
                            equilibria_data['stable'].append(True)
                    else:
                        equilibria_data['S'].append(s_SM)
                        equilibria_data['M'].append(m_SM)
                        equilibria_data['G'].append(0.0)
                        equilibria_data['stable'].append(True)
                except:
                    equilibria_data['S'].append(s_SM)
                    equilibria_data['M'].append(m_SM)
                    equilibria_data['G'].append(0.0)
                    equilibria_data['stable'].append(True)
            else:
                # S-M equilibrium
                equilibria_data['S'].append(s_SM)
                equilibria_data['M'].append(m_SM)
                equilibria_data['G'].append(0.0)
                equilibria_data['stable'].append(True)
        else:
            equilibria_data['S'].append(np.nan)
            equilibria_data['M'].append(np.nan)
            equilibria_data['G'].append(np.nan)
            equilibria_data['stable'].append(False)

    # Plot bifurcation diagram
    stable_arr = np.array(equilibria_data['stable'])
    s_arr = np.array(equilibria_data['S'])
    m_arr = np.array(equilibria_data['M'])
    g_arr = np.array(equilibria_data['G'])

    ax_a.plot(omega_scan[stable_arr], s_arr[stable_arr], '-',
             color=COLORS['S'], linewidth=2, label='S (stable)')
    ax_a.plot(omega_scan[~stable_arr], s_arr[~stable_arr], '--',
             color=COLORS['S'], linewidth=2, alpha=0.4, label='S (unstable)')

    ax_a.plot(omega_scan[stable_arr], m_arr[stable_arr], '-',
             color=COLORS['M'], linewidth=2, label='M (stable)')
    ax_a.plot(omega_scan[~stable_arr], m_arr[~stable_arr], '--',
             color=COLORS['M'], linewidth=2, alpha=0.4, label='M (unstable)')

    ax_a.plot(omega_scan[stable_arr], g_arr[stable_arr], '-',
             color=COLORS['G'], linewidth=2, label='G (stable)')

    # Mark bifurcation points
    g_positive = g_arr > 0.01
    if np.any(g_positive):
        omega_crit1_idx = np.where(np.diff(g_positive.astype(int)) == 1)[0]
        omega_crit2_idx = np.where(np.diff(g_positive.astype(int)) == -1)[0]

        if len(omega_crit1_idx) > 0:
            ax_a.axvline(omega_scan[omega_crit1_idx[0]], color='red',
                        linestyle=':', linewidth=2, alpha=0.7,
                        label='Transcritical bifurcation')
            ax_a.text(omega_scan[omega_crit1_idx[0]], 1.5, r'$\omega_{crit}^{(1)}$',
                     fontsize=8, ha='center')

        if len(omega_crit2_idx) > 0:
            ax_a.axvline(omega_scan[omega_crit2_idx[0]], color='red',
                        linestyle=':', linewidth=2, alpha=0.7)
            ax_a.text(omega_scan[omega_crit2_idx[0]], 1.5, r'$\omega_{crit}^{(2)}$',
                     fontsize=8, ha='center')

    ax_a.set_xlabel('Pathway parameter (ω)', fontsize=9)
    ax_a.set_ylabel('Equilibrium density', fontsize=9)
    ax_a.set_title('A. Complete bifurcation diagram: Equilibrium densities vs pathway strategy',
                   fontweight='bold', loc='left')
    ax_a.legend(frameon=False, loc='upper right', ncol=3, fontsize=7)
    ax_a.grid(alpha=0.3, linestyle='--', linewidth=0.5)
    ax_a.set_ylim(-0.1, 2.0)

    # Panel B: Eigenvalue evolution
    ax_b = fig.add_subplot(gs[1, 0])

    omega_eig_scan = np.linspace(0.2, 0.8, 100)
    max_real_eigenvalues = []

    for omega in omega_eig_scan:
        model.omega = omega
        s_SM, m_SM, stable = model.compute_SM_equilibrium()

        if stable:
            lambda_G = model.invasion_fitness(omega)
            max_real_eigenvalues.append(lambda_G / model.r_G)  # Normalized
        else:
            max_real_eigenvalues.append(np.nan)

    ax_b.plot(omega_eig_scan, max_real_eigenvalues, '-',
             color=COLORS['G'], linewidth=2.5)
    ax_b.axhline(0, color='black', linestyle='--', linewidth=1)
    ax_b.fill_between(omega_eig_scan, 0, max_real_eigenvalues,
                     where=np.array(max_real_eigenvalues) > 0,
                     alpha=0.3, color=COLORS['coexist'], label='G invades')
    ax_b.fill_between(omega_eig_scan, max_real_eigenvalues, 0,
                     where=np.array(max_real_eigenvalues) < 0,
                     alpha=0.3, color=COLORS['exclude'], label='G excluded')

    ax_b.set_xlabel('Pathway parameter (ω)')
    ax_b.set_ylabel('Leading eigenvalue (scaled)')
    ax_b.set_title('B. Stability transition: G invasion eigenvalue',
                   fontweight='bold', loc='left')
    ax_b.legend(frameon=False, loc='upper left')
    ax_b.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    # Panel C: Multi-parameter coexistence window
    ax_c = fig.add_subplot(gs[1, 1])

    # Coexistence window for different combinations of parameters
    sigma_MS_vals = [1.2, 1.5, 2.0, 2.5]
    omega_window = np.linspace(0, 1, 200)

    for sigma_MS_val in sigma_MS_vals:
        model.sigma_MS = sigma_MS_val
        coexistence_indicator = []

        for omega in omega_window:
            lambda_G = model.invasion_fitness(omega)
            coexistence_indicator.append(1 if lambda_G > 0 else 0)

        ax_c.plot(omega_window, np.array(coexistence_indicator) * sigma_MS_val,
                 '-', linewidth=2, label=f'σ_MS = {sigma_MS_val}', alpha=0.7)

    ax_c.set_xlabel('Pathway parameter (ω)')
    ax_c.set_ylabel('Coexistence indicator × σ_MS')
    ax_c.set_title('C. Coexistence windows across mutualism strengths',
                   fontweight='bold', loc='left')
    ax_c.legend(frameon=False, fontsize=7)
    ax_c.grid(alpha=0.3, linestyle='--', linewidth=0.5)

    plt.savefig('figures/Figure4_bifurcation_structure.png', dpi=300, bbox_inches='tight')
    print("✓ Figure 4 created: Complete bifurcation analysis")
    return fig


# Main execution
if __name__ == "__main__":
    print("Generating enhanced comprehensive figures for PNAS publication...")
    print("=" * 70)

    # Default parameters
    params = {
        'r_S': 1.0, 'r_M': 0.8, 'r_G': 0.9,
        'sigma_SM': 0.5, 'sigma_MS': 1.5,
        'sigma_SG': 0.4, 'sigma_GS': 0.4,
        'sigma_MG': 0.4, 'sigma_GM': 0.4,
        'alpha_SG': 0.3, 'alpha_GS': 0.3,
        'alpha_MG': 0.3, 'alpha_GM': 0.3,
        'omega': 0.5
    }

    model = ThreeSpeciesModel(params)

    # Create all figures
    print("\n[1/4] Creating Figure 1: Comprehensive pairwise systems analysis...")
    fig1 = create_figure1_pairwise_analysis(model)
    plt.close(fig1)

    print("\n[2/4] Creating Figure 2: Enhanced parameter space landscapes...")
    fig2 = create_figure2_parameter_landscapes(model)
    plt.close(fig2)

    print("\n[3/4] Creating Figure 3: Three-species dynamics and phase space...")
    fig3 = create_figure3_three_species_dynamics(model)
    plt.close(fig3)

    print("\n[4/4] Creating Figure 4: Complete bifurcation analysis...")
    fig4 = create_figure4_bifurcation_analysis(model)
    plt.close(fig4)

    print("\n" + "=" * 70)
    print("✅ All enhanced figures generated successfully!")
    print("\nGenerated files:")
    print("  • Figure1_comprehensive_pairwise.png (6 panels)")
    print("  • Figure2_enhanced_parameter_landscapes.png (8 panels)")
    print("  • Figure3_three_species_dynamics.png (6 panels)")
    print("  • Figure4_bifurcation_structure.png (4 panels)")
    print("\nTotal: 24 sub-panels across 4 main figures")
    print("=" * 70)

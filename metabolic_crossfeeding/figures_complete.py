"""
Comprehensive Figures for American Naturalist Manuscript

Generates publication-quality figures covering:
1. Model schematic and biological context
2. Fitness landscape and trade-offs (both species)
3. Selection gradients for all four strategies
4. Phase portraits and evolutionary trajectories
5. Stability analysis (eigenvalues, bifurcation)
6. Pairwise invasibility plots
7. Parameter sensitivity analysis
8. Basin of attraction analysis
9. Cheater invasion analysis
10. Summary and biological interpretation

Author: Jian Wang
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, gridspec
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, Wedge
from matplotlib.colors import LinearSegmentedColormap
from scipy.integrate import odeint
import os

from model_complete import CrossFeedingModel

# Publication style
plt.rcParams.update({
    'font.size': 9,
    'font.family': 'serif',
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'lines.linewidth': 1.2,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Color scheme
C = {
    'A': '#E64B35',      # Species A (red-orange)
    'B': '#4DBBD5',      # Species B (blue)
    'F1': '#00A087',     # Amino acid F1 (teal)
    'F2': '#3C5488',     # Amino acid F2 (navy)
    'stable': '#009E73', # Stable equilibrium (green)
    'unstable': '#D55E00', # Unstable (orange)
    'neutral': '#999999',  # Neutral (gray)
}


def ensure_dir():
    """Create figures directory."""
    fig_dir = os.path.join(os.path.dirname(__file__), 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    return fig_dir


# =============================================================================
# FIGURE 1: Model Schematic (Full Page)
# =============================================================================

def figure1_model_schematic(save=True):
    """
    Comprehensive model schematic showing:
    - Bioreactor setup
    - Two species with investment decisions
    - Amino acid production and consumption
    - Trade-off structure
    """
    fig = plt.figure(figsize=(7, 8))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1.2, 1, 1], hspace=0.35, wspace=0.3)

    # Panel A: Biological context
    ax = fig.add_subplot(gs[0, :])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title('(A) Metabolic cross-feeding in a bioreactor', fontweight='bold', loc='left')

    # Bioreactor vessel
    vessel = Rectangle((0.2, 0.15), 0.6, 0.7, fill=False,
                       edgecolor='black', linewidth=2)
    ax.add_patch(vessel)

    # Species A (left side)
    for x, y in [(0.32, 0.6), (0.38, 0.45), (0.28, 0.35)]:
        c = Circle((x, y), 0.05, color=C['A'], alpha=0.8)
        ax.add_patch(c)
        ax.text(x, y, 'A', ha='center', va='center', fontsize=10,
               fontweight='bold', color='white')

    # Species B (right side)
    for x, y in [(0.62, 0.55), (0.68, 0.4), (0.58, 0.7)]:
        c = Circle((x, y), 0.05, color=C['B'], alpha=0.8)
        ax.add_patch(c)
        ax.text(x, y, 'B', ha='center', va='center', fontsize=10,
               fontweight='bold', color='white')

    # Amino acid pool (center)
    pool = Rectangle((0.42, 0.35), 0.16, 0.25, fill=True,
                    facecolor='lightyellow', edgecolor='black', alpha=0.7)
    ax.add_patch(pool)
    ax.text(0.5, 0.52, 'Public', ha='center', va='center', fontsize=8)
    ax.text(0.5, 0.47, 'Amino Acid', ha='center', va='center', fontsize=8)
    ax.text(0.5, 0.42, 'Pool', ha='center', va='center', fontsize=8)
    ax.text(0.45, 0.37, '$F_1$', ha='center', fontsize=9, color=C['F1'], fontweight='bold')
    ax.text(0.55, 0.37, '$F_2$', ha='center', fontsize=9, color=C['F2'], fontweight='bold')

    # Arrows: production
    ax.annotate('', xy=(0.42, 0.48), xytext=(0.35, 0.5),
                arrowprops=dict(arrowstyle='->', color=C['F1'], lw=1.5))
    ax.annotate('', xy=(0.58, 0.48), xytext=(0.65, 0.5),
                arrowprops=dict(arrowstyle='->', color=C['F2'], lw=1.5))

    # Inflow/outflow
    ax.annotate('', xy=(0.2, 0.5), xytext=(0.08, 0.5),
                arrowprops=dict(arrowstyle='->', color='gray', lw=2))
    ax.text(0.05, 0.55, 'Glucose\ninflow', fontsize=7, ha='center')

    ax.annotate('', xy=(0.92, 0.5), xytext=(0.8, 0.5),
                arrowprops=dict(arrowstyle='->', color='gray', lw=2))
    ax.text(0.95, 0.55, 'Dilution\nD', fontsize=7, ha='center')

    # Investment labels
    ax.text(0.3, 0.78, 'Species A', fontsize=10, ha='center', color=C['A'], fontweight='bold')
    ax.text(0.3, 0.72, '$f_{A1}, f_{A2}$', fontsize=9, ha='center')
    ax.text(0.7, 0.78, 'Species B', fontsize=10, ha='center', color=C['B'], fontweight='bold')
    ax.text(0.7, 0.72, '$f_{B1}, f_{B2}$', fontsize=9, ha='center')

    # Panel B: Trade-off structure
    ax = fig.add_subplot(gs[1, 0])
    f = np.linspace(0, 0.6, 100)
    ax.fill_between(f, 0, 1-f, alpha=0.3, color='gray', label='Growth allocation')
    ax.fill_between(f, 1-f, 1, alpha=0.3, color=C['F1'], label='Investment')
    ax.plot(f, 1-f, 'k-', lw=2)
    ax.set_xlabel('Total investment ($f_1 + f_2$)')
    ax.set_ylabel('Resource allocation')
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 1)
    ax.legend(loc='upper right', fontsize=7)
    ax.set_title('(B) Investment-growth trade-off', fontweight='bold', loc='left')

    # Panel C: Public goods dynamics
    ax = fig.add_subplot(gs[1, 1])
    P = np.linspace(0.01, 1, 100)
    for alpha in [0.5, 0.75, 1.0]:
        g = P ** (2*alpha)
        ax.plot(P, g, label=f'$\\alpha = {alpha}$', lw=1.5)
    ax.set_xlabel('Total amino acid ($P_1 \\cdot P_2$)')
    ax.set_ylabel('Growth rate $g$')
    ax.legend(loc='lower right', fontsize=7)
    ax.set_title('(C) Growth depends on both amino acids', fontweight='bold', loc='left')

    # Panel D: Two equilibria
    ax = fig.add_subplot(gs[2, 0])

    # Symmetric
    theta = np.linspace(0, 2*np.pi, 100)
    r = 0.3
    # Pie chart for symmetric
    wedges = [
        Wedge((0.25, 0.5), 0.15, 0, 90, facecolor=C['F1'], alpha=0.7),
        Wedge((0.25, 0.5), 0.15, 90, 180, facecolor=C['F2'], alpha=0.7),
        Wedge((0.25, 0.5), 0.15, 180, 360, facecolor='lightgray', alpha=0.7),
    ]
    for w in wedges:
        ax.add_patch(w)
    ax.text(0.25, 0.5, 'A', ha='center', va='center', fontsize=12, fontweight='bold')
    ax.text(0.25, 0.25, 'Symmetric\n$f_{A1}=f_{A2}=0.25$', ha='center', fontsize=7)

    wedges2 = [
        Wedge((0.75, 0.5), 0.15, 0, 90, facecolor=C['F1'], alpha=0.7),
        Wedge((0.75, 0.5), 0.15, 90, 180, facecolor=C['F2'], alpha=0.7),
        Wedge((0.75, 0.5), 0.15, 180, 360, facecolor='lightgray', alpha=0.7),
    ]
    for w in wedges2:
        ax.add_patch(w)
    ax.text(0.75, 0.5, 'B', ha='center', va='center', fontsize=12, fontweight='bold')
    ax.text(0.75, 0.25, 'Symmetric\n$f_{B1}=f_{B2}=0.25$', ha='center', fontsize=7)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title('(D) Symmetric equilibrium', fontweight='bold', loc='left')

    # Panel E: Division of labor
    ax = fig.add_subplot(gs[2, 1])

    # Species A specializes on F1
    wedges3 = [
        Wedge((0.25, 0.5), 0.15, 0, 180, facecolor=C['F1'], alpha=0.9),
        Wedge((0.25, 0.5), 0.15, 180, 360, facecolor='lightgray', alpha=0.7),
    ]
    for w in wedges3:
        ax.add_patch(w)
    ax.text(0.25, 0.5, 'A', ha='center', va='center', fontsize=12, fontweight='bold')
    ax.text(0.25, 0.25, 'Specialist\n$f_{A1}=0.5, f_{A2}=0$', ha='center', fontsize=7)

    # Species B specializes on F2
    wedges4 = [
        Wedge((0.75, 0.5), 0.15, 0, 180, facecolor=C['F2'], alpha=0.9),
        Wedge((0.75, 0.5), 0.15, 180, 360, facecolor='lightgray', alpha=0.7),
    ]
    for w in wedges4:
        ax.add_patch(w)
    ax.text(0.75, 0.5, 'B', ha='center', va='center', fontsize=12, fontweight='bold')
    ax.text(0.75, 0.25, 'Specialist\n$f_{B1}=0, f_{B2}=0.5$', ha='center', fontsize=7)

    # Arrow showing exchange
    ax.annotate('', xy=(0.6, 0.55), xytext=(0.4, 0.55),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax.text(0.5, 0.62, 'Exchange', ha='center', fontsize=8)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title('(E) Division of labor equilibrium', fontweight='bold', loc='left')

    if save:
        fig_dir = ensure_dir()
        plt.savefig(os.path.join(fig_dir, 'fig1_schematic.pdf'))
        plt.savefig(os.path.join(fig_dir, 'fig1_schematic.png'))

    return fig


# =============================================================================
# FIGURE 2: Selection Gradients (All Four)
# =============================================================================

def figure2_selection_gradients(save=True):
    """
    Selection gradients for all four investment strategies.
    """
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    model = CrossFeedingModel()

    f_range = np.linspace(0.01, 0.6, 50)
    _, state_sym = model.symmetric_equilibrium()
    f_sym = state_sym[0]

    # Panel A: ∂W_A/∂f_A1
    ax = axes[0, 0]
    for f_A2_fixed in [0.0, 0.1, 0.25]:
        grads = []
        for f_A1 in f_range:
            g = model.selection_gradient_A1(f_A1, f_A2_fixed, 0.25, 0.25)
            grads.append(g)
        label = f'$f_{{A2}} = {f_A2_fixed}$'
        ax.plot(f_range, grads, label=label, lw=1.5)

    ax.axhline(0, ls='--', color='gray', lw=0.5)
    ax.axvline(f_sym, ls=':', color=C['unstable'], alpha=0.5)
    ax.set_xlabel('$f_{A1}$')
    ax.set_ylabel('$\\partial W_A / \\partial f_{A1}$')
    ax.set_title('(A) Selection on A investing in $F_1$', fontweight='bold', loc='left')
    ax.legend(loc='best', fontsize=7)

    # Panel B: ∂W_A/∂f_A2
    ax = axes[0, 1]
    for f_A1_fixed in [0.0, 0.1, 0.25]:
        grads = []
        for f_A2 in f_range:
            g = model.selection_gradient_A2(f_A1_fixed, f_A2, 0.25, 0.25)
            grads.append(g)
        label = f'$f_{{A1}} = {f_A1_fixed}$'
        ax.plot(f_range, grads, label=label, lw=1.5)

    ax.axhline(0, ls='--', color='gray', lw=0.5)
    ax.set_xlabel('$f_{A2}$')
    ax.set_ylabel('$\\partial W_A / \\partial f_{A2}$')
    ax.set_title('(B) Selection on A investing in $F_2$', fontweight='bold', loc='left')
    ax.legend(loc='best', fontsize=7)

    # Panel C: ∂W_B/∂f_B1
    ax = axes[1, 0]
    for f_B2_fixed in [0.0, 0.1, 0.25]:
        grads = []
        for f_B1 in f_range:
            g = model.selection_gradient_B1(0.25, 0.25, f_B1, f_B2_fixed)
            grads.append(g)
        label = f'$f_{{B2}} = {f_B2_fixed}$'
        ax.plot(f_range, grads, label=label, lw=1.5)

    ax.axhline(0, ls='--', color='gray', lw=0.5)
    ax.set_xlabel('$f_{B1}$')
    ax.set_ylabel('$\\partial W_B / \\partial f_{B1}$')
    ax.set_title('(C) Selection on B investing in $F_1$', fontweight='bold', loc='left')
    ax.legend(loc='best', fontsize=7)

    # Panel D: ∂W_B/∂f_B2
    ax = axes[1, 1]
    for f_B1_fixed in [0.0, 0.1, 0.25]:
        grads = []
        for f_B2 in f_range:
            g = model.selection_gradient_B2(0.25, 0.25, f_B1_fixed, f_B2)
            grads.append(g)
        label = f'$f_{{B1}} = {f_B1_fixed}$'
        ax.plot(f_range, grads, label=label, lw=1.5)

    ax.axhline(0, ls='--', color='gray', lw=0.5)
    ax.set_xlabel('$f_{B2}$')
    ax.set_ylabel('$\\partial W_B / \\partial f_{B2}$')
    ax.set_title('(D) Selection on B investing in $F_2$', fontweight='bold', loc='left')
    ax.legend(loc='best', fontsize=7)

    plt.tight_layout()

    if save:
        fig_dir = ensure_dir()
        plt.savefig(os.path.join(fig_dir, 'fig2_gradients.pdf'))
        plt.savefig(os.path.join(fig_dir, 'fig2_gradients.png'))

    return fig


# =============================================================================
# FIGURE 3: Phase Portrait and Trajectories
# =============================================================================

def figure3_phase_portrait(save=True):
    """
    Phase portrait showing evolutionary dynamics in strategy space.
    """
    fig, axes = plt.subplots(1, 3, figsize=(11, 4))
    model = CrossFeedingModel({'gamma': 1.0, 'alpha': 1.0, 'D': 0.1, 'sigma': 0.02,
                               'mu_A': 1.0, 'mu_B': 1.0, 'c_A': 1.0, 'c_B': 1.0,
                               'eta_A1': 1.0, 'eta_A2': 1.0, 'eta_B1': 1.0, 'eta_B2': 1.0,
                               'delta': 0.0})

    # Panel A: Species A strategy space
    ax = axes[0]
    n = 15
    f1_range = np.linspace(0.02, 0.48, n)
    f2_range = np.linspace(0.02, 0.48, n)
    F1, F2 = np.meshgrid(f1_range, f2_range)
    dF1 = np.zeros_like(F1)
    dF2 = np.zeros_like(F2)

    for i in range(n):
        for j in range(n):
            f_A1, f_A2 = F1[i, j], F2[i, j]
            # Assume B at symmetric
            state = [f_A1, f_A2, 0.25, 0.25]
            derivs = model.evolutionary_dynamics(state, 0)
            dF1[i, j] = derivs[0]
            dF2[i, j] = derivs[1]

    mag = np.sqrt(dF1**2 + dF2**2)
    ax.quiver(F1, F2, dF1/mag, dF2/mag, mag, cmap='viridis', alpha=0.6)
    ax.plot([0, 0.5], [0.5, 0], 'k--', lw=0.5, alpha=0.5)  # constraint
    ax.plot(0.25, 0.25, 'o', color=C['unstable'], markersize=10,
            markeredgecolor='black', label='Symmetric')
    ax.set_xlabel('$f_{A1}$')
    ax.set_ylabel('$f_{A2}$')
    ax.set_title('(A) Species A strategy space', fontweight='bold', loc='left')
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 0.5)
    ax.legend(loc='upper right', fontsize=7)

    # Panel B: Specialization coordinates
    ax = axes[1]
    n = 12
    spec_A = np.linspace(-0.35, 0.35, n)
    spec_B = np.linspace(-0.35, 0.35, n)
    SA, SB = np.meshgrid(spec_A, spec_B)
    dSA = np.zeros_like(SA)
    dSB = np.zeros_like(SB)

    f_base = 0.25
    for i in range(n):
        for j in range(n):
            s_A = SA[i, j]
            s_B = SB[i, j]
            f_A1 = np.clip(f_base + s_A/2, 0.01, 0.49)
            f_A2 = np.clip(f_base - s_A/2, 0.01, 0.49)
            f_B1 = np.clip(f_base - s_B/2, 0.01, 0.49)
            f_B2 = np.clip(f_base + s_B/2, 0.01, 0.49)

            state = [f_A1, f_A2, f_B1, f_B2]
            derivs = model.evolutionary_dynamics(state, 0)
            dSA[i, j] = derivs[0] - derivs[1]
            dSB[i, j] = derivs[3] - derivs[2]

    mag = np.sqrt(dSA**2 + dSB**2)
    mag[mag == 0] = 1
    ax.quiver(SA, SB, dSA/mag, dSB/mag, mag, cmap='viridis', alpha=0.6)

    # Trajectories
    ics = [
        [0.26, 0.24, 0.24, 0.26],
        [0.24, 0.26, 0.26, 0.24],
        [0.27, 0.23, 0.23, 0.27],
        [0.23, 0.27, 0.27, 0.23],
    ]
    for ic in ics:
        t, sol = model.simulate(ic, t_max=500, n_points=300)
        sA = sol[:, 0] - sol[:, 1]
        sB = sol[:, 3] - sol[:, 2]
        ax.plot(sA, sB, '-', color=C['A'], lw=1, alpha=0.7)
        ax.plot(sA[0], sB[0], 'o', color=C['A'], markersize=4)

    ax.plot(0, 0, 'o', color=C['unstable'], markersize=10, markeredgecolor='black')
    ax.plot(0.5, 0.5, 's', color=C['stable'], markersize=10, markeredgecolor='black')
    ax.plot(-0.5, -0.5, 's', color=C['stable'], markersize=10, markeredgecolor='black')

    ax.axhline(0, ls='--', color='gray', lw=0.5, alpha=0.5)
    ax.axvline(0, ls='--', color='gray', lw=0.5, alpha=0.5)
    ax.set_xlabel('A specialization ($f_{A1} - f_{A2}$)')
    ax.set_ylabel('B specialization ($f_{B2} - f_{B1}$)')
    ax.set_title('(B) Coevolutionary dynamics', fontweight='bold', loc='left')
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(-0.6, 0.6)

    # Panel C: Time series
    ax = axes[2]
    ic = [0.26, 0.24, 0.24, 0.26]
    t, sol = model.simulate(ic, t_max=400, n_points=400)

    ax.plot(t, sol[:, 0], '-', color=C['A'], lw=1.5, label='$f_{A1}$')
    ax.plot(t, sol[:, 1], '--', color=C['A'], lw=1.5, label='$f_{A2}$')
    ax.plot(t, sol[:, 2], '--', color=C['B'], lw=1.5, label='$f_{B1}$')
    ax.plot(t, sol[:, 3], '-', color=C['B'], lw=1.5, label='$f_{B2}$')

    ax.axhline(0.25, ls=':', color='gray', alpha=0.5)
    ax.axhline(0.5, ls=':', color=C['stable'], alpha=0.5)
    ax.text(20, 0.27, 'Symmetric', fontsize=7, color='gray')
    ax.text(20, 0.52, 'Specialized', fontsize=7, color=C['stable'])

    ax.set_xlabel('Evolutionary time')
    ax.set_ylabel('Investment fraction')
    ax.set_title('(C) Evolutionary trajectory', fontweight='bold', loc='left')
    ax.legend(loc='center right', fontsize=7, ncol=2)
    ax.set_xlim(0, 400)
    ax.set_ylim(0, 0.6)

    plt.tight_layout()

    if save:
        fig_dir = ensure_dir()
        plt.savefig(os.path.join(fig_dir, 'fig3_phase.pdf'))
        plt.savefig(os.path.join(fig_dir, 'fig3_phase.png'))

    return fig


# =============================================================================
# FIGURE 4: Stability Analysis
# =============================================================================

def figure4_stability(save=True):
    """
    Stability analysis: eigenvalues and parameter effects.
    """
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))

    # Panel A: Eigenvalues at symmetric equilibrium
    ax = axes[0, 0]
    model = CrossFeedingModel()
    _, state_sym = model.symmetric_equilibrium()
    stab = model.stability_analysis(state_sym)

    eigs = stab['eigenvalues']
    ax.scatter(np.real(eigs), np.imag(eigs), s=100, c='black', zorder=10)
    ax.axhline(0, ls='-', color='gray', lw=0.5)
    ax.axvline(0, ls='-', color='gray', lw=0.5)
    ax.fill_betweenx([-0.01, 0.01], 0, 0.01, alpha=0.2, color=C['unstable'])
    ax.fill_betweenx([-0.01, 0.01], -0.01, 0, alpha=0.2, color=C['stable'])
    ax.set_xlabel('Real part')
    ax.set_ylabel('Imaginary part')
    ax.set_title('(A) Eigenvalues at symmetric eq.', fontweight='bold', loc='left')
    ax.text(0.003, 0.005, 'Unstable', fontsize=8, color=C['unstable'])

    # Panel B: Effect of α on equilibria
    ax = axes[0, 1]
    alpha_range = np.linspace(0.3, 1.0, 30)
    f_sym_vals = []
    f_div_vals = []

    for alpha in alpha_range:
        f_sym_vals.append(alpha / (2 * (1 + alpha)))
        f_div_vals.append(alpha / (1 + alpha))

    ax.plot(alpha_range, f_sym_vals, '-', color=C['unstable'], lw=2, label='Symmetric $f^*$')
    ax.plot(alpha_range, f_div_vals, '-', color=C['stable'], lw=2, label='Division of labor $f^*_{div}$')
    ax.set_xlabel('Returns to investment ($\\alpha$)')
    ax.set_ylabel('Equilibrium investment')
    ax.legend(loc='best', fontsize=8)
    ax.set_title('(B) Equilibria vs. $\\alpha$', fontweight='bold', loc='left')

    # Panel C: Fitness comparison
    ax = axes[1, 0]
    alpha_range = np.linspace(0.3, 1.0, 50)
    W_sym = []
    W_div = []

    for alpha in alpha_range:
        m = CrossFeedingModel({'gamma': 1.0, 'alpha': alpha, 'D': 0.05,
                               'mu_A': 1.0, 'mu_B': 1.0, 'c_A': 1.0, 'c_B': 1.0,
                               'eta_A1': 1.0, 'eta_A2': 1.0, 'eta_B1': 1.0, 'eta_B2': 1.0,
                               'sigma': 0.01, 'delta': 0.0})
        f_s = alpha / (2 * (1 + alpha))
        f_d = alpha / (1 + alpha)

        # Symmetric fitness
        P_s = 2 * f_s
        g_s = (P_s ** alpha) * (P_s ** alpha)
        W_s = (1 - 2*f_s) * g_s - 0.05

        # Division of labor fitness
        P_d = f_d
        g_d = (P_d ** alpha) * (P_d ** alpha)
        W_d = (1 - f_d) * g_d - 0.05

        W_sym.append(W_s)
        W_div.append(W_d)

    ax.plot(alpha_range, W_sym, '--', color=C['unstable'], lw=2, label='Symmetric')
    ax.plot(alpha_range, W_div, '-', color=C['stable'], lw=2, label='Division of labor')
    ax.fill_between(alpha_range, W_sym, W_div,
                   where=np.array(W_div) >= np.array(W_sym),
                   alpha=0.2, color=C['stable'])
    ax.set_xlabel('Returns to investment ($\\alpha$)')
    ax.set_ylabel('Equilibrium fitness')
    ax.legend(loc='best', fontsize=8)
    ax.set_title('(C) Fitness advantage of specialization', fontweight='bold', loc='left')

    # Panel D: Basin of attraction size
    ax = axes[1, 1]

    # Simulate from many initial conditions
    n_samples = 200
    np.random.seed(42)
    converge_A1B2 = 0
    converge_A2B1 = 0

    model = CrossFeedingModel()
    for _ in range(n_samples):
        # Random initial condition
        f_A1 = np.random.uniform(0.1, 0.4)
        f_A2 = np.random.uniform(0.1, 0.4)
        f_B1 = np.random.uniform(0.1, 0.4)
        f_B2 = np.random.uniform(0.1, 0.4)

        t, sol = model.simulate([f_A1, f_A2, f_B1, f_B2], t_max=500, n_points=100)
        final = sol[-1]

        # Check which equilibrium
        if final[0] > final[1] and final[3] > final[2]:
            converge_A1B2 += 1
        elif final[1] > final[0] and final[2] > final[3]:
            converge_A2B1 += 1

    # Bar plot
    bars = ax.bar(['A→$F_1$, B→$F_2$', 'A→$F_2$, B→$F_1$'],
                  [converge_A1B2/n_samples*100, converge_A2B1/n_samples*100],
                  color=[C['stable'], C['B']], edgecolor='black')
    ax.set_ylabel('% of initial conditions')
    ax.set_title('(D) Basin of attraction', fontweight='bold', loc='left')
    ax.set_ylim(0, 60)

    plt.tight_layout()

    if save:
        fig_dir = ensure_dir()
        plt.savefig(os.path.join(fig_dir, 'fig4_stability.pdf'))
        plt.savefig(os.path.join(fig_dir, 'fig4_stability.png'))

    return fig


# =============================================================================
# FIGURE 5: Invasion Analysis
# =============================================================================

def figure5_invasion(save=True):
    """
    Invasion analysis and pairwise invasibility.
    """
    fig, axes = plt.subplots(1, 3, figsize=(11, 4))
    model = CrossFeedingModel()

    # Panel A: Invasion fitness landscape
    ax = axes[0]
    n = 40
    f_res = np.linspace(0.1, 0.5, n)
    f_mut = np.linspace(0.1, 0.5, n)
    R, M = np.meshgrid(f_res, f_mut)
    invasion = np.zeros_like(R)

    for i in range(n):
        for j in range(n):
            # Resident symmetric: f_A1 = f_A2 = f_res[j]/2
            # Mutant tries: f_A1 = f_mut[i]/2, f_A2 = f_mut[i]/2
            f_r = f_res[j]
            f_m = f_mut[i]

            # Environment set by resident
            P1_res = f_r  # f_A1 + f_B1 = f_r/2 + f_r/2
            P2_res = f_r

            g = model.growth_rate(P1_res, P2_res)
            W_mut = (1 - f_m) * g - model.params['D']
            W_res = (1 - f_r) * g - model.params['D']
            invasion[i, j] = W_mut - W_res

    im = ax.contourf(R, M, invasion, levels=20, cmap='RdBu_r', vmin=-0.05, vmax=0.05)
    ax.contour(R, M, invasion, levels=[0], colors='black', linewidths=2)
    ax.plot([0.1, 0.5], [0.1, 0.5], 'k--', lw=1)
    plt.colorbar(im, ax=ax, label='Invasion fitness')
    ax.set_xlabel('Resident investment')
    ax.set_ylabel('Mutant investment')
    ax.set_title('(A) Pairwise invasibility', fontweight='bold', loc='left')

    # Panel B: Selection gradient along evolutionary trajectory
    ax = axes[1]
    f_range = np.linspace(0.1, 0.5, 50)

    # Gradient when B is at f_B
    for f_B_total in [0.3, 0.4, 0.5]:
        grads = []
        for f_A_total in f_range:
            # A invests f_A_total/2 in each, B invests f_B_total/2 in each
            g = model.selection_gradient_A1(f_A_total/2, f_A_total/2,
                                            f_B_total/2, f_B_total/2)
            grads.append(g)
        ax.plot(f_range, grads, label=f'$f_B = {f_B_total}$', lw=1.5)

    ax.axhline(0, ls='--', color='gray', lw=0.5)
    ax.set_xlabel("A's total investment")
    ax.set_ylabel('Selection gradient $\\partial W_A / \\partial f_{A1}$')
    ax.legend(loc='best', fontsize=8)
    ax.set_title("(B) Selection depends on partner's strategy", fontweight='bold', loc='left')

    # Panel C: Cheater invasion
    ax = axes[2]

    # Cheater invasion fitness vs. resident investment
    f_range = np.linspace(0.1, 0.6, 50)
    cheater_fitness = []

    for f_total in f_range:
        # Resident state: both species at f_total/2 each
        state = [f_total/2, f_total/2, f_total/2, f_total/2]
        inv_fit, W_c, W_r = model.cheater_invasion_fitness(state)
        cheater_fitness.append(inv_fit)

    ax.plot(f_range, cheater_fitness, 'k-', lw=2)
    ax.axhline(0, ls='--', color='gray', lw=0.5)
    ax.fill_between(f_range, cheater_fitness, 0,
                   where=np.array(cheater_fitness) > 0,
                   alpha=0.3, color=C['unstable'], label='Cheater invades')
    ax.fill_between(f_range, cheater_fitness, 0,
                   where=np.array(cheater_fitness) < 0,
                   alpha=0.3, color=C['stable'], label='Cheater excluded')

    ax.axvline(0.5, ls=':', color='gray', alpha=0.5)
    ax.text(0.51, 0.02, 'Symmetric\neq.', fontsize=7)

    ax.set_xlabel('Resident total investment')
    ax.set_ylabel('Cheater invasion fitness')
    ax.legend(loc='best', fontsize=8)
    ax.set_title('(C) Cheater invasion analysis', fontweight='bold', loc='left')

    plt.tight_layout()

    if save:
        fig_dir = ensure_dir()
        plt.savefig(os.path.join(fig_dir, 'fig5_invasion.pdf'))
        plt.savefig(os.path.join(fig_dir, 'fig5_invasion.png'))

    return fig


# =============================================================================
# FIGURE 6: Robustness Analysis
# =============================================================================

def figure6_robustness(save=True):
    """
    Robustness to parameter variation and asymmetry.
    """
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))

    # Panel A: Effect of dilution rate D
    ax = axes[0, 0]
    D_range = np.linspace(0.01, 0.3, 30)
    converge_to_DOL = []

    for D in D_range:
        model = CrossFeedingModel({'gamma': 1.0, 'alpha': 1.0, 'D': D,
                                   'mu_A': 1.0, 'mu_B': 1.0, 'c_A': 1.0, 'c_B': 1.0,
                                   'eta_A1': 1.0, 'eta_A2': 1.0, 'eta_B1': 1.0, 'eta_B2': 1.0,
                                   'sigma': 0.01, 'delta': 0.0})
        t, sol = model.simulate([0.26, 0.24, 0.24, 0.26], t_max=500, n_points=100)
        final = sol[-1]

        # Check if converged to DOL
        spec_A = abs(final[0] - final[1]) / (final[0] + final[1] + 1e-10)
        spec_B = abs(final[3] - final[2]) / (final[2] + final[3] + 1e-10)
        converge_to_DOL.append((spec_A + spec_B) / 2)

    ax.plot(D_range, converge_to_DOL, 'k-', lw=2)
    ax.axhline(1, ls='--', color=C['stable'], alpha=0.5)
    ax.set_xlabel('Dilution rate $D$')
    ax.set_ylabel('Specialization index')
    ax.set_title('(A) Effect of dilution rate', fontweight='bold', loc='left')
    ax.set_ylim(0, 1.1)

    # Panel B: Effect of asymmetric growth rates
    ax = axes[0, 1]
    mu_ratio = np.linspace(0.5, 2.0, 30)
    f_A1_final = []
    f_B2_final = []

    for ratio in mu_ratio:
        model = CrossFeedingModel({'gamma': 1.0, 'alpha': 1.0, 'D': 0.1,
                                   'mu_A': 1.0, 'mu_B': ratio, 'c_A': 1.0, 'c_B': 1.0,
                                   'eta_A1': 1.0, 'eta_A2': 1.0, 'eta_B1': 1.0, 'eta_B2': 1.0,
                                   'sigma': 0.01, 'delta': 0.0})
        t, sol = model.simulate([0.26, 0.24, 0.24, 0.26], t_max=600, n_points=100)
        f_A1_final.append(sol[-1, 0])
        f_B2_final.append(sol[-1, 3])

    ax.plot(mu_ratio, f_A1_final, '-', color=C['A'], lw=2, label='$f_{A1}^*$')
    ax.plot(mu_ratio, f_B2_final, '-', color=C['B'], lw=2, label='$f_{B2}^*$')
    ax.axvline(1.0, ls=':', color='gray', alpha=0.5)
    ax.set_xlabel('Growth rate ratio $\\mu_B / \\mu_A$')
    ax.set_ylabel('Equilibrium investment')
    ax.legend(loc='best', fontsize=8)
    ax.set_title('(B) Asymmetric growth rates', fontweight='bold', loc='left')

    # Panel C: Effect of asymmetric production efficiency
    ax = axes[1, 0]
    eta_A1_range = np.linspace(0.5, 2.0, 30)
    f_A1_final = []
    f_A2_final = []

    for eta in eta_A1_range:
        model = CrossFeedingModel({'gamma': 1.0, 'alpha': 1.0, 'D': 0.1,
                                   'mu_A': 1.0, 'mu_B': 1.0, 'c_A': 1.0, 'c_B': 1.0,
                                   'eta_A1': eta, 'eta_A2': 1.0, 'eta_B1': 1.0, 'eta_B2': 1.0,
                                   'sigma': 0.01, 'delta': 0.0})
        t, sol = model.simulate([0.3, 0.2, 0.2, 0.3], t_max=600, n_points=100)
        f_A1_final.append(sol[-1, 0])
        f_A2_final.append(sol[-1, 1])

    ax.plot(eta_A1_range, f_A1_final, '-', color=C['F1'], lw=2, label='$f_{A1}^*$')
    ax.plot(eta_A1_range, f_A2_final, '-', color=C['F2'], lw=2, label='$f_{A2}^*$')
    ax.axvline(1.0, ls=':', color='gray', alpha=0.5)
    ax.set_xlabel("A's efficiency for $F_1$ ($\\eta_{A1}$)")
    ax.set_ylabel("A's equilibrium investment")
    ax.legend(loc='best', fontsize=8)
    ax.set_title('(C) Comparative advantage', fontweight='bold', loc='left')

    # Panel D: Specialization bonus effect
    ax = axes[1, 1]
    delta_range = np.linspace(0, 0.3, 30)
    spec_index = []

    for delta in delta_range:
        model = CrossFeedingModel({'gamma': 1.0, 'alpha': 1.0, 'D': 0.1,
                                   'mu_A': 1.0, 'mu_B': 1.0, 'c_A': 1.0, 'c_B': 1.0,
                                   'eta_A1': 1.0, 'eta_A2': 1.0, 'eta_B1': 1.0, 'eta_B2': 1.0,
                                   'sigma': 0.01, 'delta': delta})
        t, sol = model.simulate([0.26, 0.24, 0.24, 0.26], t_max=400, n_points=100)
        final = sol[-1]
        s_A = abs(final[0] - final[1]) / (final[0] + final[1] + 1e-10)
        s_B = abs(final[3] - final[2]) / (final[2] + final[3] + 1e-10)
        spec_index.append((s_A + s_B) / 2)

    ax.plot(delta_range, spec_index, 'k-', lw=2)
    ax.set_xlabel('Specialization bonus $\\delta$')
    ax.set_ylabel('Final specialization index')
    ax.set_title('(D) Effect of specialization bonus', fontweight='bold', loc='left')
    ax.set_ylim(0, 1.1)

    plt.tight_layout()

    if save:
        fig_dir = ensure_dir()
        plt.savefig(os.path.join(fig_dir, 'fig6_robustness.pdf'))
        plt.savefig(os.path.join(fig_dir, 'fig6_robustness.png'))

    return fig


# =============================================================================
# MAIN
# =============================================================================

def generate_all():
    """Generate all figures."""
    print("Generating comprehensive figures...")
    print()

    fig_dir = ensure_dir()
    print(f"Output: {fig_dir}")
    print()

    figs = [
        ("Figure 1: Model schematic", figure1_model_schematic),
        ("Figure 2: Selection gradients", figure2_selection_gradients),
        ("Figure 3: Phase portrait", figure3_phase_portrait),
        ("Figure 4: Stability analysis", figure4_stability),
        ("Figure 5: Invasion analysis", figure5_invasion),
        ("Figure 6: Robustness", figure6_robustness),
    ]

    for name, func in figs:
        print(f"{name}...")
        func(save=True)
        plt.close()
        print("  Done!")

    print()
    print("All figures generated!")


if __name__ == "__main__":
    generate_all()

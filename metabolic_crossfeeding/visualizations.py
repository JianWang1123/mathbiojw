"""
Visualization Module for Metabolic Cross-Feeding Model

Generates publication-quality figures demonstrating:
1. Phase portraits of evolutionary dynamics
2. Fitness landscapes
3. Time series of evolution toward division of labor
4. Stability analysis diagrams
5. Game-theoretic payoff matrices

Author: Jian Wang
Date: 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import os

# Import model
from crossfeeding_model import (
    AnalyticalCrossFeedingModel,
    FullChemostatModel
)

# Set publication-quality defaults
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.figsize': (10, 8),
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'text.usetex': False,
    'font.family': 'serif',
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
})

# Color schemes
COLORS = {
    'species_A': '#E64B35',  # Red
    'species_B': '#4DBBD5',  # Blue
    'F1': '#00A087',         # Teal
    'F2': '#3C5488',         # Dark blue
    'stable': '#00A087',     # Green
    'unstable': '#E64B35',   # Red
    'neutral': '#7E6148',    # Brown
}


def ensure_figures_dir():
    """Create figures directory if it doesn't exist."""
    fig_dir = os.path.join(os.path.dirname(__file__), 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    return fig_dir


# =============================================================================
# FIGURE 1: Model Schematic
# =============================================================================

def plot_model_schematic(save=True):
    """
    Create a schematic diagram of the cross-feeding model.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Draw bioreactor vessel
    vessel = plt.Rectangle((0.15, 0.15), 0.7, 0.6, fill=False,
                           edgecolor='black', linewidth=3)
    ax.add_patch(vessel)

    # Species A (red circles)
    for x, y in [(0.25, 0.5), (0.35, 0.35), (0.28, 0.65)]:
        circle = plt.Circle((x, y), 0.04, color=COLORS['species_A'],
                            alpha=0.8, label='Species A' if x == 0.25 else '')
        ax.add_patch(circle)
        ax.text(x, y, 'A', ha='center', va='center', fontsize=12,
               fontweight='bold', color='white')

    # Species B (blue circles)
    for x, y in [(0.65, 0.45), (0.72, 0.6), (0.58, 0.3)]:
        circle = plt.Circle((x, y), 0.04, color=COLORS['species_B'],
                            alpha=0.8, label='Species B' if x == 0.65 else '')
        ax.add_patch(circle)
        ax.text(x, y, 'B', ha='center', va='center', fontsize=12,
               fontweight='bold', color='white')

    # Amino acids F1 (small teal dots)
    for x, y in [(0.4, 0.55), (0.45, 0.4), (0.5, 0.65), (0.38, 0.3)]:
        ax.plot(x, y, 'o', markersize=8, color=COLORS['F1'], alpha=0.7)

    # Amino acids F2 (small dark blue dots)
    for x, y in [(0.55, 0.5), (0.48, 0.35), (0.6, 0.68), (0.52, 0.25)]:
        ax.plot(x, y, 's', markersize=8, color=COLORS['F2'], alpha=0.7)

    # Arrows showing production
    ax.annotate('', xy=(0.4, 0.5), xytext=(0.29, 0.5),
                arrowprops=dict(arrowstyle='->', color=COLORS['F1'], lw=2))
    ax.annotate('', xy=(0.55, 0.45), xytext=(0.61, 0.45),
                arrowprops=dict(arrowstyle='->', color=COLORS['F2'], lw=2))

    # Inflow arrow
    ax.annotate('', xy=(0.15, 0.45), xytext=(0.02, 0.45),
                arrowprops=dict(arrowstyle='->', color='gray', lw=3))
    ax.text(0.02, 0.52, 'Glucose\n(inflow)', fontsize=10, ha='left')

    # Outflow arrow
    ax.annotate('', xy=(0.98, 0.45), xytext=(0.85, 0.45),
                arrowprops=dict(arrowstyle='->', color='gray', lw=3))
    ax.text(0.92, 0.52, 'Dilution\n(outflow)', fontsize=10, ha='center')

    # Legend for amino acids
    ax.plot([], [], 'o', markersize=10, color=COLORS['F1'], label='Amino acid F1')
    ax.plot([], [], 's', markersize=10, color=COLORS['F2'], label='Amino acid F2')

    # Investment fractions annotation
    ax.text(0.25, 0.82, r'$f_{A,1}, f_{A,2}$', fontsize=12, ha='center',
           color=COLORS['species_A'])
    ax.text(0.75, 0.82, r'$f_{B,1}, f_{B,2}$', fontsize=12, ha='center',
           color=COLORS['species_B'])

    # Title
    ax.set_title('Cross-Feeding Model in Chemostat', fontsize=16, fontweight='bold')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')

    # Legend
    ax.legend(loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.02))

    if save:
        fig_dir = ensure_figures_dir()
        plt.savefig(os.path.join(fig_dir, 'fig1_model_schematic.png'))
        plt.savefig(os.path.join(fig_dir, 'fig1_model_schematic.pdf'))

    return fig, ax


# =============================================================================
# FIGURE 2: Fitness Landscape
# =============================================================================

def plot_fitness_landscape(save=True):
    """
    Plot fitness landscape showing direct costs and indirect benefits.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Parameters
    gamma = 1.0
    D = 0.1

    # Panel A: Direct cost (growth allocation)
    ax = axes[0]
    f = np.linspace(0, 0.5, 100)
    growth_alloc = 1 - 2*f  # Assuming f_A1 = f_A2 = f
    ax.plot(f, growth_alloc, 'k-', lw=2.5)
    ax.fill_between(f, 0, growth_alloc, alpha=0.3, color='gray')
    ax.set_xlabel(r'Investment fraction $f$')
    ax.set_ylabel('Growth allocation\n' + r'$(1 - f_{A,1} - f_{A,2})$')
    ax.set_title('(A) Direct Cost of Investment', fontweight='bold')
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 1.1)
    ax.axhline(0.5, ls='--', color='gray', alpha=0.5)
    ax.axvline(0.25, ls='--', color='gray', alpha=0.5)
    ax.text(0.26, 0.52, r'$f^* = 0.25$', fontsize=10)

    # Panel B: Indirect benefit (amino acid production)
    ax = axes[1]
    P = 2*f  # Total production when f_A1 = f_A2 = f_B1 = f_B2 = f
    g = gamma * P * P  # Growth rate = gamma * P1 * P2
    ax.plot(f, g, 'b-', lw=2.5, color=COLORS['F1'])
    ax.fill_between(f, 0, g, alpha=0.3, color=COLORS['F1'])
    ax.set_xlabel(r'Investment fraction $f$')
    ax.set_ylabel(r'Growth rate $g = \gamma P_1 P_2$')
    ax.set_title('(B) Indirect Benefit (Production)', fontweight='bold')
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 1.1)

    # Panel C: Net fitness
    ax = axes[2]
    W = growth_alloc * g - D  # Net fitness
    ax.plot(f, W, 'k-', lw=2.5)
    ax.fill_between(f, np.minimum(W, 0), W, where=W>0, alpha=0.3, color=COLORS['stable'])
    ax.fill_between(f, W, np.maximum(W, 0), where=W<0, alpha=0.3, color=COLORS['unstable'])

    # Find and mark optimum
    idx_max = np.argmax(W)
    f_opt = f[idx_max]
    W_opt = W[idx_max]
    ax.plot(f_opt, W_opt, 'ko', markersize=10)
    ax.annotate(f'Maximum\n$f^* = {f_opt:.2f}$', xy=(f_opt, W_opt),
               xytext=(f_opt+0.08, W_opt+0.05),
               fontsize=10, ha='left',
               arrowprops=dict(arrowstyle='->', color='black'))

    ax.axhline(0, ls='-', color='gray', alpha=0.5, lw=1)
    ax.set_xlabel(r'Investment fraction $f$')
    ax.set_ylabel(r'Net fitness $W = (1-2f) \cdot g - D$')
    ax.set_title('(C) Net Fitness Landscape', fontweight='bold')
    ax.set_xlim(0, 0.5)

    plt.tight_layout()

    if save:
        fig_dir = ensure_figures_dir()
        plt.savefig(os.path.join(fig_dir, 'fig2_fitness_landscape.png'))
        plt.savefig(os.path.join(fig_dir, 'fig2_fitness_landscape.pdf'))

    return fig, axes


# =============================================================================
# FIGURE 3: Phase Portrait of Evolutionary Dynamics
# =============================================================================

def plot_phase_portrait(save=True):
    """
    Plot phase portrait showing evolutionary trajectories.

    Reduced to 2D by fixing: f_A2 = f_B1 = 0 (complete specialization assumption)
    Variables: f_A1 (A's investment in F1) vs f_B2 (B's investment in F2)
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Model with specialization efficiency
    model = AnalyticalCrossFeedingModel(gamma=1.0, D=0.1, sigma=0.01, eta=1.2)

    # Panel A: Vector field
    ax = axes[0]

    n_grid = 20
    f_A1_range = np.linspace(0.05, 0.95, n_grid)
    f_B2_range = np.linspace(0.05, 0.95, n_grid)
    F_A1, F_B2 = np.meshgrid(f_A1_range, f_B2_range)

    # Compute vector field (assuming f_A2 = f_B1 = 0.1 for symmetry)
    df_A1 = np.zeros_like(F_A1)
    df_B2 = np.zeros_like(F_B2)

    for i in range(n_grid):
        for j in range(n_grid):
            f_A1 = F_A1[i, j]
            f_B2 = F_B2[i, j]
            f_A2 = 0.1  # Small residual investment
            f_B1 = 0.1

            grad_A1 = model.selection_gradient_A1(f_A1, f_A2, f_B1, f_B2)
            grad_B2 = model.selection_gradient_B2(f_A1, f_A2, f_B1, f_B2)

            # Replicator dynamics
            df_A1[i, j] = model.sigma * f_A1 * (1-f_A1) * grad_A1
            df_B2[i, j] = model.sigma * f_B2 * (1-f_B2) * grad_B2

    # Normalize for visualization
    magnitude = np.sqrt(df_A1**2 + df_B2**2)
    magnitude[magnitude == 0] = 1
    df_A1_norm = df_A1 / magnitude
    df_B2_norm = df_B2 / magnitude

    # Plot vector field
    ax.quiver(F_A1, F_B2, df_A1_norm, df_B2_norm, magnitude,
              cmap='viridis', alpha=0.7, scale=25)

    # Plot some trajectories
    initial_conditions = [
        (0.25, 0.25),  # Symmetric start
        (0.1, 0.1),
        (0.4, 0.1),
        (0.1, 0.4),
        (0.3, 0.3),
        (0.5, 0.5),
    ]

    colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(initial_conditions)))

    for ic, color in zip(initial_conditions, colors):
        f_A1_0, f_B2_0 = ic
        state0 = [f_A1_0, 0.1, 0.1, f_B2_0]
        t, sol = model.simulate_evolution(state0, (0, 500), n_points=500)
        ax.plot(sol[:, 0], sol[:, 3], '-', color=color, lw=1.5, alpha=0.8)
        ax.plot(sol[0, 0], sol[0, 3], 'o', color=color, markersize=6)
        ax.plot(sol[-1, 0], sol[-1, 3], 's', color=color, markersize=8)

    # Mark equilibria
    ax.plot(0.25, 0.25, 'ko', markersize=12, markerfacecolor='white',
           markeredgewidth=2, label='Symmetric eq. (unstable)')
    ax.plot(0.5, 0.5, 'k^', markersize=12, markerfacecolor=COLORS['stable'],
           markeredgewidth=2, label='Division of labor (stable)')

    ax.set_xlabel(r"Species A investment in F1 ($f_{A,1}$)", fontsize=12)
    ax.set_ylabel(r"Species B investment in F2 ($f_{B,2}$)", fontsize=12)
    ax.set_title('(A) Phase Portrait: Evolution of Specialization', fontweight='bold')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc='lower right', fontsize=9)

    # Add diagonal line
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, lw=1)
    ax.text(0.7, 0.65, 'Symmetric\nline', fontsize=9, alpha=0.5, rotation=45)

    # Panel B: Time series from symmetric start
    ax = axes[1]

    state0 = [0.26, 0.24, 0.24, 0.26]  # Small perturbation from symmetric
    t, sol = model.simulate_evolution(state0, (0, 800), n_points=1000)

    ax.plot(t, sol[:, 0], '-', color=COLORS['species_A'], lw=2,
           label=r'$f_{A,1}$ (A invests in F1)')
    ax.plot(t, sol[:, 1], '--', color=COLORS['species_A'], lw=2,
           label=r'$f_{A,2}$ (A invests in F2)')
    ax.plot(t, sol[:, 2], '--', color=COLORS['species_B'], lw=2,
           label=r'$f_{B,1}$ (B invests in F1)')
    ax.plot(t, sol[:, 3], '-', color=COLORS['species_B'], lw=2,
           label=r'$f_{B,2}$ (B invests in F2)')

    ax.axhline(0.25, ls=':', color='gray', alpha=0.5)
    ax.text(50, 0.27, 'Symmetric eq.', fontsize=9, alpha=0.7)

    ax.set_xlabel('Evolutionary time', fontsize=12)
    ax.set_ylabel('Investment fraction', fontsize=12)
    ax.set_title('(B) Evolution from Near-Symmetric State', fontweight='bold')
    ax.legend(loc='right', fontsize=9)
    ax.set_xlim(0, 800)
    ax.set_ylim(0, 1)

    plt.tight_layout()

    if save:
        fig_dir = ensure_figures_dir()
        plt.savefig(os.path.join(fig_dir, 'fig3_phase_portrait.png'))
        plt.savefig(os.path.join(fig_dir, 'fig3_phase_portrait.pdf'))

    return fig, axes


# =============================================================================
# FIGURE 4: Stability Analysis
# =============================================================================

def plot_stability_analysis(save=True):
    """
    Visualize stability of symmetric vs specialized equilibria.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: Eigenvalues as function of eta
    ax = axes[0]

    eta_range = np.linspace(1.0, 1.5, 50)
    eigenvalues_sym = []
    eigenvalues_div = []

    for eta in eta_range:
        model = AnalyticalCrossFeedingModel(gamma=1.0, D=0.1, sigma=0.01, eta=eta)

        # Symmetric equilibrium Jacobian (approximate)
        f_star = 0.25
        eps = 1e-4

        # Compute Jacobian numerically
        state_sym = [f_star, f_star, f_star, f_star]
        J_sym = np.zeros((4, 4))

        for i in range(4):
            state_plus = state_sym.copy()
            state_minus = state_sym.copy()
            state_plus[i] += eps
            state_minus[i] -= eps

            deriv_plus = model.evolutionary_dynamics(state_plus, 0)
            deriv_minus = model.evolutionary_dynamics(state_minus, 0)

            for j in range(4):
                J_sym[j, i] = (deriv_plus[j] - deriv_minus[j]) / (2*eps)

        eigs_sym = np.linalg.eigvals(J_sym)
        eigenvalues_sym.append(np.max(np.real(eigs_sym)))

        # Division of labor equilibrium
        state_div = [0.45, 0.05, 0.05, 0.45]
        J_div = np.zeros((4, 4))

        for i in range(4):
            state_plus = state_div.copy()
            state_minus = state_div.copy()
            state_plus[i] += eps
            state_minus[i] -= eps

            deriv_plus = model.evolutionary_dynamics(state_plus, 0)
            deriv_minus = model.evolutionary_dynamics(state_minus, 0)

            for j in range(4):
                J_div[j, i] = (deriv_plus[j] - deriv_minus[j]) / (2*eps)

        eigs_div = np.linalg.eigvals(J_div)
        eigenvalues_div.append(np.max(np.real(eigs_div)))

    ax.plot(eta_range, eigenvalues_sym, '-', color=COLORS['unstable'], lw=2.5,
           label='Symmetric equilibrium')
    ax.plot(eta_range, eigenvalues_div, '-', color=COLORS['stable'], lw=2.5,
           label='Division of labor')
    ax.axhline(0, ls='--', color='gray', lw=1)
    ax.fill_between(eta_range, 0, 0.02, alpha=0.2, color=COLORS['unstable'])
    ax.fill_between(eta_range, -0.02, 0, alpha=0.2, color=COLORS['stable'])

    ax.text(1.25, 0.008, 'Unstable\nregion', fontsize=9, ha='center')
    ax.text(1.25, -0.012, 'Stable\nregion', fontsize=9, ha='center')

    ax.set_xlabel(r'Specialization efficiency $\eta$', fontsize=12)
    ax.set_ylabel('Maximum eigenvalue', fontsize=12)
    ax.set_title('(A) Stability vs Specialization Efficiency', fontweight='bold')
    ax.legend(loc='best', fontsize=9)

    # Panel B: Bifurcation diagram
    ax = axes[1]

    eta_range = np.linspace(1.0, 1.5, 30)
    f_A1_eq = []
    f_A2_eq = []

    for eta in eta_range:
        model = AnalyticalCrossFeedingModel(gamma=1.0, D=0.1, sigma=0.01, eta=eta)

        # Simulate from near-symmetric initial condition
        state0 = [0.26, 0.24, 0.24, 0.26]
        t, sol = model.simulate_evolution(state0, (0, 1000), n_points=500)

        f_A1_eq.append(sol[-1, 0])
        f_A2_eq.append(sol[-1, 1])

    ax.plot(eta_range, f_A1_eq, 'o-', color=COLORS['species_A'], lw=2,
           markersize=4, label=r'$f_{A,1}^*$ (F1 investment)')
    ax.plot(eta_range, f_A2_eq, 's-', color=COLORS['species_B'], lw=2,
           markersize=4, label=r'$f_{A,2}^*$ (F2 investment)')
    ax.axhline(0.25, ls='--', color='gray', alpha=0.5)

    # Mark bifurcation point
    ax.axvline(1.0, ls=':', color='gray', alpha=0.5)
    ax.text(1.02, 0.1, 'Bifurcation\npoint', fontsize=9, alpha=0.7)

    ax.set_xlabel(r'Specialization efficiency $\eta$', fontsize=12)
    ax.set_ylabel('Equilibrium investment', fontsize=12)
    ax.set_title('(B) Bifurcation Diagram', fontweight='bold')
    ax.legend(loc='best', fontsize=9)
    ax.set_ylim(0, 0.6)

    # Panel C: Fitness at equilibrium
    ax = axes[2]

    gamma = 1.0
    D = 0.1

    eta_range = np.linspace(1.0, 1.5, 50)
    W_sym = []
    W_div = []

    for eta in eta_range:
        # Symmetric: f = 0.25 for all
        f = 0.25
        P = 0.5
        g = gamma * P * P
        W_s = (1 - 2*f) * g - D
        W_sym.append(W_s)

        # Division of labor: f_A1 = f_B2 = 0.5, f_A2 = f_B1 = 0
        f_div = 0.5
        P_div = 0.5
        g_div = gamma * P_div * P_div
        spec_bonus = (eta - 1) * f_div  # Specialization bonus
        W_d = (1 - f_div + spec_bonus) * g_div - D
        W_div.append(W_d)

    ax.plot(eta_range, W_sym, '-', color=COLORS['unstable'], lw=2.5,
           label='Symmetric')
    ax.plot(eta_range, W_div, '-', color=COLORS['stable'], lw=2.5,
           label='Division of labor')
    ax.fill_between(eta_range, W_sym, W_div, where=np.array(W_div)>np.array(W_sym),
                   alpha=0.3, color=COLORS['stable'])

    ax.set_xlabel(r'Specialization efficiency $\eta$', fontsize=12)
    ax.set_ylabel('Equilibrium fitness $W^*$', fontsize=12)
    ax.set_title('(C) Fitness Comparison', fontweight='bold')
    ax.legend(loc='best', fontsize=9)

    plt.tight_layout()

    if save:
        fig_dir = ensure_figures_dir()
        plt.savefig(os.path.join(fig_dir, 'fig4_stability_analysis.png'))
        plt.savefig(os.path.join(fig_dir, 'fig4_stability_analysis.pdf'))

    return fig, axes


# =============================================================================
# FIGURE 5: Game Theory Payoff Matrix
# =============================================================================

def plot_game_theory(save=True):
    """
    Visualize the game-theoretic structure and Nash equilibria.
    """
    fig = plt.figure(figsize=(14, 5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.2, 1, 1])

    # Parameters
    gamma = 1.0
    eta = 1.2

    def payoff(strategy_A, strategy_B):
        strats = {
            'F1': (0.5, 0.0),
            'F2': (0.0, 0.5),
            'G': (0.25, 0.25)
        }

        f_A1, f_A2 = strats[strategy_A]
        f_B1, f_B2 = strats[strategy_B]

        P1 = f_A1 + f_B1
        P2 = f_A2 + f_B2

        if P1 < 0.01 or P2 < 0.01:
            return -1, -1

        g = gamma * np.sqrt(P1) * np.sqrt(P2)

        spec_A = abs(f_A1 - f_A2) / (f_A1 + f_A2 + 0.01)
        spec_B = abs(f_B1 - f_B2) / (f_B1 + f_B2 + 0.01)

        W_A = (1 - f_A1 - f_A2 + (eta-1)*spec_A*(f_A1+f_A2)) * g
        W_B = (1 - f_B1 - f_B2 + (eta-1)*spec_B*(f_B1+f_B2)) * g

        return W_A, W_B

    strategies = ['F1', 'F2', 'G']
    strategy_labels = ['Specialize\nF1', 'Specialize\nF2', 'Generalist']

    # Panel A: Payoff matrix heatmap
    ax = fig.add_subplot(gs[0])

    payoff_A = np.zeros((3, 3))
    payoff_B = np.zeros((3, 3))

    for i, sA in enumerate(strategies):
        for j, sB in enumerate(strategies):
            pA, pB = payoff(sA, sB)
            payoff_A[i, j] = pA
            payoff_B[i, j] = pB

    im = ax.imshow(payoff_A, cmap='RdYlGn', aspect='equal', vmin=-0.5, vmax=0.5)

    # Add text annotations
    for i in range(3):
        for j in range(3):
            text_color = 'white' if abs(payoff_A[i,j]) > 0.3 else 'black'
            ax.text(j, i, f'{payoff_A[i,j]:.2f}\n({payoff_B[i,j]:.2f})',
                   ha='center', va='center', fontsize=10, color=text_color)

    # Mark Nash equilibria
    nash_eq = [(0, 1), (1, 0)]  # (F1, F2) and (F2, F1)
    for i, j in nash_eq:
        rect = plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False,
                             edgecolor='gold', linewidth=4)
        ax.add_patch(rect)

    ax.set_xticks(range(3))
    ax.set_yticks(range(3))
    ax.set_xticklabels(strategy_labels, fontsize=10)
    ax.set_yticklabels(strategy_labels, fontsize=10)
    ax.set_xlabel('Species B Strategy', fontsize=12)
    ax.set_ylabel('Species A Strategy', fontsize=12)
    ax.set_title('(A) Payoff Matrix\n(Species A, Species B)', fontweight='bold')

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Payoff', fontsize=10)

    # Panel B: Best response dynamics
    ax = fig.add_subplot(gs[1])

    # Continuous strategy space
    f_A_range = np.linspace(0, 1, 50)  # f_A1 (with f_A2 = 1 - f_A1 scaled)
    f_B_range = np.linspace(0, 1, 50)

    best_response_A = []
    best_response_B = []

    for f_B in f_B_range:
        # Find best response for A given B's strategy
        best_W_A = -np.inf
        best_f_A = 0
        for f_A in np.linspace(0, 1, 100):
            f_A1 = f_A * 0.5
            f_A2 = (1-f_A) * 0.5
            f_B1 = (1-f_B) * 0.5
            f_B2 = f_B * 0.5

            P1 = f_A1 + f_B1
            P2 = f_A2 + f_B2

            if P1 > 0.01 and P2 > 0.01:
                g = gamma * np.sqrt(P1) * np.sqrt(P2)
                spec_A = abs(f_A1 - f_A2) / (f_A1 + f_A2 + 0.01)
                W_A = (1 - f_A1 - f_A2 + (eta-1)*spec_A*(f_A1+f_A2)) * g

                if W_A > best_W_A:
                    best_W_A = W_A
                    best_f_A = f_A

        best_response_A.append(best_f_A)

    for f_A in f_A_range:
        # Find best response for B given A's strategy
        best_W_B = -np.inf
        best_f_B = 0
        for f_B in np.linspace(0, 1, 100):
            f_A1 = f_A * 0.5
            f_A2 = (1-f_A) * 0.5
            f_B1 = (1-f_B) * 0.5
            f_B2 = f_B * 0.5

            P1 = f_A1 + f_B1
            P2 = f_A2 + f_B2

            if P1 > 0.01 and P2 > 0.01:
                g = gamma * np.sqrt(P1) * np.sqrt(P2)
                spec_B = abs(f_B1 - f_B2) / (f_B1 + f_B2 + 0.01)
                W_B = (1 - f_B1 - f_B2 + (eta-1)*spec_B*(f_B1+f_B2)) * g

                if W_B > best_W_B:
                    best_W_B = W_B
                    best_f_B = f_B

        best_response_B.append(best_f_B)

    ax.plot(f_B_range, best_response_A, '-', color=COLORS['species_A'], lw=2.5,
           label="A's best response")
    ax.plot(best_response_B, f_A_range, '-', color=COLORS['species_B'], lw=2.5,
           label="B's best response")

    # Mark Nash equilibria (intersections)
    ax.plot([0, 1], [1, 0], 'ko', markersize=12, markerfacecolor='gold',
           markeredgewidth=2, zorder=10)
    ax.annotate('Nash\nequilibria', xy=(0.1, 0.9), xytext=(0.3, 0.95),
               fontsize=10, ha='center',
               arrowprops=dict(arrowstyle='->', color='gray'))

    ax.set_xlabel(r'B specialization on F2 ($f_B$)', fontsize=12)
    ax.set_ylabel(r'A specialization on F1 ($f_A$)', fontsize=12)
    ax.set_title('(B) Best Response Functions', fontweight='bold')
    ax.legend(loc='center', fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Panel C: Coordination diagram
    ax = fig.add_subplot(gs[2])

    # Draw coordination game outcome space
    outcomes = {
        'Both\nGeneralist': (0.5, 0.5, 0.3, 'gray'),
        'A: F1\nB: F2': (0.15, 0.85, 0.45, COLORS['stable']),
        'A: F2\nB: F1': (0.85, 0.15, 0.45, COLORS['stable']),
        'Both\nF1': (0.15, 0.15, 0.1, COLORS['unstable']),
        'Both\nF2': (0.85, 0.85, 0.1, COLORS['unstable']),
    }

    for label, (x, y, size, color) in outcomes.items():
        circle = plt.Circle((x, y), size/3, color=color, alpha=0.6)
        ax.add_patch(circle)
        ax.text(x, y, label, ha='center', va='center', fontsize=9,
               fontweight='bold' if 'Nash' in label or size > 0.4 else 'normal')

    # Arrows showing evolutionary flow
    ax.annotate('', xy=(0.2, 0.8), xytext=(0.45, 0.55),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
    ax.annotate('', xy=(0.8, 0.2), xytext=(0.55, 0.45),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

    ax.text(0.5, 0.02, 'Evolutionary\nattractors', fontsize=10, ha='center',
           style='italic')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('(C) Coordination Game Outcome', fontweight='bold')

    plt.tight_layout()

    if save:
        fig_dir = ensure_figures_dir()
        plt.savefig(os.path.join(fig_dir, 'fig5_game_theory.png'))
        plt.savefig(os.path.join(fig_dir, 'fig5_game_theory.pdf'))

    return fig


# =============================================================================
# FIGURE 6: Full Eco-Evolutionary Dynamics
# =============================================================================

def plot_eco_evo_dynamics(save=True):
    """
    Simulate and plot full ecological-evolutionary dynamics.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Initialize model
    model = FullChemostatModel()
    model.params['sigma'] = 0.005  # Slower evolution for clarity
    model.params['eta'] = 1.2

    # Initial state: near-symmetric with small perturbation
    # [N_A, N_B, F1, F2, S, f_A1, f_A2, f_B1, f_B2]
    initial_state = [0.1, 0.1, 10.0, 10.0, 1000.0, 0.26, 0.24, 0.24, 0.26]

    # Simulate
    t, sol = model.simulate_eco_evo(initial_state, (0, 2000), n_points=2000)

    N_A, N_B = sol[:, 0], sol[:, 1]
    F1, F2 = sol[:, 2], sol[:, 3]
    S = sol[:, 4]
    f_A1, f_A2 = sol[:, 5], sol[:, 6]
    f_B1, f_B2 = sol[:, 7], sol[:, 8]

    # Panel A: Population dynamics
    ax = axes[0, 0]
    ax.plot(t, N_A, '-', color=COLORS['species_A'], lw=2, label='Species A')
    ax.plot(t, N_B, '-', color=COLORS['species_B'], lw=2, label='Species B')
    ax.set_xlabel('Time (hours)', fontsize=12)
    ax.set_ylabel('Population density (OD)', fontsize=12)
    ax.set_title('(A) Population Dynamics', fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.set_xlim(0, 2000)

    # Panel B: Resource dynamics
    ax = axes[0, 1]
    ax.plot(t, F1, '-', color=COLORS['F1'], lw=2, label='Amino acid F1')
    ax.plot(t, F2, '-', color=COLORS['F2'], lw=2, label='Amino acid F2')
    ax.set_xlabel('Time (hours)', fontsize=12)
    ax.set_ylabel('Concentration (µM)', fontsize=12)
    ax.set_title('(B) Amino Acid Dynamics', fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.set_xlim(0, 2000)

    # Panel C: Investment strategy evolution
    ax = axes[1, 0]
    ax.plot(t, f_A1, '-', color=COLORS['species_A'], lw=2, label=r'$f_{A,1}$')
    ax.plot(t, f_A2, '--', color=COLORS['species_A'], lw=2, label=r'$f_{A,2}$')
    ax.plot(t, f_B1, '--', color=COLORS['species_B'], lw=2, label=r'$f_{B,1}$')
    ax.plot(t, f_B2, '-', color=COLORS['species_B'], lw=2, label=r'$f_{B,2}$')

    ax.axhline(0.25, ls=':', color='gray', alpha=0.5, lw=1)
    ax.text(100, 0.27, 'Symmetric', fontsize=9, alpha=0.7)

    ax.set_xlabel('Time (hours)', fontsize=12)
    ax.set_ylabel('Investment fraction', fontsize=12)
    ax.set_title('(C) Evolution of Investment Strategies', fontweight='bold')
    ax.legend(loc='right', fontsize=9, ncol=2)
    ax.set_xlim(0, 2000)
    ax.set_ylim(0, 0.6)

    # Panel D: Specialization index over time
    ax = axes[1, 1]

    spec_A = np.abs(f_A1 - f_A2) / (f_A1 + f_A2 + 1e-10)
    spec_B = np.abs(f_B1 - f_B2) / (f_B1 + f_B2 + 1e-10)

    ax.plot(t, spec_A, '-', color=COLORS['species_A'], lw=2,
           label='Species A specialization')
    ax.plot(t, spec_B, '-', color=COLORS['species_B'], lw=2,
           label='Species B specialization')

    ax.axhline(1.0, ls='--', color='gray', alpha=0.5, lw=1)
    ax.text(100, 1.02, 'Complete\nspecialization', fontsize=9, alpha=0.7)

    ax.set_xlabel('Time (hours)', fontsize=12)
    ax.set_ylabel('Specialization index\n' + r'$|f_1 - f_2| / (f_1 + f_2)$', fontsize=12)
    ax.set_title('(D) Division of Labor Index', fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.set_xlim(0, 2000)
    ax.set_ylim(0, 1.1)

    plt.tight_layout()

    if save:
        fig_dir = ensure_figures_dir()
        plt.savefig(os.path.join(fig_dir, 'fig6_eco_evo_dynamics.png'))
        plt.savefig(os.path.join(fig_dir, 'fig6_eco_evo_dynamics.pdf'))

    return fig, axes


# =============================================================================
# FIGURE 7: Summary Diagram - Why Division of Labor Evolves
# =============================================================================

def plot_summary_diagram(save=True):
    """
    Create a summary figure explaining why division of labor evolves.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # Panel A: Trade-off between investment and growth
    ax = axes[0]

    f = np.linspace(0, 0.5, 100)

    # Generalist: invests equally in both
    growth_gen = 1 - 2*f
    production_gen = 2*f

    # Specialist: invests all in one
    growth_spec = 1 - f
    production_spec = f

    ax.plot(production_gen, growth_gen, '-', color='gray', lw=3,
           label='Generalist path')
    ax.plot(production_spec, growth_spec, '-', color=COLORS['stable'], lw=3,
           label='Specialist path')

    # Mark key points
    ax.plot(0.5, 0.5, 'ko', markersize=12, label='Symmetric eq.')
    ax.plot(0.5, 0.5, 'k^', markersize=12, markerfacecolor=COLORS['stable'],
           label='Division of labor eq.')

    ax.set_xlabel('Total amino acid investment', fontsize=12)
    ax.set_ylabel('Growth allocation', fontsize=12)
    ax.set_title('(A) Investment-Growth Trade-off', fontweight='bold')
    ax.legend(loc='best', fontsize=9)

    # Panel B: Public goods dilemma
    ax = axes[1]

    # Illustrate the public goods nature
    theta = np.linspace(0, 2*np.pi, 100)
    r = 0.3

    # Species A circle
    ax.add_patch(plt.Circle((0.3, 0.5), r, color=COLORS['species_A'], alpha=0.3))
    ax.text(0.3, 0.5, 'A', fontsize=20, ha='center', va='center',
           fontweight='bold', color=COLORS['species_A'])

    # Species B circle
    ax.add_patch(plt.Circle((0.7, 0.5), r, color=COLORS['species_B'], alpha=0.3))
    ax.text(0.7, 0.5, 'B', fontsize=20, ha='center', va='center',
           fontweight='bold', color=COLORS['species_B'])

    # Overlap region (public goods pool)
    overlap = plt.Circle((0.5, 0.5), 0.15, color='gold', alpha=0.5)
    ax.add_patch(overlap)
    ax.text(0.5, 0.5, 'F1+F2', fontsize=10, ha='center', va='center')

    # Arrows showing production
    ax.annotate('', xy=(0.42, 0.55), xytext=(0.32, 0.55),
                arrowprops=dict(arrowstyle='->', color=COLORS['F1'], lw=2))
    ax.annotate('', xy=(0.58, 0.45), xytext=(0.68, 0.45),
                arrowprops=dict(arrowstyle='->', color=COLORS['F2'], lw=2))

    ax.text(0.25, 0.7, 'Produces F1', fontsize=10, ha='center', color=COLORS['F1'])
    ax.text(0.75, 0.3, 'Produces F2', fontsize=10, ha='center', color=COLORS['F2'])

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('(B) Mutualistic Cross-Feeding', fontweight='bold')

    # Panel C: Evolutionary advantage of specialization
    ax = axes[2]

    categories = ['Generalist\n(baseline)', 'Specialist\n(with partner)']
    width = 0.35

    # Fitness components
    growth_alloc = [0.5, 0.5]  # Same for both when total investment equal
    efficiency = [1.0, 1.2]   # Specialization bonus
    total_fitness = [g * e for g, e in zip(growth_alloc, efficiency)]

    x = np.arange(len(categories))

    bars1 = ax.bar(x - width/2, growth_alloc, width, label='Growth allocation',
                   color='lightgray', edgecolor='black')
    bars2 = ax.bar(x + width/2, total_fitness, width, label='Effective fitness',
                   color=COLORS['stable'], edgecolor='black')

    ax.set_ylabel('Relative value', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=11)
    ax.set_title('(C) Fitness Advantage of Specialization', fontweight='bold')
    ax.legend(loc='upper left', fontsize=9)
    ax.set_ylim(0, 0.8)

    # Add annotation
    ax.annotate('', xy=(1.17, 0.6), xytext=(0.17, 0.5),
                arrowprops=dict(arrowstyle='->', color='black', lw=2))
    ax.text(0.67, 0.65, '+20%', fontsize=12, fontweight='bold', color=COLORS['stable'])

    plt.tight_layout()

    if save:
        fig_dir = ensure_figures_dir()
        plt.savefig(os.path.join(fig_dir, 'fig7_summary.png'))
        plt.savefig(os.path.join(fig_dir, 'fig7_summary.pdf'))

    return fig, axes


# =============================================================================
# MAIN FUNCTION: Generate All Figures
# =============================================================================

def generate_all_figures():
    """
    Generate all figures for the paper.
    """
    print("="*60)
    print("Generating all figures for Cross-Feeding Model")
    print("="*60)

    fig_dir = ensure_figures_dir()
    print(f"\nFigures will be saved to: {fig_dir}\n")

    print("Figure 1: Model schematic...")
    plot_model_schematic(save=True)
    plt.close()
    print("  Done!")

    print("Figure 2: Fitness landscape...")
    plot_fitness_landscape(save=True)
    plt.close()
    print("  Done!")

    print("Figure 3: Phase portrait...")
    plot_phase_portrait(save=True)
    plt.close()
    print("  Done!")

    print("Figure 4: Stability analysis...")
    plot_stability_analysis(save=True)
    plt.close()
    print("  Done!")

    print("Figure 5: Game theory...")
    plot_game_theory(save=True)
    plt.close()
    print("  Done!")

    print("Figure 6: Eco-evolutionary dynamics...")
    plot_eco_evo_dynamics(save=True)
    plt.close()
    print("  Done!")

    print("Figure 7: Summary diagram...")
    plot_summary_diagram(save=True)
    plt.close()
    print("  Done!")

    print("\n" + "="*60)
    print("All figures generated successfully!")
    print("="*60)


if __name__ == "__main__":
    generate_all_figures()

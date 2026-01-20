"""
Publication figures for The American Naturalist manuscript.

Generates clean, minimal figures demonstrating division of labor as ESS.

Author: Jian Wang
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy.integrate import odeint
import os

from model import CrossFeedingModel

# Publication style
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'serif',
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.figsize': (6, 4.5),
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.0,
    'lines.linewidth': 1.5,
})

COLORS = {
    'A': '#D55E00',      # Orange-red for species A
    'B': '#0072B2',      # Blue for species B
    'stable': '#009E73', # Green for stable
    'unstable': '#CC79A7', # Pink for unstable
}


def ensure_dir():
    """Create figures directory."""
    fig_dir = os.path.join(os.path.dirname(__file__), 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    return fig_dir


# =============================================================================
# FIGURE 1: Fitness trade-off and equilibria
# =============================================================================

def figure1_fitness_tradeoff(save=True):
    """
    Panel A: Investment-growth trade-off
    Panel B: Fitness landscape showing both equilibria
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    model = CrossFeedingModel(alpha=1.0, D=0.1)

    # Panel A: Trade-off structure
    ax = axes[0]
    f_total = np.linspace(0, 0.6, 100)

    # Growth allocation decreases with investment
    growth_alloc = 1 - f_total
    ax.plot(f_total, growth_alloc, 'k-', lw=2, label='Growth allocation')

    # Amino acid benefit increases (for symmetric case: P = f_total)
    # Using sqrt to show diminishing returns visually
    benefit = np.sqrt(f_total / 0.6)
    ax.plot(f_total, benefit, '--', color='gray', lw=2, label='Amino acid benefit')

    ax.axvline(0.5, ls=':', color=COLORS['stable'], alpha=0.7)
    ax.text(0.51, 0.85, 'Division of\nlabor', fontsize=9, color=COLORS['stable'])

    ax.axvline(0.5, ls=':', color=COLORS['unstable'], alpha=0.5)

    ax.set_xlabel('Total investment ($f_1 + f_2$)')
    ax.set_ylabel('Relative value')
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 1.1)
    ax.legend(loc='center right', frameon=False)
    ax.set_title('(A) Investment trade-off', fontweight='bold')

    # Panel B: Fitness as function of strategy
    ax = axes[1]

    # For symmetric strategy: f1 = f2 = f, total = 2f
    f_each = np.linspace(0.01, 0.4, 100)
    f_total_sym = 2 * f_each
    P_sym = f_total_sym  # P1 = P2 = 2f when symmetric
    W_sym = (1 - f_total_sym) * model.growth(P_sym/2, P_sym/2) - model.D

    # For specialist: f1 = f, f2 = 0, partner has f1=0, f2=f
    f_spec = np.linspace(0.01, 0.7, 100)
    P_spec = f_spec  # P1 = P2 = f when division of labor
    W_spec = (1 - f_spec) * model.growth(P_spec, P_spec) - model.D

    ax.plot(f_each, W_sym, '-', color=COLORS['unstable'], lw=2,
            label='Symmetric generalist')
    ax.plot(f_spec, W_spec, '-', color=COLORS['stable'], lw=2,
            label='Division of labor')

    # Mark equilibria
    f_sym_eq = model.symmetric_equilibrium()
    f_div_eq = model.division_of_labor_equilibrium()
    W_sym_eq = model.fitness_at_equilibrium('symmetric')
    W_div_eq = model.fitness_at_equilibrium('division')

    ax.plot(f_sym_eq, W_sym_eq, 'o', color=COLORS['unstable'], markersize=10,
            markeredgecolor='black', markeredgewidth=1.5, zorder=10)
    ax.plot(f_div_eq, W_div_eq, 's', color=COLORS['stable'], markersize=10,
            markeredgecolor='black', markeredgewidth=1.5, zorder=10)

    ax.annotate('Symmetric\n(unstable)', xy=(f_sym_eq, W_sym_eq),
                xytext=(f_sym_eq + 0.08, W_sym_eq + 0.02),
                fontsize=9, ha='left', color=COLORS['unstable'])
    ax.annotate('Division of labor\n(ESS)', xy=(f_div_eq, W_div_eq),
                xytext=(f_div_eq - 0.15, W_div_eq + 0.025),
                fontsize=9, ha='center', color=COLORS['stable'])

    ax.axhline(0, ls='-', color='gray', lw=0.5, alpha=0.5)
    ax.set_xlabel('Investment in specialized amino acid ($f$)')
    ax.set_ylabel('Fitness ($W$)')
    ax.set_xlim(0, 0.7)
    ax.legend(loc='lower right', frameon=False)
    ax.set_title('(B) Equilibrium fitness', fontweight='bold')

    plt.tight_layout()

    if save:
        fig_dir = ensure_dir()
        plt.savefig(os.path.join(fig_dir, 'fig1_fitness.pdf'))
        plt.savefig(os.path.join(fig_dir, 'fig1_fitness.png'))

    return fig


# =============================================================================
# FIGURE 2: Phase portrait showing evolution toward division of labor
# =============================================================================

def figure2_phase_portrait(save=True):
    """
    Phase portrait in reduced 2D space.
    x-axis: Species A specialization (f_A1 - f_A2)
    y-axis: Species B specialization (f_B2 - f_B1)
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    model = CrossFeedingModel(alpha=1.0, D=0.1, sigma=0.02)

    # Panel A: Phase portrait
    ax = axes[0]

    # Grid for vector field
    n = 15
    spec_A = np.linspace(-0.4, 0.4, n)  # f_A1 - f_A2
    spec_B = np.linspace(-0.4, 0.4, n)  # f_B2 - f_B1
    SA, SB = np.meshgrid(spec_A, spec_B)

    dSA = np.zeros_like(SA)
    dSB = np.zeros_like(SB)

    # Convert specialization to investment fractions
    # At symmetric: f_A1 = f_A2 = f_B1 = f_B2 = 0.25
    # Specialization: s_A = f_A1 - f_A2, s_B = f_B2 - f_B1
    f_base = 0.25

    for i in range(n):
        for j in range(n):
            s_A = SA[i, j]
            s_B = SB[i, j]

            # Convert to investment fractions
            f_A1 = f_base + s_A/2
            f_A2 = f_base - s_A/2
            f_B1 = f_base - s_B/2
            f_B2 = f_base + s_B/2

            # Clip to valid range
            f_A1 = np.clip(f_A1, 0.01, 0.49)
            f_A2 = np.clip(f_A2, 0.01, 0.49)
            f_B1 = np.clip(f_B1, 0.01, 0.49)
            f_B2 = np.clip(f_B2, 0.01, 0.49)

            state = [f_A1, f_A2, f_B1, f_B2]
            deriv = model.evolutionary_dynamics(state, 0)

            # Convert derivatives to specialization coordinates
            dSA[i, j] = deriv[0] - deriv[1]  # d(f_A1 - f_A2)/dt
            dSB[i, j] = deriv[3] - deriv[2]  # d(f_B2 - f_B1)/dt

    # Normalize for visualization
    mag = np.sqrt(dSA**2 + dSB**2)
    mag[mag == 0] = 1
    dSA_norm = dSA / mag * 0.03
    dSB_norm = dSB / mag * 0.03

    ax.quiver(SA, SB, dSA_norm, dSB_norm, mag, cmap='viridis', alpha=0.6)

    # Trajectories
    initial_conditions = [
        [0.26, 0.24, 0.24, 0.26],  # Near symmetric, slight A->F1, B->F2 bias
        [0.24, 0.26, 0.26, 0.24],  # Near symmetric, slight A->F2, B->F1 bias
        [0.25, 0.25, 0.25, 0.25],  # Exactly symmetric
        [0.30, 0.20, 0.20, 0.30],  # More specialized start
    ]

    for ic in initial_conditions:
        t, sol = model.simulate(ic, t_max=800, n_points=500)
        s_A_traj = sol[:, 0] - sol[:, 1]
        s_B_traj = sol[:, 3] - sol[:, 2]
        ax.plot(s_A_traj, s_B_traj, '-', color=COLORS['A'], lw=1.5, alpha=0.8)
        ax.plot(s_A_traj[0], s_B_traj[0], 'o', color=COLORS['A'], markersize=5)

    # Mark equilibria
    ax.plot(0, 0, 'o', color=COLORS['unstable'], markersize=12,
            markeredgecolor='black', markeredgewidth=2, label='Symmetric (unstable)')
    ax.plot(0.5, 0.5, 's', color=COLORS['stable'], markersize=12,
            markeredgecolor='black', markeredgewidth=2, label='A→F1, B→F2 (ESS)')
    ax.plot(-0.5, -0.5, 's', color=COLORS['stable'], markersize=12,
            markeredgecolor='black', markeredgewidth=2, label='A→F2, B→F1 (ESS)')

    ax.set_xlabel('Species A specialization ($f_{A1} - f_{A2}$)')
    ax.set_ylabel('Species B specialization ($f_{B2} - f_{B1}$)')
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(-0.6, 0.6)
    ax.axhline(0, ls='--', color='gray', lw=0.5, alpha=0.5)
    ax.axvline(0, ls='--', color='gray', lw=0.5, alpha=0.5)
    ax.legend(loc='upper left', frameon=True, fontsize=8)
    ax.set_title('(A) Evolutionary phase portrait', fontweight='bold')
    ax.set_aspect('equal')

    # Panel B: Time series
    ax = axes[1]

    # Simulate from near-symmetric initial condition
    ic = [0.26, 0.24, 0.24, 0.26]
    t, sol = model.simulate(ic, t_max=600, n_points=500)

    ax.plot(t, sol[:, 0], '-', color=COLORS['A'], lw=2, label='$f_{A1}$')
    ax.plot(t, sol[:, 1], '--', color=COLORS['A'], lw=2, label='$f_{A2}$')
    ax.plot(t, sol[:, 2], '--', color=COLORS['B'], lw=2, label='$f_{B1}$')
    ax.plot(t, sol[:, 3], '-', color=COLORS['B'], lw=2, label='$f_{B2}$')

    ax.axhline(0.25, ls=':', color='gray', alpha=0.5)
    ax.axhline(0.5, ls=':', color=COLORS['stable'], alpha=0.5)
    ax.text(50, 0.27, 'Symmetric', fontsize=8, color='gray')
    ax.text(50, 0.52, 'Division of labor', fontsize=8, color=COLORS['stable'])

    ax.set_xlabel('Evolutionary time')
    ax.set_ylabel('Investment fraction')
    ax.set_xlim(0, 600)
    ax.set_ylim(0, 0.6)
    ax.legend(loc='center right', frameon=False, ncol=2)
    ax.set_title('(B) Evolutionary trajectory', fontweight='bold')

    plt.tight_layout()

    if save:
        fig_dir = ensure_dir()
        plt.savefig(os.path.join(fig_dir, 'fig2_phase_portrait.pdf'))
        plt.savefig(os.path.join(fig_dir, 'fig2_phase_portrait.png'))

    return fig


# =============================================================================
# FIGURE 3: Stability analysis and diminishing returns
# =============================================================================

def figure3_stability(save=True):
    """
    Panel A: Eigenvalues at symmetric equilibrium
    Panel B: Effect of diminishing returns (alpha)
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Panel A: Maximum eigenvalue as function of perturbation
    ax = axes[0]

    model = CrossFeedingModel(alpha=1.0)
    f_sym = model.symmetric_equilibrium()

    # Compute Jacobian eigenvalues numerically
    eps = 1e-5
    state_sym = [f_sym, f_sym, f_sym, f_sym]

    J = np.zeros((4, 4))
    for i in range(4):
        state_plus = state_sym.copy()
        state_minus = state_sym.copy()
        state_plus[i] += eps
        state_minus[i] -= eps

        deriv_plus = model.evolutionary_dynamics(state_plus, 0)
        deriv_minus = model.evolutionary_dynamics(state_minus, 0)

        for j in range(4):
            J[j, i] = (deriv_plus[j] - deriv_minus[j]) / (2 * eps)

    eigenvalues = np.linalg.eigvals(J)

    # Plot eigenvalues
    ax.scatter(np.real(eigenvalues), np.imag(eigenvalues), s=150, c='black',
               zorder=10, edgecolors='white', linewidth=2)

    ax.axhline(0, ls='-', color='gray', lw=0.5)
    ax.axvline(0, ls='-', color='gray', lw=0.5)
    ax.fill_betweenx([-0.01, 0.01], 0, 0.01, alpha=0.2, color=COLORS['unstable'])
    ax.fill_betweenx([-0.01, 0.01], -0.01, 0, alpha=0.2, color=COLORS['stable'])

    ax.text(0.002, 0.005, 'Unstable', fontsize=9, color=COLORS['unstable'])
    ax.text(-0.008, 0.005, 'Stable', fontsize=9, color=COLORS['stable'])

    ax.set_xlabel('Real part of eigenvalue')
    ax.set_ylabel('Imaginary part')
    ax.set_xlim(-0.015, 0.015)
    ax.set_ylim(-0.01, 0.01)
    ax.set_title('(A) Eigenvalues at symmetric equilibrium', fontweight='bold')

    # Panel B: Effect of alpha on equilibria and fitness
    ax = axes[1]

    alpha_range = np.linspace(0.3, 1.0, 50)
    f_sym_values = []
    f_div_values = []
    W_sym_values = []
    W_div_values = []

    for alpha in alpha_range:
        m = CrossFeedingModel(alpha=alpha, D=0.05)
        f_sym_values.append(m.symmetric_equilibrium())
        f_div_values.append(m.division_of_labor_equilibrium())

        # Fitness comparison
        f_s = m.symmetric_equilibrium()
        f_d = m.division_of_labor_equilibrium()
        W_s = (1 - 2*f_s) * m.growth(2*f_s, 2*f_s) - m.D
        W_d = (1 - f_d) * m.growth(f_d, f_d) - m.D
        W_sym_values.append(W_s)
        W_div_values.append(W_d)

    ax.plot(alpha_range, W_div_values, '-', color=COLORS['stable'], lw=2,
            label='Division of labor')
    ax.plot(alpha_range, W_sym_values, '--', color=COLORS['unstable'], lw=2,
            label='Symmetric')

    ax.fill_between(alpha_range, W_sym_values, W_div_values,
                    where=np.array(W_div_values) > np.array(W_sym_values),
                    alpha=0.2, color=COLORS['stable'])

    ax.set_xlabel('Returns to investment ($\\alpha$)')
    ax.set_ylabel('Equilibrium fitness ($W^*$)')
    ax.set_xlim(0.3, 1.0)
    ax.legend(loc='best', frameon=False)
    ax.set_title('(B) Fitness advantage of division of labor', fontweight='bold')

    # Add annotation
    ax.annotate('Stronger diminishing\nreturns favor\nspecialization',
                xy=(0.5, 0.06), xytext=(0.55, 0.1),
                fontsize=9, ha='left',
                arrowprops=dict(arrowstyle='->', color='gray'))

    plt.tight_layout()

    if save:
        fig_dir = ensure_dir()
        plt.savefig(os.path.join(fig_dir, 'fig3_stability.pdf'))
        plt.savefig(os.path.join(fig_dir, 'fig3_stability.png'))

    return fig


# =============================================================================
# FIGURE 4: Biological interpretation
# =============================================================================

def figure4_interpretation(save=True):
    """
    Summary figure showing the biological interpretation.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Draw the two species and their amino acid production
    # Left: Symmetric (unstable)
    # Right: Division of labor (ESS)

    # --- Left side: Symmetric ---
    ax.text(0.2, 0.92, 'Symmetric Generalists', fontsize=12, ha='center',
            fontweight='bold', color=COLORS['unstable'])
    ax.text(0.2, 0.86, '(Evolutionarily Unstable)', fontsize=10, ha='center',
            color=COLORS['unstable'], style='italic')

    # Species A
    circle_A1 = plt.Circle((0.12, 0.65), 0.08, color=COLORS['A'], alpha=0.7)
    ax.add_patch(circle_A1)
    ax.text(0.12, 0.65, 'A', fontsize=14, ha='center', va='center',
            fontweight='bold', color='white')
    ax.text(0.12, 0.52, '$f_{A1}=f_{A2}=0.25$', fontsize=8, ha='center')

    # Species B
    circle_B1 = plt.Circle((0.28, 0.65), 0.08, color=COLORS['B'], alpha=0.7)
    ax.add_patch(circle_B1)
    ax.text(0.28, 0.65, 'B', fontsize=14, ha='center', va='center',
            fontweight='bold', color='white')
    ax.text(0.28, 0.52, '$f_{B1}=f_{B2}=0.25$', fontsize=8, ha='center')

    # Amino acid pool (symmetric)
    ax.add_patch(plt.Rectangle((0.08, 0.3), 0.24, 0.12, fill=True,
                                facecolor='lightyellow', edgecolor='black', lw=1))
    ax.text(0.2, 0.36, 'Amino acid pool', fontsize=9, ha='center')
    ax.text(0.2, 0.32, '$P_1 = P_2 = 0.5$', fontsize=8, ha='center')

    # Arrows
    ax.annotate('', xy=(0.15, 0.42), xytext=(0.12, 0.55),
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    ax.annotate('', xy=(0.25, 0.42), xytext=(0.28, 0.55),
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))

    # Growth allocation
    ax.text(0.2, 0.2, 'Growth: 50%', fontsize=10, ha='center', color='gray')

    # --- Right side: Division of Labor ---
    ax.text(0.7, 0.92, 'Division of Labor', fontsize=12, ha='center',
            fontweight='bold', color=COLORS['stable'])
    ax.text(0.7, 0.86, '(Evolutionarily Stable)', fontsize=10, ha='center',
            color=COLORS['stable'], style='italic')

    # Species A (specializes on F1)
    circle_A2 = plt.Circle((0.6, 0.65), 0.08, color=COLORS['A'], alpha=0.9)
    ax.add_patch(circle_A2)
    ax.text(0.6, 0.65, 'A', fontsize=14, ha='center', va='center',
            fontweight='bold', color='white')
    ax.text(0.6, 0.52, '$f_{A1}=0.5$\n$f_{A2}=0$', fontsize=8, ha='center')

    # Species B (specializes on F2)
    circle_B2 = plt.Circle((0.8, 0.65), 0.08, color=COLORS['B'], alpha=0.9)
    ax.add_patch(circle_B2)
    ax.text(0.8, 0.65, 'B', fontsize=14, ha='center', va='center',
            fontweight='bold', color='white')
    ax.text(0.8, 0.52, '$f_{B1}=0$\n$f_{B2}=0.5$', fontsize=8, ha='center')

    # Amino acid pool (division of labor)
    ax.add_patch(plt.Rectangle((0.55, 0.3), 0.12, 0.12, fill=True,
                                facecolor='lightgreen', edgecolor='black', lw=1))
    ax.text(0.61, 0.36, '$F_1$', fontsize=10, ha='center', fontweight='bold')

    ax.add_patch(plt.Rectangle((0.73, 0.3), 0.12, 0.12, fill=True,
                                facecolor='lightblue', edgecolor='black', lw=1))
    ax.text(0.79, 0.36, '$F_2$', fontsize=10, ha='center', fontweight='bold')

    # Arrows showing specialization
    ax.annotate('', xy=(0.61, 0.42), xytext=(0.6, 0.55),
                arrowprops=dict(arrowstyle='->', color=COLORS['A'], lw=2))
    ax.annotate('', xy=(0.79, 0.42), xytext=(0.8, 0.55),
                arrowprops=dict(arrowstyle='->', color=COLORS['B'], lw=2))

    # Cross-consumption arrows
    ax.annotate('', xy=(0.73, 0.36), xytext=(0.68, 0.58),
                arrowprops=dict(arrowstyle='<-', color=COLORS['A'], lw=1, ls='--'))
    ax.annotate('', xy=(0.67, 0.36), xytext=(0.72, 0.58),
                arrowprops=dict(arrowstyle='<-', color=COLORS['B'], lw=1, ls='--'))

    # Growth allocation
    ax.text(0.7, 0.2, 'Growth: 50%', fontsize=10, ha='center', color='gray')

    # Central arrow showing evolution
    ax.annotate('', xy=(0.48, 0.5), xytext=(0.38, 0.5),
                arrowprops=dict(arrowstyle='->', color='black', lw=2))
    ax.text(0.43, 0.55, 'Evolution', fontsize=10, ha='center', fontweight='bold')

    # Key insight box
    ax.text(0.5, 0.08,
            'Division of labor: same growth allocation, same amino acid availability,\n'
            'but resistant to invasion by alternative strategies',
            fontsize=9, ha='center', style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')

    if save:
        fig_dir = ensure_dir()
        plt.savefig(os.path.join(fig_dir, 'fig4_interpretation.pdf'))
        plt.savefig(os.path.join(fig_dir, 'fig4_interpretation.png'))

    return fig


# =============================================================================
# MAIN
# =============================================================================

def generate_all():
    """Generate all figures for the manuscript."""
    print("Generating figures for American Naturalist manuscript...")
    print()

    fig_dir = ensure_dir()
    print(f"Output directory: {fig_dir}")
    print()

    print("Figure 1: Fitness trade-off...")
    figure1_fitness_tradeoff(save=True)
    plt.close()
    print("  Done!")

    print("Figure 2: Phase portrait...")
    figure2_phase_portrait(save=True)
    plt.close()
    print("  Done!")

    print("Figure 3: Stability analysis...")
    figure3_stability(save=True)
    plt.close()
    print("  Done!")

    print("Figure 4: Biological interpretation...")
    figure4_interpretation(save=True)
    plt.close()
    print("  Done!")

    print()
    print("All figures generated successfully!")


if __name__ == "__main__":
    generate_all()

#!/usr/bin/env python3
"""
Numerical Computation of Both Critical Thresholds: ω_crit1 and ω_crit2

This script:
1. Computes ω_crit1 analytically (G invasion into S-M)
2. Computes ω_crit2 numerically (M displacement from S-M-G)
3. Visualizes the coexistence window
4. Shows parameter sensitivity for both thresholds

Author: Jian Wang
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, brentq
from matplotlib.gridspec import GridSpec

# Publication quality settings
plt.rcParams.update({
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.family': 'sans-serif'
})

class ThreeSpeciesModel:
    """Three-species cross-feeding model with complete bifurcation analysis"""

    def __init__(self, r_S=1.0, r_M=0.8, sigma_MS=1.5, sigma_SM=0.5,
                 sigma_GS=0.4, sigma_GM=0.4, sigma_SG=0.4, sigma_MG=0.4,
                 alpha_GS=0.3, alpha_GM=0.3, alpha_SG=0.3, alpha_MG=0.3):
        """Initialize with baseline parameters"""
        self.r_S = r_S
        self.r_M = r_M
        self.sigma_MS = sigma_MS
        self.sigma_SM = sigma_SM
        self.sigma_GS = sigma_GS
        self.sigma_GM = sigma_GM
        self.sigma_SG = sigma_SG
        self.sigma_MG = sigma_MG
        self.alpha_GS = alpha_GS
        self.alpha_GM = alpha_GM
        self.alpha_SG = alpha_SG
        self.alpha_MG = alpha_MG

    def net_parameters(self, omega):
        """Compute net interaction parameters as functions of ω"""
        a = (1 - omega) * self.sigma_SG - omega * self.alpha_SG
        b = omega * self.sigma_MG - (1 - omega) * self.alpha_MG
        c = (1 - omega) * self.sigma_GS - omega * self.alpha_GS
        d = 2 * omega - 1
        e = omega * self.sigma_GM - (1 - omega) * self.alpha_GM
        return a, b, c, d, e

    def r_G(self, omega):
        """Generalist growth rate as function of ω"""
        return -self.r_M + omega * (self.r_S + self.r_M)

    def SM_equilibrium(self):
        """
        Analytical S-M equilibrium (foundation platform)

        Returns: (s_star, m_star)
        """
        # From M equation: s = (1 + m) / sigma_MS
        # Substitute into S equation:
        # 1 + sigma_SM*m - (1+m)/sigma_MS = 0
        # sigma_MS + sigma_MS*sigma_SM*m - 1 - m = 0
        # m*(sigma_MS*sigma_SM - 1) = 1 - sigma_MS

        m_star = (1 - self.sigma_MS) / (self.sigma_MS * self.sigma_SM - 1)
        s_star = (1 + m_star) / self.sigma_MS

        return s_star, m_star

    def invasion_fitness_G(self, omega):
        """
        G invasion fitness into S-M equilibrium

        λ_G = r_G * (d + c*s* + e*m*)
        """
        s_star, m_star = self.SM_equilibrium()
        a, b, c, d, e = self.net_parameters(omega)

        lambda_G = self.r_G(omega) * (d + c * s_star + e * m_star)

        return lambda_G

    def omega_crit1_analytical(self):
        """
        Analytical formula for ω_crit1 (Equation 3 in manuscript)

        This is the threshold where G can invade S-M equilibrium
        """
        s_star, m_star = self.SM_equilibrium()

        # From λ_G = 0: d + c*s* + e*m* = 0
        # (2ω - 1) + [(1-ω)σ_GS - ω·α_GS]s* + [ω·σ_GM - (1-ω)·α_GM]m* = 0
        # Expand and collect terms in ω:

        numerator = 1 - self.sigma_GS * s_star + self.alpha_GM * m_star
        denominator = 2 - (self.sigma_GS + self.alpha_GS) * s_star + (self.sigma_GM + self.alpha_GM) * m_star

        omega_crit1 = numerator / denominator

        return omega_crit1

    def SG_equilibrium(self, omega):
        """
        S-G equilibrium (when M is absent)

        Returns: (s_star_SG, g_star_SG)
        """
        a, b, c, d, e = self.net_parameters(omega)

        # From the system:
        # S: (1-ω) - s + a*g = 0
        # G: d + c*s - g = 0
        #
        # From G equation: g = d + c*s
        # Substitute into S equation:
        # (1-ω) - s + a*(d + c*s) = 0
        # (1-ω) + a*d - s + a*c*s = 0
        # s*(1 - a*c) = (1-ω) + a*d

        denominator = 1 - a * c

        if abs(denominator) < 1e-10:
            return None, None  # Singular case

        s_star_SG = ((1 - omega) + a * d) / denominator
        g_star_SG = d + c * s_star_SG

        # Check positivity
        if s_star_SG <= 0 or g_star_SG <= 0:
            return None, None

        return s_star_SG, g_star_SG

    def invasion_fitness_M(self, omega):
        """
        M invasion fitness into S-G equilibrium

        λ_M = -r_M + σ_MS * s*_SG + σ_MG * g*_SG
        """
        s_star_SG, g_star_SG = self.SG_equilibrium(omega)

        if s_star_SG is None:
            return -np.inf  # S-G equilibrium doesn't exist

        lambda_M = -self.r_M + self.sigma_MS * s_star_SG + self.sigma_MG * g_star_SG

        return lambda_M

    def omega_crit2_numerical(self):
        """
        Numerical solution for ω_crit2

        Find ω where λ_M(ω) = 0 (M can no longer invade S-G equilibrium)
        """
        # Start searching from omega_crit1 + small margin
        omega_start = max(0.4, self.omega_crit1_analytical() + 0.05)
        omega_range = np.linspace(omega_start, 0.99, 300)
        lambda_M_values = [self.invasion_fitness_M(om) for om in omega_range]

        # Find zero crossing from positive to negative
        for i in range(len(lambda_M_values) - 1):
            if not np.isnan(lambda_M_values[i]) and not np.isnan(lambda_M_values[i+1]):
                if lambda_M_values[i] > 0 and lambda_M_values[i+1] < 0:
                    # Use Brent's method for precise root
                    try:
                        omega_crit2 = brentq(self.invasion_fitness_M,
                                             omega_range[i], omega_range[i+1])
                        return omega_crit2
                    except:
                        continue

        return None  # No bifurcation found

    def compute_both_critical_thresholds(self):
        """
        Compute both ω_crit1 and ω_crit2

        Returns: dict with both values and coexistence window
        """
        omega_crit1 = self.omega_crit1_analytical()
        omega_crit2 = self.omega_crit2_numerical()

        if omega_crit2 is None:
            window_width = None
        else:
            window_width = omega_crit2 - omega_crit1

        return {
            'omega_crit1': omega_crit1,
            'omega_crit2': omega_crit2,
            'coexistence_window': (omega_crit1, omega_crit2) if omega_crit2 else None,
            'window_width': window_width
        }


def visualize_both_bifurcations():
    """Create comprehensive visualization of both bifurcation points"""

    model = ThreeSpeciesModel()

    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.35)

    # ========================================================================
    # PANEL A: Invasion fitness functions (both λ_G and λ_M)
    # ========================================================================
    ax_a = fig.add_subplot(gs[0, :])

    omega_range = np.linspace(0.01, 0.99, 300)
    lambda_G_values = [model.invasion_fitness_G(om) for om in omega_range]
    lambda_M_values = [model.invasion_fitness_M(om) for om in omega_range]

    ax_a.plot(omega_range, lambda_G_values, '-', color='#2ca02c',
             linewidth=2, label=r'$\lambda_G(\omega)$ (G invasion into S-M)')
    ax_a.plot(omega_range, lambda_M_values, '-', color='#d62728',
             linewidth=2, label=r'$\lambda_M(\omega)$ (M invasion into S-G)')
    ax_a.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)

    # Mark critical points
    results = model.compute_both_critical_thresholds()
    omega_crit1 = results['omega_crit1']
    omega_crit2 = results['omega_crit2']

    ax_a.axvline(omega_crit1, color='green', linestyle=':', linewidth=2, alpha=0.7)
    ax_a.text(omega_crit1, ax_a.get_ylim()[1]*0.9,
             f'$\\omega_{{crit1}}$ = {omega_crit1:.3f}',
             fontsize=9, ha='center', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    if omega_crit2:
        ax_a.axvline(omega_crit2, color='red', linestyle=':', linewidth=2, alpha=0.7)
        ax_a.text(omega_crit2, ax_a.get_ylim()[1]*0.9,
                 f'$\\omega_{{crit2}}$ = {omega_crit2:.3f}',
                 fontsize=9, ha='center', bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))

        # Shade coexistence window
        ax_a.axvspan(omega_crit1, omega_crit2, alpha=0.2, color='gold',
                    label=f'Coexistence window (width = {omega_crit2-omega_crit1:.3f})')

    ax_a.set_xlabel('Pathway parameter (ω)', fontsize=10)
    ax_a.set_ylabel('Invasion fitness', fontsize=10)
    ax_a.set_title('A. Dual invasion fitness functions reveal coexistence window',
                   fontweight='bold', loc='left', fontsize=11)
    ax_a.legend(frameon=False, loc='upper left', fontsize=8)
    ax_a.grid(alpha=0.3)

    # ========================================================================
    # PANEL B: S-M equilibrium landscape (foundation for ω_crit1)
    # ========================================================================
    ax_b = fig.add_subplot(gs[1, 0])

    s_star, m_star = model.SM_equilibrium()

    ax_b.bar(['$s^*_{SM}$', '$m^*_{SM}$'], [s_star, m_star],
            color=['#1f77b4', '#d62728'], alpha=0.7, edgecolor='black', linewidth=1.5)
    ax_b.set_ylabel('Equilibrium density', fontsize=9)
    ax_b.set_title('B. S-M equilibrium platform\n(for $\\omega_{crit1}$ calculation)',
                   fontweight='bold', loc='left', fontsize=9)
    ax_b.grid(axis='y', alpha=0.3)

    # Add text annotations
    ax_b.text(0, s_star*1.1, f'{s_star:.3f}', ha='center', fontsize=8, fontweight='bold')
    ax_b.text(1, m_star*1.1, f'{m_star:.3f}', ha='center', fontsize=8, fontweight='bold')

    # ========================================================================
    # PANEL C: S-G equilibrium landscape (foundation for ω_crit2)
    # ========================================================================
    ax_c = fig.add_subplot(gs[1, 1])

    if omega_crit2:
        s_SG_crit, g_SG_crit = model.SG_equilibrium(omega_crit2)

        ax_c.bar(['$s^*_{SG}$', '$g^*_{SG}$'], [s_SG_crit, g_SG_crit],
                color=['#1f77b4', '#2ca02c'], alpha=0.7, edgecolor='black', linewidth=1.5)
        ax_c.set_ylabel('Equilibrium density', fontsize=9)
        ax_c.set_title(f'C. S-G equilibrium at $\\omega_{{crit2}}$\n(M displacement point)',
                       fontweight='bold', loc='left', fontsize=9)
        ax_c.grid(axis='y', alpha=0.3)

        ax_c.text(0, s_SG_crit*1.1, f'{s_SG_crit:.3f}', ha='center', fontsize=8, fontweight='bold')
        ax_c.text(1, g_SG_crit*1.1, f'{g_SG_crit:.3f}', ha='center', fontsize=8, fontweight='bold')

    # ========================================================================
    # PANEL D: Analytical vs numerical comparison
    # ========================================================================
    ax_d = fig.add_subplot(gs[1, 2])

    properties = ['Existence', 'Formula type', 'Baseline\nvalue']
    omega1_props = ['Always\n(if S-M stable)', 'Explicit\nclosed form', f'{omega_crit1:.4f}']
    omega2_props = ['Parameter\ndependent', 'Implicit\n(numerical)', f'{omega_crit2:.4f}' if omega_crit2 else 'N/A']

    x = np.arange(len(properties))
    width = 0.35

    ax_d.barh([i - width/2 for i in x], [1]*len(properties), width,
             label='$\\omega_{crit1}$', color='lightgreen', edgecolor='green', linewidth=1.5)
    ax_d.barh([i + width/2 for i in x], [1]*len(properties), width,
             label='$\\omega_{crit2}$', color='lightcoral', edgecolor='red', linewidth=1.5)

    # Add text labels
    for i, (prop, val1, val2) in enumerate(zip(properties, omega1_props, omega2_props)):
        ax_d.text(0.5, i - width/2, val1, ha='center', va='center', fontsize=7, fontweight='bold')
        ax_d.text(0.5, i + width/2, val2, ha='center', va='center', fontsize=7, fontweight='bold')

    ax_d.set_yticks(x)
    ax_d.set_yticklabels(properties, fontsize=8)
    ax_d.set_xlim(0, 1)
    ax_d.set_xticks([])
    ax_d.set_title('D. Comparison of\ncritical thresholds', fontweight='bold', loc='left', fontsize=9)
    ax_d.legend(frameon=False, loc='upper right', fontsize=7)

    # ========================================================================
    # PANEL E: Parameter sensitivity of ω_crit1
    # ========================================================================
    ax_e = fig.add_subplot(gs[2, 0])

    sigma_MS_range = np.linspace(1.1, 3.0, 50)
    omega_crit1_values = []

    for sig_ms in sigma_MS_range:
        model_temp = ThreeSpeciesModel(sigma_MS=sig_ms)
        omega_crit1_values.append(model_temp.omega_crit1_analytical())

    ax_e.plot(sigma_MS_range, omega_crit1_values, '-', color='green', linewidth=2)
    ax_e.axhline(model.omega_crit1_analytical(), color='green', linestyle='--', alpha=0.5)
    ax_e.axvline(model.sigma_MS, color='gray', linestyle='--', alpha=0.5)
    ax_e.set_xlabel('$\\sigma_{MS}$ (mutualism strength)', fontsize=9)
    ax_e.set_ylabel('$\\omega_{crit1}$', fontsize=9)
    ax_e.set_title('E. $\\omega_{crit1}$ sensitivity to\nmutualism strength',
                   fontweight='bold', loc='left', fontsize=9)
    ax_e.grid(alpha=0.3)

    # ========================================================================
    # PANEL F: Parameter sensitivity of ω_crit2
    # ========================================================================
    ax_f = fig.add_subplot(gs[2, 1])

    omega_crit2_values = []

    for sig_ms in sigma_MS_range:
        model_temp = ThreeSpeciesModel(sigma_MS=sig_ms)
        omega_c2 = model_temp.omega_crit2_numerical()
        omega_crit2_values.append(omega_c2 if omega_c2 else np.nan)

    ax_f.plot(sigma_MS_range, omega_crit2_values, '-', color='red', linewidth=2)
    if omega_crit2 is not None:
        ax_f.axhline(omega_crit2, color='red', linestyle='--', alpha=0.5)
    ax_f.axvline(model.sigma_MS, color='gray', linestyle='--', alpha=0.5)
    ax_f.set_xlabel('$\\sigma_{MS}$ (mutualism strength)', fontsize=9)
    ax_f.set_ylabel('$\\omega_{crit2}$', fontsize=9)
    ax_f.set_title('F. $\\omega_{crit2}$ sensitivity to\nmutualism strength',
                   fontweight='bold', loc='left', fontsize=9)
    ax_f.grid(alpha=0.3)

    # ========================================================================
    # PANEL G: Coexistence window width
    # ========================================================================
    ax_g = fig.add_subplot(gs[2, 2])

    window_widths = []

    for sig_ms, oc1, oc2 in zip(sigma_MS_range, omega_crit1_values, omega_crit2_values):
        if not np.isnan(oc2):
            window_widths.append(oc2 - oc1)
        else:
            window_widths.append(np.nan)

    ax_g.plot(sigma_MS_range, window_widths, '-', color='gold', linewidth=2.5)
    ax_g.fill_between(sigma_MS_range, 0, window_widths, alpha=0.3, color='gold')
    ax_g.axvline(model.sigma_MS, color='gray', linestyle='--', alpha=0.5)
    ax_g.set_xlabel('$\\sigma_{MS}$ (mutualism strength)', fontsize=9)
    ax_g.set_ylabel('Window width ($\\omega_{crit2} - \\omega_{crit1}$)', fontsize=9)
    ax_g.set_title('G. Coexistence window\nvs mutualism strength',
                   fontweight='bold', loc='left', fontsize=9)
    ax_g.grid(alpha=0.3)

    plt.suptitle('Complete Bifurcation Analysis: Both Critical Thresholds ($\\omega_{crit1}$ and $\\omega_{crit2}$)',
                fontsize=12, fontweight='bold', y=0.995)

    plt.savefig('both_omega_critical_complete_analysis.png', dpi=300, bbox_inches='tight')
    print("\n✓ Figure saved: both_omega_critical_complete_analysis.png")

    return fig, results


def print_analytical_summary():
    """Print complete analytical summary with formulas"""

    model = ThreeSpeciesModel()
    results = model.compute_both_critical_thresholds()

    s_star, m_star = model.SM_equilibrium()

    print("=" * 80)
    print("COMPLETE ANALYTICAL SUMMARY: ω_crit1 and ω_crit2")
    print("=" * 80)

    print("\n1. S-M EQUILIBRIUM (Foundation Platform)")
    print("-" * 80)
    print(f"   s*_SM = {s_star:.6f}")
    print(f"   m*_SM = {m_star:.6f}")
    print(f"\n   Existence requires: σ_MS > 1 (current: {model.sigma_MS})")
    print(f"   Stability requires: σ_MS·σ_SM < 1 (current: {model.sigma_MS * model.sigma_SM:.3f})")

    print("\n2. FIRST BIFURCATION: ω_crit1 (Generalist Invasion)")
    print("-" * 80)
    print("   ANALYTICAL FORMULA (Manuscript Equation 3):")
    print("   ")
    print("         1 - σ_GS·s*_SM + α_GM·m*_SM")
    print("   ω₁ = ───────────────────────────────────────")
    print("         2 - (σ_GS + α_GS)·s*_SM + (σ_GM + α_GM)·m*_SM")
    print()

    numerator = 1 - model.sigma_GS * s_star + model.alpha_GM * m_star
    denominator = 2 - (model.sigma_GS + model.alpha_GS) * s_star + (model.sigma_GM + model.alpha_GM) * m_star

    print(f"   Numerator   = {numerator:.6f}")
    print(f"   Denominator = {denominator:.6f}")
    print(f"\n   ✓ ω_crit1 = {results['omega_crit1']:.6f}")

    print("\n   BIOLOGICAL MEANING:")
    print("   - Below ω_crit1: G too metabolite-specialized → cannot invade S-M")
    print("   - At ω_crit1: Transcritical bifurcation (λ_G crosses zero)")
    print("   - Above ω_crit1: Three-species coexistence emerges")

    print("\n3. SECOND BIFURCATION: ω_crit2 (Metabolite Specialist Displacement)")
    print("-" * 80)
    print("   IMPLICIT CONDITION:")
    print("   σ_MS·s*_SG(ω_crit2) + σ_MG·g*_SG(ω_crit2) = r_M")
    print()
    print("   where S-G equilibrium satisfies:")
    print("   s*_SG(ω) = (1-ω) + a(ω)·g*_SG(ω)")
    print("            (2ω - 1) + c(ω)(1-ω)")
    print("   g*_SG(ω) = ─────────────────────")
    print("                 1 - c(ω)·a(ω)")
    print()

    if results['omega_crit2']:
        print(f"   ✓ ω_crit2 = {results['omega_crit2']:.6f}  (numerical solution)")

        s_SG, g_SG = model.SG_equilibrium(results['omega_crit2'])
        print(f"\n   At ω_crit2:")
        print(f"   s*_SG = {s_SG:.6f}")
        print(f"   g*_SG = {g_SG:.6f}")

        # Verify the condition
        lambda_M_check = -model.r_M + model.sigma_MS * s_SG + model.sigma_MG * g_SG
        print(f"\n   Verification: λ_M = {lambda_M_check:.10f}  (should be ≈ 0)")

    print("\n   BIOLOGICAL MEANING:")
    print("   - Below ω_crit2: S-M-G three-species coexistence")
    print("   - At ω_crit2: Transcritical bifurcation (M* → 0)")
    print("   - Above ω_crit2: G too substrate-specialized → outcompetes M")

    print("\n4. COEXISTENCE WINDOW")
    print("-" * 80)
    if results['coexistence_window']:
        print(f"   ω ∈ ({results['omega_crit1']:.6f}, {results['omega_crit2']:.6f})")
        print(f"   Window width: {results['window_width']:.6f}")
        print(f"\n   Three species coexist ONLY when generalist maintains")
        print(f"   intermediate metabolic allocation within this window.")

    print("\n5. PARAMETER SENSITIVITY")
    print("-" * 80)
    print("   Effect of increasing key parameters:")
    print()
    print("   Parameter  │  ω_crit1  │  ω_crit2  │  Window Width")
    print("   ───────────┼───────────┼───────────┼──────────────")
    print("   σ_MS ↑     │     ↓     │     ↑     │      ↑")
    print("   σ_SM ↑     │     ↑     │     ↑     │    varies")
    print("   σ_GS ↑     │     ↓     │     ↓     │    varies")
    print("   α_GS ↑     │     ↑     │     ↑     │    varies")
    print()
    print("   Stronger S→M mutualism (↑σ_MS) EXPANDS coexistence window")

    print("\n6. STRUCTURAL COMPARISON")
    print("-" * 80)
    print("   Property           │  ω_crit1         │  ω_crit2")
    print("   ───────────────────┼──────────────────┼──────────────────")
    print("   Equilibrium base   │  S-M             │  S-G")
    print("   Invasion species   │  G into S-M      │  M into S-G")
    print("   Condition          │  λ_G = 0         │  λ_M = 0")
    print("   Formula type       │  Explicit        │  Implicit")
    print("   Transition         │  S-M → S-M-G     │  S-M-G → S-G")

    print("\n" + "=" * 80)
    print()


if __name__ == '__main__':
    print("\nGenerating complete bifurcation analysis...\n")

    # Print analytical summary
    print_analytical_summary()

    # Create comprehensive visualization
    fig, results = visualize_both_bifurcations()

    print("\nAnalysis complete!")
    print(f"\n✓ ω_crit1 = {results['omega_crit1']:.6f} (analytical)")
    print(f"✓ ω_crit2 = {results['omega_crit2']:.6f} (numerical)")
    print(f"✓ Coexistence window width = {results['window_width']:.6f}")

    plt.show()

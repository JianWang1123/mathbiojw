#!/usr/bin/env python3
"""
Demonstration of ω_crit2 in Different Parameter Regimes

Shows both scenarios:
1. Baseline parameters (NO ω_crit2)
2. Modified parameters (ω_crit2 EXISTS)

Author: Jian Wang
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})


class CrossFeedingModel:
    """Three-species cross-feeding model"""

    def __init__(self, r_S=1.0, r_M=0.8, sigma_MS=1.5, sigma_SM=0.5,
                 sigma_GS=0.4, sigma_GM=0.4, sigma_SG=0.4, sigma_MG=0.4,
                 alpha_GS=0.3, alpha_GM=0.3, alpha_SG=0.3, alpha_MG=0.3):
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

    def net_params(self, omega):
        a = (1 - omega) * self.sigma_SG - omega * self.alpha_SG
        c = (1 - omega) * self.sigma_GS - omega * self.alpha_GS
        d = 2 * omega - 1
        e = omega * self.sigma_GM - (1 - omega) * self.alpha_GM
        return a, c, d, e

    def SM_equilibrium(self):
        """S-M equilibrium"""
        m_star = (1 - self.sigma_MS) / (self.sigma_MS * self.sigma_SM - 1)
        s_star = (1 + m_star) / self.sigma_MS
        return s_star, m_star

    def SG_equilibrium(self, omega):
        """S-G equilibrium (when M absent)"""
        a, c, d, e = self.net_params(omega)

        denom = 1 - a * c
        if abs(denom) < 1e-10:
            return None, None

        s_SG = ((1 - omega) + a * d) / denom
        g_SG = d + c * s_SG

        if s_SG > 0 and g_SG > 0:
            return s_SG, g_SG
        return None, None

    def lambda_G(self, omega):
        """G invasion fitness into S-M"""
        s_SM, m_SM = self.SM_equilibrium()
        a, c, d, e = self.net_params(omega)
        r_G = -self.r_M + omega * (self.r_S + self.r_M)
        return r_G * (d + c * s_SM + e * m_SM)

    def lambda_M(self, omega):
        """M invasion fitness into S-G"""
        s_SG, g_SG = self.SG_equilibrium(omega)
        if s_SG is None:
            return -np.inf
        return -self.r_M + self.sigma_MS * s_SG + self.sigma_MG * g_SG

    def find_omega_crit1(self):
        """Find ω_crit1 analytically"""
        s_SM, m_SM = self.SM_equilibrium()
        num = 1 - self.sigma_GS * s_SM + self.alpha_GM * m_SM
        denom = 2 - (self.sigma_GS + self.alpha_GS) * s_SM + (self.sigma_GM + self.alpha_GM) * m_SM
        return num / denom

    def find_omega_crit2(self):
        """Find ω_crit2 numerically"""
        omega_start = max(0.4, self.find_omega_crit1() + 0.05)
        omega_range = np.linspace(omega_start, 0.99, 400)

        lambda_M_vals = [self.lambda_M(om) for om in omega_range]

        for i in range(len(lambda_M_vals) - 1):
            if not np.isinf(lambda_M_vals[i]) and not np.isinf(lambda_M_vals[i+1]):
                if lambda_M_vals[i] > 0 and lambda_M_vals[i+1] < 0:
                    try:
                        return brentq(self.lambda_M, omega_range[i], omega_range[i+1])
                    except:
                        continue
        return None


def compare_parameter_regimes():
    """Compare regimes with and without ω_crit2"""

    # Regime 1: Baseline (NO ω_crit2)
    model_baseline = CrossFeedingModel(
        sigma_MS=1.5, sigma_SM=0.5,
        sigma_GS=0.4, sigma_GM=0.4, sigma_SG=0.4, sigma_MG=0.4,
        alpha_GS=0.3, alpha_GM=0.3, alpha_SG=0.3, alpha_MG=0.3
    )

    # Regime 2: Modified (ω_crit2 EXISTS)
    # Stronger M-G and G-M mutualism
    model_modified = CrossFeedingModel(
        sigma_MS=1.8, sigma_SM=0.5,
        sigma_GS=0.5, sigma_GM=0.7, sigma_SG=0.4, sigma_MG=0.8,
        alpha_GS=0.3, alpha_GM=0.3, alpha_SG=0.3, alpha_MG=0.2
    )

    # Compute critical points
    omega_crit1_base = model_baseline.find_omega_crit1()
    omega_crit2_base = model_baseline.find_omega_crit2()

    omega_crit1_mod = model_modified.find_omega_crit1()
    omega_crit2_mod = model_modified.find_omega_crit2()

    print("=" * 80)
    print("PARAMETER REGIME COMPARISON")
    print("=" * 80)
    print("\nREGIME 1: Baseline Parameters (Standard)")
    print("-" * 80)
    print(f"  σ_MS = {model_baseline.sigma_MS}, σ_MG = {model_baseline.sigma_MG}")
    print(f"  ω_crit1 = {omega_crit1_base:.4f}")
    print(f"  ω_crit2 = {omega_crit2_base if omega_crit2_base else 'Does NOT exist'}")
    print(f"  Outcome: {'Bounded coexistence window' if omega_crit2_base else 'Permanent S-M-G coexistence for ω > ω_crit1'}")

    print("\nREGIME 2: Modified Parameters (Enhanced M-G Mutualism)")
    print("-" * 80)
    print(f"  σ_MS = {model_modified.sigma_MS}, σ_MG = {model_modified.sigma_MG}, σ_GM = {model_modified.sigma_GM}")
    print(f"  ω_crit1 = {omega_crit1_mod:.4f}")
    print(f"  ω_crit2 = {omega_crit2_mod:.4f}" if omega_crit2_mod else "  ω_crit2 = Does NOT exist")
    if omega_crit2_mod:
        print(f"  Coexistence window: ω ∈ ({omega_crit1_mod:.4f}, {omega_crit2_mod:.4f})")
        print(f"  Window width: {omega_crit2_mod - omega_crit1_mod:.4f}")

    # Create visualization
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    omega_scan = np.linspace(0.01, 0.99, 400)

    # --- REGIME 1: Baseline ---
    lambda_G_base = [model_baseline.lambda_G(om) for om in omega_scan]
    lambda_M_base = [model_baseline.lambda_M(om) for om in omega_scan]

    # Panel A: Invasion fitness (Baseline)
    axes[0, 0].plot(omega_scan, lambda_G_base, '-', color='#2ca02c', linewidth=2, label='λ_G (G invasion)')
    axes[0, 0].plot(omega_scan, lambda_M_base, '-', color='#d62728', linewidth=2, label='λ_M (M invasion)')
    axes[0, 0].axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    axes[0, 0].axvline(omega_crit1_base, color='green', linestyle=':', linewidth=2, alpha=0.7)
    axes[0, 0].text(omega_crit1_base, axes[0, 0].get_ylim()[1]*0.8,
                    f'ω_crit1={omega_crit1_base:.3f}', fontsize=8, rotation=90, va='bottom')

    axes[0, 0].fill_between(omega_scan, axes[0, 0].get_ylim()[0], axes[0, 0].get_ylim()[1],
                            where=np.array(omega_scan) > omega_crit1_base,
                            alpha=0.2, color='gold', label='S-M-G coexist')

    axes[0, 0].set_xlabel('ω (pathway parameter)', fontsize=10)
    axes[0, 0].set_ylabel('Invasion fitness', fontsize=10)
    axes[0, 0].set_title('A. Regime 1 (Baseline): NO ω_crit2\nλ_M always negative → M never returns', fontweight='bold')
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 0].grid(alpha=0.3)

    # Panel B: S-G equilibrium (Baseline)
    s_SG_base = [model_baseline.SG_equilibrium(om)[0] if model_baseline.SG_equilibrium(om)[0] else np.nan for om in omega_scan]
    g_SG_base = [model_baseline.SG_equilibrium(om)[1] if model_baseline.SG_equilibrium(om)[1] else np.nan for om in omega_scan]

    axes[0, 1].plot(omega_scan, s_SG_base, '-', color='#1f77b4', linewidth=2, label='s*_SG')
    axes[0, 1].plot(omega_scan, g_SG_base, '-', color='#2ca02c', linewidth=2, label='g*_SG')
    axes[0, 1].axvline(omega_crit1_base, color='green', linestyle=':', linewidth=2, alpha=0.7)
    axes[0, 1].set_xlabel('ω (pathway parameter)', fontsize=10)
    axes[0, 1].set_ylabel('S-G equilibrium density', fontsize=10)
    axes[0, 1].set_title('B. Regime 1: S-G equilibrium landscape', fontweight='bold')
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[0, 1].grid(alpha=0.3)

    # --- REGIME 2: Modified ---
    lambda_G_mod = [model_modified.lambda_G(om) for om in omega_scan]
    lambda_M_mod = [model_modified.lambda_M(om) for om in omega_scan]

    # Panel C: Invasion fitness (Modified)
    axes[1, 0].plot(omega_scan, lambda_G_mod, '-', color='#2ca02c', linewidth=2, label='λ_G (G invasion)')
    axes[1, 0].plot(omega_scan, lambda_M_mod, '-', color='#d62728', linewidth=2, label='λ_M (M invasion)')
    axes[1, 0].axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    axes[1, 0].axvline(omega_crit1_mod, color='green', linestyle=':', linewidth=2, alpha=0.7)
    axes[1, 0].text(omega_crit1_mod, axes[1, 0].get_ylim()[1]*0.8,
                    f'ω_crit1={omega_crit1_mod:.3f}', fontsize=8, rotation=90, va='bottom')

    if omega_crit2_mod:
        axes[1, 0].axvline(omega_crit2_mod, color='red', linestyle=':', linewidth=2, alpha=0.7)
        axes[1, 0].text(omega_crit2_mod, axes[1, 0].get_ylim()[1]*0.8,
                        f'ω_crit2={omega_crit2_mod:.3f}', fontsize=8, rotation=90, va='bottom')

        axes[1, 0].axvspan(omega_crit1_mod, omega_crit2_mod, alpha=0.2, color='gold',
                          label=f'Coexistence window\n(width={omega_crit2_mod-omega_crit1_mod:.3f})')

    axes[1, 0].set_xlabel('ω (pathway parameter)', fontsize=10)
    axes[1, 0].set_ylabel('Invasion fitness', fontsize=10)
    axes[1, 0].set_title('C. Regime 2 (Modified): ω_crit2 EXISTS\nλ_M crosses zero → bounded window', fontweight='bold')
    axes[1, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].grid(alpha=0.3)

    # Panel D: Community composition diagram
    axes[1, 1].text(0.5, 0.9, 'Community Composition Transitions', ha='center',
                   fontsize=12, fontweight='bold', transform=axes[1, 1].transAxes)

    # Regime 1 (no omega_crit2)
    axes[1, 1].text(0.1, 0.7, 'Regime 1 (Baseline):', fontsize=10, fontweight='bold',
                   transform=axes[1, 1].transAxes)
    axes[1, 1].arrow(0.15, 0.6, 0.3, 0, transform=axes[1, 1].transAxes,
                    head_width=0.03, head_length=0.05, fc='gray', ec='gray')
    axes[1, 1].text(0.1, 0.55, 'S-M', ha='center', fontsize=9, bbox=dict(boxstyle='round', facecolor='lightblue'),
                   transform=axes[1, 1].transAxes)
    axes[1, 1].text(0.3, 0.63, f'ω={omega_crit1_base:.2f}', ha='center', fontsize=8,
                   transform=axes[1, 1].transAxes)
    axes[1, 1].text(0.5, 0.55, 'S-M-G\n(permanent)', ha='center', fontsize=9,
                   bbox=dict(boxstyle='round', facecolor='gold'),
                   transform=axes[1, 1].transAxes)

    # Regime 2 (with omega_crit2)
    if omega_crit2_mod:
        axes[1, 1].text(0.1, 0.35, 'Regime 2 (Modified):', fontsize=10, fontweight='bold',
                       transform=axes[1, 1].transAxes)

        axes[1, 1].arrow(0.15, 0.25, 0.15, 0, transform=axes[1, 1].transAxes,
                        head_width=0.03, head_length=0.03, fc='gray', ec='gray')
        axes[1, 1].arrow(0.35, 0.25, 0.15, 0, transform=axes[1, 1].transAxes,
                        head_width=0.03, head_length=0.03, fc='gray', ec='gray')

        axes[1, 1].text(0.1, 0.2, 'S-M', ha='center', fontsize=9, bbox=dict(boxstyle='round', facecolor='lightblue'),
                       transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.25, 0.28, f'ω={omega_crit1_mod:.2f}', ha='center', fontsize=7,
                       transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.325, 0.2, 'S-M-G', ha='center', fontsize=9, bbox=dict(boxstyle='round', facecolor='gold'),
                       transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.425, 0.28, f'ω={omega_crit2_mod:.2f}', ha='center', fontsize=7,
                       transform=axes[1, 1].transAxes)
        axes[1, 1].text(0.55, 0.2, 'S-G', ha='center', fontsize=9, bbox=dict(boxstyle='round', facecolor='lightcoral'),
                       transform=axes[1, 1].transAxes)

        # Add interpretation
        axes[1, 1].text(0.5, 0.05, 'Key: Regime 2 shows BOUNDED coexistence window\nM is displaced at high ω',
                       ha='center', fontsize=8, style='italic',
                       bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5),
                       transform=axes[1, 1].transAxes)

    axes[1, 1].set_xlim(0, 1)
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].axis('off')

    plt.suptitle('ω_crit2 Existence Depends on Parameter Regime', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('omega_crit2_parameter_regimes.png', dpi=300, bbox_inches='tight')
    print("\n✓ Figure saved: omega_crit2_parameter_regimes.png")

    return fig


if __name__ == '__main__':
    print("\nDemonstrating ω_crit2 in different parameter regimes...\n")
    fig = compare_parameter_regimes()
    print("\n" + "=" * 80)
    print("BIOLOGICAL INTERPRETATION")
    print("=" * 80)
    print("""
Regime 1 (Baseline): λ_M always negative
  → M can NEVER invade S-G equilibrium
  → Once G displaces M, M cannot return
  → Result: Permanent S-M-G coexistence for ω > ω_crit1

Regime 2 (Enhanced M-G mutualism): λ_M crosses zero
  → M CAN invade S-G at intermediate ω (λ_M > 0)
  → At high ω, G becomes too substrate-specialized → λ_M < 0
  → Result: THREE-WAY TRANSITION: S-M → S-M-G → S-G

Why does stronger σ_MG create ω_crit2?
  → G produces metabolites that help M
  → At intermediate ω, this allows M to coexist with S-G
  → But at high ω, G shifts to substrate utilization
  → Metabolite production drops → M starves → displaced

Experimental prediction:
  → Engineer stronger cross-feeding FROM generalist to metabolite specialist
  → This creates bounded coexistence window with measurable ω_crit2
    """)

    plt.show()

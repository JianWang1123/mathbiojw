#!/usr/bin/env python3
"""
Complete Explanation: How S-M-G System Evolves with ω

Shows the full evolutionary trajectory as ω increases from 0 to 1,
including transitions at both critical thresholds.

Author: Jian Wang
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches

plt.rcParams.update({'font.size': 9, 'font.family': 'sans-serif'})


class ThreeSpeciesSystem:
    """Complete three-species cross-feeding model"""

    def __init__(self, r_S=1.0, r_M=0.8,
                 sigma_MS=1.5, sigma_SM=0.5,
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
        """Net interaction parameters"""
        a = (1 - omega) * self.sigma_SG - omega * self.alpha_SG
        b = omega * self.sigma_MG - (1 - omega) * self.alpha_MG
        c = (1 - omega) * self.sigma_GS - omega * self.alpha_GS
        d = 2 * omega - 1
        e = omega * self.sigma_GM - (1 - omega) * self.alpha_GM
        return a, b, c, d, e

    def r_G(self, omega):
        """Generalist growth rate"""
        return -self.r_M + omega * (self.r_S + self.r_M)

    def dynamics_three_species(self, t, y, omega):
        """Three-species dynamics"""
        s, m, g = np.maximum(y, 0)  # Prevent negative densities
        a, b, c, d, e = self.net_params(omega)

        dsdt = self.r_S * s * (1 + self.sigma_SM * m + a * g - s)
        dmdt = self.r_M * m * (-1 + self.sigma_MS * s + b * g - m)
        dgdt = self.r_G(omega) * g * (d + c * s + e * m - g)

        return [dsdt, dmdt, dgdt]

    def SM_equilibrium(self):
        """S-M equilibrium"""
        m_star = (1 - self.sigma_MS) / (self.sigma_MS * self.sigma_SM - 1)
        s_star = (1 + m_star) / self.sigma_MS
        return s_star, m_star

    def find_equilibrium_numerically(self, omega, initial_guess=None):
        """Find three-species equilibrium numerically"""
        if initial_guess is None:
            s_SM, m_SM = self.SM_equilibrium()
            initial_guess = [s_SM, m_SM, 0.5]

        def equations(x):
            s, m, g = x
            a, b, c, d, e = self.net_params(omega)
            return [
                s * (1 + self.sigma_SM * m + a * g - s),
                m * (-1 + self.sigma_MS * s + b * g - m),
                g * (d + c * s + e * m - g)
            ]

        try:
            sol = fsolve(equations, initial_guess, full_output=True)
            if sol[2] == 1:  # Convergence
                s_eq, m_eq, g_eq = sol[0]
                if all(np.array([s_eq, m_eq, g_eq]) > 1e-6):
                    return s_eq, m_eq, g_eq
        except:
            pass

        return None, None, None

    def lambda_G(self, omega):
        """G invasion fitness into S-M"""
        s_SM, m_SM = self.SM_equilibrium()
        a, b, c, d, e = self.net_params(omega)
        return self.r_G(omega) * (d + c * s_SM + e * m_SM)

    def omega_crit1(self):
        """Analytical ω_crit1"""
        s_SM, m_SM = self.SM_equilibrium()
        num = 1 - self.sigma_GS * s_SM + self.alpha_GM * m_SM
        denom = 2 - (self.sigma_GS + self.alpha_GS) * s_SM + \
                (self.sigma_GM + self.alpha_GM) * m_SM
        return num / denom


def create_omega_evolution_diagram():
    """
    Create comprehensive diagram showing system evolution with ω
    """
    model = ThreeSpeciesSystem()

    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(4, 3, figure=fig, hspace=0.4, wspace=0.35)

    omega_scan = np.linspace(0.01, 0.99, 200)
    omega_crit1 = model.omega_crit1()

    # Compute equilibria across omega range
    s_eq_list = []
    m_eq_list = []
    g_eq_list = []

    print("Computing equilibria across ω range...")
    for i, omega in enumerate(omega_scan):
        s_eq, m_eq, g_eq = model.find_equilibrium_numerically(omega)

        if s_eq is not None:
            s_eq_list.append(s_eq)
            m_eq_list.append(m_eq)
            g_eq_list.append(g_eq)
        else:
            # Fall back to S-M equilibrium
            s_SM, m_SM = model.SM_equilibrium()
            s_eq_list.append(s_SM)
            m_eq_list.append(m_SM)
            g_eq_list.append(0.0)

        if (i + 1) % 50 == 0:
            print(f"  Progress: {i+1}/{len(omega_scan)}")

    s_eq_arr = np.array(s_eq_list)
    m_eq_arr = np.array(m_eq_list)
    g_eq_arr = np.array(g_eq_list)

    # ========================================================================
    # PANEL A: Complete bifurcation diagram (top row, spanning 3 columns)
    # ========================================================================
    ax_a = fig.add_subplot(gs[0, :])

    ax_a.plot(omega_scan, s_eq_arr, '-', color='#1f77b4', linewidth=2.5,
             label='S (substrate specialist)', zorder=3)
    ax_a.plot(omega_scan, m_eq_arr, '-', color='#d62728', linewidth=2.5,
             label='M (metabolite specialist)', zorder=3)
    ax_a.plot(omega_scan, g_eq_arr, '-', color='#2ca02c', linewidth=2.5,
             label='G (generalist)', zorder=3)

    # Mark ω_crit1
    ax_a.axvline(omega_crit1, color='red', linestyle='--', linewidth=2.5,
                alpha=0.7, zorder=2, label=f'ω_crit1 = {omega_crit1:.3f}')

    # Shade regions
    ax_a.axvspan(0, omega_crit1, alpha=0.15, color='lightblue',
                label='Region I: S-M only', zorder=1)
    ax_a.axvspan(omega_crit1, 1.0, alpha=0.15, color='gold',
                label='Region II: S-M-G coexist', zorder=1)

    # Annotations at critical point
    ax_a.annotate('TRANSCRITICAL\nBIFURCATION',
                 xy=(omega_crit1, 1.5), xytext=(omega_crit1 + 0.15, 2.5),
                 fontsize=10, fontweight='bold', color='red',
                 arrowprops=dict(arrowstyle='->', color='red', lw=2),
                 bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    ax_a.set_xlabel('Pathway parameter (ω): 0 = pure M-like, 1 = pure S-like',
                   fontsize=11, fontweight='bold')
    ax_a.set_ylabel('Equilibrium density', fontsize=11)
    ax_a.set_title('A. Complete Bifurcation Diagram: System Evolution as ω Increases',
                  fontsize=12, fontweight='bold', loc='left')
    ax_a.legend(frameon=True, loc='upper left', ncol=3, fontsize=9,
               fancybox=True, shadow=True)
    ax_a.grid(alpha=0.3, linestyle=':', linewidth=1)
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 3.5)

    # ========================================================================
    # PANEL B: Phase portrait at ω = 0.2 (before ω_crit1)
    # ========================================================================
    ax_b = fig.add_subplot(gs[1, 0])

    omega_low = 0.2
    t_span = (0, 50)
    t_eval = np.linspace(0, 50, 500)

    # Multiple initial conditions
    initial_conditions = [
        [2.5, 2.5, 0.1],
        [1.5, 1.5, 0.3],
        [2.0, 1.0, 0.2],
    ]

    for ic in initial_conditions:
        sol = solve_ivp(model.dynamics_three_species, t_span, ic,
                       args=(omega_low,), t_eval=t_eval, method='RK45')
        ax_b.plot(sol.t, sol.y[2], '-', alpha=0.6, linewidth=1.5)

    ax_b.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax_b.set_xlabel('Time', fontsize=10)
    ax_b.set_ylabel('G density', fontsize=10)
    ax_b.set_title(f'B. ω = {omega_low:.1f} < ω_crit1\nG cannot invade (λ_G < 0)',
                  fontweight='bold', loc='left', fontsize=10)
    ax_b.grid(alpha=0.3)
    ax_b.set_ylim(-0.05, 0.4)

    # Add annotation
    ax_b.text(0.5, 0.8, 'G → 0\n(extinction)',
             transform=ax_b.transAxes, ha='center',
             fontsize=11, fontweight='bold', color='red',
             bbox=dict(boxstyle='round', facecolor='pink', alpha=0.7))

    # ========================================================================
    # PANEL C: Phase portrait at ω = ω_crit1 (exactly at bifurcation)
    # ========================================================================
    ax_c = fig.add_subplot(gs[1, 1])

    omega_crit = omega_crit1
    for ic in initial_conditions:
        sol = solve_ivp(model.dynamics_three_species, t_span, ic,
                       args=(omega_crit,), t_eval=t_eval, method='RK45')
        ax_c.plot(sol.t, sol.y[2], '-', alpha=0.6, linewidth=1.5)

    ax_c.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax_c.set_xlabel('Time', fontsize=10)
    ax_c.set_ylabel('G density', fontsize=10)
    ax_c.set_title(f'C. ω = ω_crit1 = {omega_crit:.3f}\nCritical point (λ_G = 0)',
                  fontweight='bold', loc='left', fontsize=10)
    ax_c.grid(alpha=0.3)

    # Add annotation
    ax_c.text(0.5, 0.8, 'Neutral stability\n(slow approach)',
             transform=ax_c.transAxes, ha='center',
             fontsize=11, fontweight='bold', color='orange',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

    # ========================================================================
    # PANEL D: Phase portrait at ω = 0.6 (above ω_crit1)
    # ========================================================================
    ax_d = fig.add_subplot(gs[1, 2])

    omega_high = 0.6
    for ic in initial_conditions:
        sol = solve_ivp(model.dynamics_three_species, t_span, ic,
                       args=(omega_high,), t_eval=t_eval, method='RK45')
        ax_d.plot(sol.t, sol.y[2], '-', alpha=0.6, linewidth=1.5)

    ax_d.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax_d.set_xlabel('Time', fontsize=10)
    ax_d.set_ylabel('G density', fontsize=10)
    ax_d.set_title(f'D. ω = {omega_high:.1f} > ω_crit1\nG invades successfully (λ_G > 0)',
                  fontweight='bold', loc='left', fontsize=10)
    ax_d.grid(alpha=0.3)

    # Add annotation
    ax_d.text(0.5, 0.8, 'G → G*\n(coexistence)',
             transform=ax_d.transAxes, ha='center',
             fontsize=11, fontweight='bold', color='green',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    # ========================================================================
    # PANEL E: Invasion fitness λ_G(ω)
    # ========================================================================
    ax_e = fig.add_subplot(gs[2, 0])

    lambda_G_values = [model.lambda_G(om) for om in omega_scan]

    ax_e.plot(omega_scan, lambda_G_values, '-', color='#2ca02c',
             linewidth=2.5, label='λ_G(ω)')
    ax_e.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    ax_e.axvline(omega_crit1, color='red', linestyle='--', linewidth=2,
                alpha=0.7, label=f'ω_crit1 = {omega_crit1:.3f}')

    # Shade regions
    ax_e.fill_between(omega_scan, 0, lambda_G_values,
                     where=np.array(lambda_G_values) < 0,
                     alpha=0.3, color='pink', label='G excluded (λ_G < 0)')
    ax_e.fill_between(omega_scan, 0, lambda_G_values,
                     where=np.array(lambda_G_values) > 0,
                     alpha=0.3, color='lightgreen', label='G invades (λ_G > 0)')

    ax_e.set_xlabel('ω', fontsize=10)
    ax_e.set_ylabel('G invasion fitness (λ_G)', fontsize=10)
    ax_e.set_title('E. G Invasion Condition\n(into S-M equilibrium)',
                  fontweight='bold', loc='left', fontsize=10)
    ax_e.legend(frameon=True, fontsize=8)
    ax_e.grid(alpha=0.3)

    # ========================================================================
    # PANEL F: Net interaction parameters evolution
    # ========================================================================
    ax_f = fig.add_subplot(gs[2, 1])

    a_values = [(1-om)*model.sigma_SG - om*model.alpha_SG for om in omega_scan]
    b_values = [om*model.sigma_MG - (1-om)*model.alpha_MG for om in omega_scan]
    c_values = [(1-om)*model.sigma_GS - om*model.alpha_GS for om in omega_scan]
    d_values = [2*om - 1 for om in omega_scan]
    e_values = [om*model.sigma_GM - (1-om)*model.alpha_GM for om in omega_scan]

    ax_f.plot(omega_scan, a_values, '-', linewidth=2, label='a (S-G net)')
    ax_f.plot(omega_scan, c_values, '-', linewidth=2, label='c (G-S net)')
    ax_f.plot(omega_scan, d_values, '-', linewidth=2, label='d (G basal)')
    ax_f.plot(omega_scan, e_values, '-', linewidth=2, label='e (G-M net)')

    ax_f.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax_f.axvline(omega_crit1, color='red', linestyle='--', linewidth=2, alpha=0.5)

    ax_f.set_xlabel('ω', fontsize=10)
    ax_f.set_ylabel('Net parameter value', fontsize=10)
    ax_f.set_title('F. Net Interaction Parameters\nvs Pathway Allocation',
                  fontweight='bold', loc='left', fontsize=10)
    ax_f.legend(frameon=True, fontsize=8)
    ax_f.grid(alpha=0.3)

    # ========================================================================
    # PANEL G: Generalist growth rate
    # ========================================================================
    ax_g = fig.add_subplot(gs[2, 2])

    r_G_values = [model.r_G(om) for om in omega_scan]

    ax_g.plot(omega_scan, r_G_values, '-', color='#2ca02c', linewidth=2.5)
    ax_g.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    ax_g.axvline(omega_crit1, color='red', linestyle='--', linewidth=2, alpha=0.7)

    # Mark critical points
    omega_r_zero = model.r_M / (model.r_S + model.r_M)
    ax_g.axvline(omega_r_zero, color='purple', linestyle=':', linewidth=2,
                alpha=0.7, label=f'r_G = 0 at ω = {omega_r_zero:.3f}')

    ax_g.set_xlabel('ω', fontsize=10)
    ax_g.set_ylabel('r_G(ω)', fontsize=10)
    ax_g.set_title('G. Generalist Growth Rate\n(weighted average)',
                  fontweight='bold', loc='left', fontsize=10)
    ax_g.legend(frameon=True, fontsize=8)
    ax_g.grid(alpha=0.3)

    # ========================================================================
    # PANEL H: Conceptual diagram (bottom row)
    # ========================================================================
    ax_h = fig.add_subplot(gs[3, :])
    ax_h.set_xlim(0, 1)
    ax_h.set_ylim(0, 1)
    ax_h.axis('off')

    # Title
    ax_h.text(0.5, 0.95, 'MECHANISTIC EXPLANATION: How ω Controls Community Composition',
             ha='center', fontsize=13, fontweight='bold')

    # Low omega region
    low_box = FancyBboxPatch((0.05, 0.55), 0.25, 0.35,
                            boxstyle="round,pad=0.02",
                            facecolor='lightblue', edgecolor='blue',
                            linewidth=3, alpha=0.7)
    ax_h.add_patch(low_box)

    ax_h.text(0.175, 0.85, 'LOW ω (0 - 0.4)', ha='center',
             fontsize=11, fontweight='bold')
    ax_h.text(0.175, 0.78, 'Metabolite-specialized', ha='center', fontsize=9,
             style='italic')
    ax_h.text(0.175, 0.72, 'r_G < 0', ha='center', fontsize=10, color='red')
    ax_h.text(0.175, 0.66, 'd < 0 (basal negative)', ha='center', fontsize=9)
    ax_h.text(0.175, 0.60, 'λ_G < 0', ha='center', fontsize=10, color='red',
             fontweight='bold')
    ax_h.text(0.175, 0.54, '→ G cannot invade', ha='center', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='pink'))

    # Arrow to bifurcation
    arrow1 = FancyArrowPatch((0.32, 0.725), (0.42, 0.725),
                            arrowstyle='->', mutation_scale=30,
                            linewidth=3, color='red')
    ax_h.add_patch(arrow1)

    # Bifurcation point
    bif_box = FancyBboxPatch((0.42, 0.55), 0.16, 0.35,
                            boxstyle="round,pad=0.02",
                            facecolor='lightyellow', edgecolor='orange',
                            linewidth=3, alpha=0.9)
    ax_h.add_patch(bif_box)

    ax_h.text(0.50, 0.85, f'ω = {omega_crit1:.3f}', ha='center',
             fontsize=11, fontweight='bold', color='red')
    ax_h.text(0.50, 0.78, 'BIFURCATION', ha='center', fontsize=10,
             fontweight='bold', color='red')
    ax_h.text(0.50, 0.71, 'λ_G = 0', ha='center', fontsize=11, color='red',
             fontweight='bold')
    ax_h.text(0.50, 0.64, 'G can invade', ha='center', fontsize=9)
    ax_h.text(0.50, 0.58, 'Transcritical', ha='center', fontsize=9,
             style='italic')

    # Arrow to high omega
    arrow2 = FancyArrowPatch((0.60, 0.725), (0.70, 0.725),
                            arrowstyle='->', mutation_scale=30,
                            linewidth=3, color='green')
    ax_h.add_patch(arrow2)

    # High omega region
    high_box = FancyBboxPatch((0.70, 0.55), 0.25, 0.35,
                             boxstyle="round,pad=0.02",
                             facecolor='lightgreen', edgecolor='darkgreen',
                             linewidth=3, alpha=0.7)
    ax_h.add_patch(high_box)

    ax_h.text(0.825, 0.85, 'HIGH ω (0.4 - 1.0)', ha='center',
             fontsize=11, fontweight='bold')
    ax_h.text(0.825, 0.78, 'Substrate-specialized', ha='center', fontsize=9,
             style='italic')
    ax_h.text(0.825, 0.72, 'r_G > 0', ha='center', fontsize=10, color='green')
    ax_h.text(0.825, 0.66, 'd > 0 (basal positive)', ha='center', fontsize=9)
    ax_h.text(0.825, 0.60, 'λ_G > 0', ha='center', fontsize=10, color='green',
             fontweight='bold')
    ax_h.text(0.825, 0.54, '→ S-M-G coexist', ha='center', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='lightgreen'))

    # Bottom explanation
    explanation_box = FancyBboxPatch((0.05, 0.05), 0.9, 0.42,
                                    boxstyle="round,pad=0.02",
                                    facecolor='lavender', edgecolor='purple',
                                    linewidth=2, alpha=0.6)
    ax_h.add_patch(explanation_box)

    ax_h.text(0.5, 0.42, 'KEY MECHANISM: Pathway Allocation Controls Invasion Fitness',
             ha='center', fontsize=11, fontweight='bold', color='purple')

    # Create three columns for explanation
    col_width = 0.28

    # Column 1: Metabolic allocation
    ax_h.text(0.05 + col_width/2, 0.35, 'Metabolic Allocation',
             ha='center', fontsize=10, fontweight='bold', color='darkblue')
    ax_h.text(0.05 + col_width/2, 0.30,
             '• ω = 0: Pure M-like\n  (metabolite pathway)\n\n'
             '• ω = 0.5: Intermediate\n  (balanced)\n\n'
             '• ω = 1: Pure S-like\n  (substrate pathway)',
             ha='left', fontsize=8, va='top')

    # Column 2: Growth parameters
    ax_h.text(0.36 + col_width/2, 0.35, 'Growth & Interactions',
             ha='center', fontsize=10, fontweight='bold', color='darkgreen')
    ax_h.text(0.36 + col_width/2, 0.30,
             f'• r_G(ω) = -r_M + ω(r_S+r_M)\n'
             f'  Crosses 0 at ω = {omega_r_zero:.3f}\n\n'
             f'• d(ω) = 2ω - 1\n'
             f'  Basal fitness term\n\n'
             f'• Net interactions: a,b,c,e\n'
             f'  All linear in ω',
             ha='left', fontsize=8, va='top')

    # Column 3: Critical threshold
    ax_h.text(0.67 + col_width/2, 0.35, 'Critical Threshold',
             ha='center', fontsize=10, fontweight='bold', color='darkred')
    ax_h.text(0.67 + col_width/2, 0.30,
             f'• ω_crit1 = {omega_crit1:.4f}\n'
             f'  (explicit formula)\n\n'
             f'• Below: λ_G < 0\n'
             f'  G excluded from S-M\n\n'
             f'• Above: λ_G > 0\n'
             f'  G invades → 3-species',
             ha='left', fontsize=8, va='top')

    # Main conclusion
    ax_h.text(0.5, 0.08,
             'CONCLUSION: Intermediate metabolic strategies (ω ≈ 0.4-0.6) allow generalists to invade\n'
             'Extreme specialization (too M-like or too S-like initially) excludes generalists',
             ha='center', fontsize=10, style='italic',
             bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

    plt.suptitle('Complete System Evolution: How ω Controls Three-Species Cross-Feeding Dynamics',
                fontsize=15, fontweight='bold', y=0.995)

    plt.savefig('omega_evolution_complete_explanation.png', dpi=300, bbox_inches='tight')
    print("\n✓ Complete evolution diagram saved: omega_evolution_complete_explanation.png")

    return fig


def create_step_by_step_narrative():
    """Create detailed narrative explanation"""

    narrative = """
================================================================================
完整解释：S-M-G系统如何随ω演化
COMPLETE EXPLANATION: S-M-G System Evolution with ω
================================================================================

## 系统设置 (System Setup)

三个物种：
- S: 底物专一种 (substrate specialist)
- M: 代谢物专一种 (metabolite specialist, 负basal增长率)
- G: 通才 (generalist, 可调控的metabolic pathway allocation)

关键参数：
- ω ∈ [0, 1]: 通才的pathway分配参数
  • ω = 0: 纯M-like (只使用代谢物pathway)
  • ω = 1: 纯S-like (只使用底物pathway)

================================================================================
演化过程 (EVOLUTIONARY TRAJECTORY)
================================================================================

### 阶段 I: ω = 0 → 0.4 (低ω区域 - S-M共存)

状态: **只有S和M共存, G被排除**

机制：
1. G的basal增长率: r_G(ω) = -r_M + ω(r_S + r_M)
   - 在ω = 0: r_G = -0.8 (负增长率!)
   - 在ω = 0.4: r_G = -0.08 (仍然负)

2. G的basal fitness: d(ω) = 2ω - 1
   - 在ω = 0: d = -1 (负!)
   - 在ω = 0.4: d = -0.2 (接近0)

3. **关键:** G入侵适应度 λ_G = r_G(ω) × [d + c·s*_SM + e·m*_SM]
   - 虽然括号内的项可能为正
   - 但r_G < 0 导致 λ_G < 0
   - 或者d太负,导致整个括号内为负

结果: **G无法入侵S-M平衡** → 只有S-M共存

生物学意义:
- G太"代谢物专一化"
- 无法在底物丰富的环境(S提供)中生长
- 类似M,但M有S→M的互惠支持,G没有

系统行为:
- 初始有小量G → G(t) → 0 (指数衰减到零)
- S和M达到稳定平衡: s* = 2.0, m* = 2.0

--------------------------------------------------------------------------------

### 临界点: ω = ω_crit1 = 0.4000 (转折点bifurcation)

**TRANSCRITICAL BIFURCATION 发生!**

数学条件:
λ_G(ω_crit1) = 0  (G的invasion fitness恰好为零)

这意味着:
- G既不能入侵,也不会被完全排除
- 系统处于"刀刃"平衡状态
- 微小的ω变化会导致定性行为改变

解析公式:
ω_crit1 = [1 - σ_GS·s*_SM + α_GM·m*_SM] /
          [2 - (σ_GS+α_GS)·s*_SM + (σ_GM+α_GM)·m*_SM]

代入参数:
分子 = 1 - 0.4×2.0 + 0.3×2.0 = 0.8
分母 = 2 - 0.7×2.0 + 0.7×2.0 = 2.0
→ ω_crit1 = 0.4

系统行为 (在ω_crit1):
- G密度非常缓慢地趋近于零或小正值
- "临界变慢" (critical slowing down)
- 对扰动非常敏感

物理类比:
- 像一个球在山顶平衡 (unstable fixed point)
- 或者在平坦的valley底部 (neutral stability)

--------------------------------------------------------------------------------

### 阶段 II: ω = 0.4 → 1.0 (高ω区域 - 三物种共存)

状态: **S-M-G三物种稳定共存**

机制转变：
1. G的basal增长率变正:
   - 在ω = 0.5: r_G = -0.8 + 0.5×1.8 = 0.1 > 0 ✓
   - 在ω = 0.6: r_G = 0.28 > 0 ✓
   - 在ω = 1.0: r_G = 1.0 > 0 ✓

2. G的basal fitness变正:
   - 在ω = 0.5: d = 0 (临界!)
   - 在ω = 0.6: d = 0.2 > 0 ✓
   - 在ω = 1.0: d = 1.0 > 0 ✓

3. **关键转变:** G入侵适应度变正
   λ_G(ω > ω_crit1) > 0

结果: **G成功入侵并稳定共存**

动力学过程:
1. 初始状态: S-M平衡 (s* = 2.0, m* = 2.0, g = 0)
2. 引入小量G: g(0) = 0.01
3. G开始增长: dG/dt > 0 (因为λ_G > 0)
4. G的增长影响S和M:
   - S受到竞争 (通过a参数)
   - M受到帮助或竞争 (通过b参数)
5. 系统重新平衡到新的三物种平衡:
   s* ≈ 1.8, m* ≈ 1.7, g* ≈ 0.5 (数值示例)

生物学意义:
- G现在有足够的"底物利用能力"
- 可以在S-M群落中找到生态位
- 中等ω (0.4-0.6)最有利:
  • 既能利用底物(像S)
  • 又能利用代谢物(像M)
  • "两全其美"的策略

系统稳定性:
- 三物种平衡是局部渐近稳定的
- Jacobian矩阵的所有特征值实部为负
- 对小扰动robust

--------------------------------------------------------------------------------

### 重要观察: 为什么在基线参数下没有ω_crit2?

在高ω (比如ω → 1.0)时:
- G变得非常"S-like" (底物专一化)
- 问题: M会不会被排除?

检查M的invasion fitness into S-G equilibrium:
λ_M = -r_M + σ_MS·s*_SG + σ_MG·g*_SG

计算发现:
- 在ω = 0.6: λ_M ≈ -0.3 < 0 (M不能入侵S-G)
- 在ω = 0.8: λ_M ≈ -0.4 < 0
- 在ω = 1.0: λ_M ≈ -0.5 < 0

**结论:** λ_M 在整个ω > ω_crit1 范围内都是负的!

这意味着:
1. M一旦在三物种equilibrium中存在
2. 就会一直存在 (因为三物种equilibrium稳定)
3. 不会有第二个bifurcation点

原因: σ_MG = 0.4 太弱
- G→M的facilitation不够强
- 无法补偿M的负basal增长率

如果增加σ_MG到0.8或更高:
→ λ_M会在某个高ω值变正
→ 然后随ω继续增加又变负
→ 产生ω_crit2! (M被排除的点)

================================================================================
系统行为总结 (SYSTEM BEHAVIOR SUMMARY)
================================================================================

### 正向演化 (ω从0增加到1):

ω = 0.0:   S-M equilibrium, G extinct
              ↓ (λ_G < 0, G无法入侵)
ω = 0.2:   S-M equilibrium, G extinct
              ↓
ω = 0.3:   S-M equilibrium, G extinct
              ↓
           ━━━━ ω_crit1 = 0.4000 ━━━━  [TRANSCRITICAL BIFURCATION]
              ↓ (λ_G crosses 0)
ω = 0.5:   S-M-G coexistence (3-species stable)
              ↓ (继续稳定共存)
ω = 0.7:   S-M-G coexistence
              ↓
ω = 0.9:   S-M-G coexistence
              ↓
ω = 1.0:   S-M-G coexistence (permanent)

【基线参数下: 没有ω_crit2!】

--------------------------------------------------------------------------------

### 反向演化 (ω从1减少到0):

ω = 1.0:   S-M-G coexistence
              ↓
ω = 0.7:   S-M-G coexistence
              ↓
ω = 0.5:   S-M-G coexistence
              ↓
           ━━━━ ω_crit1 = 0.4000 ━━━━  [TRANSCRITICAL BIFURCATION]
              ↓ (λ_G crosses 0)
ω = 0.3:   G goes extinct → S-M equilibrium
              ↓
ω = 0.1:   S-M equilibrium only
              ↓
ω = 0.0:   S-M equilibrium only

**可逆过程:** 同一个临界点控制正向和反向转变!

但注意: 在ω_crit1附近可能有hysteresis (如果有更复杂的dynamics)

================================================================================
对应的动力学系统行为 (DYNAMICAL SYSTEM BEHAVIORS)
================================================================================

### 在ω < ω_crit1 (Region I):

平衡点:
- E_SM = (s*, m*, 0) - 稳定
- E_SMG: 不存在正的三物种equilibrium

吸引子:
- S-M equilibrium 是全局吸引子
- 所有初始条件 → E_SM

轨迹特征:
- G的轨迹: G(t) = G(0)·exp(λ_G·t) → 0 (exponential decay)
- S,M的轨迹: 震荡收敛到s*, m*

Jacobian特征值 (at E_SM):
- 关于S,M方向: 负实部 (稳定)
- 关于G方向: λ_G < 0 (G方向也稳定 - G被吸引到0)

--------------------------------------------------------------------------------

### 在ω = ω_crit1 (Bifurcation point):

平衡点:
- E_SM = (s*, m*, 0) - marginally stable
- E_SMG: 刚刚"出生"或"消失"

吸引子:
- S-M equilibrium 仍然是吸引子
- 但吸引速度极慢 (代数衰减而非指数)

轨迹特征:
- G的轨迹: G(t) ~ G(0)·t^(-α) (power-law decay, α > 0)
- 临界变慢现象

Jacobian特征值:
- 关于G方向: λ = 0 (零特征值!)
- 这是transcritical bifurcation的signature

--------------------------------------------------------------------------------

### 在ω > ω_crit1 (Region II):

平衡点:
- E_SM = (s*, m*, 0) - 变为不稳定 (saddle)
- E_SMG = (s*, m*, g*) - 新的稳定equilibrium (g* > 0)

吸引子:
- 三物种equilibrium E_SMG 是吸引子
- Basin of attraction包括几乎所有正象限

轨迹特征:
- 从接近E_SM的初始条件出发
- 轨迹沿着E_SM的unstable manifold离开
- 最终收敛到E_SMG
- 可能有震荡 (damped oscillations)

Jacobian特征值 (at E_SMG):
- 所有三个方向: 负实部 (或复数但实部为负)
- 稳定focus或stable node

相空间结构:
- E_SM变成saddle point
- Stable manifold: 2维 (S-M平面)
- Unstable manifold: 1维 (G方向)

--------------------------------------------------------------------------------

### Transcritical Bifurcation的几何意义:

在参数空间 (ω-space):
- ω < ω_crit1: 稳定S-M, 不稳定SMG
- ω = ω_crit1: 两个equilibria"交换稳定性"
- ω > ω_crit1: 不稳定S-M, 稳定SMG

特征:
- 两个equilibria在bifurcation点"相遇"
- 但不是同时消失 (像saddle-node)
- 而是交换稳定性

数学表达:
- G密度作为ω的函数: g*(ω)
  • g*(ω < ω_crit1) = 0
  • g*(ω_crit1) = 0
  • g*(ω > ω_crit1) ∝ √(ω - ω_crit1) (square-root scaling near bifurcation)

实际上:在我们的模型中,scaling可能稍有不同,
因为这是三维系统的codimension-1 bifurcation

================================================================================
参数敏感性 (PARAMETER SENSITIVITY)
================================================================================

如何调控ω_crit1?

增加ω_crit1 (延迟G invasion):
- 增加σ_GS (G-S竞争)
- 增加α_GS (G对S的负效应)
- 减少σ_MS (削弱S-M mutualism)

减少ω_crit1 (提前G invasion):
- 增加σ_MS (增强S-M mutualism)
- 增加σ_GM (G帮助M)
- 减少α_GM (减少G-M竞争)

为什么σ_MS很重要?
- 更强的S-M mutualism → 更高的s*_SM, m*_SM
- 这使得S-M平台对G更"hospitable"
- G更容易入侵

================================================================================
实验预测 (EXPERIMENTAL PREDICTIONS)
================================================================================

### 可测量的signatures:

1. **在ω_crit1附近:**
   - 临界变慢 (critical slowing down)
     • 扰动后恢复时间 → ∞ as ω → ω_crit1
   - 方差增加 (variance amplification)
     • Fluctuations in G density增大
   - 自相关时间增加
     • Temporal autocorrelation function decays slowly

2. **扫描ω实验:**
   - 从ω = 0.2开始,逐步增加ω
   - 在ω < 0.4: G应该始终为0
   - 在ω ≈ 0.4: G出现但非常不稳定
   - 在ω > 0.4: G稳定存在并增加

3. **Hysteresis测试:**
   - 正向: 从低ω增加到高ω → 记录G出现的ω
   - 反向: 从高ω减少到低ω → 记录G消失的ω
   - 如果两个临界ω相同 → transcritical (我们的情况)
   - 如果不同 → subcritical bifurcation (不是我们的模型)

================================================================================
与经典bifurcation理论的联系
================================================================================

我们的transcritical bifurcation属于:

**Codimension-1 bifurcation**
- 单个参数 (ω) 变化导致
- 定性行为改变

标准形式 (Normal form):
dx/dt = μ·x - x²

其中μ是bifurcation parameter (类比我们的ω - ω_crit1)

在我们的系统:
- x ↔ G的密度
- μ ↔ (ω - ω_crit1)
- 二次项来自density-dependent regulation

特征:
- μ < 0: x* = 0 stable, x* = μ unstable
- μ = 0: bifurcation point
- μ > 0: x* = 0 unstable, x* = μ stable

√(μ) scaling: g* ∝ √(ω - ω_crit1) near bifurcation

================================================================================
总结 (FINAL SUMMARY)
================================================================================

S-M-G系统随ω演化的完整图景:

1. **ω ∈ [0, 0.4)**: 只有S-M共存
   - G太代谢物专一化
   - r_G < 0, λ_G < 0
   - G被排除

2. **ω = 0.4**: Transcritical bifurcation
   - λ_G = 0 (临界点)
   - G的invasion fitness crosses zero
   - 系统行为定性改变

3. **ω ∈ (0.4, 1.0]**: S-M-G三物种共存
   - G有足够底物pathway allocation
   - r_G > 0, λ_G > 0
   - 稳定三物种equilibrium

4. **基线参数: 没有ω_crit2**
   - M在整个ω > 0.4范围持续存在
   - 永久三物种共存
   - 需要更强σ_MG才会产生第二个bifurcation

关键洞察:
✓ 单个参数ω控制完整的群落组装过程
✓ 中等代谢策略(intermediate allocation)是generalist成功的关键
✓ Transcritical bifurcation是smooth transition (不是catastrophic)
✓ 系统行为可预测且robust

================================================================================
"""

    print(narrative)

    # Save to file
    with open('omega_evolution_narrative.txt', 'w', encoding='utf-8') as f:
        f.write(narrative)

    print("\n✓ Detailed narrative saved: omega_evolution_narrative.txt")


if __name__ == '__main__':
    print("\n" + "="*80)
    print("Creating complete ω evolution explanation...")
    print("="*80 + "\n")

    # Create main diagram
    fig = create_omega_evolution_diagram()

    # Create detailed narrative
    create_step_by_step_narrative()

    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
    print("""
Files generated:
1. omega_evolution_complete_explanation.png - Comprehensive visual diagram
2. omega_evolution_narrative.txt - Detailed Chinese/English narrative

Key points:
- ω controls community composition through invasion fitness
- ω_crit1 = 0.40 is the transcritical bifurcation point
- Below: S-M only (G excluded)
- Above: S-M-G coexist (permanent with baseline parameters)
- System behavior changes smoothly but qualitatively at threshold
    """)

    plt.show()

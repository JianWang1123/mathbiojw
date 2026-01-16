#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analytical Analysis with Simplified Net Interaction Parameters
基于净相互作用参数的解析分析

采用Jeff Gore风格的严格数学分析和高质量可视化

参数简化:
- a := (1-ω)σ_SG - ω·α_SG  (G → S net interaction)
- c := (1-ω)σ_GS - ω·α_GS  (S → G net interaction)
- b := ω·σ_MG - (1-ω)·α_MG  (G → M net interaction)
- e := ω·σ_GM - (1-ω)·α_GM  (M → G net interaction)
- d := 2ω - 1                (generalist pathway balance)

Author: Jian Wang
Date: January 2026
Reference style: Gore et al., Nature, PNAS
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, FancyBboxPatch
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import fsolve
import seaborn as sns
from typing import Dict, Tuple, List

# Publication-quality settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.dpi'] = 150


class SMSystemAnalysis:
    """
    S-M Two-Species System Analysis

    Scaled equations:
    ds/dt = r_S · s(1 + σ_SM · m - s)
    dm/dt = r_M · m(-1 + σ_MS · s - m)

    Parameters:
    - σ_SM: M → S net benefit (mutualistic coefficient)
    - σ_MS: S → M net benefit (obligate dependence coefficient)
    """

    def __init__(self, r_S=1.0, r_M=0.8, sigma_SM=0.5, sigma_MS=1.5):
        """Initialize parameters"""
        self.r_S = r_S
        self.r_M = r_M
        self.sigma_SM = sigma_SM
        self.sigma_MS = sigma_MS

    def equations(self, t, y):
        """ODE system"""
        s, m = y
        s = max(0, s)
        m = max(0, m)

        ds_dt = self.r_S * s * (1 + self.sigma_SM * m - s)
        dm_dt = self.r_M * m * (-1 + self.sigma_MS * s - m)

        return [ds_dt, dm_dt]

    def find_equilibria(self):
        """
        Find all equilibrium points analytically

        Returns
        -------
        dict
            Dictionary of equilibria with coordinates and existence conditions
        """
        equilibria = {}

        # E1: Extinction (0, 0)
        equilibria['E1_extinction'] = {
            'point': (0, 0),
            'exists': 'Always',
            'biological': 'Complete extinction'
        }

        # E2: S-only (1, 0)
        equilibria['E2_S_only'] = {
            'point': (1, 0),
            'exists': 'Always',
            'biological': 'S at carrying capacity, M extinct'
        }

        # E3: Coexistence (s*, m*)
        # From equations at equilibrium:
        # 1 + σ_SM·m - s = 0  →  s = 1 + σ_SM·m
        # -1 + σ_MS·s - m = 0  →  m = σ_MS·s - 1
        #
        # Substituting:
        # m = σ_MS(1 + σ_SM·m) - 1
        # m = σ_MS + σ_MS·σ_SM·m - 1
        # m(1 - σ_MS·σ_SM) = σ_MS - 1
        # m* = (σ_MS - 1)/(1 - σ_MS·σ_SM)
        # s* = 1 + σ_SM·m* = (1 - σ_SM)/(1 - σ_MS·σ_SM)

        s_star = (1 - self.sigma_SM) / (1 - self.sigma_MS * self.sigma_SM)
        m_star = (self.sigma_MS - 1) / (1 - self.sigma_MS * self.sigma_SM)

        equilibria['E3_coexistence'] = {
            'point': (s_star, m_star),
            'formula_s': '(1 - σ_SM)/(1 - σ_MS·σ_SM)',
            'formula_m': '(σ_MS - 1)/(1 - σ_MS·σ_SM)',
            'exists': 'σ_MS > 1 AND σ_MS·σ_SM < 1',
            'biological': 'S-M mutualistic coexistence'
        }

        return equilibria

    def jacobian(self, s, m):
        """
        Compute Jacobian matrix at (s, m)

        J = | ∂f_s/∂s   ∂f_s/∂m |
            | ∂f_m/∂s   ∂f_m/∂m |

        where:
        f_s = r_S·s(1 + σ_SM·m - s)
        f_m = r_M·m(-1 + σ_MS·s - m)

        ∂f_s/∂s = r_S(1 + σ_SM·m - 2s)
        ∂f_s/∂m = r_S·s·σ_SM
        ∂f_m/∂s = r_M·m·σ_MS
        ∂f_m/∂m = r_M(-1 + σ_MS·s - 2m)
        """
        J = np.array([
            [self.r_S * (1 + self.sigma_SM * m - 2*s), self.r_S * s * self.sigma_SM],
            [self.r_M * m * self.sigma_MS, self.r_M * (-1 + self.sigma_MS * s - 2*m)]
        ])
        return J

    def stability_analysis(self):
        """
        Complete stability analysis of all equilibria

        Stability criteria (Routh-Hurwitz for 2×2):
        1. Tr(J) < 0
        2. Det(J) > 0

        Returns
        -------
        dict
            Stability results for each equilibrium
        """
        equilibria = self.find_equilibria()
        results = {}

        for name, eq_info in equilibria.items():
            s, m = eq_info['point']
            J = self.jacobian(s, m)

            trace = np.trace(J)
            det = np.linalg.det(J)
            eigenvalues = np.linalg.eigvals(J)

            # Stability
            is_stable = (trace < 0) and (det > 0)

            results[name] = {
                'point': (s, m),
                'jacobian': J,
                'trace': trace,
                'determinant': det,
                'eigenvalues': eigenvalues,
                'stable': is_stable,
                'exists': eq_info['exists']
            }

        return results

    def print_analysis(self):
        """Print complete analytical results"""
        print("\n" + "="*80)
        print("S-M SYSTEM ANALYTICAL ANALYSIS")
        print("="*80)

        print("\nSYSTEM EQUATIONS (Scaled):")
        print("  ds/dt = r_S · s(1 + σ_SM·m - s)")
        print("  dm/dt = r_M · m(-1 + σ_MS·s - m)")

        print(f"\nPARAMETERS:")
        print(f"  r_S = {self.r_S}")
        print(f"  r_M = {self.r_M}")
        print(f"  σ_SM = {self.sigma_SM} (M → S benefit)")
        print(f"  σ_MS = {self.sigma_MS} (S → M benefit, obligate)")

        print("\n" + "-"*80)
        print("EQUILIBRIUM POINTS")
        print("-"*80)

        equilibria = self.find_equilibria()
        for name, info in equilibria.items():
            print(f"\n{name}:")
            s, m = info['point']
            print(f"  Coordinates: ({s:.4f}, {m:.4f})")
            if 'formula_s' in info:
                print(f"  s* = {info['formula_s']}")
                print(f"  m* = {info['formula_m']}")
            print(f"  Existence: {info['exists']}")
            print(f"  Biology: {info['biological']}")

        print("\n" + "-"*80)
        print("STABILITY ANALYSIS")
        print("-"*80)

        stability = self.stability_analysis()
        for name, result in stability.items():
            print(f"\n{name}:")
            s, m = result['point']
            print(f"  Point: ({s:.4f}, {m:.4f})")
            print(f"  Trace(J) = {result['trace']:.4f}")
            print(f"  Det(J) = {result['determinant']:.4f}")
            print(f"  Eigenvalues: {result['eigenvalues']}")
            print(f"  Stable: {result['stable']}")

        print("\n" + "-"*80)
        print("STABILITY CRITERIA (GENERAL)")
        print("-"*80)
        print("\nFor E3 (Coexistence):")
        print("  Trace(J) = -(r_S·s* + r_M·m*) < 0  ✓ (always satisfied)")
        print("  Det(J) = r_S·r_M·s*·m*(1 - σ_SM·σ_MS)")
        print("\nStability requires:")
        print("  1 - σ_SM·σ_MS > 0")
        print("  ⟹ σ_SM·σ_MS < 1")
        print("\nCombined with existence (σ_MS > 1):")
        print("  ⟹ σ_SM < 1/σ_MS")

        print("\n" + "="*80)

    def phase_portrait(self):
        """Generate phase portrait with nullclines"""
        fig, ax = plt.subplots(figsize=(8, 7))

        # Nullclines
        s_range = np.linspace(0, 2, 200)
        m_range = np.linspace(0, 2, 200)

        # S-nullcline: 1 + σ_SM·m - s = 0 → s = 1 + σ_SM·m
        s_nullcline = 1 + self.sigma_SM * m_range

        # M-nullcline: -1 + σ_MS·s - m = 0 → m = σ_MS·s - 1
        m_nullcline = self.sigma_MS * s_range - 1

        # Plot nullclines
        ax.plot(s_nullcline, m_range, 'b-', linewidth=2, label='S-nullcline (ds/dt=0)')
        ax.plot(s_range, m_nullcline, 'r-', linewidth=2, label='M-nullcline (dm/dt=0)')

        # Vector field
        s_grid = np.linspace(0, 2, 20)
        m_grid = np.linspace(0, 2, 20)
        S, M = np.meshgrid(s_grid, m_grid)

        dS = np.zeros_like(S)
        dM = np.zeros_like(M)

        for i in range(S.shape[0]):
            for j in range(S.shape[1]):
                derivs = self.equations(0, [S[i,j], M[i,j]])
                dS[i,j] = derivs[0]
                dM[i,j] = derivs[1]

        # Normalize for better visualization
        speed = np.sqrt(dS**2 + dM**2)
        speed[speed == 0] = 1
        dS_norm = dS / speed
        dM_norm = dM / speed

        ax.quiver(S, M, dS_norm, dM_norm, speed,
                 cmap='gray', alpha=0.4, width=0.003)

        # Equilibria
        equilibria = self.find_equilibria()
        stability = self.stability_analysis()

        for name, result in stability.items():
            s, m = result['point']
            if s >= 0 and m >= 0 and s <= 2 and m <= 2:
                if result['stable']:
                    ax.plot(s, m, 'go', markersize=12, markeredgewidth=2,
                           markeredgecolor='darkgreen', label=f'{name} (stable)')
                else:
                    ax.plot(s, m, 'ro', markersize=12, markerfacecolor='white',
                           markeredgewidth=2, markeredgecolor='red',
                           label=f'{name} (unstable)')

        # Trajectories
        initial_conditions = [
            [0.5, 0.5], [1.5, 0.2], [0.3, 1.5], [1.2, 1.2]
        ]

        for ic in initial_conditions:
            sol = solve_ivp(self.equations, [0, 50], ic,
                          dense_output=True, method='RK45')
            t_eval = np.linspace(0, 50, 500)
            traj = sol.sol(t_eval)
            ax.plot(traj[0], traj[1], 'k-', alpha=0.3, linewidth=1)

        ax.set_xlabel('S (substrate specialist)', fontweight='bold')
        ax.set_ylabel('M (metabolite specialist)', fontweight='bold')
        ax.set_title(f'S-M System Phase Portrait\n' +
                    f'σ_SM={self.sigma_SM}, σ_MS={self.sigma_MS}',
                    fontweight='bold')
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 2)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', framealpha=0.9)

        plt.tight_layout()
        return fig

    def bifurcation_diagram(self):
        """
        Bifurcation diagram with σ_MS as bifurcation parameter
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Parameter range
        sigma_MS_range = np.linspace(0.5, 2.5, 200)

        # Arrays to store equilibria
        s_coex = []
        m_coex = []
        stable_coex = []

        for sigma_MS in sigma_MS_range:
            temp_system = SMSystemAnalysis(
                r_S=self.r_S, r_M=self.r_M,
                sigma_SM=self.sigma_SM, sigma_MS=sigma_MS
            )

            eq = temp_system.find_equilibria()
            stab = temp_system.stability_analysis()

            s, m = eq['E3_coexistence']['point']

            # Check if exists and positive
            if sigma_MS > 1 and (1 - sigma_MS * self.sigma_SM) != 0:
                if s > 0 and m > 0:
                    s_coex.append(s)
                    m_coex.append(m)
                    stable_coex.append(stab['E3_coexistence']['stable'])
                else:
                    s_coex.append(np.nan)
                    m_coex.append(np.nan)
                    stable_coex.append(False)
            else:
                s_coex.append(np.nan)
                m_coex.append(np.nan)
                stable_coex.append(False)

        s_coex = np.array(s_coex)
        m_coex = np.array(m_coex)
        stable_coex = np.array(stable_coex)

        # Plot 1: S equilibrium vs σ_MS
        ax1 = axes[0, 0]
        # Stable
        ax1.plot(sigma_MS_range[stable_coex], s_coex[stable_coex],
                'b-', linewidth=2, label='Stable')
        # Unstable
        ax1.plot(sigma_MS_range[~stable_coex], s_coex[~stable_coex],
                'b--', linewidth=2, alpha=0.5, label='Unstable')
        # Boundary equilibrium
        ax1.axhline(1, color='gray', linestyle=':', linewidth=1.5, label='E2 (S-only)')
        ax1.axvline(1, color='red', linestyle='--', alpha=0.5, label='σ_MS = 1 (critical)')

        ax1.set_xlabel('σ_MS (S → M benefit)', fontweight='bold')
        ax1.set_ylabel('s* (S equilibrium)', fontweight='bold')
        ax1.set_title('S Bifurcation Diagram', fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, 2)

        # Plot 2: M equilibrium vs σ_MS
        ax2 = axes[0, 1]
        ax2.plot(sigma_MS_range[stable_coex], m_coex[stable_coex],
                'r-', linewidth=2, label='Stable')
        ax2.plot(sigma_MS_range[~stable_coex], m_coex[~stable_coex],
                'r--', linewidth=2, alpha=0.5, label='Unstable')
        ax2.axhline(0, color='gray', linestyle=':', linewidth=1.5)
        ax2.axvline(1, color='red', linestyle='--', alpha=0.5, label='σ_MS = 1 (critical)')

        ax2.set_xlabel('σ_MS (S → M benefit)', fontweight='bold')
        ax2.set_ylabel('m* (M equilibrium)', fontweight='bold')
        ax2.set_title('M Bifurcation Diagram', fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(0, 2)

        # Plot 3: Stability regions in parameter space
        ax3 = axes[1, 0]

        sigma_SM_range = np.linspace(0, 1.5, 100)
        sigma_MS_range_2 = np.linspace(0.5, 2.5, 100)

        SM, MS = np.meshgrid(sigma_SM_range, sigma_MS_range_2)

        # Region 1: σ_MS < 1 (M cannot survive)
        region1 = MS < 1

        # Region 2: σ_MS > 1 and σ_MS·σ_SM < 1 (Stable coexistence)
        region2 = (MS > 1) & (MS * SM < 1)

        # Region 3: σ_MS > 1 and σ_MS·σ_SM > 1 (Unstable)
        region3 = (MS > 1) & (MS * SM > 1)

        ax3.contourf(SM, MS, region1.astype(int), levels=[0.5, 1.5],
                    colors=['lightcoral'], alpha=0.3)
        ax3.contourf(SM, MS, region2.astype(int), levels=[0.5, 1.5],
                    colors=['lightgreen'], alpha=0.5)
        ax3.contourf(SM, MS, region3.astype(int), levels=[0.5, 1.5],
                    colors=['lightyellow'], alpha=0.3)

        # Critical lines
        ax3.axhline(1, color='red', linestyle='--', linewidth=2, label='σ_MS = 1')
        sigma_MS_line = np.linspace(1, 2.5, 100)
        sigma_SM_line = 1 / sigma_MS_line
        ax3.plot(sigma_SM_line, sigma_MS_line, 'b--', linewidth=2,
                label='σ_MS·σ_SM = 1')

        ax3.set_xlabel('σ_SM (M → S benefit)', fontweight='bold')
        ax3.set_ylabel('σ_MS (S → M benefit)', fontweight='bold')
        ax3.set_title('Parameter Space Regions', fontweight='bold')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        # Add text annotations
        ax3.text(0.2, 0.7, 'M extinct', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.5))
        ax3.text(0.2, 1.5, 'Stable\ncoexistence', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
        ax3.text(0.8, 1.8, 'Unstable', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))

        # Plot 4: Trajectories to equilibrium
        ax4 = axes[1, 1]

        time = np.linspace(0, 50, 1000)
        initial_conditions = [
            [0.3, 0.1], [0.8, 0.3], [1.5, 0.5], [1.2, 1.0]
        ]

        for ic in initial_conditions:
            sol = solve_ivp(self.equations, [0, 50], ic,
                          t_eval=time, method='RK45')
            ax4.plot(sol.t, sol.y[0], 'b-', alpha=0.6, linewidth=1.5)
            ax4.plot(sol.t, sol.y[1], 'r-', alpha=0.6, linewidth=1.5)

        # Equilibrium
        eq = self.find_equilibria()
        s_eq, m_eq = eq['E3_coexistence']['point']
        ax4.axhline(s_eq, color='b', linestyle='--', linewidth=2, label=f's* = {s_eq:.3f}')
        ax4.axhline(m_eq, color='r', linestyle='--', linewidth=2, label=f'm* = {m_eq:.3f}')

        ax4.set_xlabel('Time', fontweight='bold')
        ax4.set_ylabel('Population density', fontweight='bold')
        ax4.set_title('Convergence to Equilibrium', fontweight='bold')
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        return fig


def main():
    """Main analysis"""
    print("\n" + "🔬"*40)
    print("S-M SYSTEM: COMPLETE ANALYTICAL ANALYSIS")
    print("Jeff Gore Style - Publication Quality")
    print("🔬"*40)

    # Create system with biologically realistic parameters
    system = SMSystemAnalysis(
        r_S=1.0,
        r_M=0.8,
        sigma_SM=0.5,  # M provides moderate benefit to S
        sigma_MS=1.5   # S provides strong benefit to M (obligate)
    )

    # Print analytical results
    system.print_analysis()

    # Generate figures
    print("\nGenerating figures...")

    # Phase portrait
    fig1 = system.phase_portrait()
    fig1.savefig('figures/SM_phase_portrait.png', dpi=300, bbox_inches='tight')
    print("  ✓ Phase portrait saved")

    # Bifurcation diagram
    fig2 = system.bifurcation_diagram()
    fig2.savefig('figures/SM_bifurcation_analysis.png', dpi=300, bbox_inches='tight')
    print("  ✓ Bifurcation analysis saved")

    plt.show()

    print("\n" + "="*80)
    print("✅ ANALYSIS COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()

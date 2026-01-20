#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Complete Analysis for Publication
三物种交叉喂养系统的完整分析

生成所有论文需要的图片：
1. Section 1: 成对相互作用 (S-M, S-G, M-G)
2. Section 2: 三物种系统分析

Author: Jian Wang
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import seaborn as sns

# Publication settings
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.dpi'] = 300

# Color palette
COLORS = {
    'S': '#1f77b4',  # blue
    'M': '#ff7f0e',  # orange
    'G': '#2ca02c',  # green
    'stable': '#2ca02c',
    'unstable': '#d62728'
}


# ============================================================================
# SECTION 1.1: S-M SYSTEM
# ============================================================================

class SMSystem:
    """S-M two-species system"""

    def __init__(self, r_S=1.0, r_M=0.8, sigma_SM=0.5, sigma_MS=1.5):
        self.r_S = r_S
        self.r_M = r_M
        self.sigma_SM = sigma_SM
        self.sigma_MS = sigma_MS

    def equations(self, t, y):
        s, m = np.maximum(y, 0)
        ds = self.r_S * s * (1 + self.sigma_SM * m - s)
        dm = self.r_M * m * (-1 + self.sigma_MS * s - m)
        return [ds, dm]

    def equilibria(self):
        """Return all equilibria"""
        eq = {}
        eq['E1'] = (0, 0)
        eq['E2'] = (1, 0)

        # Coexistence
        denom = 1 - self.sigma_MS * self.sigma_SM
        if abs(denom) > 1e-10:
            s_star = (1 - self.sigma_SM) / denom
            m_star = (self.sigma_MS - 1) / denom
            eq['E3'] = (s_star, m_star)

        return eq

    def jacobian(self, s, m):
        J = np.array([
            [self.r_S * (1 + self.sigma_SM * m - 2*s), self.r_S * s * self.sigma_SM],
            [self.r_M * m * self.sigma_MS, self.r_M * (-1 + self.sigma_MS * s - 2*m)]
        ])
        return J

    def is_stable(self, s, m):
        J = self.jacobian(s, m)
        eigs = np.linalg.eigvals(J)
        return np.all(np.real(eigs) < 0)


# ============================================================================
# SECTION 1.2: S-G SYSTEM
# ============================================================================

class SGSystem:
    """S-G two-species system with net interaction parameters"""

    def __init__(self, r_S=1.0, r_G=0.9, a=0.2, c=0.3, d=0.0):
        self.r_S = r_S
        self.r_G = r_G
        self.a = a  # G -> S net effect
        self.c = c  # S -> G net effect
        self.d = d  # 2ω - 1

    def equations(self, t, y):
        s, g = np.maximum(y, 0)
        ds = self.r_S * s * (1 + self.a * g - s)
        dg = self.r_G * g * (self.d + self.c * s - g)
        return [ds, dg]

    def equilibria(self):
        eq = {}
        eq['E1'] = (0, 0)
        eq['E2'] = (1, 0)

        # G-only (if d > 0)
        if self.d > 0:
            eq['E3_G'] = (0, self.d)

        # Coexistence
        denom = 1 - self.a * self.c
        if abs(denom) > 1e-10:
            g_star = (self.d + self.c) / denom
            s_star = (1 + self.a * self.d) / denom
            eq['E4_SG'] = (s_star, g_star)

        return eq

    def jacobian(self, s, g):
        J = np.array([
            [self.r_S * (1 + self.a * g - 2*s), self.r_S * s * self.a],
            [self.r_G * g * self.c, self.r_G * (self.d + self.c * s - 2*g)]
        ])
        return J

    def is_stable(self, s, g):
        J = self.jacobian(s, g)
        eigs = np.linalg.eigvals(J)
        return np.all(np.real(eigs) < 0)


# ============================================================================
# SECTION 1.3: M-G SYSTEM
# ============================================================================

class MGSystem:
    """M-G two-species system"""

    def __init__(self, r_M=0.8, r_G=0.9, b=0.3, e=0.2, d=0.0):
        self.r_M = r_M
        self.r_G = r_G
        self.b = b  # G -> M net effect
        self.e = e  # M -> G net effect
        self.d = d  # 2ω - 1

    def equations(self, t, y):
        m, g = np.maximum(y, 0)
        dm = self.r_M * m * (-1 + self.b * g - m)
        dg = self.r_G * g * (self.d + self.e * m - g)
        return [dm, dg]

    def equilibria(self):
        eq = {}
        eq['E1'] = (0, 0)

        # M cannot exist alone (base rate = -1)
        # G-only (if d > 0)
        if self.d > 0:
            eq['E2_G'] = (0, self.d)

        # Coexistence (if M can survive with G's help)
        denom = 1 - self.b * self.e
        if abs(denom) > 1e-10:
            m_star = (self.b * self.d - 1) / denom
            g_star = (self.d + self.e) / denom
            eq['E3_MG'] = (m_star, g_star)

        return eq

    def jacobian(self, m, g):
        J = np.array([
            [self.r_M * (-1 + self.b * g - 2*m), self.r_M * m * self.b],
            [self.r_G * g * self.e, self.r_G * (self.d + self.e * m - 2*g)]
        ])
        return J

    def is_stable(self, m, g):
        J = self.jacobian(m, g)
        eigs = np.linalg.eigvals(J)
        return np.all(np.real(eigs) < 0)


# ============================================================================
# SECTION 2: THREE-SPECIES SYSTEM
# ============================================================================

class ThreeSpeciesSystem:
    """Complete three-species system"""

    def __init__(self, r_S=1.0, r_M=0.8, r_G=0.9,
                 sigma_SM=0.5, sigma_MS=1.5,
                 a=0.2, b=0.3, c=0.3, d=0.0, e=0.2):
        self.r_S = r_S
        self.r_M = r_M
        self.r_G = r_G
        self.sigma_SM = sigma_SM
        self.sigma_MS = sigma_MS
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e

    def equations(self, t, y):
        s, m, g = np.maximum(y, 0)

        ds = self.r_S * s * (1 + self.sigma_SM * m + self.a * g - s)
        dm = self.r_M * m * (-1 + self.sigma_MS * s + self.b * g - m)
        dg = self.r_G * g * (self.d + self.c * s + self.e * m - g)

        return [ds, dm, dg]

    def find_interior_equilibrium(self):
        """Find three-species coexistence equilibrium numerically"""
        def system(X):
            s, m, g = X
            return [
                1 + self.sigma_SM * m + self.a * g - s,
                -1 + self.sigma_MS * s + self.b * g - m,
                self.d + self.c * s + self.e * m - g
            ]

        # Try multiple initial guesses
        for guess in [[0.5, 0.5, 0.5], [1.0, 1.0, 1.0], [0.8, 0.8, 0.8]]:
            try:
                sol = fsolve(system, guess)
                if np.all(sol > 0) and np.linalg.norm(system(sol)) < 1e-6:
                    return sol
            except:
                continue

        return None

    def G_invasion_rate_SM(self):
        """Compute G's invasion rate at S-M equilibrium"""
        # S-M equilibrium
        denom = 1 - self.sigma_MS * self.sigma_SM
        s_SM = (1 - self.sigma_SM) / denom
        m_SM = (self.sigma_MS - 1) / denom

        # G's growth rate at (s_SM, m_SM, 0)
        invasion_rate = self.d + self.c * s_SM + self.e * m_SM

        return invasion_rate, s_SM, m_SM

    def jacobian(self, s, m, g):
        J = np.array([
            [self.r_S * (1 + self.sigma_SM*m + self.a*g - 2*s),
             self.r_S * s * self.sigma_SM,
             self.r_S * s * self.a],
            [self.r_M * m * self.sigma_MS,
             self.r_M * (-1 + self.sigma_MS*s + self.b*g - 2*m),
             self.r_M * m * self.b],
            [self.r_G * g * self.c,
             self.r_G * g * self.e,
             self.r_G * (self.d + self.c*s + self.e*m - 2*g)]
        ])
        return J


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def plot_pairwise_systems():
    """Generate Figure 1: All pairwise interactions"""

    fig = plt.figure(figsize=(15, 5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.3)

    # ========================================================================
    # Panel A: S-M System
    # ========================================================================
    ax1 = fig.add_subplot(gs[0, 0])

    sm_system = SMSystem(sigma_SM=0.5, sigma_MS=1.5)
    eq_sm = sm_system.equilibria()

    # Phase portrait
    s_range = np.linspace(0, 2.5, 300)
    m_range = np.linspace(0, 2.5, 300)

    # S-nullcline: s = 1 + σ_SM·m
    s_nullcline = 1 + sm_system.sigma_SM * m_range
    # M-nullcline: m = σ_MS·s - 1
    m_nullcline = sm_system.sigma_MS * s_range - 1

    ax1.plot(s_nullcline, m_range, 'b-', linewidth=2, label='S-nullcline', alpha=0.7)
    ax1.plot(s_range, m_nullcline, 'orange', linewidth=2, label='M-nullcline', alpha=0.7)

    # Vector field
    s_grid = np.linspace(0, 2.5, 15)
    m_grid = np.linspace(0, 2.5, 15)
    S, M = np.meshgrid(s_grid, m_grid)
    dS = np.zeros_like(S)
    dM = np.zeros_like(M)

    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            deriv = sm_system.equations(0, [S[i,j], M[i,j]])
            dS[i,j] = deriv[0]
            dM[i,j] = deriv[1]

    speed = np.sqrt(dS**2 + dM**2)
    speed[speed == 0] = 1
    ax1.quiver(S, M, dS/speed, dM/speed, speed,
               cmap='Greys', alpha=0.3, width=0.003)

    # Equilibria
    for name, (s, m) in eq_sm.items():
        if 0 <= s <= 2.5 and 0 <= m <= 2.5:
            stable = sm_system.is_stable(s, m)
            if stable:
                ax1.plot(s, m, 'o', color=COLORS['stable'],
                        markersize=10, markeredgewidth=2,
                        markeredgecolor='darkgreen', zorder=10)
            else:
                ax1.plot(s, m, 'o', color='white',
                        markersize=10, markeredgewidth=2,
                        markeredgecolor=COLORS['unstable'], zorder=10)

    # Trajectories
    for ic in [[0.5, 0.2], [1.5, 1.5], [0.3, 1.0]]:
        sol = solve_ivp(sm_system.equations, [0, 50], ic,
                       dense_output=True, method='RK45')
        t_eval = np.linspace(0, 50, 500)
        traj = sol.sol(t_eval)
        ax1.plot(traj[0], traj[1], 'k-', alpha=0.3, linewidth=0.8)

    ax1.set_xlabel('$s$ (S density)', fontweight='bold')
    ax1.set_ylabel('$m$ (M density)', fontweight='bold')
    ax1.set_title('(A) S-M System', fontweight='bold', loc='left')
    ax1.set_xlim(0, 2.5)
    ax1.set_ylim(0, 2.5)
    ax1.grid(True, alpha=0.2)
    ax1.legend(loc='upper right', framealpha=0.9)

    # ========================================================================
    # Panel B: S-G System
    # ========================================================================
    ax2 = fig.add_subplot(gs[0, 1])

    sg_system = SGSystem(a=0.2, c=0.3, d=0.0)
    eq_sg = sg_system.equilibria()

    # Nullclines
    s_nullcline_sg = 1 + sg_system.a * m_range
    g_nullcline_sg = (sg_system.d + sg_system.c * s_range)

    ax2.plot(s_nullcline_sg, m_range, 'b-', linewidth=2, label='S-nullcline', alpha=0.7)
    ax2.plot(s_range, g_nullcline_sg, color=COLORS['G'], linewidth=2,
            label='G-nullcline', alpha=0.7)

    # Vector field
    G_grid = np.linspace(0, 2.5, 15)
    S_grid, G_grid_mesh = np.meshgrid(s_grid, G_grid)
    dS_sg = np.zeros_like(S_grid)
    dG_sg = np.zeros_like(G_grid_mesh)

    for i in range(S_grid.shape[0]):
        for j in range(S_grid.shape[1]):
            deriv = sg_system.equations(0, [S_grid[i,j], G_grid_mesh[i,j]])
            dS_sg[i,j] = deriv[0]
            dG_sg[i,j] = deriv[1]

    speed_sg = np.sqrt(dS_sg**2 + dG_sg**2)
    speed_sg[speed_sg == 0] = 1
    ax2.quiver(S_grid, G_grid_mesh, dS_sg/speed_sg, dG_sg/speed_sg, speed_sg,
               cmap='Greys', alpha=0.3, width=0.003)

    # Equilibria
    for name, point in eq_sg.items():
        if len(point) == 2:
            s, g = point
            if 0 <= s <= 2.5 and 0 <= g <= 2.5:
                stable = sg_system.is_stable(s, g)
                if stable:
                    ax2.plot(s, g, 'o', color=COLORS['stable'],
                            markersize=10, markeredgewidth=2,
                            markeredgecolor='darkgreen', zorder=10)
                else:
                    ax2.plot(s, g, 'o', color='white',
                            markersize=10, markeredgewidth=2,
                            markeredgecolor=COLORS['unstable'], zorder=10)

    # Trajectories
    for ic in [[0.5, 0.5], [1.5, 1.0], [0.8, 1.5]]:
        sol = solve_ivp(sg_system.equations, [0, 50], ic,
                       dense_output=True, method='RK45')
        t_eval = np.linspace(0, 50, 500)
        traj = sol.sol(t_eval)
        ax2.plot(traj[0], traj[1], 'k-', alpha=0.3, linewidth=0.8)

    ax2.set_xlabel('$s$ (S density)', fontweight='bold')
    ax2.set_ylabel('$g$ (G density)', fontweight='bold')
    ax2.set_title('(B) S-G System', fontweight='bold', loc='left')
    ax2.set_xlim(0, 2.5)
    ax2.set_ylim(0, 2.5)
    ax2.grid(True, alpha=0.2)
    ax2.legend(loc='upper right', framealpha=0.9)

    # ========================================================================
    # Panel C: M-G System
    # ========================================================================
    ax3 = fig.add_subplot(gs[0, 2])

    mg_system = MGSystem(b=0.3, e=0.2, d=0.0)
    eq_mg = mg_system.equilibria()

    # Nullclines
    m_nullcline_mg = mg_system.b * m_range - 1
    g_nullcline_mg = mg_system.d + mg_system.e * s_range

    ax3.plot(m_nullcline_mg, m_range, color='orange', linewidth=2,
            label='M-nullcline', alpha=0.7)
    ax3.plot(s_range, g_nullcline_mg, color=COLORS['G'], linewidth=2,
            label='G-nullcline', alpha=0.7)

    # Vector field
    M_grid = np.linspace(0, 2.5, 15)
    M_grid_mesh, G_grid_mesh = np.meshgrid(M_grid, G_grid)
    dM_mg = np.zeros_like(M_grid_mesh)
    dG_mg = np.zeros_like(G_grid_mesh)

    for i in range(M_grid_mesh.shape[0]):
        for j in range(M_grid_mesh.shape[1]):
            deriv = mg_system.equations(0, [M_grid_mesh[i,j], G_grid_mesh[i,j]])
            dM_mg[i,j] = deriv[0]
            dG_mg[i,j] = deriv[1]

    speed_mg = np.sqrt(dM_mg**2 + dG_mg**2)
    speed_mg[speed_mg == 0] = 1
    ax3.quiver(M_grid_mesh, G_grid_mesh, dM_mg/speed_mg, dG_mg/speed_mg, speed_mg,
               cmap='Greys', alpha=0.3, width=0.003)

    # Equilibria
    for name, point in eq_mg.items():
        if len(point) == 2:
            m, g = point
            if 0 <= m <= 2.5 and 0 <= g <= 2.5:
                stable = mg_system.is_stable(m, g)
                if stable:
                    ax3.plot(m, g, 'o', color=COLORS['stable'],
                            markersize=10, markeredgewidth=2,
                            markeredgecolor='darkgreen', zorder=10)
                else:
                    ax3.plot(m, g, 'o', color='white',
                            markersize=10, markeredgewidth=2,
                            markeredgecolor=COLORS['unstable'], zorder=10)

    # Trajectories
    for ic in [[0.5, 0.5], [1.0, 1.5], [0.3, 1.0]]:
        sol = solve_ivp(mg_system.equations, [0, 50], ic,
                       dense_output=True, method='RK45')
        t_eval = np.linspace(0, 50, 500)
        traj = sol.sol(t_eval)
        ax3.plot(traj[0], traj[1], 'k-', alpha=0.3, linewidth=0.8)

    ax3.set_xlabel('$m$ (M density)', fontweight='bold')
    ax3.set_ylabel('$g$ (G density)', fontweight='bold')
    ax3.set_title('(C) M-G System', fontweight='bold', loc='left')
    ax3.set_xlim(0, 2.5)
    ax3.set_ylim(0, 2.5)
    ax3.grid(True, alpha=0.2)
    ax3.legend(loc='upper right', framealpha=0.9)

    plt.savefig('figures/Figure1_pairwise_systems.png',
                dpi=300, bbox_inches='tight')
    print("✓ Figure 1 saved: figures/Figure1_pairwise_systems.png")

    return fig


def plot_G_invasion_analysis():
    """Generate Figure 2: G invasion of S-M consortia"""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Fixed parameters
    sigma_SM = 0.5
    sigma_MS = 1.5

    # S-M equilibrium
    denom = 1 - sigma_MS * sigma_SM
    s_SM = (1 - sigma_SM) / denom
    m_SM = (sigma_MS - 1) / denom

    # ========================================================================
    # Panel A: G invasion rate as function of ω
    # ========================================================================
    ax1 = axes[0, 0]

    omega_range = np.linspace(0, 1, 200)
    invasion_rates = []

    # Example parameters
    sigma_GS, sigma_GM = 0.4, 0.4
    alpha_GS, alpha_GM = 0.3, 0.3

    for omega in omega_range:
        c = (1 - omega) * sigma_GS - omega * alpha_GS
        e = omega * sigma_GM - (1 - omega) * alpha_GM
        d = 2 * omega - 1

        invasion_rate = d + c * s_SM + e * m_SM
        invasion_rates.append(invasion_rate)

    invasion_rates = np.array(invasion_rates)

    # Plot invasion rate
    ax1.plot(omega_range, invasion_rates, 'g-', linewidth=2.5)
    ax1.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.5)
    ax1.fill_between(omega_range, 0, invasion_rates,
                     where=invasion_rates>0, alpha=0.3, color='lightgreen',
                     label='G can invade')
    ax1.fill_between(omega_range, invasion_rates, 0,
                     where=invasion_rates<0, alpha=0.3, color='lightcoral',
                     label='G excluded')

    # Find critical omega
    zero_crossings = np.where(np.diff(np.sign(invasion_rates)))[0]
    if len(zero_crossings) > 0:
        omega_crit = omega_range[zero_crossings[0]]
        ax1.axvline(omega_crit, color='red', linestyle='--', linewidth=2,
                   label=f'$\\omega_{{crit}}$ = {omega_crit:.3f}')

    ax1.set_xlabel('$\\omega$ (pathway parameter)', fontweight='bold')
    ax1.set_ylabel('G invasion rate', fontweight='bold')
    ax1.set_title('(A) G Invasion Rate vs $\\omega$', fontweight='bold', loc='left')
    ax1.legend(framealpha=0.9)
    ax1.grid(True, alpha=0.3)

    # ========================================================================
    # Panel B: Parameter components
    # ========================================================================
    ax2 = axes[0, 1]

    d_vals = 2 * omega_range - 1
    c_vals = (1 - omega_range) * sigma_GS - omega_range * alpha_GS
    e_vals = omega_range * sigma_GM - (1 - omega_range) * alpha_GM

    ax2.plot(omega_range, d_vals, 'b-', linewidth=2, label='$d = 2\\omega - 1$')
    ax2.plot(omega_range, c_vals * s_SM, 'r-', linewidth=2, label='$c \\cdot s^*_{SM}$')
    ax2.plot(omega_range, e_vals * m_SM, 'g-', linewidth=2, label='$e \\cdot m^*_{SM}$')
    ax2.plot(omega_range, invasion_rates, 'k--', linewidth=2.5,
            label='Total = $d + cs^* + em^*$')

    ax2.axhline(0, color='k', linestyle=':', linewidth=1, alpha=0.5)
    ax2.set_xlabel('$\\omega$', fontweight='bold')
    ax2.set_ylabel('Parameter values', fontweight='bold')
    ax2.set_title('(B) Components of Invasion Rate', fontweight='bold', loc='left')
    ax2.legend(framealpha=0.9)
    ax2.grid(True, alpha=0.3)

    # ========================================================================
    # Panel C: 2D parameter space (ω, σ_MS)
    # ========================================================================
    ax3 = axes[1, 0]

    omega_2d = np.linspace(0, 1, 100)
    sigma_MS_2d = np.linspace(1.0, 2.5, 100)

    Omega, Sigma_MS = np.meshgrid(omega_2d, sigma_MS_2d)
    InvasionRate = np.zeros_like(Omega)

    for i in range(len(sigma_MS_2d)):
        for j in range(len(omega_2d)):
            omega_val = omega_2d[j]
            sigma_MS_val = sigma_MS_2d[i]

            # Compute S-M equilibrium
            denom_temp = 1 - sigma_MS_val * sigma_SM
            s_SM_temp = (1 - sigma_SM) / denom_temp
            m_SM_temp = (sigma_MS_val - 1) / denom_temp

            # Compute invasion rate
            c_temp = (1 - omega_val) * sigma_GS - omega_val * alpha_GS
            e_temp = omega_val * sigma_GM - (1 - omega_val) * alpha_GM
            d_temp = 2 * omega_val - 1

            InvasionRate[i, j] = d_temp + c_temp * s_SM_temp + e_temp * m_SM_temp

    # Contour plot
    levels = [-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0]
    cs = ax3.contourf(Omega, Sigma_MS, InvasionRate, levels=levels,
                     cmap='RdYlGn', alpha=0.7)
    ct = ax3.contour(Omega, Sigma_MS, InvasionRate, levels=[0],
                    colors='red', linewidths=2.5)

    plt.colorbar(cs, ax=ax3, label='Invasion rate')
    ax3.clabel(ct, inline=True, fontsize=9, fmt='%0.2f')

    ax3.set_xlabel('$\\omega$', fontweight='bold')
    ax3.set_ylabel('$\\sigma_{MS}$', fontweight='bold')
    ax3.set_title('(C) Parameter Space: G Invasion', fontweight='bold', loc='left')

    # ========================================================================
    # Panel D: Time series showing invasion
    # ========================================================================
    ax4 = axes[1, 1]

    # Simulate with ω where G can invade
    omega_invade = 0.6
    c_inv = (1 - omega_invade) * sigma_GS - omega_invade * alpha_GS
    e_inv = omega_invade * sigma_GM - (1 - omega_invade) * alpha_GM
    d_inv = 2 * omega_invade - 1

    system = ThreeSpeciesSystem(sigma_SM=sigma_SM, sigma_MS=sigma_MS,
                                a=0.2, b=0.3, c=c_inv, d=d_inv, e=e_inv)

    # Start near S-M equilibrium, small G perturbation
    y0 = [s_SM, m_SM, 0.01]
    sol = solve_ivp(system.equations, [0, 100], y0,
                   dense_output=True, method='RK45')
    t_eval = np.linspace(0, 100, 1000)
    traj = sol.sol(t_eval)

    ax4.plot(t_eval, traj[0], color=COLORS['S'], linewidth=2, label='S')
    ax4.plot(t_eval, traj[1], color=COLORS['M'], linewidth=2, label='M')
    ax4.plot(t_eval, traj[2], color=COLORS['G'], linewidth=2, label='G')

    ax4.set_xlabel('Time', fontweight='bold')
    ax4.set_ylabel('Population density', fontweight='bold')
    ax4.set_title(f'(D) G Invasion Dynamics ($\\omega$ = {omega_invade})',
                 fontweight='bold', loc='left')
    ax4.legend(framealpha=0.9)
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('figures/Figure2_G_invasion.png', dpi=300, bbox_inches='tight')
    print("✓ Figure 2 saved: figures/Figure2_G_invasion.png")

    return fig


def plot_bifurcation_diagram():
    """Generate Figure 3: Complete bifurcation diagram"""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Parameters
    sigma_SM = 0.5
    sigma_MS = 1.5
    sigma_GS, sigma_GM = 0.4, 0.4
    alpha_GS, alpha_GM = 0.3, 0.3

    omega_range = np.linspace(0.0, 1.0, 150)

    # Store results
    s_vals, m_vals, g_vals = [], [], []
    stable_vals = []

    for omega in omega_range:
        # Compute parameters
        c = (1 - omega) * sigma_GS - omega * alpha_GS
        e = omega * sigma_GM - (1 - omega) * alpha_GM
        d = 2 * omega - 1
        a = 0.2
        b = 0.3

        system = ThreeSpeciesSystem(sigma_SM=sigma_SM, sigma_MS=sigma_MS,
                                    a=a, b=b, c=c, d=d, e=e)

        # Try to find equilibrium
        eq = system.find_interior_equilibrium()

        if eq is not None and np.all(eq > 0.01):
            s_vals.append(eq[0])
            m_vals.append(eq[1])
            g_vals.append(eq[2])

            # Check stability
            J = system.jacobian(*eq)
            eigs = np.linalg.eigvals(J)
            stable = np.all(np.real(eigs) < 0)
            stable_vals.append(stable)
        else:
            s_vals.append(np.nan)
            m_vals.append(np.nan)
            g_vals.append(np.nan)
            stable_vals.append(False)

    s_vals = np.array(s_vals)
    m_vals = np.array(m_vals)
    g_vals = np.array(g_vals)
    stable_vals = np.array(stable_vals)

    # ========================================================================
    # Panel A: S population
    # ========================================================================
    ax1 = axes[0, 0]

    # Plot stable/unstable branches
    ax1.plot(omega_range[stable_vals], s_vals[stable_vals],
            'o', color=COLORS['S'], markersize=3, label='Stable')
    ax1.plot(omega_range[~stable_vals], s_vals[~stable_vals],
            'x', color=COLORS['S'], markersize=3, alpha=0.3, label='Unstable')

    # S-M equilibrium (no G)
    denom = 1 - sigma_MS * sigma_SM
    s_SM = (1 - sigma_SM) / denom
    ax1.axhline(s_SM, color='gray', linestyle=':', linewidth=1.5,
               label='S-M equilibrium')

    ax1.set_xlabel('$\\omega$', fontweight='bold')
    ax1.set_ylabel('$s^*$ (S density)', fontweight='bold')
    ax1.set_title('(A) S Population at Equilibrium', fontweight='bold', loc='left')
    ax1.legend(framealpha=0.9)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 3)

    # ========================================================================
    # Panel B: M population
    # ========================================================================
    ax2 = axes[0, 1]

    ax2.plot(omega_range[stable_vals], m_vals[stable_vals],
            'o', color=COLORS['M'], markersize=3, label='Stable')
    ax2.plot(omega_range[~stable_vals], m_vals[~stable_vals],
            'x', color=COLORS['M'], markersize=3, alpha=0.3, label='Unstable')

    # S-M equilibrium
    m_SM = (sigma_MS - 1) / denom
    ax2.axhline(m_SM, color='gray', linestyle=':', linewidth=1.5,
               label='S-M equilibrium')
    ax2.axhline(0, color='k', linestyle='-', linewidth=1, alpha=0.3)

    ax2.set_xlabel('$\\omega$', fontweight='bold')
    ax2.set_ylabel('$m^*$ (M density)', fontweight='bold')
    ax2.set_title('(B) M Population at Equilibrium', fontweight='bold', loc='left')
    ax2.legend(framealpha=0.9)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 3)

    # ========================================================================
    # Panel C: G population
    # ========================================================================
    ax3 = axes[1, 0]

    ax3.plot(omega_range[stable_vals], g_vals[stable_vals],
            'o', color=COLORS['G'], markersize=3, label='Stable')
    ax3.plot(omega_range[~stable_vals], g_vals[~stable_vals],
            'x', color=COLORS['G'], markersize=3, alpha=0.3, label='Unstable')

    ax3.axhline(0, color='gray', linestyle=':', linewidth=1.5,
               label='S-M equilibrium (G=0)')

    ax3.set_xlabel('$\\omega$', fontweight='bold')
    ax3.set_ylabel('$g^*$ (G density)', fontweight='bold')
    ax3.set_title('(C) G Population at Equilibrium', fontweight='bold', loc='left')
    ax3.legend(framealpha=0.9)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(-0.2, 3)

    # ========================================================================
    # Panel D: All species together
    # ========================================================================
    ax4 = axes[1, 1]

    ax4.plot(omega_range, s_vals, '-', color=COLORS['S'], linewidth=2,
            alpha=0.7, label='S')
    ax4.plot(omega_range, m_vals, '-', color=COLORS['M'], linewidth=2,
            alpha=0.7, label='M')
    ax4.plot(omega_range, g_vals, '-', color=COLORS['G'], linewidth=2,
            alpha=0.7, label='G')

    # Mark coexistence region
    coexist = np.where(~np.isnan(g_vals) & (g_vals > 0.1))[0]
    if len(coexist) > 0:
        omega_min = omega_range[coexist[0]]
        omega_max = omega_range[coexist[-1]]
        ax4.axvspan(omega_min, omega_max, alpha=0.2, color='gold',
                   label='Three-species coexistence')
        ax4.axvline(omega_min, color='red', linestyle='--', linewidth=1.5,
                   alpha=0.7)
        ax4.axvline(omega_max, color='red', linestyle='--', linewidth=1.5,
                   alpha=0.7)

    ax4.set_xlabel('$\\omega$', fontweight='bold')
    ax4.set_ylabel('Population density', fontweight='bold')
    ax4.set_title('(D) Complete Bifurcation Diagram', fontweight='bold', loc='left')
    ax4.legend(framealpha=0.9)
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(0, 3)

    plt.tight_layout()
    plt.savefig('figures/Figure3_bifurcation_diagram.png',
                dpi=300, bbox_inches='tight')
    print("✓ Figure 3 saved: figures/Figure3_bifurcation_diagram.png")

    return fig


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Generate all figures for the paper"""

    print("\n" + "="*70)
    print("GENERATING ALL FIGURES FOR PUBLICATION")
    print("="*70 + "\n")

    # Section 1: Pairwise interactions
    print("Section 1: Pairwise interactions...")
    plot_pairwise_systems()

    # Section 2.2: G invasion
    print("\nSection 2.2: G invasion analysis...")
    plot_G_invasion_analysis()

    # Section 2.3: Bifurcation analysis
    print("\nSection 2.3: Bifurcation diagram...")
    plot_bifurcation_diagram()

    print("\n" + "="*70)
    print("✅ ALL FIGURES GENERATED SUCCESSFULLY")
    print("="*70)
    print("\nGenerated files:")
    print("  - figures/Figure1_pairwise_systems.png")
    print("  - figures/Figure2_G_invasion.png")
    print("  - figures/Figure3_bifurcation_diagram.png")
    print("\n" + "="*70)


if __name__ == "__main__":
    main()

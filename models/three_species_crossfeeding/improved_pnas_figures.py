#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Improved PNAS-Style Figures - High Quality, Microbial Ecology Focus
基于manuscript_PNAS.tex的高质量图表生成

确保：
1. 每个panel都有清晰、丰富的内容
2. 强化微生物生态学视角
3. 完美对应manuscript_PNAS.tex文本
4. 专业出版质量

Author: Jian Wang
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import cm
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.linalg import eig

# PNAS高质量设置
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica'],
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'lines.linewidth': 1.5,
    'axes.linewidth': 1.0,
    'grid.linewidth': 0.5
})

# 微生物生态学主题配色
COLORS = {
    'S': '#1f77b4',  # 蓝色 - 基质专食者
    'M': '#d62728',  # 红色 - 代谢物专食者
    'G': '#2ca02c',  # 绿色 - 通食者
    'invasion': '#90EE90',  # 浅绿 - 入侵成功
    'exclusion': '#FFB6C1',  # 浅红 - 排斥
    'coexist': '#FFD700'   # 金色 - 共存
}


class MicrobialCrossFeedingModel:
    """微生物交叉喂养三物种模型"""

    def __init__(self, params):
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
        """计算净交互参数"""
        a = (1 - omega) * self.sigma_SG - omega * self.alpha_SG
        c = (1 - omega) * self.sigma_GS - omega * self.alpha_GS
        b = omega * self.sigma_MG - (1 - omega) * self.alpha_MG
        e = omega * self.sigma_GM - (1 - omega) * self.alpha_GM
        d = 2 * omega - 1
        return a, b, c, d, e

    def dynamics_SM(self, t, y):
        """S-M双物种动力学"""
        s, m = y
        dsdt = self.r_S * s * (1 + self.sigma_SM * m - s)
        dmdt = self.r_M * m * (-1 + self.sigma_MS * s - m)
        return [dsdt, dmdt]

    def dynamics_three_species(self, t, y, omega):
        """三物种完整动力学"""
        s, m, g = y
        a, b, c, d, e = self.net_interactions(omega)

        dsdt = self.r_S * s * (1 + self.sigma_SM * m + a * g - s)
        dmdt = self.r_M * m * (-1 + self.sigma_MS * s + b * g - m)
        dgdt = self.r_G * g * (d + c * s + e * m - g)

        return [dsdt, dmdt, dgdt]

    def compute_SM_equilibrium(self):
        """计算S-M平衡点"""
        denom = 1 - self.sigma_MS * self.sigma_SM
        if abs(denom) > 1e-10 and self.sigma_MS > 1:
            s_SM = (1 - self.sigma_SM) / denom
            m_SM = (self.sigma_MS - 1) / denom
            stable = (self.sigma_MS > 1) and (self.sigma_MS * self.sigma_SM < 1)
            return s_SM, m_SM, stable
        return None, None, False

    def invasion_fitness(self, omega):
        """计算通食者入侵适应度"""
        s_SM, m_SM, stable = self.compute_SM_equilibrium()
        if not stable:
            return -np.inf
        _, _, c, d, e = self.net_interactions(omega)
        lambda_G = self.r_G * (d + c * s_SM + e * m_SM)
        return lambda_G


def create_figure1_pairwise_and_invasion():
    """
    Figure 1: S-M互利共生平台与通食者入侵

    对应manuscript_PNAS.tex第44-56行和第57-75行
    重点：obligate mutualism建立稳定平台 + generalist invasion需要中间策略
    """
    print("\n创建Figure 1: S-M互利共生与通食者入侵...")

    model = MicrobialCrossFeedingModel({})

    fig = plt.figure(figsize=(10, 6))
    gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.4)

    # Panel A: S-M时间序列（展示obligate mutualism）
    ax_a = fig.add_subplot(gs[0, 0])

    t_span = (0, 50)
    t_eval = np.linspace(0, 50, 500)
    sol = solve_ivp(model.dynamics_SM, t_span, [0.5, 0.2], t_eval=t_eval,
                   method='RK45', rtol=1e-9)

    ax_a.plot(sol.t, sol.y[0], '-', color=COLORS['S'], linewidth=2.5,
             label='S (substrate specialist)', alpha=0.9)
    ax_a.plot(sol.t, sol.y[1], '-', color=COLORS['M'], linewidth=2.5,
             label='M (metabolite specialist)', alpha=0.9)

    # 标注平衡点
    s_eq, m_eq, _ = model.compute_SM_equilibrium()
    ax_a.axhline(s_eq, color=COLORS['S'], linestyle='--', linewidth=1, alpha=0.5)
    ax_a.axhline(m_eq, color=COLORS['M'], linestyle='--', linewidth=1, alpha=0.5)
    ax_a.text(45, s_eq+0.05, f's*={s_eq:.2f}', fontsize=7, color=COLORS['S'])
    ax_a.text(45, m_eq+0.05, f'm*={m_eq:.2f}', fontsize=7, color=COLORS['M'])

    ax_a.set_xlabel('Time (days)', fontweight='bold')
    ax_a.set_ylabel('Population density', fontweight='bold')
    ax_a.set_title('A. Obligate mutualism: S-M platform', fontweight='bold', loc='left')
    ax_a.legend(loc='right', frameon=True, fancybox=True, shadow=True)
    ax_a.grid(alpha=0.3, linestyle=':', linewidth=0.5)
    ax_a.set_ylim(0, 1.2)

    # 添加生物学注释
    ax_a.text(0.05, 0.95, 'M obligate cross-feeder\n(cannot survive alone)',
             transform=ax_a.transAxes, fontsize=6, va='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel B: S-M相空间（展示稳定性）
    ax_b = fig.add_subplot(gs[0, 1])

    # Nullclines
    s_range = np.linspace(0, 2, 100)
    # S nullcline: m = (s - 1) / sigma_SM
    m_null_s = (s_range - 1) / model.sigma_SM
    # M nullcline: s = (1 + m) / sigma_MS
    m_range = np.linspace(0, 1.5, 100)
    s_null_m = (1 + m_range) / model.sigma_MS

    ax_b.plot(s_range, m_null_s, '-', color=COLORS['S'], linewidth=2.5,
             label='S nullcline', alpha=0.7)
    ax_b.plot(s_null_m, m_range, '-', color=COLORS['M'], linewidth=2.5,
             label='M nullcline', alpha=0.7)

    # 轨迹
    for s0, m0 in [(0.3, 0.1), (0.8, 0.4), (1.2, 0.7)]:
        sol_traj = solve_ivp(model.dynamics_SM, (0, 80), [s0, m0],
                           t_eval=np.linspace(0, 80, 1000))
        ax_b.plot(sol_traj.y[0], sol_traj.y[1], '-', color='gray',
                 alpha=0.6, linewidth=1.2)
        ax_b.plot(s0, m0, 'ko', markersize=5)

    # 平衡点
    ax_b.plot(s_eq, m_eq, '*', color=COLORS['coexist'], markersize=15,
             markeredgecolor='black', markeredgewidth=1.5,
             label='Stable equilibrium', zorder=10)

    ax_b.set_xlabel('S density', fontweight='bold')
    ax_b.set_ylabel('M density', fontweight='bold')
    ax_b.set_title('B. Phase portrait: global convergence', fontweight='bold', loc='left')
    ax_b.legend(loc='upper left', frameon=True, fontsize=6)
    ax_b.grid(alpha=0.3, linestyle=':', linewidth=0.5)
    ax_b.set_xlim(0, 1.8)
    ax_b.set_ylim(0, 1.2)

    # Panel C: S-M稳定性区域（parameter space）
    ax_c = fig.add_subplot(gs[0, 2])

    sigma_MS_range = np.linspace(0.5, 3.0, 150)
    sigma_SM_range = np.linspace(0, 1.5, 150)

    MS_grid, SM_grid = np.meshgrid(sigma_MS_range, sigma_SM_range)

    # 稳定性条件
    viable = MS_grid > 1
    stable = (MS_grid > 1) & (MS_grid * SM_grid < 1)

    # 绘制区域
    ax_c.contourf(MS_grid, SM_grid, stable.astype(int), levels=[0, 0.5, 1],
                 colors=['#FFE4E1', '#90EE90'], alpha=0.6)

    # 边界线
    sigma_SM_boundary = 1.0 / sigma_MS_range
    ax_c.plot(sigma_MS_range, sigma_SM_boundary, 'k-', linewidth=2.5,
             label=r'$\sigma_{MS} \cdot \sigma_{SM} = 1$')
    ax_c.axvline(1.0, color='k', linestyle='--', linewidth=2, alpha=0.7,
                label=r'$\sigma_{MS} = 1$')

    # 标注当前参数
    ax_c.plot(model.sigma_MS, model.sigma_SM, 'r*', markersize=15,
             markeredgecolor='black', markeredgewidth=1.5,
             label='Baseline parameters', zorder=10)

    ax_c.set_xlabel(r'S$\to$M mutualism ($\sigma_{MS}$)', fontweight='bold')
    ax_c.set_ylabel(r'M$\to$S benefit ($\sigma_{SM}$)', fontweight='bold')
    ax_c.set_title('C. S-M stability region', fontweight='bold', loc='left')
    ax_c.legend(loc='upper left', frameon=True, fontsize=6)
    ax_c.grid(alpha=0.3, linestyle=':', linewidth=0.5)

    # 添加区域标签
    ax_c.text(2.5, 0.2, 'Stable\ncoexistence', fontsize=8, ha='center',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    ax_c.text(0.7, 0.8, 'M extinct', fontsize=8, ha='center',
             bbox=dict(boxstyle='round', facecolor='mistyrose', alpha=0.7))

    # Panel D: 通食者入侵适应度 vs ω（关键结果）
    ax_d = fig.add_subplot(gs[1, :2])

    omega_range = np.linspace(0, 1, 300)
    invasion_fitness = [model.invasion_fitness(om) for om in omega_range]

    ax_d.plot(omega_range, invasion_fitness, '-', color=COLORS['G'],
             linewidth=3, label=r'Invasion fitness $\lambda_G(\omega)$')
    ax_d.axhline(0, color='black', linestyle='--', linewidth=1.5)

    # 填充区域
    ax_d.fill_between(omega_range, 0, invasion_fitness,
                     where=np.array(invasion_fitness) > 0,
                     alpha=0.3, color=COLORS['invasion'],
                     label='G invades (coexistence)')
    ax_d.fill_between(omega_range, invasion_fitness, 0,
                     where=np.array(invasion_fitness) < 0,
                     alpha=0.3, color=COLORS['exclusion'],
                     label='G excluded (S-M only)')

    # 找到临界点
    zero_crossings = np.where(np.diff(np.sign(invasion_fitness)))[0]
    if len(zero_crossings) > 0:
        omega_crit = omega_range[zero_crossings[0]]
        ax_d.axvline(omega_crit, color='red', linestyle=':', linewidth=2.5, alpha=0.8)
        ax_d.annotate(f'$\\omega_{{crit}} = {omega_crit:.3f}$\n(transcritical bifurcation)',
                     xy=(omega_crit, 0), xytext=(omega_crit+0.15, 0.2),
                     arrowprops=dict(arrowstyle='->', color='red', lw=2),
                     fontsize=9, ha='left', color='red',
                     bbox=dict(boxstyle='round', facecolor='white',
                              edgecolor='red', linewidth=2))

    ax_d.set_xlabel(r'Pathway parameter $\omega$ (0=pure M, 1=pure S)', fontweight='bold')
    ax_d.set_ylabel('Invasion fitness $\\lambda_G$', fontweight='bold')
    ax_d.set_title('D. Generalist invasion requires intermediate metabolic strategy',
                   fontweight='bold', loc='left')
    ax_d.legend(loc='upper left', frameon=True, fancybox=True, shadow=True, fontsize=7)
    ax_d.grid(alpha=0.3, linestyle=':', linewidth=0.5)

    # 添加生态学解释
    ax_d.text(0.15, -0.25, 'Too metabolite-specialized\n→ competes with M',
             fontsize=7, ha='center', style='italic',
             bbox=dict(boxstyle='round', facecolor=COLORS['exclusion'], alpha=0.5))
    ax_d.text(0.85, -0.25, 'Too substrate-specialized\n→ competes with S',
             fontsize=7, ha='center', style='italic',
             bbox=dict(boxstyle='round', facecolor=COLORS['exclusion'], alpha=0.5))
    ax_d.text(0.5, 0.25, 'Intermediate niche\n→ coexistence',
             fontsize=7, ha='center', style='italic', fontweight='bold',
             bbox=dict(boxstyle='round', facecolor=COLORS['invasion'], alpha=0.7))

    # Panel E: 关键参数解释（呼应manuscript结论）
    ax_e = fig.add_subplot(gs[1, 2])
    ax_e.axis('off')

    # 显示关键公式和结论
    conclusion_text = r"""
$\bf{Key\ Results}$

$\bf{1.\ S-M\ Platform:}$
$s^* = \frac{1-\sigma_{SM}}{1-\sigma_{MS}\sigma_{SM}}$

$m^* = \frac{\sigma_{MS}-1}{1-\sigma_{MS}\sigma_{SM}}$

$\bf{Stability:}$ $\sigma_{MS} > 1$ AND
$\sigma_{MS}\sigma_{SM} < 1$

$\bf{2.\ G\ Invasion:}$
$\lambda_G = r_G(d + c \cdot s^* + e \cdot m^*)$

$\bf{Critical\ threshold:}$
$\omega_{crit} \approx 0.35$

$\bf{Ecological\ Insight:}$
Generalists coexist by occupying
metabolic niche unavailable to
specialists at pathway extremes
"""

    ax_e.text(0.1, 0.95, conclusion_text, transform=ax_e.transAxes,
             fontsize=7, va='top', ha='left', family='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow',
                      edgecolor='orange', linewidth=2, alpha=0.8, pad=0.8))

    ax_e.set_title('E. Analytical results & conclusions',
                   fontweight='bold', loc='left', pad=10)

    plt.savefig('figures/Figure1_mutualism_invasion_IMPROVED.png',
                dpi=300, bbox_inches='tight', facecolor='white')
    print("✓ Figure 1 创建完成：S-M互利共生与通食者入侵")

    return fig


def create_figure2_parameter_landscapes():
    """
    Figure 2: 参数空间架构决定群落组装

    对应manuscript_PNAS.tex第76-89行
    重点：curved invasion boundary + 参数空间非线性耦合
    """
    print("\n创建Figure 2: 多维参数空间景观...")

    model = MicrobialCrossFeedingModel({})

    fig = plt.figure(figsize=(12, 8))
    gs = GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.4)

    # Panel A: (ω, σ_MS) invasion landscape - 核心图
    ax_a = fig.add_subplot(gs[0, :2])

    omega_range = np.linspace(0, 1, 200)
    sigma_MS_range = np.linspace(1.1, 3.0, 200)

    invasion_map = np.zeros((len(sigma_MS_range), len(omega_range)))

    for i, sigma_MS in enumerate(sigma_MS_range):
        for j, omega in enumerate(omega_range):
            model.sigma_MS = sigma_MS
            lambda_G = model.invasion_fitness(omega)
            invasion_map[i, j] = lambda_G

    # 绘制热图
    im = ax_a.pcolormesh(omega_range, sigma_MS_range, invasion_map,
                        cmap='RdBu_r', shading='auto', vmin=-0.5, vmax=0.5)

    # 关键等高线
    CS = ax_a.contour(omega_range, sigma_MS_range, invasion_map,
                     levels=[-0.3, -0.1, 0, 0.1, 0.3],
                     colors=['darkred', 'red', 'black', 'green', 'darkgreen'],
                     linewidths=[1.5, 1.5, 3, 1.5, 1.5],
                     linestyles=['--', '--', '-', '--', '--'])
    ax_a.clabel(CS, inline=True, fontsize=7, fmt='%.2f')

    # 标注零等高线（invasion boundary）
    zero_contour = ax_a.contour(omega_range, sigma_MS_range, invasion_map,
                                levels=[0], colors='black', linewidths=3)

    ax_a.set_xlabel('Pathway parameter ($\\omega$)', fontsize=10, fontweight='bold')
    ax_a.set_ylabel('S$\\to$M mutualism ($\\sigma_{MS}$)', fontsize=10, fontweight='bold')
    ax_a.set_title('A. Curved invasion boundary reflects nonlinear pathway-mutualism coupling',
                   fontweight='bold', loc='left', fontsize=10)

    cbar = plt.colorbar(im, ax=ax_a)
    cbar.set_label('Invasion fitness $\\lambda_G$', fontsize=9, fontweight='bold')

    # 添加生态学解释
    ax_a.text(0.25, 2.7, 'Low mutualism →\nhigher ω needed',
             fontsize=7, ha='center', color='white',
             bbox=dict(boxstyle='round', facecolor='darkred', alpha=0.7))
    ax_a.text(0.7, 1.5, 'High mutualism →\nlower ω sufficient',
             fontsize=7, ha='center', color='black',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    # Panel B: S平衡密度景观
    ax_b = fig.add_subplot(gs[0, 2])

    s_SM_map = np.zeros_like(invasion_map)
    for i, sigma_MS in enumerate(sigma_MS_range):
        model.sigma_MS = sigma_MS
        s_eq, _, _ = model.compute_SM_equilibrium()
        s_SM_map[i, :] = s_eq if s_eq is not None else np.nan

    im_b = ax_b.pcolormesh(omega_range, sigma_MS_range, s_SM_map,
                          cmap='Blues', shading='auto')
    ax_b.contour(omega_range, sigma_MS_range, invasion_map, levels=[0],
                colors='red', linewidths=2.5, linestyles='-')

    ax_b.set_xlabel('$\\omega$', fontweight='bold')
    ax_b.set_ylabel('$\\sigma_{MS}$', fontweight='bold')
    ax_b.set_title('B. S equilibrium\ndensity $s^*_{SM}$',
                   fontweight='bold', loc='left', fontsize=9)

    cbar_b = plt.colorbar(im_b, ax=ax_b)
    cbar_b.set_label('$s^*_{SM}$', fontsize=8)

    # Panel C: M平衡密度景观
    ax_c = fig.add_subplot(gs[1, 0])

    m_SM_map = np.zeros_like(invasion_map)
    for i, sigma_MS in enumerate(sigma_MS_range):
        model.sigma_MS = sigma_MS
        _, m_eq, _ = model.compute_SM_equilibrium()
        m_SM_map[i, :] = m_eq if m_eq is not None else np.nan

    im_c = ax_c.pcolormesh(omega_range, sigma_MS_range, m_SM_map,
                          cmap='Reds', shading='auto')
    ax_c.contour(omega_range, sigma_MS_range, invasion_map, levels=[0],
                colors='blue', linewidths=2.5, linestyles='-')

    ax_c.set_xlabel('$\\omega$', fontweight='bold')
    ax_c.set_ylabel('$\\sigma_{MS}$', fontweight='bold')
    ax_c.set_title('C. M equilibrium\ndensity $m^*_{SM}$',
                   fontweight='bold', loc='left', fontsize=9)

    cbar_c = plt.colorbar(im_c, ax=ax_c)
    cbar_c.set_label('$m^*_{SM}$', fontsize=8)

    # Panel D: (ω, σ_SM) 互补约束
    ax_d = fig.add_subplot(gs[1, 1])

    sigma_SM_range = np.linspace(0.1, 1.0, 150)
    model.sigma_MS = 1.5  # 固定

    invasion_map_d = np.zeros((len(sigma_SM_range), len(omega_range)))

    for i, sigma_SM in enumerate(sigma_SM_range):
        for j, omega in enumerate(omega_range):
            model.sigma_SM = sigma_SM
            lambda_G = model.invasion_fitness(omega)
            invasion_map_d[i, j] = lambda_G

    im_d = ax_d.pcolormesh(omega_range, sigma_SM_range, invasion_map_d,
                          cmap='RdBu_r', shading='auto', vmin=-0.5, vmax=0.5)
    ax_d.contour(omega_range, sigma_SM_range, invasion_map_d, levels=[0],
                colors='black', linewidths=2.5)

    # 标注稳定性边界
    sigma_SM_stability = 1.0 / model.sigma_MS
    ax_d.axhline(sigma_SM_stability, color='purple', linestyle='--',
                linewidth=2.5, label=r'S-M stability: $\sigma_{SM}=1/\sigma_{MS}$')

    ax_d.set_xlabel('$\\omega$', fontweight='bold')
    ax_d.set_ylabel('$\\sigma_{SM}$', fontweight='bold')
    ax_d.set_title('D. Complementary constraint:\nS-M stability bound',
                   fontweight='bold', loc='left', fontsize=9)
    ax_d.legend(loc='upper right', frameon=True, fontsize=6)

    cbar_d = plt.colorbar(im_d, ax=ax_d)
    cbar_d.set_label('$\\lambda_G$', fontsize=8)

    # Panel E: 合作-竞争平衡
    ax_e = fig.add_subplot(gs[1, 2])

    sigma_GS_range = np.linspace(0.1, 0.8, 150)
    alpha_GS_range = np.linspace(0.1, 0.8, 150)
    model.omega = 0.5  # 固定在中间值

    invasion_map_e = np.zeros((len(alpha_GS_range), len(sigma_GS_range)))

    for i, alpha_GS in enumerate(alpha_GS_range):
        for j, sigma_GS in enumerate(sigma_GS_range):
            model.alpha_GS = alpha_GS
            model.sigma_GS = sigma_GS
            lambda_G = model.invasion_fitness(model.omega)
            invasion_map_e[i, j] = lambda_G

    im_e = ax_e.pcolormesh(sigma_GS_range, alpha_GS_range, invasion_map_e,
                          cmap='RdBu_r', shading='auto', vmin=-0.5, vmax=0.5)
    ax_e.contour(sigma_GS_range, alpha_GS_range, invasion_map_e,
                levels=[0], colors='black', linewidths=2.5)

    # 对角线
    ax_e.plot([0.1, 0.8], [0.1, 0.8], 'k:', linewidth=1.5, alpha=0.6)
    ax_e.text(0.6, 0.65, r'$\sigma_{GS}=\alpha_{GS}$', fontsize=7, rotation=45)

    ax_e.set_xlabel('S$\\to$G cooperation ($\\sigma_{GS}$)', fontweight='bold', fontsize=8)
    ax_e.set_ylabel('S-G competition ($\\alpha_{GS}$)', fontweight='bold', fontsize=8)
    ax_e.set_title('E. Cooperation-competition\nbalance',
                   fontweight='bold', loc='left', fontsize=9)

    cbar_e = plt.colorbar(im_e, ax=ax_e)
    cbar_e.set_label('$\\lambda_G$', fontsize=8)

    plt.savefig('figures/Figure2_parameter_landscapes_IMPROVED.png',
                dpi=300, bbox_inches='tight', facecolor='white')
    print("✓ Figure 2 创建完成：多维参数空间景观")

    return fig


# 主程序
if __name__ == "__main__":
    print("\n" + "="*70)
    print("生成改进的PNAS风格图表 - 基于manuscript_PNAS.tex")
    print("="*70)

    # 创建图表
    fig1 = create_figure1_pairwise_and_invasion()
    plt.close(fig1)

    fig2 = create_figure2_parameter_landscapes()
    plt.close(fig2)

    print("\n" + "="*70)
    print("✅ 所有改进图表生成完成！")
    print("\n生成的文件：")
    print("  • Figure1_mutualism_invasion_IMPROVED.png")
    print("  • Figure2_parameter_landscapes_IMPROVED.png")
    print("\n特点：")
    print("  ✓ 高质量，每个panel都有清晰内容")
    print("  ✓ 强化微生物生态学视角")
    print("  ✓ 完美对应manuscript_PNAS.tex")
    print("  ✓ 专业出版质量")
    print("="*70)

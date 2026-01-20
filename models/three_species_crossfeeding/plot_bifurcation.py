#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
跨临界分岔图 - ω作为分岔参数
Transcritical Bifurcation Diagram - ω as Bifurcation Parameter

展示三物种模型中ω变化时平衡点的演化和稳定性变化
Shows evolution of equilibria and stability changes as ω varies

Author: Jian Wang
Date: January 2026
"""

import sys
sys.path.append('src')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from three_species_model import ThreeSpeciesModel

# 设置中文字体（如果可用）
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.figsize'] = (14, 10)
plt.rcParams['font.size'] = 11


def compute_bifurcation_diagram(omega_values, model_params, N0):
    """
    计算分岔图数据

    Parameters
    ----------
    omega_values : array
        ω值的数组
    model_params : dict
        模型参数
    N0 : array
        初始条件

    Returns
    -------
    dict
        包含分岔数据的字典
    """
    n_omega = len(omega_values)

    # 存储结果
    results = {
        'omega': omega_values,
        'N_S_sim': np.zeros(n_omega),
        'N_M_sim': np.zeros(n_omega),
        'N_G_sim': np.zeros(n_omega),
        'N_S_stable': [],
        'N_M_stable': [],
        'N_G_stable': [],
        'N_S_unstable': [],
        'N_M_unstable': [],
        'N_G_unstable': [],
        'omega_stable': [],
        'omega_unstable': [],
        'coexistence': np.zeros(n_omega, dtype=bool)
    }

    for i, omega in enumerate(omega_values):
        # 更新模型参数
        model = ThreeSpeciesModel(model_params.copy())
        model.params['omega'] = omega

        # 模拟到平衡态
        sol = model.simulate(N0, (0, 500))
        results['N_S_sim'][i] = sol['N_S'][-1]
        results['N_M_sim'][i] = sol['N_M'][-1]
        results['N_G_sim'][i] = sol['N_G'][-1]

        # 检查是否三物种共存
        final_state = np.array([sol['N_S'][-1], sol['N_M'][-1], sol['N_G'][-1]])
        results['coexistence'][i] = np.all(final_state > 1.0)

        # 寻找平衡点并分析稳定性
        equilibria = model.find_equilibria(n_attempts=30)

        for eq in equilibria:
            # 只关注非零平衡点
            if np.sum(eq > 0.1) > 0:
                stability = model.stability_analysis(eq)

                if stability['stable']:
                    results['N_S_stable'].append(eq[0])
                    results['N_M_stable'].append(eq[1])
                    results['N_G_stable'].append(eq[2])
                    results['omega_stable'].append(omega)
                else:
                    results['N_S_unstable'].append(eq[0])
                    results['N_M_unstable'].append(eq[1])
                    results['N_G_unstable'].append(eq[2])
                    results['omega_unstable'].append(omega)

    return results


def plot_bifurcation_diagram_comprehensive():
    """绘制完整的跨临界分岔图"""

    print("="*70)
    print("绘制跨临界分岔图...")
    print("="*70)

    # 参数设置
    params = {
        'r_S': 1.0, 'r_M': 0.8, 'r_G': 0.9,
        'K_S': 100.0, 'K_M': 100.0, 'K_G': 100.0,
        'sigma_SM': 0.5, 'sigma_MS': 0.6,
        'sigma_SG': 0.3, 'sigma_GS': 0.4,
        'sigma_MG': 0.3, 'sigma_GM': 0.4,
        'alpha_SG': 0.4, 'alpha_MG': 0.4,
        'alpha_GS': 0.3, 'alpha_GM': 0.3,
    }

    # ω值范围
    omega_values = np.linspace(0.0, 1.0, 100)
    N0 = np.array([50.0, 50.0, 50.0])

    print("\n计算分岔数据...")
    results = compute_bifurcation_diagram(omega_values, params, N0)

    # 创建图形
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.25)

    # ========================================================================
    # 图1: 经典分岔图（所有物种）
    # ========================================================================
    ax1 = fig.add_subplot(gs[0, :])

    # 绘制模拟轨迹
    ax1.plot(omega_values, results['N_S_sim'], 'b-', linewidth=2, alpha=0.6, label='S-specialist')
    ax1.plot(omega_values, results['N_M_sim'], 'r-', linewidth=2, alpha=0.6, label='M-specialist')
    ax1.plot(omega_values, results['N_G_sim'], 'g-', linewidth=2, alpha=0.6, label='Generalist')

    # 标记共存区域
    coexist_omega = omega_values[results['coexistence']]
    if len(coexist_omega) > 0:
        omega_min_coexist = coexist_omega.min()
        omega_max_coexist = coexist_omega.max()
        ax1.axvspan(omega_min_coexist, omega_max_coexist,
                   alpha=0.15, color='gold', zorder=0,
                   label='Three-species coexistence')

        # 标注临界点
        ax1.axvline(omega_min_coexist, color='orange', linestyle='--', linewidth=1.5, alpha=0.7)
        ax1.axvline(omega_max_coexist, color='orange', linestyle='--', linewidth=1.5, alpha=0.7)

        ax1.text(omega_min_coexist, 105, r'$\omega_{crit1}$',
                fontsize=12, ha='center', color='orange', fontweight='bold')
        ax1.text(omega_max_coexist, 105, r'$\omega_{crit2}$',
                fontsize=12, ha='center', color='orange', fontweight='bold')

    ax1.set_xlabel(r'Pathway weighting parameter $\omega$', fontsize=13, fontweight='bold')
    ax1.set_ylabel('Population density at equilibrium', fontsize=13, fontweight='bold')
    ax1.set_title('Transcritical Bifurcation Diagram: Complete View',
                 fontsize=15, fontweight='bold', pad=15)
    ax1.legend(loc='upper left', fontsize=11, framealpha=0.9)
    ax1.grid(True, alpha=0.3, linestyle=':')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(-5, 110)

    # ========================================================================
    # 图2: S-specialist 分岔图（详细）
    # ========================================================================
    ax2 = fig.add_subplot(gs[1, 0])

    ax2.plot(omega_values, results['N_S_sim'], 'b-', linewidth=2.5, label='Trajectory')

    # 稳定和不稳定平衡点
    if len(results['omega_stable']) > 0:
        ax2.scatter(results['omega_stable'], results['N_S_stable'],
                   c='darkblue', s=30, marker='o', alpha=0.6,
                   label='Stable equilibria', zorder=5)
    if len(results['omega_unstable']) > 0:
        ax2.scatter(results['omega_unstable'], results['N_S_unstable'],
                   c='lightblue', s=30, marker='x', alpha=0.6,
                   label='Unstable equilibria', zorder=5)

    ax2.set_xlabel(r'$\omega$', fontsize=12, fontweight='bold')
    ax2.set_ylabel(r'$N_S$ (S-specialist)', fontsize=12, fontweight='bold')
    ax2.set_title('S-specialist Bifurcation', fontsize=13, fontweight='bold')
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 1)

    # ========================================================================
    # 图3: M-specialist 分岔图（详细）
    # ========================================================================
    ax3 = fig.add_subplot(gs[1, 1])

    ax3.plot(omega_values, results['N_M_sim'], 'r-', linewidth=2.5, label='Trajectory')

    if len(results['omega_stable']) > 0:
        ax3.scatter(results['omega_stable'], results['N_M_stable'],
                   c='darkred', s=30, marker='o', alpha=0.6,
                   label='Stable equilibria', zorder=5)
    if len(results['omega_unstable']) > 0:
        ax3.scatter(results['omega_unstable'], results['N_M_unstable'],
                   c='lightcoral', s=30, marker='x', alpha=0.6,
                   label='Unstable equilibria', zorder=5)

    # 标记M灭绝区域
    M_extinct = results['N_M_sim'] < 1.0
    if np.any(M_extinct):
        extinct_regions = np.where(M_extinct)[0]
        if len(extinct_regions) > 0:
            ax3.axhline(1.0, color='red', linestyle=':', alpha=0.5, label='Extinction threshold')

    ax3.set_xlabel(r'$\omega$', fontsize=12, fontweight='bold')
    ax3.set_ylabel(r'$N_M$ (M-specialist)', fontsize=12, fontweight='bold')
    ax3.set_title('M-specialist Bifurcation (Obligate Cross-feeder)',
                 fontsize=13, fontweight='bold')
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 1)

    # ========================================================================
    # 图4: Generalist 分岔图（详细）
    # ========================================================================
    ax4 = fig.add_subplot(gs[2, 0])

    ax4.plot(omega_values, results['N_G_sim'], 'g-', linewidth=2.5, label='Trajectory')

    if len(results['omega_stable']) > 0:
        ax4.scatter(results['omega_stable'], results['N_G_stable'],
                   c='darkgreen', s=30, marker='o', alpha=0.6,
                   label='Stable equilibria', zorder=5)
    if len(results['omega_unstable']) > 0:
        ax4.scatter(results['omega_unstable'], results['N_G_unstable'],
                   c='lightgreen', s=30, marker='x', alpha=0.6,
                   label='Unstable equilibria', zorder=5)

    ax4.set_xlabel(r'$\omega$', fontsize=12, fontweight='bold')
    ax4.set_ylabel(r'$N_G$ (Generalist)', fontsize=12, fontweight='bold')
    ax4.set_title('Generalist Bifurcation', fontsize=13, fontweight='bold')
    ax4.legend(loc='best', fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(0, 1)

    # ========================================================================
    # 图5: 相图（状态空间区域）
    # ========================================================================
    ax5 = fig.add_subplot(gs[2, 1])

    # 计算物种组成
    species_composition = np.zeros((len(omega_values), 4))
    for i in range(len(omega_values)):
        N_S = results['N_S_sim'][i]
        N_M = results['N_M_sim'][i]
        N_G = results['N_G_sim'][i]

        threshold = 1.0
        if N_S > threshold and N_M > threshold and N_G > threshold:
            species_composition[i, 0] = 1  # 三物种
        elif N_S > threshold and N_M > threshold:
            species_composition[i, 1] = 1  # S-M
        elif N_S > threshold and N_G > threshold:
            species_composition[i, 2] = 1  # S-G
        else:
            species_composition[i, 3] = 1  # 其他

    # 堆叠区域图
    ax5.fill_between(omega_values, 0, species_composition[:, 0],
                     alpha=0.7, color='gold', label='S-M-G coexistence')
    ax5.fill_between(omega_values, 0, species_composition[:, 1],
                     alpha=0.7, color='purple', label='S-M coexistence')
    ax5.fill_between(omega_values, 0, species_composition[:, 2],
                     alpha=0.7, color='cyan', label='S-G coexistence')
    ax5.fill_between(omega_values, 0, species_composition[:, 3],
                     alpha=0.7, color='gray', label='Other')

    ax5.set_xlabel(r'$\omega$', fontsize=12, fontweight='bold')
    ax5.set_ylabel('Community State', fontsize=12, fontweight='bold')
    ax5.set_title('Community Composition Phase Diagram', fontsize=13, fontweight='bold')
    ax5.set_yticks([0, 1])
    ax5.set_yticklabels(['Absent', 'Present'])
    ax5.legend(loc='upper right', fontsize=10)
    ax5.grid(True, alpha=0.3, axis='x')
    ax5.set_xlim(0, 1)

    # 添加总标题
    fig.suptitle('Transcritical Bifurcation Analysis: Three-Species Cross-Feeding Model',
                fontsize=17, fontweight='bold', y=0.995)

    # 保存
    plt.savefig('figures/transcritical_bifurcation_diagram.png',
                dpi=300, bbox_inches='tight')
    print("\n✓ 分岔图已保存: figures/transcritical_bifurcation_diagram.png")

    plt.show()

    # 输出分析结果
    print("\n" + "="*70)
    print("分岔分析结果:")
    print("="*70)

    if len(coexist_omega) > 0:
        print(f"\n三物种共存区间: ω ∈ [{omega_min_coexist:.3f}, {omega_max_coexist:.3f}]")
        print(f"共存窗口宽度: Δω = {omega_max_coexist - omega_min_coexist:.3f}")
        print(f"参数空间占比: {100*len(coexist_omega)/len(omega_values):.1f}%")

        print(f"\n临界点:")
        print(f"  ω_crit1 ≈ {omega_min_coexist:.3f} (G入侵S-M平衡点)")
        print(f"  ω_crit2 ≈ {omega_max_coexist:.3f} (M无法继续生存)")
    else:
        print("\n⚠ 当前参数下无三物种共存")
        print("建议: 增大 σ_MS 或调整其他合作系数")

    print("\n" + "="*70)


def plot_bifurcation_schematic():
    """绘制分岔示意图（理论图）"""

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # ========================================================================
    # 左图: 跨临界分岔示意图
    # ========================================================================
    ax1 = axes[0]

    omega = np.linspace(0, 1, 200)

    # 理论曲线
    # S-M平衡点的N_G分支（ω < ω_crit1时为0）
    omega_crit1 = 0.3
    omega_crit2 = 0.7

    # S-M平衡点的分支
    N_G_SM = np.where(omega < omega_crit1, 0, np.nan)

    # 三物种平衡点的N_G分支（从ω_crit1开始增长）
    N_G_3sp = np.where((omega >= omega_crit1) & (omega <= omega_crit2),
                       80*(omega - omega_crit1)/(omega_crit2 - omega_crit1),
                       np.nan)

    # S-G平衡点的分支
    N_G_SG = np.where(omega > omega_crit2,
                     80 + 15*(omega - omega_crit2)/(1 - omega_crit2),
                     np.nan)

    # 绘制
    ax1.plot(omega, N_G_SM, 'g-', linewidth=3, label='S-M equilibrium (N_G=0)')
    ax1.plot(omega, N_G_3sp, 'g-', linewidth=3, label='Three-species equilibrium')
    ax1.plot(omega, N_G_SG, 'g-', linewidth=3, label='S-G equilibrium')

    # 不稳定分支
    N_G_unstable1 = np.where(omega >= omega_crit1, 0, np.nan)
    ax1.plot(omega, N_G_unstable1, 'g--', linewidth=2, alpha=0.5, label='Unstable branch')

    # 标记分岔点
    ax1.plot(omega_crit1, 0, 'ro', markersize=12, label='Bifurcation point')
    ax1.plot(omega_crit2, 80*(omega_crit2-omega_crit1)/(omega_crit2-omega_crit1),
            'ro', markersize=12)

    # 注释
    ax1.annotate(r'$\omega_{crit1}$', xy=(omega_crit1, 0),
                xytext=(omega_crit1-0.1, 20),
                fontsize=13, fontweight='bold', color='red',
                arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

    ax1.annotate(r'$\omega_{crit2}$', xy=(omega_crit2, 80),
                xytext=(omega_crit2+0.1, 60),
                fontsize=13, fontweight='bold', color='red',
                arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

    # 区域标注
    ax1.text(0.15, 50, 'S-M\nstable', fontsize=12, ha='center',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax1.text(0.5, 50, 'Three-species\ncoexistence', fontsize=12, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    ax1.text(0.85, 50, 'S-G\nstable', fontsize=12, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

    ax1.set_xlabel(r'Bifurcation parameter $\omega$', fontsize=13, fontweight='bold')
    ax1.set_ylabel(r'$N_G$ (Generalist density)', fontsize=13, fontweight='bold')
    ax1.set_title('Transcritical Bifurcation Schematic', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle=':')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(-10, 110)
    ax1.legend(loc='upper left', fontsize=10)

    # ========================================================================
    # 右图: 分岔机制示意图
    # ========================================================================
    ax2 = axes[1]
    ax2.axis('off')

    # 绘制流程图
    box_props = dict(boxstyle='round,pad=0.5', facecolor='lightblue',
                    edgecolor='black', linewidth=2)

    # ω < ω_crit1
    ax2.text(0.5, 0.85, r'$\omega < \omega_{crit1}$',
            transform=ax2.transAxes, fontsize=14, fontweight='bold',
            ha='center', bbox=box_props)
    ax2.text(0.5, 0.75, 'S-M coexistence\nG excluded\n(G cannot invade)',
            transform=ax2.transAxes, fontsize=11, ha='center',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    # 箭头
    ax2.annotate('', xy=(0.5, 0.65), xytext=(0.5, 0.72),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=2, color='red'))
    ax2.text(0.52, 0.685, r'$\omega$ increases', transform=ax2.transAxes,
            fontsize=10, color='red', fontweight='bold')

    # ω = ω_crit1
    ax2.text(0.5, 0.60, r'$\omega = \omega_{crit1}$',
            transform=ax2.transAxes, fontsize=14, fontweight='bold',
            ha='center', bbox=dict(boxstyle='round', facecolor='orange',
                                  edgecolor='red', linewidth=2))
    ax2.text(0.5, 0.52, 'Transcritical bifurcation\nG invades successfully',
            transform=ax2.transAxes, fontsize=11, ha='center',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

    # 箭头
    ax2.annotate('', xy=(0.5, 0.42), xytext=(0.5, 0.49),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=2, color='red'))

    # ω_crit1 < ω < ω_crit2
    ax2.text(0.5, 0.37, r'$\omega_{crit1} < \omega < \omega_{crit2}$',
            transform=ax2.transAxes, fontsize=14, fontweight='bold',
            ha='center', bbox=box_props)
    ax2.text(0.5, 0.27, 'THREE-SPECIES COEXISTENCE\nAll species present',
            transform=ax2.transAxes, fontsize=11, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    # 箭头
    ax2.annotate('', xy=(0.5, 0.17), xytext=(0.5, 0.24),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', lw=2, color='red'))

    # ω > ω_crit2
    ax2.text(0.5, 0.12, r'$\omega > \omega_{crit2}$',
            transform=ax2.transAxes, fontsize=14, fontweight='bold',
            ha='center', bbox=box_props)
    ax2.text(0.5, 0.02, 'S-G coexistence\nM excluded\n(M cannot survive)',
            transform=ax2.transAxes, fontsize=11, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))

    ax2.set_title('Bifurcation Mechanism', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig('figures/bifurcation_schematic.png', dpi=300, bbox_inches='tight')
    print("\n✓ 分岔示意图已保存: figures/bifurcation_schematic.png")
    plt.show()


if __name__ == "__main__":
    print("\n" + "🔬"*35)
    print("跨临界分岔分析 - Transcritical Bifurcation Analysis")
    print("🔬"*35 + "\n")

    # 绘制完整分岔图
    plot_bifurcation_diagram_comprehensive()

    print("\n" + "-"*70)

    # 绘制理论示意图
    print("\n绘制分岔机制示意图...")
    plot_bifurcation_schematic()

    print("\n" + "="*70)
    print("✅ 所有分岔图已生成完成！")
    print("="*70)
    print("\n生成的图片:")
    print("  1. figures/transcritical_bifurcation_diagram.png - 完整分岔图")
    print("  2. figures/bifurcation_schematic.png - 分岔机制示意图")
    print("\n" + "="*70)

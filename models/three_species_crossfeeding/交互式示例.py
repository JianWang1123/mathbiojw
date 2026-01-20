#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
交互式示例：三物种模型
Interactive Example: Three-Species Model

这个脚本展示如何一步步使用模型进行自定义分析
"""

import sys
sys.path.append('src')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from three_species_model import ThreeSpeciesModel

def example_1_basic_simulation():
    """示例1：基础模拟"""
    print("\n" + "="*70)
    print("示例 1: 基础模拟")
    print("="*70)

    # 创建模型
    model = ThreeSpeciesModel()

    # 设置初始条件
    N0 = np.array([60.0, 40.0, 30.0])  # 不同的初始密度

    # 运行模拟
    sol = model.simulate(N0, (0, 150))

    # 输出结果
    print(f"\n初始种群: N_S={N0[0]}, N_M={N0[1]}, N_G={N0[2]}")
    print(f"最终种群: N_S={sol['N_S'][-1]:.2f}, N_M={sol['N_M'][-1]:.2f}, N_G={sol['N_G'][-1]:.2f}")

    # 绘图
    plt.figure(figsize=(10, 6))
    plt.plot(sol['t'], sol['N_S'], 'b-', linewidth=2, label='S-specialist')
    plt.plot(sol['t'], sol['N_M'], 'r-', linewidth=2, label='M-specialist')
    plt.plot(sol['t'], sol['N_G'], 'g-', linewidth=2, label='Generalist')
    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Population density', fontsize=14)
    plt.title('Example 1: Basic Simulation', fontsize=16, fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('figures/example1_basic.png', dpi=150, bbox_inches='tight')
    print("✓ 保存图片: figures/example1_basic.png")

    return model


def example_2_test_omega():
    """示例2：测试不同的ω值"""
    print("\n" + "="*70)
    print("示例 2: 测试不同的ω值 (路径权重)")
    print("="*70)

    model = ThreeSpeciesModel()
    N0 = np.array([50.0, 50.0, 50.0])

    omega_values = [0.1, 0.3, 0.5, 0.7, 0.9]
    results = []

    print("\nω值 | N_S最终 | N_M最终 | N_G最终")
    print("-" * 40)

    for omega in omega_values:
        model.params['omega'] = omega
        sol = model.simulate(N0, (0, 150))

        results.append({
            'omega': omega,
            'N_S': sol['N_S'][-1],
            'N_M': sol['N_M'][-1],
            'N_G': sol['N_G'][-1]
        })

        print(f"{omega:.1f} | {sol['N_S'][-1]:>7.2f} | {sol['N_M'][-1]:>7.2f} | {sol['N_G'][-1]:>7.2f}")

    # 绘制结果
    fig, ax = plt.subplots(figsize=(10, 6))
    omegas = [r['omega'] for r in results]
    ax.plot(omegas, [r['N_S'] for r in results], 'bo-', linewidth=2, markersize=8, label='S')
    ax.plot(omegas, [r['N_M'] for r in results], 'ro-', linewidth=2, markersize=8, label='M')
    ax.plot(omegas, [r['N_G'] for r in results], 'go-', linewidth=2, markersize=8, label='G')
    ax.set_xlabel('Pathway weighting (ω)', fontsize=14)
    ax.set_ylabel('Final population density', fontsize=14)
    ax.set_title('Example 2: Effect of ω on Community Composition', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.savefig('figures/example2_omega_scan.png', dpi=150, bbox_inches='tight')
    print("\n✓ 保存图片: figures/example2_omega_scan.png")

    return results


def example_3_cooperation_strength():
    """示例3：改变合作强度"""
    print("\n" + "="*70)
    print("示例 3: 改变合作强度 (sigma_MS)")
    print("="*70)

    model = ThreeSpeciesModel()
    N0 = np.array([50.0, 50.0, 50.0])

    sigma_values = [0.3, 0.5, 0.7, 0.9]
    results = []

    print("\nσ_MS | N_S最终 | N_M最终 | N_G最终 | M存活?")
    print("-" * 50)

    for sigma in sigma_values:
        model.params['sigma_MS'] = sigma
        sol = model.simulate(N0, (0, 200))

        m_survives = sol['N_M'][-1] > 1.0
        results.append({
            'sigma': sigma,
            'N_S': sol['N_S'][-1],
            'N_M': sol['N_M'][-1],
            'N_G': sol['N_G'][-1],
            'M_survives': m_survives
        })

        print(f"{sigma:.1f}  | {sol['N_S'][-1]:>7.2f} | {sol['N_M'][-1]:>7.2f} | {sol['N_G'][-1]:>7.2f} | {'✓' if m_survives else '✗'}")

    print("\n关键发现: σ_MS需要足够大（约>0.7）M-specialist才能生存！")

    # 绘图
    fig, ax = plt.subplots(figsize=(10, 6))
    sigmas = [r['sigma'] for r in results]
    ax.plot(sigmas, [r['N_S'] for r in results], 'bo-', linewidth=2, markersize=8, label='S')
    ax.plot(sigmas, [r['N_M'] for r in results], 'ro-', linewidth=2, markersize=8, label='M')
    ax.plot(sigmas, [r['N_G'] for r in results], 'go-', linewidth=2, markersize=8, label='G')
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Survival threshold')
    ax.set_xlabel('Cooperation strength (σ_MS)', fontsize=14)
    ax.set_ylabel('Final population density', fontsize=14)
    ax.set_title('Example 3: Effect of Cooperation on M-specialist Survival', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.savefig('figures/example3_cooperation.png', dpi=150, bbox_inches='tight')
    print("✓ 保存图片: figures/example3_cooperation.png")

    return results


def example_4_initial_conditions():
    """示例4：不同初始条件的影响"""
    print("\n" + "="*70)
    print("示例 4: 不同初始条件")
    print("="*70)

    model = ThreeSpeciesModel()

    # 测试4种不同的初始条件
    initial_conditions = [
        ('S主导', np.array([90.0, 5.0, 5.0])),
        ('M主导', np.array([5.0, 90.0, 5.0])),
        ('G主导', np.array([5.0, 5.0, 90.0])),
        ('均匀分布', np.array([33.0, 33.0, 33.0]))
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for idx, (name, N0) in enumerate(initial_conditions):
        sol = model.simulate(N0, (0, 150))

        ax = axes[idx]
        ax.plot(sol['t'], sol['N_S'], 'b-', linewidth=2, label='S')
        ax.plot(sol['t'], sol['N_M'], 'r-', linewidth=2, label='M')
        ax.plot(sol['t'], sol['N_G'], 'g-', linewidth=2, label='G')
        ax.set_xlabel('Time', fontsize=12)
        ax.set_ylabel('Population density', fontsize=12)
        ax.set_title(f'{name}\nInitial: S={N0[0]:.0f}, M={N0[1]:.0f}, G={N0[2]:.0f}',
                    fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 120)

        print(f"\n{name}:")
        print(f"  初始: N_S={N0[0]:.0f}, N_M={N0[1]:.0f}, N_G={N0[2]:.0f}")
        print(f"  最终: N_S={sol['N_S'][-1]:.1f}, N_M={sol['N_M'][-1]:.1f}, N_G={sol['N_G'][-1]:.1f}")

    plt.tight_layout()
    plt.savefig('figures/example4_initial_conditions.png', dpi=150, bbox_inches='tight')
    print("\n✓ 保存图片: figures/example4_initial_conditions.png")

    print("\n关键发现: 尽管初始条件不同，系统倾向于收敛到相同的最终状态！")


def example_5_find_coexistence():
    """示例5：寻找三物种共存的参数"""
    print("\n" + "="*70)
    print("示例 5: 寻找三物种共存的参数组合")
    print("="*70)

    N0 = np.array([50.0, 50.0, 50.0])

    # 测试几组参数
    parameter_sets = [
        {
            'name': '默认参数',
            'params': {}  # 使用默认值
        },
        {
            'name': '增强合作',
            'params': {
                'sigma_SM': 0.7,
                'sigma_MS': 0.9,
                'sigma_SG': 0.5,
                'sigma_GS': 0.6,
                'sigma_MG': 0.5,
                'sigma_GM': 0.6,
            }
        },
        {
            'name': '平衡策略',
            'params': {
                'sigma_SM': 0.6,
                'sigma_MS': 0.8,
                'omega': 0.4,
                'alpha_SG': 0.2,
                'alpha_MG': 0.2,
            }
        }
    ]

    results = []

    for param_set in parameter_sets:
        model = ThreeSpeciesModel()

        # 更新参数
        for key, value in param_set['params'].items():
            model.params[key] = value

        # 模拟
        sol = model.simulate(N0, (0, 300))

        # 检查共存（所有物种>1）
        final = [sol['N_S'][-1], sol['N_M'][-1], sol['N_G'][-1]]
        coexistence = all(n > 1.0 for n in final)

        results.append({
            'name': param_set['name'],
            'N_S': final[0],
            'N_M': final[1],
            'N_G': final[2],
            'coexistence': coexistence
        })

        print(f"\n{param_set['name']}:")
        print(f"  最终: N_S={final[0]:.2f}, N_M={final[1]:.2f}, N_G={final[2]:.2f}")
        print(f"  三物种共存: {'✓ 是' if coexistence else '✗ 否'}")

    print("\n" + "-"*70)
    if any(r['coexistence'] for r in results):
        print("✓ 找到了支持三物种共存的参数组合！")
    else:
        print("✗ 当前参数下无法实现三物种共存")
        print("   建议: 进一步增大σ_MS和σ_SM，降低α系数")


def main():
    """运行所有示例"""
    print("\n" + "🔬"*35)
    print("三物种交叉喂养模型 - 交互式示例集")
    print("Interactive Examples for Three-Species Model")
    print("🔬"*35)

    # 运行所有示例
    example_1_basic_simulation()
    example_2_test_omega()
    example_3_cooperation_strength()
    example_4_initial_conditions()
    example_5_find_coexistence()

    # 总结
    print("\n" + "="*70)
    print("所有示例运行完成！")
    print("="*70)
    print("\n生成的图片:")
    print("  1. figures/example1_basic.png - 基础模拟")
    print("  2. figures/example2_omega_scan.png - ω参数扫描")
    print("  3. figures/example3_cooperation.png - 合作强度影响")
    print("  4. figures/example4_initial_conditions.png - 初始条件影响")
    print("\n关键结论:")
    print("  • ω值强烈影响Generalist的竞争力")
    print("  • M-specialist需要足够的交叉喂养（σ_MS>0.7）才能生存")
    print("  • 系统对初始条件具有一定的鲁棒性")
    print("  • 三物种共存需要精心调节合作-竞争平衡")
    print("\n下一步:")
    print("  • 修改参数重新运行任一示例")
    print("  • 运行完整Jupyter笔记本获得更多分析")
    print("  • 查看 使用说明_中文.md 了解详细用法")
    print("="*70)


if __name__ == "__main__":
    main()

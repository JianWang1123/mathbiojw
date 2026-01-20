#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
快速入门：三物种交叉喂养模型
Quick Start Guide for Three-Species Cross-Feeding Model

这个脚本演示如何运行模型的最基本功能
This script demonstrates the most basic functionality

运行方法 (How to run):
    python3 quick_start.py
"""

import sys
sys.path.append('src')

import numpy as np
import matplotlib
matplotlib.use('Agg')  # 不需要显示窗口
import matplotlib.pyplot as plt

print("="*70)
print("三物种交叉喂养模型 - 快速入门")
print("Three-Species Cross-Feeding Model - Quick Start")
print("="*70)

# ============================================================================
# 步骤 1: 导入模型
# Step 1: Import the model
# ============================================================================
print("\n步骤 1: 导入模型...")
print("Step 1: Importing model...")

try:
    from three_species_model import ThreeSpeciesModel
    print("✓ 模型导入成功！")
    print("✓ Model imported successfully!")
except ImportError as e:
    print(f"✗ 导入失败: {e}")
    print("请先安装依赖: pip install numpy scipy matplotlib")
    print("Please install dependencies: pip install numpy scipy matplotlib")
    sys.exit(1)

# ============================================================================
# 步骤 2: 创建模型实例
# Step 2: Create model instance
# ============================================================================
print("\n步骤 2: 创建模型实例...")
print("Step 2: Creating model instance...")

model = ThreeSpeciesModel()
print("✓ 模型创建成功！")
print("✓ Model created successfully!")

# 显示参数
print("\n当前参数 (Current parameters):")
print(f"  生长率 (Growth rates): r_S={model.params['r_S']}, r_M={model.params['r_M']}, r_G={model.params['r_G']}")
print(f"  路径权重 (Pathway weight): ω = {model.params['omega']}")

# ============================================================================
# 步骤 3: 设置初始条件
# Step 3: Set initial conditions
# ============================================================================
print("\n步骤 3: 设置初始条件...")
print("Step 3: Setting initial conditions...")

N0 = np.array([50.0, 50.0, 50.0])  # 三个物种的初始密度相同
print(f"  初始种群 (Initial populations): N_S={N0[0]}, N_M={N0[1]}, N_G={N0[2]}")

# ============================================================================
# 步骤 4: 运行模拟
# Step 4: Run simulation
# ============================================================================
print("\n步骤 4: 运行模拟 (t = 0 到 100)...")
print("Step 4: Running simulation (t = 0 to 100)...")

sol = model.simulate(N0, (0, 100))

if sol['success']:
    print("✓ 模拟成功！")
    print("✓ Simulation successful!")
    print(f"\n最终种群密度 (Final populations):")
    print(f"  N_S = {sol['N_S'][-1]:.2f}")
    print(f"  N_M = {sol['N_M'][-1]:.2f}")
    print(f"  N_G = {sol['N_G'][-1]:.2f}")
else:
    print(f"✗ 模拟失败: {sol['message']}")
    sys.exit(1)

# ============================================================================
# 步骤 5: 绘制时间序列图
# Step 5: Plot time series
# ============================================================================
print("\n步骤 5: 绘制结果...")
print("Step 5: Plotting results...")

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(sol['t'], sol['N_S'], 'b-', linewidth=2, label='S (底物专家 / Substrate specialist)')
ax.plot(sol['t'], sol['N_M'], 'r-', linewidth=2, label='M (代谢物专家 / Metabolite specialist)')
ax.plot(sol['t'], sol['N_G'], 'g-', linewidth=2, label='G (通才 / Generalist)')

ax.set_xlabel('时间 (Time)', fontsize=14)
ax.set_ylabel('种群密度 (Population density)', fontsize=14)
ax.set_title('三物种动力学 (Three-Species Dynamics)', fontsize=16, fontweight='bold')
ax.legend(loc='best', fontsize=11)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/quick_start_timeseries.png', dpi=150, bbox_inches='tight')
print("✓ 图片已保存到: figures/quick_start_timeseries.png")
print("✓ Figure saved to: figures/quick_start_timeseries.png")

# ============================================================================
# 步骤 6: 寻找平衡点
# Step 6: Find equilibria
# ============================================================================
print("\n步骤 6: 寻找平衡点...")
print("Step 6: Finding equilibria...")

equilibria = model.find_equilibria(n_attempts=30)
print(f"✓ 找到 {len(equilibria)} 个平衡点")
print(f"✓ Found {len(equilibria)} equilibria")

for i, eq in enumerate(equilibria):
    eco_type = model.classify_equilibrium_ecology(eq)
    stability = model.stability_analysis(eq)

    print(f"\n  平衡点 {i+1} (Equilibrium {i+1}):")
    print(f"    种群 (Populations): N_S={eq[0]:.2f}, N_M={eq[1]:.2f}, N_G={eq[2]:.2f}")
    print(f"    类型 (Type): {eco_type}")
    print(f"    稳定性 (Stability): {'稳定 (Stable)' if stability['stable'] else '不稳定 (Unstable)'}")

# ============================================================================
# 步骤 7: 测试不同的 ω 值
# Step 7: Test different ω values
# ============================================================================
print("\n步骤 7: 测试不同的路径权重 (ω)...")
print("Step 7: Testing different pathway weights (ω)...")

omega_values = [0.2, 0.5, 0.8]
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for idx, omega in enumerate(omega_values):
    model.params['omega'] = omega
    sol = model.simulate(N0, (0, 100))

    ax = axes[idx]
    ax.plot(sol['t'], sol['N_S'], 'b-', linewidth=2, label='S')
    ax.plot(sol['t'], sol['N_M'], 'r-', linewidth=2, label='M')
    ax.plot(sol['t'], sol['N_G'], 'g-', linewidth=2, label='G')

    ax.set_xlabel('时间 (Time)', fontsize=12)
    ax.set_ylabel('种群密度', fontsize=12)
    ax.set_title(f'ω = {omega}', fontsize=14, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    print(f"  ω = {omega}: 最终种群 = N_S={sol['N_S'][-1]:.1f}, N_M={sol['N_M'][-1]:.1f}, N_G={sol['N_G'][-1]:.1f}")

plt.tight_layout()
plt.savefig('figures/quick_start_omega_comparison.png', dpi=150, bbox_inches='tight')
print("\n✓ 图片已保存到: figures/quick_start_omega_comparison.png")
print("✓ Figure saved to: figures/quick_start_omega_comparison.png")

# ============================================================================
# 总结
# Summary
# ============================================================================
print("\n" + "="*70)
print("运行完成！")
print("Run completed successfully!")
print("="*70)
print("\n生成的文件 (Generated files):")
print("  1. figures/quick_start_timeseries.png - 时间序列图")
print("  2. figures/quick_start_omega_comparison.png - ω 比较图")
print("\n下一步 (Next steps):")
print("  • 查看生成的图片了解结果")
print("  • 运行完整分析: jupyter notebook notebooks/three_species_phase_analysis.ipynb")
print("  • 修改参数重新运行此脚本")
print("="*70)

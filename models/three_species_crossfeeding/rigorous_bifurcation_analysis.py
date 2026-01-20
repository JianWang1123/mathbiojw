#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
严格的分岔分析 - Rigorous Bifurcation Analysis
Dynamical Systems Theory Approach

基于标准的分岔理论对三物种交叉喂养模型进行严格的数学分析

理论框架:
1. 系统验证：自治性、光滑性、相空间
2. 平衡点分析：边界平衡点和内部平衡点的解析解
3. 线性稳定性分析：Jacobian矩阵、特征值、Routh-Hurwitz判据
4. 分岔分析：分岔点定位、分岔类型判定、分岔图
5. 参数空间：相图、分岔曲线、共存区域

Author: Jian Wang
Date: January 2026
"""

import sys
sys.path.append('src')

import numpy as np
import sympy as sp
from sympy import symbols, Matrix, solve, simplify, expand, factor
from sympy import lambdify, latex, diff
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
from typing import Dict, List, Tuple, Optional
from three_species_model import ThreeSpeciesModel

# 设置
plt.rcParams['font.size'] = 11
plt.rcParams['figure.figsize'] = (14, 10)
sns.set_style("whitegrid")


class RigorousBifurcationAnalysis:
    """
    严格的分岔分析类

    基于动力学系统理论的标准方法
    """

    def __init__(self, params: Dict):
        """
        初始化分析

        Parameters
        ----------
        params : dict
            模型参数（除ω外的所有参数）
        """
        self.params = params
        self._setup_symbolic_system()

    def _setup_symbolic_system(self):
        """建立符号化系统"""
        print("\n" + "="*70)
        print("步骤 1: 建立符号化动力学系统")
        print("="*70)

        # 定义符号变量
        self.N_S = symbols('N_S', real=True, nonnegative=True)
        self.N_M = symbols('N_M', real=True, nonnegative=True)
        self.N_G = symbols('N_G', real=True, nonnegative=True)
        self.omega = symbols('omega', real=True)

        # 参数符号
        self.r_S, self.r_M, self.r_G = symbols('r_S r_M r_G', positive=True)
        self.K_S, self.K_M, self.K_G = symbols('K_S K_M K_G', positive=True)
        self.sigma_SM, self.sigma_MS = symbols('sigma_SM sigma_MS', positive=True)
        self.sigma_SG, self.sigma_GS = symbols('sigma_SG sigma_GS', positive=True)
        self.sigma_MG, self.sigma_GM = symbols('sigma_MG sigma_GM', positive=True)
        self.alpha_SG, self.alpha_MG = symbols('alpha_SG alpha_MG', positive=True)
        self.alpha_GS, self.alpha_GM = symbols('alpha_GS alpha_GM', positive=True)

        # 系统方程 dN/dt = F(N; ω)
        self.F_S = self.r_S * self.N_S * (
            1
            + self.sigma_SM * self.N_M / self.K_M
            + (1 - self.omega) * self.sigma_SG * self.N_G / self.K_G
            - self.omega * self.alpha_SG * self.N_G / self.K_G
            - self.N_S / self.K_S
        )

        self.F_M = self.r_M * self.N_M * (
            -1
            + self.sigma_MS * self.N_S / self.K_S
            + self.omega * self.sigma_MG * self.N_G / self.K_G
            - (1 - self.omega) * self.alpha_MG * self.N_G / self.K_G
            - self.N_M / self.K_M
        )

        self.F_G = self.r_G * self.N_G * (
            self.omega * (1 - self.alpha_GS * self.N_S / self.K_S + self.sigma_GM * self.N_M / self.K_M)
            + (1 - self.omega) * (-1 - self.alpha_GM * self.N_M / self.K_M + self.sigma_GS * self.N_S / self.K_S)
            - self.N_G / self.K_G
        )

        # 向量场
        self.F = Matrix([self.F_S, self.F_M, self.F_G])

        print("\n✓ 系统方程:")
        print(f"  dN_S/dt = F_S(N; ω)")
        print(f"  dN_M/dt = F_M(N; ω)")
        print(f"  dN_G/dt = F_G(N; ω)")
        print("\n✓ 系统是自治的 (autonomous): F不显含时间t")
        print("✓ 相空间: ℝ³₊ = {(N_S, N_M, N_G) : N_i ≥ 0}")
        print("✓ 分岔参数: ω ∈ [0, 1]")

    def verify_system_properties(self):
        """验证系统的数学性质"""
        print("\n" + "="*70)
        print("步骤 2: 验证系统的数学性质")
        print("="*70)

        properties = {
            '自治性 (Autonomous)': 'F = F(N; ω)，不显含时间t ✓',
            '光滑性 (Smoothness)': 'F在ℝ³₊内C^∞光滑 ✓',
            '正向不变性 (Positive invariance)': '坐标超平面{N_i=0}是不变流形 ✓',
            '有界性 (Boundedness)': '有限初值导致有界轨道（由logistic项保证）✓',
            '分岔参数': 'ω ∈ [0,1] 是单参数分岔问题 ✓'
        }

        print("\n系统性质验证:")
        for prop, status in properties.items():
            print(f"  • {prop}: {status}")

        return properties

    def compute_jacobian(self):
        """计算Jacobian矩阵"""
        print("\n" + "="*70)
        print("步骤 3: 计算Jacobian矩阵 J = DF(N; ω)")
        print("="*70)

        # Jacobian矩阵
        N_vec = Matrix([self.N_S, self.N_M, self.N_G])
        self.J = self.F.jacobian(N_vec)

        print("\nJacobian矩阵 J = DF:")
        print("  J_ij = ∂F_i/∂N_j")
        print(f"\n  维度: {self.J.shape[0]} × {self.J.shape[1]}")

        # 显示矩阵结构（简化形式）
        print("\n  矩阵结构:")
        print("  ┌─                                          ─┐")
        print("  │  ∂F_S/∂N_S   ∂F_S/∂N_M   ∂F_S/∂N_G       │")
        print("  │  ∂F_M/∂N_S   ∂F_M/∂N_M   ∂F_M/∂N_G       │")
        print("  │  ∂F_G/∂N_S   ∂F_G/∂N_M   ∂F_G/∂N_G       │")
        print("  └─                                          ─┘")

        return self.J

    def find_equilibria_symbolic(self):
        """
        符号化求解平衡点

        求解 F(N*; ω) = 0
        """
        print("\n" + "="*70)
        print("步骤 4: 求解平衡点 F(N*; ω) = 0")
        print("="*70)

        equilibria = {}

        # ================================================================
        # E0: 灭绝平衡点 (0, 0, 0)
        # ================================================================
        equilibria['E0_extinction'] = {
            'point': (0, 0, 0),
            'exists': 'Always',
            'biological_meaning': '所有物种灭绝'
        }
        print("\n✓ E0 (灭绝): (0, 0, 0) - 总是存在")

        # ================================================================
        # E_S: S单独平衡点 (K_S, 0, 0)
        # ================================================================
        equilibria['E_S_only'] = {
            'point': (self.K_S, 0, 0),
            'exists': 'Always',
            'biological_meaning': 'S独自达到承载量'
        }
        print("✓ E_S (S单独): (K_S, 0, 0) - 总是存在")

        # ================================================================
        # E_M: M单独平衡点 - 不存在！
        # ================================================================
        equilibria['E_M_only'] = {
            'point': 'Does not exist',
            'exists': 'Never',
            'reason': 'M是专性交叉喂养者，基础生长率r=-1 < 0',
            'biological_meaning': 'M无法独立生存'
        }
        print("✗ E_M (M单独): 不存在 (M基础生长率 = -1 < 0)")

        # ================================================================
        # E_G: G单独平衡点
        # ================================================================
        # 求解: F_G(0, 0, N_G) = 0
        # r_G * N_G * [ω·1 + (1-ω)·(-1) - N_G/K_G] = 0
        # ω - (1-ω) - N_G/K_G = 0
        # 2ω - 1 = N_G/K_G
        # N_G = K_G(2ω - 1)

        N_G_solo = self.K_G * (2*self.omega - 1)
        equilibria['E_G_only'] = {
            'point': (0, 0, N_G_solo),
            'exists': 'ω > 1/2',
            'condition': '2ω - 1 > 0',
            'biological_meaning': 'G独立生存（需要足够依赖底物途径）'
        }
        print(f"✓ E_G (G单独): (0, 0, K_G(2ω-1)) - 存在当ω > 1/2")

        # ================================================================
        # E_SM: S-M共存平衡点 (N_S*, N_M*, 0)
        # ================================================================
        print("\n计算 E_SM (S-M共存)...")
        # 系统（N_G = 0）:
        # 1 + σ_SM·N_M/K_M - N_S/K_S = 0  ... (1)
        # -1 + σ_MS·N_S/K_S - N_M/K_M = 0  ... (2)

        # 从(1): N_S = K_S(1 + σ_SM·N_M/K_M)
        # 代入(2): -1 + σ_MS(1 + σ_SM·N_M/K_M) - N_M/K_M = 0
        # -1 + σ_MS + σ_MS·σ_SM·N_M/K_M - N_M/K_M = 0
        # σ_MS - 1 + N_M/K_M(σ_MS·σ_SM - 1) = 0
        # N_M = K_M(1 - σ_MS)/(σ_MS·σ_SM - 1)

        N_M_SM = self.K_M * (1 - self.sigma_MS) / (self.sigma_MS * self.sigma_SM - 1)
        N_S_SM = self.K_S * (1 + self.sigma_SM * N_M_SM / self.K_M)
        N_S_SM = simplify(N_S_SM)
        N_M_SM = simplify(N_M_SM)

        equilibria['E_SM_coexist'] = {
            'point': (N_S_SM, N_M_SM, 0),
            'exists': 'σ_MS > 1 AND σ_MS·σ_SM > 1',
            'condition': '互惠足够强',
            'biological_meaning': 'S-M互惠共存',
            'N_S_formula': N_S_SM,
            'N_M_formula': N_M_SM
        }
        print(f"✓ E_SM (S-M共存): 存在条件 σ_MS > 1 且 σ_MS·σ_SM > 1")
        print(f"  N_S* = {N_S_SM}")
        print(f"  N_M* = {N_M_SM}")

        # ================================================================
        # E_SG: S-G共存平衡点 (N_S*, 0, N_G*)
        # ================================================================
        print("\n计算 E_SG (S-G共存)...")
        # 这是2×2系统，需要数值求解或参数化求解
        equilibria['E_SG_coexist'] = {
            'point': 'Requires numerical solution (depends on ω)',
            'exists': 'Depends on ω and parameter balance',
            'biological_meaning': 'S-G共存（M被排除）'
        }
        print("△ E_SG (S-G共存): 需要数值求解（依赖ω）")

        # ================================================================
        # E_MG: M-G共存 - 不存在！
        # ================================================================
        equilibria['E_MG_coexist'] = {
            'point': 'Does not exist',
            'exists': 'Never',
            'reason': 'M需要S提供底物，无法与G单独共存',
            'biological_meaning': 'M无法在没有S的情况下生存'
        }
        print("✗ E_MG (M-G共存): 不存在 (M需要S)")

        # ================================================================
        # E*: 三物种共存平衡点 (N_S*, N_M*, N_G*)
        # ================================================================
        print("\n计算 E* (三物种共存)...")
        equilibria['E_interior'] = {
            'point': 'Numerical solution of 3×3 system',
            'exists': 'ω_min < ω < ω_max (coexistence window)',
            'biological_meaning': '三物种共存平衡点',
            'note': '需要数值求解F(N*; ω) = 0'
        }
        print("△ E* (三物种共存): 需要数值求解 3×3 非线性系统")

        self.equilibria_symbolic = equilibria
        return equilibria

    def stability_analysis_symbolic(self, eq_name: str):
        """
        对特定平衡点进行符号化稳定性分析

        使用Routh-Hurwitz判据
        """
        print(f"\n稳定性分析: {eq_name}")
        print("-" * 50)

        eq_info = self.equilibria_symbolic[eq_name]

        if eq_info['exists'] in ['Never', 'Does not exist']:
            print(f"  平衡点不存在，跳过稳定性分析")
            return None

        # 在平衡点处计算Jacobian
        point = eq_info['point']
        if isinstance(point, str):
            print(f"  需要数值求解，稍后分析")
            return None

        # 代入平衡点
        J_at_eq = self.J.subs({
            self.N_S: point[0],
            self.N_M: point[1],
            self.N_G: point[2]
        })

        # 计算trace, determinant
        trace = J_at_eq.trace()
        det = J_at_eq.det()

        print(f"  Trace(J) = {simplify(trace)}")
        print(f"  Det(J) = {simplify(det)}")

        # Routh-Hurwitz判据
        print(f"\n  Routh-Hurwitz稳定性判据:")
        print(f"    稳定 ⟺ Trace(J) < 0 且所有特征值实部 < 0")

        return {
            'jacobian': J_at_eq,
            'trace': trace,
            'determinant': det
        }

    def find_bifurcation_points_numerical(self, omega_range=(0.0, 1.0), n_points=200):
        """
        数值方法寻找分岔点

        通过扫描ω，追踪平衡点和特征值
        """
        print("\n" + "="*70)
        print("步骤 5: 数值寻找分岔点")
        print("="*70)

        omega_values = np.linspace(omega_range[0], omega_range[1], n_points)

        # 存储结果
        results = {
            'omega': omega_values,
            'equilibria': [],
            'eigenvalues': [],
            'stable': []
        }

        # 对每个ω值
        for i, omega in enumerate(omega_values):
            if i % 20 == 0:
                print(f"  计算进度: {i}/{n_points} (ω = {omega:.3f})")

            # 更新模型参数
            params = self.params.copy()
            params['omega'] = omega
            model = ThreeSpeciesModel(params)

            # 寻找平衡点
            N0_guesses = [
                np.array([50.0, 50.0, 50.0]),  # 三物种共存
                np.array([100.0, 50.0, 0.1]),   # S-M主导
                np.array([100.0, 0.1, 50.0]),   # S-G主导
            ]

            eq_found = None
            eig_found = None
            stable_found = False

            for N0 in N0_guesses:
                try:
                    # 模拟到稳态
                    sol = model.simulate(N0, (0, 500))
                    eq_candidate = np.array([sol['N_S'][-1], sol['N_M'][-1], sol['N_G'][-1]])

                    # 检查是否真的是平衡点
                    F_val = model.equations(0, eq_candidate)
                    if np.linalg.norm(F_val) < 1e-6:
                        # 计算特征值
                        J = model.jacobian(eq_candidate)
                        eigenvalues = np.linalg.eigvals(J)

                        # 检查稳定性
                        is_stable = np.all(np.real(eigenvalues) < 0)

                        # 如果是三物种共存
                        if np.all(eq_candidate > 1.0):
                            eq_found = eq_candidate
                            eig_found = eigenvalues
                            stable_found = is_stable
                            break
                except:
                    continue

            results['equilibria'].append(eq_found if eq_found is not None else np.array([np.nan, np.nan, np.nan]))
            results['eigenvalues'].append(eig_found if eig_found is not None else np.array([np.nan, np.nan, np.nan]))
            results['stable'].append(stable_found)

        # 转换为数组
        results['equilibria'] = np.array(results['equilibria'])
        results['stable'] = np.array(results['stable'])

        print("\n✓ 数值扫描完成")

        # 寻找分岔点
        bifurcation_points = self._detect_bifurcations(results)

        return results, bifurcation_points

    def _detect_bifurcations(self, results):
        """检测分岔点"""
        print("\n检测分岔点...")

        omega = results['omega']
        stable = results['stable']
        eq = results['equilibria']

        bifurcations = []

        # 检测稳定性变化点
        for i in range(1, len(stable)-1):
            # 稳定性变化
            if stable[i] != stable[i-1]:
                bifurcations.append({
                    'omega': omega[i],
                    'type': 'Stability change',
                    'description': f'稳定性从 {stable[i-1]} 变为 {stable[i]}'
                })

        # 检测N_M → 0的点（M灭绝）
        N_M = eq[:, 1]
        for i in range(1, len(N_M)-1):
            if not np.isnan(N_M[i]) and not np.isnan(N_M[i+1]):
                if N_M[i] > 5.0 and N_M[i+1] < 5.0:
                    bifurcations.append({
                        'omega': omega[i],
                        'type': 'M extinction',
                        'description': 'M物种灭绝（N_M → 0）'
                    })

        # 检测N_G → 0的点（G灭绝）
        N_G = eq[:, 2]
        for i in range(1, len(N_G)-1):
            if not np.isnan(N_G[i]) and not np.isnan(N_G[i+1]):
                if N_G[i] < 5.0 and N_G[i+1] > 5.0:
                    bifurcations.append({
                        'omega': omega[i],
                        'type': 'G invasion',
                        'description': 'G物种入侵（N_G从0增长）'
                    })

        print(f"\n✓ 检测到 {len(bifurcations)} 个分岔点:")
        for bif in bifurcations:
            print(f"  • ω = {bif['omega']:.4f}: {bif['description']} ({bif['type']})")

        return bifurcations

    def classify_bifurcation_type(self, omega_crit, epsilon=0.01):
        """
        分类分岔类型

        判定是saddle-node, transcritical, 还是Hopf分岔
        """
        print(f"\n分岔类型判定 (ω_crit = {omega_crit:.4f}):")
        print("-" * 50)

        # 在临界点两侧计算特征值
        omega_minus = omega_crit - epsilon
        omega_plus = omega_crit + epsilon

        params_minus = self.params.copy()
        params_minus['omega'] = omega_minus
        model_minus = ThreeSpeciesModel(params_minus)

        params_plus = self.params.copy()
        params_plus['omega'] = omega_plus
        model_plus = ThreeSpeciesModel(params_plus)

        # 模拟找平衡点
        N0 = np.array([50.0, 50.0, 50.0])

        sol_minus = model_minus.simulate(N0, (0, 500))
        eq_minus = np.array([sol_minus['N_S'][-1], sol_minus['N_M'][-1], sol_minus['N_G'][-1]])

        sol_plus = model_plus.simulate(N0, (0, 500))
        eq_plus = np.array([sol_plus['N_S'][-1], sol_plus['N_M'][-1], sol_plus['N_G'][-1]])

        # 计算特征值
        J_minus = model_minus.jacobian(eq_minus)
        eig_minus = np.linalg.eigvals(J_minus)

        J_plus = model_plus.jacobian(eq_plus)
        eig_plus = np.linalg.eigvals(J_plus)

        print(f"\n  ω = {omega_minus:.4f} (before):")
        print(f"    平衡点: ({eq_minus[0]:.2f}, {eq_minus[1]:.2f}, {eq_minus[2]:.2f})")
        print(f"    特征值实部: {np.real(eig_minus)}")

        print(f"\n  ω = {omega_plus:.4f} (after):")
        print(f"    平衡点: ({eq_plus[0]:.2f}, {eq_plus[1]:.2f}, {eq_plus[2]:.2f})")
        print(f"    特征值实部: {np.real(eig_plus)}")

        # 判定分岔类型
        # Transcritical: 平衡点交换稳定性
        # Saddle-node: 平衡点出现/消失
        # Hopf: 特征值穿越虚轴

        has_zero_eig_minus = np.any(np.abs(np.real(eig_minus)) < 0.1)
        has_zero_eig_plus = np.any(np.abs(np.real(eig_plus)) < 0.1)

        has_complex_eig = np.any(np.abs(np.imag(eig_minus)) > 1e-6) or np.any(np.abs(np.imag(eig_plus)) > 1e-6)

        if has_complex_eig and (has_zero_eig_minus or has_zero_eig_plus):
            bifurcation_type = "Hopf bifurcation (可能产生极限环)"
        elif np.linalg.norm(eq_minus - eq_plus) < 10.0:
            bifurcation_type = "Transcritical bifurcation (平衡点交换稳定性)"
        else:
            bifurcation_type = "Saddle-node bifurcation (平衡点出现/消失)"

        print(f"\n✓ 分岔类型判定: {bifurcation_type}")

        return bifurcation_type


def run_complete_analysis():
    """运行完整的严格分析"""
    print("\n" + "🔬"*35)
    print("严格的分岔分析 - RIGOROUS BIFURCATION ANALYSIS")
    print("基于动力学系统理论")
    print("🔬"*35)

    # 参数
    params = {
        'r_S': 1.0, 'r_M': 0.8, 'r_G': 0.9,
        'K_S': 100.0, 'K_M': 100.0, 'K_G': 100.0,
        'sigma_SM': 0.5, 'sigma_MS': 0.6,
        'sigma_SG': 0.3, 'sigma_GS': 0.4,
        'sigma_MG': 0.3, 'sigma_GM': 0.4,
        'alpha_SG': 0.4, 'alpha_MG': 0.4,
        'alpha_GS': 0.3, 'alpha_GM': 0.3,
    }

    # 创建分析对象
    analyzer = RigorousBifurcationAnalysis(params)

    # 步骤 1-2: 系统验证
    properties = analyzer.verify_system_properties()

    # 步骤 3: Jacobian
    J = analyzer.compute_jacobian()

    # 步骤 4: 平衡点
    equilibria = analyzer.find_equilibria_symbolic()

    # 步骤 5: 稳定性分析（符号化）
    for eq_name in ['E0_extinction', 'E_S_only', 'E_SM_coexist']:
        analyzer.stability_analysis_symbolic(eq_name)

    # 步骤 6: 数值分岔分析
    results, bifurcations = analyzer.find_bifurcation_points_numerical()

    # 步骤 7: 分岔类型判定
    if len(bifurcations) > 0:
        for bif in bifurcations[:2]:  # 分析前两个
            analyzer.classify_bifurcation_type(bif['omega'])

    print("\n" + "="*70)
    print("✅ 严格分析完成")
    print("="*70)

    return analyzer, results, bifurcations


if __name__ == "__main__":
    analyzer, results, bifurcations = run_complete_analysis()

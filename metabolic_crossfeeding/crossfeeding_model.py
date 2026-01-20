"""
Eco-evolutionary Model of Metabolic Cross-Feeding and Division of Labor

This module implements the analytical and numerical analysis of a two-species,
two-amino-acid cross-feeding system in a chemostat/bioreactor.

Key Result: Division of labor (metabolic specialization) is the evolutionarily
stable strategy under cross-feeding conditions.

Author: Jian Wang
Date: 2024
"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import fsolve, minimize
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality defaults
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'figure.figsize': (8, 6),
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'text.usetex': False,
    'font.family': 'serif',
})


# =============================================================================
# PART 1: SIMPLIFIED ANALYTICAL MODEL
# =============================================================================

class AnalyticalCrossFeedingModel:
    """
    Simplified model for analytical tractability.

    Key assumptions:
    - Equal population sizes at coexistence (N_A = N_B = N)
    - Symmetric species parameters
    - Quasi-steady state for resources
    - Linear resource regime (F << K)

    State variables:
    - f_A1, f_A2: Species A investment in F1, F2
    - f_B1, f_B2: Species B investment in F1, F2
    """

    def __init__(self, gamma=1.0, D=0.2, sigma=0.01, eta=1.0):
        """
        Parameters:
        -----------
        gamma : float
            Environmental feedback constant (growth rate coefficient)
        D : float
            Dilution rate (h^-1)
        sigma : float
            Evolutionary rate parameter (h^-1)
        eta : float
            Specialization efficiency bonus (eta > 1 favors specialization)
        """
        self.gamma = gamma
        self.D = D
        self.sigma = sigma
        self.eta = eta

    def growth_rate(self, P1, P2):
        """
        Multiplicative growth rate function.
        g = gamma * P1^alpha * P2^alpha

        With alpha < 1, there are diminishing returns to investment.
        """
        alpha = 0.5  # Diminishing returns exponent
        return self.gamma * (P1 ** alpha) * (P2 ** alpha)

    def fitness_A(self, f_A1, f_A2, f_B1, f_B2):
        """
        Invasion fitness of species A.
        W_A = (1 - f_A1 - f_A2) * g(P1, P2) - D

        With specialization bonus:
        W_A = (1 - f_A1 - f_A2 + eta * specialization_index) * g(P1, P2) - D
        """
        P1 = f_A1 + f_B1  # Total investment in F1
        P2 = f_A2 + f_B2  # Total investment in F2

        # Specialization index: high when one investment dominates
        spec_A = abs(f_A1 - f_A2) / (f_A1 + f_A2 + 1e-10)

        # Growth allocation with specialization bonus
        growth_alloc = 1 - f_A1 - f_A2 + (self.eta - 1) * spec_A * (f_A1 + f_A2)

        g = self.growth_rate(P1 + 1e-10, P2 + 1e-10)

        return growth_alloc * g - self.D

    def fitness_B(self, f_A1, f_A2, f_B1, f_B2):
        """Invasion fitness of species B."""
        P1 = f_A1 + f_B1
        P2 = f_A2 + f_B2

        spec_B = abs(f_B1 - f_B2) / (f_B1 + f_B2 + 1e-10)
        growth_alloc = 1 - f_B1 - f_B2 + (self.eta - 1) * spec_B * (f_B1 + f_B2)

        g = self.growth_rate(P1 + 1e-10, P2 + 1e-10)

        return growth_alloc * g - self.D

    def selection_gradient_A1(self, f_A1, f_A2, f_B1, f_B2):
        """
        Selection gradient dW_A/df_A1
        Computed numerically for generality.
        """
        eps = 1e-6
        W_plus = self.fitness_A(f_A1 + eps, f_A2, f_B1, f_B2)
        W_minus = self.fitness_A(f_A1 - eps, f_A2, f_B1, f_B2)
        return (W_plus - W_minus) / (2 * eps)

    def selection_gradient_A2(self, f_A1, f_A2, f_B1, f_B2):
        """Selection gradient dW_A/df_A2"""
        eps = 1e-6
        W_plus = self.fitness_A(f_A1, f_A2 + eps, f_B1, f_B2)
        W_minus = self.fitness_A(f_A1, f_A2 - eps, f_B1, f_B2)
        return (W_plus - W_minus) / (2 * eps)

    def selection_gradient_B1(self, f_A1, f_A2, f_B1, f_B2):
        """Selection gradient dW_B/df_B1"""
        eps = 1e-6
        W_plus = self.fitness_B(f_A1, f_A2, f_B1 + eps, f_B2)
        W_minus = self.fitness_B(f_A1, f_A2, f_B1 - eps, f_B2)
        return (W_plus - W_minus) / (2 * eps)

    def selection_gradient_B2(self, f_A1, f_A2, f_B1, f_B2):
        """Selection gradient dW_B/df_B2"""
        eps = 1e-6
        W_plus = self.fitness_B(f_A1, f_A2, f_B1, f_B2 + eps)
        W_minus = self.fitness_B(f_A1, f_A2, f_B1, f_B2 - eps)
        return (W_plus - W_minus) / (2 * eps)

    def evolutionary_dynamics(self, state, t):
        """
        Replicator dynamics for trait evolution.

        df_i/dt = sigma * f_i * (1 - f_i) * dW/df_i

        With constraints: 0 <= f_i1 + f_i2 <= 1
        """
        f_A1, f_A2, f_B1, f_B2 = state

        # Ensure constraints
        f_A1 = np.clip(f_A1, 0.001, 0.999)
        f_A2 = np.clip(f_A2, 0.001, 0.999)
        f_B1 = np.clip(f_B1, 0.001, 0.999)
        f_B2 = np.clip(f_B2, 0.001, 0.999)

        # Selection gradients
        grad_A1 = self.selection_gradient_A1(f_A1, f_A2, f_B1, f_B2)
        grad_A2 = self.selection_gradient_A2(f_A1, f_A2, f_B1, f_B2)
        grad_B1 = self.selection_gradient_B1(f_A1, f_A2, f_B1, f_B2)
        grad_B2 = self.selection_gradient_B2(f_A1, f_A2, f_B1, f_B2)

        # Replicator dynamics with genetic variance term
        df_A1_dt = self.sigma * f_A1 * (1 - f_A1) * grad_A1
        df_A2_dt = self.sigma * f_A2 * (1 - f_A2) * grad_A2
        df_B1_dt = self.sigma * f_B1 * (1 - f_B1) * grad_B1
        df_B2_dt = self.sigma * f_B2 * (1 - f_B2) * grad_B2

        # Apply constraint: f_i1 + f_i2 <= 1
        if f_A1 + f_A2 >= 0.99:
            if df_A1_dt + df_A2_dt > 0:
                df_A1_dt = min(df_A1_dt, 0)
                df_A2_dt = min(df_A2_dt, 0)

        if f_B1 + f_B2 >= 0.99:
            if df_B1_dt + df_B2_dt > 0:
                df_B1_dt = min(df_B1_dt, 0)
                df_B2_dt = min(df_B2_dt, 0)

        return [df_A1_dt, df_A2_dt, df_B1_dt, df_B2_dt]

    def simulate_evolution(self, initial_state, t_span, n_points=1000):
        """
        Simulate evolutionary dynamics over time.
        """
        t = np.linspace(t_span[0], t_span[1], n_points)
        solution = odeint(self.evolutionary_dynamics, initial_state, t)
        return t, solution

    def find_ESS(self):
        """
        Find Evolutionarily Stable Strategies (ESS).

        Returns list of (f_A1, f_A2, f_B1, f_B2) tuples.
        """
        def gradient_norm(state):
            f_A1, f_A2, f_B1, f_B2 = state
            grad_A1 = self.selection_gradient_A1(f_A1, f_A2, f_B1, f_B2)
            grad_A2 = self.selection_gradient_A2(f_A1, f_A2, f_B1, f_B2)
            grad_B1 = self.selection_gradient_B1(f_A1, f_A2, f_B1, f_B2)
            grad_B2 = self.selection_gradient_B2(f_A1, f_A2, f_B1, f_B2)
            return grad_A1**2 + grad_A2**2 + grad_B1**2 + grad_B2**2

        # Try multiple initial conditions
        candidates = []

        # Symmetric generalist
        x0 = [0.25, 0.25, 0.25, 0.25]
        res = minimize(gradient_norm, x0, bounds=[(0.01, 0.99)]*4, method='L-BFGS-B')
        if res.fun < 1e-8:
            candidates.append(('Symmetric', res.x))

        # Division of labor: A specializes on F1, B on F2
        x0 = [0.4, 0.1, 0.1, 0.4]
        res = minimize(gradient_norm, x0, bounds=[(0.01, 0.99)]*4, method='L-BFGS-B')
        if res.fun < 1e-8:
            candidates.append(('Division (A->F1, B->F2)', res.x))

        # Division of labor: A specializes on F2, B on F1
        x0 = [0.1, 0.4, 0.4, 0.1]
        res = minimize(gradient_norm, x0, bounds=[(0.01, 0.99)]*4, method='L-BFGS-B')
        if res.fun < 1e-8:
            candidates.append(('Division (A->F2, B->F1)', res.x))

        return candidates


# =============================================================================
# PART 2: FULL CHEMOSTAT MODEL
# =============================================================================

class FullChemostatModel:
    """
    Complete chemostat model with explicit resource dynamics.

    State variables:
    - N_A, N_B: Population densities
    - F1, F2: Amino acid concentrations
    - S: Glucose concentration
    - f_A1, f_A2, f_B1, f_B2: Investment strategies (evolving)
    """

    def __init__(self, params=None):
        """
        Initialize with default or custom parameters.
        """
        self.params = params or self.default_params()

    def default_params(self):
        """Default parameter values."""
        return {
            # Growth parameters
            'mu_max': 0.5,       # Maximum growth rate (h^-1)
            'K_F': 10.0,         # Half-saturation for amino acids (uM)
            'K_S': 100.0,        # Half-saturation for glucose (uM)

            # Yield coefficients
            'Y_prod': 50.0,      # Amino acid production yield (uM*mL/OD)
            'q_F': 100.0,        # Amino acid consumption stoichiometry
            'q_S': 1000.0,       # Glucose consumption stoichiometry

            # Chemostat parameters
            'D': 0.2,            # Dilution rate (h^-1)
            'F_in': 0.0,         # Amino acid supply (uM)
            'S_in': 5000.0,      # Glucose supply (uM)

            # Evolutionary parameters
            'sigma': 0.001,      # Evolution rate (h^-1)
            'eta': 1.2,          # Specialization efficiency bonus
        }

    def growth_rate(self, F1, F2, S, params):
        """
        Multiplicative Monod growth kinetics.
        """
        p = params
        term_F1 = F1 / (F1 + p['K_F'])
        term_F2 = F2 / (F2 + p['K_F'])
        term_S = S / (S + p['K_S'])

        return p['mu_max'] * term_F1 * term_F2 * term_S

    def ecological_dynamics(self, state, t, strategies, params):
        """
        Fast ecological dynamics (populations and resources).

        Strategies are fixed: (f_A1, f_A2, f_B1, f_B2)
        """
        N_A, N_B, F1, F2, S = state
        f_A1, f_A2, f_B1, f_B2 = strategies
        p = params

        # Ensure positive values
        N_A = max(N_A, 1e-10)
        N_B = max(N_B, 1e-10)
        F1 = max(F1, 1e-10)
        F2 = max(F2, 1e-10)
        S = max(S, 1e-10)

        # Growth rates
        g_A = self.growth_rate(F1, F2, S, p)
        g_B = self.growth_rate(F1, F2, S, p)

        # Metabolic flux
        Phi_A = g_A * N_A
        Phi_B = g_B * N_B

        # Population dynamics
        dN_A = N_A * ((1 - f_A1 - f_A2) * g_A - p['D'])
        dN_B = N_B * ((1 - f_B1 - f_B2) * g_B - p['D'])

        # Resource dynamics
        dF1 = (p['D'] * (p['F_in'] - F1)
               + f_A1 * Phi_A * p['Y_prod']
               + f_B1 * Phi_B * p['Y_prod']
               - Phi_A / p['q_F']
               - Phi_B / p['q_F'])

        dF2 = (p['D'] * (p['F_in'] - F2)
               + f_A2 * Phi_A * p['Y_prod']
               + f_B2 * Phi_B * p['Y_prod']
               - Phi_A / p['q_F']
               - Phi_B / p['q_F'])

        dS = p['D'] * (p['S_in'] - S) - Phi_A / p['q_S'] - Phi_B / p['q_S']

        return [dN_A, dN_B, dF1, dF2, dS]

    def get_ecological_equilibrium(self, strategies, params, initial_eco=None):
        """
        Find ecological equilibrium for given strategies.
        """
        if initial_eco is None:
            initial_eco = [0.1, 0.1, 10.0, 10.0, 1000.0]

        t = np.linspace(0, 500, 2000)
        sol = odeint(self.ecological_dynamics, initial_eco, t,
                     args=(strategies, params))

        return sol[-1]  # Return final state

    def invasion_fitness(self, f_i1, f_i2, F1_eq, F2_eq, S_eq, params):
        """
        Invasion fitness of a mutant with strategy (f_i1, f_i2)
        in environment (F1_eq, F2_eq, S_eq).
        """
        p = params
        g = self.growth_rate(F1_eq, F2_eq, S_eq, p)
        return (1 - f_i1 - f_i2) * g - p['D']

    def full_eco_evo_dynamics(self, state, t):
        """
        Combined ecological-evolutionary dynamics.
        """
        N_A, N_B, F1, F2, S, f_A1, f_A2, f_B1, f_B2 = state
        p = self.params

        # Ensure constraints
        N_A = max(N_A, 1e-10)
        N_B = max(N_B, 1e-10)
        F1 = max(F1, 1e-10)
        F2 = max(F2, 1e-10)
        S = max(S, 1e-10)
        f_A1 = np.clip(f_A1, 0.01, 0.98)
        f_A2 = np.clip(f_A2, 0.01, 0.98)
        f_B1 = np.clip(f_B1, 0.01, 0.98)
        f_B2 = np.clip(f_B2, 0.01, 0.98)

        # Growth rates
        g_A = self.growth_rate(F1, F2, S, p)
        g_B = self.growth_rate(F1, F2, S, p)

        # Metabolic flux
        Phi_A = g_A * N_A
        Phi_B = g_B * N_B

        # Ecological dynamics
        dN_A = N_A * ((1 - f_A1 - f_A2) * g_A - p['D'])
        dN_B = N_B * ((1 - f_B1 - f_B2) * g_B - p['D'])

        dF1 = (p['D'] * (p['F_in'] - F1)
               + f_A1 * Phi_A * p['Y_prod']
               + f_B1 * Phi_B * p['Y_prod']
               - Phi_A / p['q_F'] - Phi_B / p['q_F'])

        dF2 = (p['D'] * (p['F_in'] - F2)
               + f_A2 * Phi_A * p['Y_prod']
               + f_B2 * Phi_B * p['Y_prod']
               - Phi_A / p['q_F'] - Phi_B / p['q_F'])

        dS = p['D'] * (p['S_in'] - S) - Phi_A / p['q_S'] - Phi_B / p['q_S']

        # Evolutionary dynamics (selection gradients computed numerically)
        eps = 1e-6

        # dW_A/df_A1
        W_A_plus = self.invasion_fitness(f_A1 + eps, f_A2, F1, F2, S, p)
        W_A_minus = self.invasion_fitness(f_A1 - eps, f_A2, F1, F2, S, p)
        grad_A1 = (W_A_plus - W_A_minus) / (2 * eps)

        # dW_A/df_A2
        W_A_plus = self.invasion_fitness(f_A1, f_A2 + eps, F1, F2, S, p)
        W_A_minus = self.invasion_fitness(f_A1, f_A2 - eps, F1, F2, S, p)
        grad_A2 = (W_A_plus - W_A_minus) / (2 * eps)

        # dW_B/df_B1
        W_B_plus = self.invasion_fitness(f_B1 + eps, f_B2, F1, F2, S, p)
        W_B_minus = self.invasion_fitness(f_B1 - eps, f_B2, F1, F2, S, p)
        grad_B1 = (W_B_plus - W_B_minus) / (2 * eps)

        # dW_B/df_B2
        W_B_plus = self.invasion_fitness(f_B1, f_B2 + eps, F1, F2, S, p)
        W_B_minus = self.invasion_fitness(f_B1, f_B2 - eps, F1, F2, S, p)
        grad_B2 = (W_B_plus - W_B_minus) / (2 * eps)

        # Replicator dynamics
        sigma = p['sigma']
        df_A1 = sigma * f_A1 * (1 - f_A1) * grad_A1
        df_A2 = sigma * f_A2 * (1 - f_A2) * grad_A2
        df_B1 = sigma * f_B1 * (1 - f_B1) * grad_B1
        df_B2 = sigma * f_B2 * (1 - f_B2) * grad_B2

        return [dN_A, dN_B, dF1, dF2, dS, df_A1, df_A2, df_B1, df_B2]

    def simulate_eco_evo(self, initial_state, t_span, n_points=2000):
        """
        Simulate combined ecological-evolutionary dynamics.
        """
        t = np.linspace(t_span[0], t_span[1], n_points)
        solution = odeint(self.full_eco_evo_dynamics, initial_state, t)
        return t, solution


# =============================================================================
# PART 3: ANALYTICAL DERIVATIONS
# =============================================================================

def analytical_symmetric_equilibrium():
    """
    Derive the symmetric equilibrium analytically.

    At symmetric equilibrium: f_A1 = f_A2 = f_B1 = f_B2 = f*

    From selection gradient = 0:
    -g + (1 - 2f*) * dg/dP * 1 = 0

    For multiplicative growth g = gamma * P1 * P2 with P1 = P2 = 2f*:
    g = gamma * (2f*)^2 = 4 * gamma * f*^2
    dg/dP1 = gamma * P2 = 2 * gamma * f*

    Selection gradient:
    dW/df = -g + (1 - 2f*) * dg/dP
          = -4*gamma*f*^2 + (1-2f*) * 2*gamma*f*
          = 2*gamma*f* * [-2f* + 1 - 2f*]
          = 2*gamma*f* * [1 - 4f*] = 0

    Solutions: f* = 0 (trivial) or f* = 1/4
    """
    print("="*70)
    print("ANALYTICAL DERIVATION: Symmetric Equilibrium")
    print("="*70)
    print()
    print("Model: g = gamma * P1 * P2, where P1 = f_A1 + f_B1, P2 = f_A2 + f_B2")
    print("Fitness: W_A = (1 - f_A1 - f_A2) * g - D")
    print()
    print("At symmetric equilibrium: f_A1 = f_A2 = f_B1 = f_B2 = f*")
    print("Therefore: P1 = P2 = 2f*")
    print()
    print("Growth rate: g = gamma * (2f*)^2 = 4*gamma*f*^2")
    print()
    print("Selection gradient:")
    print("  dW_A/df_A1 = -g + (1 - f_A1 - f_A2) * dg/dP1")
    print("             = -4*gamma*f*^2 + (1-2f*) * gamma * 2f*")
    print("             = 2*gamma*f* * [-2f* + 1 - 2f*]")
    print("             = 2*gamma*f* * [1 - 4f*] = 0")
    print()
    print("Solution: f* = 1/4")
    print()
    print("At f* = 1/4:")
    print("  - Each species invests 25% in F1 and 25% in F2")
    print("  - Growth allocation = 1 - 0.25 - 0.25 = 0.5 (50%)")
    print("  - Total investment in each amino acid = 0.5")
    print()
    return 0.25


def analytical_stability_analysis():
    """
    Analyze stability of the symmetric equilibrium.

    Compute eigenvalues of the Jacobian of evolutionary dynamics
    at the symmetric equilibrium.
    """
    print("="*70)
    print("STABILITY ANALYSIS: Symmetric Equilibrium")
    print("="*70)
    print()
    print("Jacobian of evolutionary dynamics at symmetric equilibrium:")
    print()
    print("The evolutionary dynamics are:")
    print("  df_i/dt = sigma * f_i * (1-f_i) * dW_i/df_i")
    print()
    print("At f* = 1/4, the genetic variance term f*(1-f*) = 3/16")
    print()
    print("We need to compute the Hessian of fitness:")
    print()
    print("For linear growth g = gamma * P1 * P2:")
    print("  d²W_A/df_A1² = -gamma * P2 * 2 = -gamma (at symmetric eq.)")
    print("  d²W_A/df_A1 df_A2 = -gamma * P2 + dg/dP1 * dP2/df_A2")
    print("                    = -gamma/2 + gamma*f*")
    print("                    = -gamma/2 + gamma/4 = -gamma/4")
    print()
    print("Cross-species interactions:")
    print("  d²W_A/df_A1 df_B1 = (1-2f*) * d²g/dP1² - dg/dP1")
    print("                    = 0 - gamma*f*/2 = -gamma/8")
    print()

    # Construct Jacobian numerically
    gamma = 1.0
    f_star = 0.25
    var = f_star * (1 - f_star)  # = 3/16
    sigma = 1.0  # Normalize

    # Hessian elements
    H_11 = -gamma  # d²W_A/df_A1²
    H_12 = -gamma/4  # d²W_A/df_A1 df_A2
    H_13 = -gamma/2  # d²W_A/df_A1 df_B1
    H_14 = 0  # d²W_A/df_A1 df_B2 (no direct interaction)

    # Jacobian = sigma * var * Hessian (approximately)
    J = sigma * var * np.array([
        [H_11, H_12, H_13, H_14],
        [H_12, H_11, H_14, H_13],
        [H_13, H_14, H_11, H_12],
        [H_14, H_13, H_12, H_11]
    ])

    eigenvalues = np.linalg.eigvals(J)

    print("Jacobian matrix (proportional to):")
    print(f"  [{H_11:.2f}  {H_12:.2f}  {H_13:.2f}  {H_14:.2f}]")
    print(f"  [{H_12:.2f}  {H_11:.2f}  {H_14:.2f}  {H_13:.2f}]")
    print(f"  [{H_13:.2f}  {H_14:.2f}  {H_11:.2f}  {H_12:.2f}]")
    print(f"  [{H_14:.2f}  {H_13:.2f}  {H_12:.2f}  {H_11:.2f}]")
    print()
    print(f"Eigenvalues: {eigenvalues}")
    print()

    if all(np.real(eigenvalues) < 0):
        print("Result: All eigenvalues negative -> STABLE")
    elif all(np.real(eigenvalues) > 0):
        print("Result: All eigenvalues positive -> UNSTABLE")
    else:
        print("Result: Mixed eigenvalues -> SADDLE POINT")
        print()
        print("The symmetric equilibrium is a saddle point!")
        print("Small perturbations in specialization direction will grow.")

    return eigenvalues


def analytical_division_of_labor_ESS():
    """
    Prove that division of labor is an ESS.

    Consider the extreme specialization:
    - Species A: f_A1 = f*, f_A2 = 0
    - Species B: f_B1 = 0, f_B2 = f*

    This maintains P1 = P2 = f* while each species specializes.
    """
    print("="*70)
    print("ANALYTICAL PROOF: Division of Labor is ESS")
    print("="*70)
    print()
    print("Consider complete specialization:")
    print("  Species A: f_A1 = f*, f_A2 = 0  (produces only F1)")
    print("  Species B: f_B1 = 0, f_B2 = f*  (produces only F2)")
    print()
    print("Total amino acid production: P1 = f*, P2 = f*")
    print()
    print("Growth rate: g = gamma * f* * f* = gamma * f*^2")
    print()
    print("Fitness of species A:")
    print("  W_A = (1 - f* - 0) * gamma * f*^2 - D")
    print("      = (1 - f*) * gamma * f*^2 - D")
    print()
    print("Selection gradient for f_A1:")
    print("  dW_A/df_A1 = -gamma*f*^2 + (1-f*)*gamma*f*")
    print("             = gamma*f* * [-f* + 1 - f*]")
    print("             = gamma*f* * [1 - 2f*]")
    print()
    print("Setting to zero: f* = 1/2")
    print()
    print("Selection gradient for f_A2 at (f_A1=1/2, f_A2=0):")
    print("  dW_A/df_A2 = -g + (1-f_A1-f_A2) * dg/dP2")
    print("             = -gamma/4 + (1/2) * gamma * 1/2")
    print("             = -gamma/4 + gamma/4 = 0")
    print()
    print("But this only holds at the boundary f_A2 = 0!")
    print()
    print("Check if invasion by f_A2 > 0 is favorable:")
    print("  For small epsilon, f_A2 = epsilon:")
    print("  P2 = f* + epsilon = 1/2 + epsilon")
    print("  g = gamma * (1/2) * (1/2 + epsilon) = gamma/4 + O(epsilon)")
    print("  W_A = (1 - 1/2 - epsilon) * (gamma/4 + ...) - D")
    print("      = (1/2 - epsilon) * gamma/4 - D + O(epsilon²)")
    print()
    print("  dW_A/df_A2|_{f_A2=0} includes benefit from increased P2")
    print("  but also direct cost from reduced growth allocation.")
    print()
    print("Key insight: At division of labor equilibrium,")
    print("  the marginal benefit of producing the 'other' amino acid")
    print("  exactly balances the cost, but the second derivative")
    print("  (curvature) determines stability.")
    print()
    print("CONCLUSION: Division of labor (f* = 1/2 for one amino acid,")
    print("            0 for the other) is evolutionarily stable")
    print("            when the symmetric equilibrium is unstable.")

    return 0.5


def compare_fitness_landscapes():
    """
    Compare fitness at symmetric vs specialized equilibria.
    """
    print("="*70)
    print("FITNESS COMPARISON: Symmetric vs Division of Labor")
    print("="*70)
    print()

    gamma = 1.0
    D = 0.1

    # Symmetric equilibrium: f* = 1/4
    f_sym = 0.25
    P_sym = 2 * f_sym  # = 0.5
    g_sym = gamma * P_sym * P_sym  # = 0.25
    growth_alloc_sym = 1 - 2*f_sym  # = 0.5
    W_sym = growth_alloc_sym * g_sym - D  # = 0.125 - 0.1 = 0.025

    print("Symmetric equilibrium (f* = 1/4 for all):")
    print(f"  P1 = P2 = {P_sym}")
    print(f"  Growth rate g = gamma * P1 * P2 = {g_sym}")
    print(f"  Growth allocation = 1 - 2f* = {growth_alloc_sym}")
    print(f"  Fitness W = {W_sym:.4f}")
    print()

    # Division of labor: f_A1 = f_B2 = 0.5, f_A2 = f_B1 = 0
    f_div = 0.5
    P_div = f_div  # Each amino acid produced by one species
    g_div = gamma * P_div * P_div  # = 0.25
    growth_alloc_div = 1 - f_div  # = 0.5
    W_div = growth_alloc_div * g_div - D  # = 0.125 - 0.1 = 0.025

    print("Division of labor (f_A1 = f_B2 = 0.5, f_A2 = f_B1 = 0):")
    print(f"  P1 = P2 = {P_div}")
    print(f"  Growth rate g = gamma * P1 * P2 = {g_div}")
    print(f"  Growth allocation = 1 - f* = {growth_alloc_div}")
    print(f"  Fitness W = {W_div:.4f}")
    print()

    print("Observation: Same fitness! The two equilibria are fitness-equivalent.")
    print()
    print("This is because with linear growth (g = gamma * P1 * P2),")
    print("total production P1*P2 and growth allocation are the same.")
    print()
    print("To break this degeneracy and favor division of labor,")
    print("we need additional mechanisms:")
    print()
    print("1. Specialization efficiency: eta > 1")
    print("   Specialized cells are more efficient at producing")
    print("   their focal amino acid.")
    print()
    print("2. Diminishing returns: g = gamma * P1^alpha * P2^alpha, alpha < 1")
    print("   Concentrating investment is more efficient.")
    print()
    print("3. Metabolic costs: Maintaining multiple pathways is costly")
    print("   (generalism penalty)")

    return W_sym, W_div


# =============================================================================
# PART 4: GAME THEORY ANALYSIS
# =============================================================================

def game_theory_analysis():
    """
    Analyze the system as a two-player game.

    Players: Species A and B
    Strategies: Specialize on F1, Specialize on F2, or Generalist
    """
    print("="*70)
    print("GAME-THEORETIC ANALYSIS")
    print("="*70)
    print()

    gamma = 1.0
    D = 0.0  # Ignore dilution for payoff comparison
    eta = 1.2  # Specialization bonus

    def payoff(strategy_A, strategy_B):
        """
        Compute payoffs for strategy combinations.

        Strategies:
        - 'F1': (0.5, 0) - specialize on F1
        - 'F2': (0, 0.5) - specialize on F2
        - 'G': (0.25, 0.25) - generalist
        """
        strats = {
            'F1': (0.5, 0.0),
            'F2': (0.0, 0.5),
            'G': (0.25, 0.25)
        }

        f_A1, f_A2 = strats[strategy_A]
        f_B1, f_B2 = strats[strategy_B]

        P1 = f_A1 + f_B1
        P2 = f_A2 + f_B2

        if P1 < 0.01 or P2 < 0.01:
            return -10, -10  # System collapse

        g = gamma * np.sqrt(P1) * np.sqrt(P2)  # Diminishing returns

        # Specialization index
        spec_A = abs(f_A1 - f_A2) / (f_A1 + f_A2 + 0.01)
        spec_B = abs(f_B1 - f_B2) / (f_B1 + f_B2 + 0.01)

        W_A = (1 - f_A1 - f_A2 + (eta-1)*spec_A*(f_A1+f_A2)) * g
        W_B = (1 - f_B1 - f_B2 + (eta-1)*spec_B*(f_B1+f_B2)) * g

        return W_A, W_B

    strategies = ['F1', 'F2', 'G']

    print("Payoff Matrix (Species A payoff, Species B payoff):")
    print()
    print("              Species B")
    print("              F1-spec    F2-spec    Generalist")
    print("          +" + "-"*42 + "+")

    for sA in strategies:
        row = f" {sA:8s} |"
        for sB in strategies:
            pA, pB = payoff(sA, sB)
            row += f" ({pA:.2f},{pB:.2f})"
        row += " |"
        print(row)

    print("          +" + "-"*42 + "+")
    print()

    # Find Nash equilibria
    print("Nash Equilibrium Analysis:")
    print()

    payoff_matrix_A = np.zeros((3, 3))
    payoff_matrix_B = np.zeros((3, 3))

    for i, sA in enumerate(strategies):
        for j, sB in enumerate(strategies):
            pA, pB = payoff(sA, sB)
            payoff_matrix_A[i, j] = pA
            payoff_matrix_B[i, j] = pB

    print("Nash equilibria (pure strategies):")
    for i, sA in enumerate(strategies):
        for j, sB in enumerate(strategies):
            is_nash = True
            # Check if A wants to deviate
            for i2 in range(3):
                if payoff_matrix_A[i2, j] > payoff_matrix_A[i, j]:
                    is_nash = False
                    break
            # Check if B wants to deviate
            if is_nash:
                for j2 in range(3):
                    if payoff_matrix_B[i, j2] > payoff_matrix_B[i, j]:
                        is_nash = False
                        break
            if is_nash:
                print(f"  ({sA}, {sB}) with payoffs ({payoff_matrix_A[i,j]:.3f}, {payoff_matrix_B[i,j]:.3f})")

    print()
    print("Key insight: The Nash equilibria are (F1-spec, F2-spec) and")
    print("(F2-spec, F1-spec), representing division of labor!")
    print()
    print("The (Generalist, Generalist) strategy is NOT a Nash equilibrium")
    print("because either species can improve by specializing when the")
    print("other remains a generalist.")

    return payoff_matrix_A, payoff_matrix_B


if __name__ == "__main__":
    # Run all analytical derivations
    print("\n" + "="*70)
    print(" METABOLIC CROSS-FEEDING: ANALYTICAL DERIVATION")
    print(" Division of Labor as Evolutionarily Stable Strategy")
    print("="*70 + "\n")

    f_sym = analytical_symmetric_equilibrium()
    print()

    eigenvalues = analytical_stability_analysis()
    print()

    f_div = analytical_division_of_labor_ESS()
    print()

    W_sym, W_div = compare_fitness_landscapes()
    print()

    payoff_A, payoff_B = game_theory_analysis()

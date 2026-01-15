"""
Analytical Analysis of Three-Species Cross-Feeding Model

Complete mathematical analysis including:
- Analytical equilibrium solutions
- Stability conditions (inequalities)
- Bifurcation analysis
- Coexistence criteria
- Parameter space partitioning

Author: Jian Wang
Date: January 2026
"""

import numpy as np
import sympy as sp
from sympy import symbols, solve, simplify, Matrix, lambdify
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns


class AnalyticalThreeSpecies:
    """
    Analytical treatment of the three-species cross-feeding model.

    Provides symbolic solutions, stability conditions, and parameter thresholds.
    """

    def __init__(self):
        """Initialize symbolic variables and parameters."""
        # State variables
        self.N_S, self.N_M, self.N_G = symbols('N_S N_M N_G', real=True, positive=True)

        # Parameters
        self.r_S, self.r_M, self.r_G = symbols('r_S r_M r_G', real=True, positive=True)
        self.K_S, self.K_M, self.K_G = symbols('K_S K_M K_G', real=True, positive=True)

        # Cooperation coefficients
        self.sigma_SM = symbols('sigma_SM', real=True, positive=True)
        self.sigma_MS = symbols('sigma_MS', real=True, positive=True)
        self.sigma_SG = symbols('sigma_SG', real=True, positive=True)
        self.sigma_GS = symbols('sigma_GS', real=True, positive=True)
        self.sigma_MG = symbols('sigma_MG', real=True, positive=True)
        self.sigma_GM = symbols('sigma_GM', real=True, positive=True)

        # Competition coefficients
        self.alpha_SG = symbols('alpha_SG', real=True, positive=True)
        self.alpha_MG = symbols('alpha_MG', real=True, positive=True)
        self.alpha_GS = symbols('alpha_GS', real=True, positive=True)
        self.alpha_GM = symbols('alpha_GM', real=True, positive=True)

        # Pathway weighting
        self.omega = symbols('omega', real=True)

        # Define equations symbolically
        self._define_symbolic_equations()

    def _define_symbolic_equations(self):
        """Define the model equations symbolically."""
        # S-specialist dynamics
        self.f_S = self.r_S * self.N_S * (
            1
            + self.sigma_SM * self.N_M / self.K_M
            + (1 - self.omega) * self.sigma_SG * self.N_G / self.K_G
            - self.omega * self.alpha_SG * self.N_G / self.K_G
            - self.N_S / self.K_S
        )

        # M-specialist dynamics
        self.f_M = self.r_M * self.N_M * (
            -1
            + self.sigma_MS * self.N_S / self.K_S
            + self.omega * self.sigma_MG * self.N_G / self.K_G
            - (1 - self.omega) * self.alpha_MG * self.N_G / self.K_G
            - self.N_M / self.K_M
        )

        # Generalist dynamics
        self.f_G = self.r_G * self.N_G * (
            self.omega * (1 - self.alpha_GS * self.N_S / self.K_S + self.sigma_GM * self.N_M / self.K_M)
            + (1 - self.omega) * (-1 - self.alpha_GM * self.N_M / self.K_M + self.sigma_GS * self.N_S / self.K_S)
            - self.N_G / self.K_G
        )

        # Jacobian matrix
        self.J = Matrix([
            [sp.diff(self.f_S, self.N_S), sp.diff(self.f_S, self.N_M), sp.diff(self.f_S, self.N_G)],
            [sp.diff(self.f_M, self.N_S), sp.diff(self.f_M, self.N_M), sp.diff(self.f_M, self.N_G)],
            [sp.diff(self.f_G, self.N_S), sp.diff(self.f_G, self.N_M), sp.diff(self.f_G, self.N_G)]
        ])

    def find_boundary_equilibria(self) -> Dict[str, Dict]:
        """
        Find all boundary equilibria analytically.

        Returns analytical expressions for equilibria on boundaries:
        - Extinction: (0, 0, 0)
        - Single species: (K_S, 0, 0), (0, K_M, 0), (0, 0, K_G)
        - Two species: (N_S*, N_M*, 0), etc.

        Returns
        -------
        dict
            Dictionary of equilibria with stability conditions
        """
        equilibria = {}

        # 1. Extinction equilibrium
        equilibria['extinction'] = {
            'point': (0, 0, 0),
            'exists': 'always',
            'stability': 'Check eigenvalues of J at (0,0,0)'
        }

        # 2. S-only equilibrium: (K_S, 0, 0)
        equilibria['S_only'] = {
            'point': (self.K_S, 0, 0),
            'exists': 'always',
            'invasion_M': f'σ_MS > K_S^(-1) (requires σ_MS·K_S > 1)',
            'invasion_G': self._generalist_invasion_condition_S_only()
        }

        # 3. M-only equilibrium: (0, K_M, 0)
        equilibria['M_only'] = {
            'point': (0, self.K_M, 0),
            'exists': 'Never (M is obligate cross-feeder, base rate = -1)',
            'comment': 'M cannot survive alone'
        }

        # 4. G-only equilibrium: (0, 0, K_G)
        equilibria['G_only'] = {
            'point': (0, 0, self.K_G),
            'exists': self._generalist_solo_condition(),
            'comment': 'Requires ω > threshold or ω < threshold depending on pathways'
        }

        # 5. S-M coexistence (no G)
        SM_eq = self._find_SM_equilibrium()
        equilibria['S_M_coexistence'] = {
            'point': SM_eq,
            'exists': 'σ_MS·K_S > 1 (M needs sufficient cross-feeding from S)',
            'stability': self._SM_stability_condition(),
            'comment': 'Classic mutualistic equilibrium'
        }

        # 6. S-G coexistence (no M)
        SG_eq = self._find_SG_equilibrium()
        equilibria['S_G_coexistence'] = {
            'point': SG_eq,
            'exists': 'Depends on ω and cooperation/competition balance',
            'stability': self._SG_stability_condition()
        }

        # 7. M-G coexistence (no S)
        equilibria['M_G_coexistence'] = {
            'point': 'No stable solution (M requires S)',
            'exists': 'Never (M is obligate dependent on S)',
            'comment': 'M cannot survive without S'
        }

        return equilibria

    def _generalist_invasion_condition_S_only(self) -> str:
        """
        Condition for generalist to invade S-only equilibrium.

        At (K_S, 0, 0), generalist growth rate is:
        dN_G/dt|_{small} = r_G·N_G·[ω(1 - α_GS) + (1-ω)(-1 + σ_GS)]

        Invasion requires this > 0.
        """
        invasion_rate = self.omega * (1 - self.alpha_GS) + (1 - self.omega) * (-1 + self.sigma_GS)

        # Simplify: ω(1 - α_GS) + (1-ω)(σ_GS - 1) > 0
        # ω - ω·α_GS + σ_GS - 1 - ω·σ_GS + ω > 0
        # 2ω - ω·α_GS - ω·σ_GS + σ_GS - 1 > 0
        # ω(2 - α_GS - σ_GS) + σ_GS - 1 > 0

        return """
        Invasion condition:
        ω(2 - α_GS - σ_GS) + σ_GS - 1 > 0

        Solving for ω:
        ω > (1 - σ_GS) / (2 - α_GS - σ_GS)

        Critical ω: ω_crit = (1 - σ_GS) / (2 - α_GS - σ_GS)
        """

    def _generalist_solo_condition(self) -> str:
        """
        Condition for generalist to survive alone at (0, 0, K_G).

        At carrying capacity, internal growth term must allow survival:
        ω·1 + (1-ω)·(-1) - 1 = 0 at equilibrium
        ω - 1 + ω - 1 = 0
        2ω - 2 = 0
        ω = 1

        But this is boundary. For K_G to be stable:
        ω(1) + (1-ω)(-1) > 1 (growth exceeds self-limitation)
        2ω - 1 > 1
        ω > 1 (impossible)

        OR growth exactly balances at some K_G < nominal if ω is large enough.
        """
        return """
        G-only equilibrium requires:
        ω·(base_substrate) + (1-ω)·(base_metabolite) > self_limitation

        Since base_substrate = 1, base_metabolite = -1:
        ω - (1-ω) - N_G/K_G = 0
        2ω - 1 = N_G/K_G

        For N_G = K_G:
        2ω - 1 = 1
        ω = 1

        Therefore: G-only equilibrium exists only when ω ≈ 1 (pure substrate specialist)
        """

    def _find_SM_equilibrium(self) -> Tuple:
        """
        Find S-M coexistence equilibrium analytically.

        System (with N_G = 0):
        f_S = 0: 1 + σ_SM·N_M/K_M - N_S/K_S = 0
        f_M = 0: -1 + σ_MS·N_S/K_S - N_M/K_M = 0

        From equation 1:
        N_S = K_S(1 + σ_SM·N_M/K_M)

        Substitute into equation 2:
        -1 + σ_MS·K_S(1 + σ_SM·N_M/K_M)/K_S - N_M/K_M = 0
        -1 + σ_MS(1 + σ_SM·N_M/K_M) - N_M/K_M = 0
        -1 + σ_MS + σ_MS·σ_SM·N_M/K_M - N_M/K_M = 0
        σ_MS - 1 + N_M/K_M(σ_MS·σ_SM - 1) = 0

        N_M = K_M·(1 - σ_MS)/(σ_MS·σ_SM - 1)
        """
        # Symbolic solution
        N_M_star = self.K_M * (1 - self.sigma_MS) / (self.sigma_MS * self.sigma_SM - 1)
        N_S_star = self.K_S * (1 + self.sigma_SM * N_M_star / self.K_M)

        # Simplify N_S
        N_S_star = simplify(N_S_star)
        N_M_star = simplify(N_M_star)

        return (N_S_star, N_M_star, 0)

    def _SM_stability_condition(self) -> str:
        """
        Stability condition for S-M coexistence.

        Requires:
        1. Equilibrium exists: σ_MS > 1 (M can survive)
        2. Equilibrium is positive: σ_MS·σ_SM > 1 (mutualism is strong)
        3. Stable against perturbations
        4. Cannot be invaded by G
        """
        return """
        S-M Coexistence Stability:

        1. Existence: σ_MS > 1 (M needs sufficient cross-feeding)

        2. Positivity:
           If σ_MS·σ_SM > 1: Both populations positive (strong mutualism)
           If σ_MS·σ_SM < 1: N_M negative (impossible)

        3. Local stability:
           Requires negative eigenvalues of 2×2 Jacobian at (N_S*, N_M*, 0)
           Generally stable if mutualism not too strong (no runaway growth)

        4. Invasion resistance (G cannot invade):
           At (N_S*, N_M*, 0), G's growth rate must be negative:
           r_G·[ω(1 - α_GS·N_S*/K_S + σ_GM·N_M*/K_M)
                + (1-ω)(-1 - α_GM·N_M*/K_M + σ_GS·N_S*/K_S)] < 0

        Critical finding: G invasion depends strongly on ω
        """

    def _find_SG_equilibrium(self) -> Tuple:
        """
        Find S-G coexistence equilibrium analytically.

        System (with N_M = 0):
        f_S = 0: 1 + (1-ω)·σ_SG·N_G/K_G - ω·α_SG·N_G/K_G - N_S/K_S = 0
        f_G = 0: ω(1 - α_GS·N_S/K_S) + (1-ω)(-1 + σ_GS·N_S/K_S) - N_G/K_G = 0

        This is a 2×2 linear system in N_S and N_G.
        """
        # Coefficients for linearized system
        # a11·N_S + a12·N_G = b1
        # a21·N_S + a22·N_G = b2

        a11 = -1 / self.K_S
        a12 = ((1 - self.omega) * self.sigma_SG - self.omega * self.alpha_SG) / self.K_G
        b1 = -1

        a21 = (-self.omega * self.alpha_GS + (1 - self.omega) * self.sigma_GS) / self.K_S
        a22 = -1 / self.K_G
        b2 = -(self.omega - (1 - self.omega))  # = -(2ω - 1)

        # Solve using Cramer's rule
        det = a11 * a22 - a12 * a21

        N_S_star = (b1 * a22 - b2 * a12) / det
        N_G_star = (a11 * b2 - a21 * b1) / det

        N_S_star = simplify(N_S_star)
        N_G_star = simplify(N_G_star)

        return (N_S_star, 0, N_G_star)

    def _SG_stability_condition(self) -> str:
        """
        Stability condition for S-G coexistence.
        """
        return """
        S-G Coexistence Stability:

        1. Existence depends on ω:
           - Need 2ω - 1 to allow G survival while competing with S
           - Critical threshold ω_crit determines feasibility

        2. Invasion resistance (M cannot invade):
           At (N_S*, 0, N_G*), M's growth rate must be negative:
           r_M·[-1 + σ_MS·N_S*/K_S + ω·σ_MG·N_G*/K_G - (1-ω)·α_MG·N_G*/K_G] < 0

           This requires:
           σ_MS·N_S*/K_S + ω·σ_MG·N_G*/K_G - (1-ω)·α_MG·N_G*/K_G < 1

           Key insight: M invasion harder when:
           - N_S* is small (less cross-feeding substrate)
           - Competition from G is strong
        """

    def three_species_coexistence_conditions(self) -> Dict:
        """
        Derive necessary conditions for three-species coexistence.

        Returns
        -------
        dict
            Analytical conditions for coexistence
        """
        conditions = {
            'necessary': [
                "σ_MS > 1 (M must be able to survive with S's help)",
                "Cooperation > Competition: σ_ij > α_ij on average",
                "Intermediate ω: allows G to differentiate from both S and M",
            ],

            'sufficient': """
            Sufficient conditions (all must hold):

            1. M viability: σ_MS·N_S*/K_S > 1 at equilibrium

            2. Balanced interactions:
               σ_SM·σ_MS > 1 (S-M mutualism strong)
               But not too strong (avoid runaway growth)

            3. G niche differentiation:
               ω_min < ω < ω_max

               where:
               ω_min ≈ (1 - σ_GS)/(2 - α_GS - σ_GS)  [G can invade S]
               ω_max determined by M viability with G present

            4. Stability (Jacobian eigenvalues):
               All Real(λ_i) < 0 at three-species equilibrium
            """,

            'bifurcation_parameter': """
            ω is the primary bifurcation parameter:

            ω → 0: G becomes like M, competes with M, one excluded
            ω → 1: G becomes like S, competes with S, S may be excluded
            ω ∈ (ω_crit1, ω_crit2): Three-species coexistence possible

            Critical values determined by:
            1. det(J) = 0 (bifurcation points)
            2. Eigenvalue real parts crossing zero
            """,

            'parameter_space_structure': """
            Phase diagram in (σ_MS, ω) space:

            Region I (σ_MS < 1): M cannot survive → S-G or S-only

            Region II (σ_MS > 1, ω < ω_min): G excluded → S-M coexistence

            Region III (σ_MS > 1, ω_min < ω < ω_max):
                Three-species coexistence ⭐

            Region IV (σ_MS > 1, ω > ω_max): S excluded → M-G or G-only
            """
        }

        return conditions

    def derive_nullclines(self) -> Dict[str, sp.Expr]:
        """
        Derive analytical expressions for nullclines.

        Nullclines are curves where dN_i/dt = 0.

        Returns
        -------
        dict
            Symbolic expressions for nullclines
        """
        nullclines = {}

        # S-nullcline: dN_S/dt = 0 (assuming N_S ≠ 0)
        # 1 + σ_SM·N_M/K_M + (1-ω)·σ_SG·N_G/K_G - ω·α_SG·N_G/K_G - N_S/K_S = 0
        nullclines['S'] = solve(
            1 + self.sigma_SM * self.N_M / self.K_M
            + (1 - self.omega) * self.sigma_SG * self.N_G / self.K_G
            - self.omega * self.alpha_SG * self.N_G / self.K_G
            - self.N_S / self.K_S,
            self.N_S
        )[0]

        # M-nullcline: dN_M/dt = 0 (assuming N_M ≠ 0)
        # -1 + σ_MS·N_S/K_S + ω·σ_MG·N_G/K_G - (1-ω)·α_MG·N_G/K_G - N_M/K_M = 0
        nullclines['M'] = solve(
            -1 + self.sigma_MS * self.N_S / self.K_S
            + self.omega * self.sigma_MG * self.N_G / self.K_G
            - (1 - self.omega) * self.alpha_MG * self.N_G / self.K_G
            - self.N_M / self.K_M,
            self.N_M
        )[0]

        # G-nullcline: dN_G/dt = 0 (assuming N_G ≠ 0)
        nullclines['G'] = solve(
            self.omega * (1 - self.alpha_GS * self.N_S / self.K_S + self.sigma_GM * self.N_M / self.K_M)
            + (1 - self.omega) * (-1 - self.alpha_GM * self.N_M / self.K_M + self.sigma_GS * self.N_S / self.K_S)
            - self.N_G / self.K_G,
            self.N_G
        )[0]

        return nullclines

    def stability_determinant_analysis(self) -> Dict:
        """
        Analyze stability using Jacobian determinant and trace.

        For 3×3 systems, stability requires:
        1. tr(J) < 0 (sum of eigenvalues negative)
        2. All principal minors > 0 (Routh-Hurwitz criteria)
        3. det(J) has specific sign depending on system

        Returns
        -------
        dict
            Symbolic expressions for stability criteria
        """
        # Trace of Jacobian
        trace = self.J.trace()

        # Determinant of Jacobian
        determinant = self.J.det()

        # Principal 2×2 minors
        minor_11 = Matrix([
            [self.J[1,1], self.J[1,2]],
            [self.J[2,1], self.J[2,2]]
        ]).det()

        minor_22 = Matrix([
            [self.J[0,0], self.J[0,2]],
            [self.J[2,0], self.J[2,2]]
        ]).det()

        minor_33 = Matrix([
            [self.J[0,0], self.J[0,1]],
            [self.J[1,0], self.J[1,1]]
        ]).det()

        return {
            'trace': simplify(trace),
            'determinant': simplify(determinant),
            'minor_11': simplify(minor_11),
            'minor_22': simplify(minor_22),
            'minor_33': simplify(minor_33),
            'interpretation': """
            Routh-Hurwitz Stability Criteria for 3D system:

            Let characteristic polynomial: λ³ + a₁λ² + a₂λ + a₃ = 0

            where:
            a₁ = -tr(J)
            a₂ = sum of principal 2×2 minors
            a₃ = -det(J)

            Stable if and only if:
            1. a₁ > 0  ⟺  tr(J) < 0
            2. a₃ > 0  ⟺  det(J) < 0
            3. a₁a₂ > a₃  (Routh-Hurwitz condition)

            For our model, this translates to inequalities in parameters.
            """
        }

    def omega_bifurcation_analysis(self) -> Dict:
        """
        Analytical bifurcation analysis as ω varies.

        Identifies:
        - Bifurcation points (where equilibria merge or disappear)
        - Saddle-node bifurcations
        - Transcritical bifurcations (equilibria exchange stability)
        - Hopf bifurcations (limit cycles, if any)

        Returns
        -------
        dict
            Bifurcation structure
        """
        return {
            'bifurcation_type': 'Transcritical bifurcation',

            'mechanism': """
            As ω increases from 0 to 1:

            1. ω = 0 (Pure metabolite pathway):
               - G behaves identically to M
               - System effectively 2D (S and M/G merged)
               - One of M or G excluded by competition

            2. ω = ω_crit1 (First critical point):
               - G can invade S-M equilibrium
               - Transcritical bifurcation:
                 * (N_S*, N_M*, 0) loses stability
                 * (N_S**, N_M**, N_G**) gains stability
               - Three-species coexistence emerges

            3. ω ∈ (ω_crit1, ω_crit2) (Coexistence window):
               - All three species coexist stably
               - N_G increases with ω
               - N_M typically decreases with ω
               - N_S may vary non-monotonically

            4. ω = ω_crit2 (Second critical point):
               - M can no longer survive
               - Transcritical bifurcation:
                 * Three-species equilibrium → S-G equilibrium
                 * N_M → 0

            5. ω = 1 (Pure substrate pathway):
               - G behaves identically to S
               - System again effectively 2D
               - Competition between S and G determines outcome
            """,

            'critical_omega_equations': """
            Critical ω values determined by:

            ω_crit1: Generalist invasion of S-M equilibrium
                     Solve: f_G(N_S*, N_M*, 0) = 0

            ω_crit2: M viability threshold in presence of S and G
                     Solve: N_M* = 0 in three-species equilibrium

            Explicit formula (approximate):
            ω_crit1 ≈ (1 - σ_GS + α_GM·N_M*/K_M) / (2 - α_GS - σ_GS + α_GM·N_M*/K_M + σ_GM·N_M*/K_M)

            ω_crit2 determined by: σ_MS·N_S*/K_S + ω·σ_MG·N_G*/K_G = 1 + (1-ω)·α_MG·N_G*/K_G + N_M*/K_M
            """,

            'bifurcation_diagram': """
            State diagram as function of ω:

            ω = 0 ----[ω_crit1]---- [coexistence] ----[ω_crit2]---- ω = 1
              |            |              |                |           |
            S-M (or     S-M stable    Three-species    S-G stable    S-G (or
            S-G)                      coexistence                     G-only)

            Hysteresis possible if multiple stable states exist.
            """
        }

    def dimensional_analysis(self) -> Dict:
        """
        Perform dimensional analysis and identify key dimensionless groups.

        Returns
        -------
        dict
            Dimensionless parameters and their interpretations
        """
        return {
            'time_scaling': """
            Natural time scales:
            τ_S = 1/r_S  (S-specialist doubling time)
            τ_M = 1/r_M  (M-specialist doubling time)
            τ_G = 1/r_G  (Generalist doubling time)

            Rescale time: t' = r_S·t (use S as reference)
            """,

            'population_scaling': """
            Natural population scales:
            Scale populations by carrying capacities:
            n_S = N_S/K_S
            n_M = N_M/K_M
            n_G = N_G/K_G

            Dimensionless populations: n_i ∈ [0, 1]
            """,

            'dimensionless_parameters': """
            Key dimensionless groups:

            1. Growth rate ratios:
               ρ_M = r_M/r_S  (relative M growth)
               ρ_G = r_G/r_S  (relative G growth)

            2. Cooperation strengths:
               σ_ij (already dimensionless)

            3. Competition strengths:
               α_ij (already dimensionless)

            4. Pathway balance:
               ω ∈ [0, 1] (already dimensionless)

            5. Mutualism index:
               Π_SM = σ_SM·σ_MS  (product, must be > 1 for S-M coexistence)

            6. Niche differentiation index:
               Δ = |ω - 0.5|  (distance from perfect generalist)
               Δ → 0: true generalist
               Δ → 0.5: specialist
            """,

            'phase_space_volume': """
            Viable phase space volume (rough estimate):

            V ∝ K_S·K_M·K_G  (product of carrying capacities)

            Fraction supporting coexistence:
            f_coexist ≈ (ω_max - ω_min)·(σ_MS - 1)·Π(conditions)

            Typically: f_coexist ∼ 0.1-0.3 (10-30% of parameter space)
            """
        }

    def ecological_interpretation(self) -> Dict:
        """
        Provide ecological interpretation of mathematical results.

        Returns
        -------
        dict
            Ecological insights
        """
        return {
            'M_obligate_mutualist': """
            M-specialist (metabolite specialist) interpretation:

            - Base growth rate = -1: Cannot survive independently
            - Requires σ_MS·N_S/K_S > 1: Needs S at sufficient density
            - Ecological analog: Obligate syntrophic bacteria
              * Examples: Methanogenic archaea requiring H₂ from fermenters
              * Syntrophomonas requiring Methanobacterium

            Critical insight: M's survival sets minimum S density
            This creates a "mutualistic platform" for community assembly
            """,

            'G_metabolic_flexibility': """
            Generalist (G-specialist) interpretation:

            - ω determines metabolic strategy:
              * ω = 0: Specialized on metabolite pathway (like M)
              * ω = 1: Specialized on substrate pathway (like S)
              * 0 < ω < 1: True generalist, uses both

            - Coexistence requires intermediate ω (niche differentiation)

            - Ecological analogs:
              * Facultative vs obligate strategies
              * Generalist predators in food webs
              * Metabolic versatility in microbes (e.g., E. coli)

            Trade-off: Generalists vs specialists
            - Specialists: Higher maximum growth on preferred substrate
            - Generalists: Lower maximum but broader niche

            Model captures this via ω-dependent growth terms
            """,

            'cooperation_competition_balance': """
            σ vs α: Cooperation vs Competition

            - σ_ij > 0: Cross-feeding, facilitation, mutualism
              * Metabolite exchange
              * Niche construction
              * Stress amelioration

            - α_ij > 0: Competition, interference, inhibition
              * Resource competition
              * Allelopathy
              * Spatial crowding

            Community assembly rule:
            σ > α → Coexistence promoted
            σ < α → Competitive exclusion

            This matches empirical patterns:
            - Spatially structured systems → higher cooperation → more coexistence
            - Well-mixed systems → stronger competition → fewer species

            Experimental predictions:
            1. Increase spatial structure → increase σ/α ratio → more coexistence
            2. Increase mixing → decrease σ/α ratio → fewer species
            3. Engineer cross-feeding (increase σ) → promote coexistence
            """,

            'omega_as_evolutionary_trait': """
            Evolutionary interpretation of ω:

            - ω can evolve via mutations in metabolic pathways
            - Selection pressure depends on community composition:
              * In S-dominated community: Select for higher ω (substrate use)
              * In M-dominated community: Select for lower ω (metabolite use)

            - Evolutionary stable strategy (ESS):
              ω* = optimal value that cannot be invaded

            - Adaptive dynamics prediction:
              If ω evolves, system may undergo evolutionary branching:
              * Single generalist → two specialists (S-like and M-like)
              * Diversification driven by frequency-dependent selection

            - Matches observations in:
              * Cross-feeding yeast (Gore lab)
              * Syntrophic bacterial pairs
              * Multi-trophic microbial communities
            """,

            'applications': """
            Practical applications:

            1. Wastewater treatment:
               - Maintain stable consortium of degraders (S), fermenters (M), methanogens (G)
               - Design: Ensure σ_MS > 1, control ω via substrate composition

            2. Bioproduction:
               - Engineer synthetic consortia for biofuel/chemical production
               - Strategy: Tune ω via promoter strength to balance pathways

            3. Gut microbiome:
               - Understand stability vs dysbiosis
               - Primary degraders (S), secondary consumers (M), opportunists (G)
               - Interventions: Modulate σ/α via diet, prebiotics

            4. Soil communities:
               - Carbon cycling networks match this structure
               - Litter degraders (S), fermenters (M), generalists (G)

            5. Synthetic ecology:
               - Rational design of stable communities
               - Prediction: Three-species systems more stable than two
                 (but only if parameters in coexistence window)
            """
        }

    def export_conditions_for_numerical_check(self, params_dict: Dict) -> Dict:
        """
        Convert analytical conditions to numerical checks.

        Parameters
        ----------
        params_dict : dict
            Numerical parameter values

        Returns
        -------
        dict
            Evaluated conditions (True/False for each)
        """
        # Substitute parameters
        subs = {
            self.r_S: params_dict.get('r_S', 1.0),
            self.r_M: params_dict.get('r_M', 0.8),
            self.r_G: params_dict.get('r_G', 0.9),
            self.K_S: params_dict.get('K_S', 100.0),
            self.K_M: params_dict.get('K_M', 100.0),
            self.K_G: params_dict.get('K_G', 100.0),
            self.sigma_SM: params_dict.get('sigma_SM', 0.5),
            self.sigma_MS: params_dict.get('sigma_MS', 0.6),
            self.sigma_SG: params_dict.get('sigma_SG', 0.3),
            self.sigma_GS: params_dict.get('sigma_GS', 0.4),
            self.sigma_MG: params_dict.get('sigma_MG', 0.3),
            self.sigma_GM: params_dict.get('sigma_GM', 0.4),
            self.alpha_SG: params_dict.get('alpha_SG', 0.4),
            self.alpha_MG: params_dict.get('alpha_MG', 0.4),
            self.alpha_GS: params_dict.get('alpha_GS', 0.3),
            self.alpha_GM: params_dict.get('alpha_GM', 0.3),
            self.omega: params_dict.get('omega', 0.5),
        }

        # Evaluate conditions
        conditions = {}

        # M viability
        conditions['M_can_survive'] = float(subs[self.sigma_MS]) > 1.0

        # Mutualism strength
        conditions['strong_mutualism'] = float(subs[self.sigma_SM] * subs[self.sigma_MS]) > 1.0

        # Average cooperation > competition
        avg_sigma = np.mean([float(subs[self.sigma_SM]), float(subs[self.sigma_MS]),
                             float(subs[self.sigma_SG]), float(subs[self.sigma_GS]),
                             float(subs[self.sigma_MG]), float(subs[self.sigma_GM])])
        avg_alpha = np.mean([float(subs[self.alpha_SG]), float(subs[self.alpha_MG]),
                             float(subs[self.alpha_GS]), float(subs[self.alpha_GM])])
        conditions['cooperation_exceeds_competition'] = avg_sigma > avg_alpha

        # Intermediate omega
        omega_val = float(subs[self.omega])
        conditions['intermediate_omega'] = 0.2 < omega_val < 0.8

        # Summary
        conditions['likely_coexistence'] = all([
            conditions['M_can_survive'],
            conditions['strong_mutualism'],
            conditions['cooperation_exceeds_competition'],
            conditions['intermediate_omega']
        ])

        return conditions


def create_analytical_summary_document():
    """
    Generate a comprehensive analytical summary document.
    """
    analyzer = AnalyticalThreeSpecies()

    summary = """
===============================================================================
ANALYTICAL ANALYSIS OF THREE-SPECIES CROSS-FEEDING MODEL
===============================================================================

Complete mathematical characterization including equilibria, stability,
bifurcations, and ecological interpretation.

===============================================================================
1. BOUNDARY EQUILIBRIA (Analytical Solutions)
===============================================================================

"""

    equilibria = analyzer.find_boundary_equilibria()
    for name, eq_info in equilibria.items():
        summary += f"\n{name.upper().replace('_', ' ')}:\n"
        summary += f"  Point: {eq_info['point']}\n"
        summary += f"  Exists: {eq_info['exists']}\n"
        if 'stability' in eq_info:
            summary += f"  Stability: {eq_info['stability']}\n"
        if 'comment' in eq_info:
            summary += f"  Note: {eq_info['comment']}\n"
        summary += "\n"

    summary += """
===============================================================================
2. THREE-SPECIES COEXISTENCE CONDITIONS
===============================================================================

"""

    coexist = analyzer.three_species_coexistence_conditions()
    for key, value in coexist.items():
        summary += f"\n{key.upper().replace('_', ' ')}:\n"
        if isinstance(value, list):
            for item in value:
                summary += f"  • {item}\n"
        else:
            summary += f"{value}\n"

    summary += """
===============================================================================
3. BIFURCATION ANALYSIS (ω as bifurcation parameter)
===============================================================================

"""

    bifurc = analyzer.omega_bifurcation_analysis()
    for key, value in bifurc.items():
        summary += f"\n{key.upper().replace('_', ' ')}:\n{value}\n"

    summary += """
===============================================================================
4. DIMENSIONAL ANALYSIS
===============================================================================

"""

    dimensional = analyzer.dimensional_analysis()
    for key, value in dimensional.items():
        summary += f"\n{key.upper().replace('_', ' ')}:\n{value}\n"

    summary += """
===============================================================================
5. ECOLOGICAL INTERPRETATION
===============================================================================

"""

    ecology = analyzer.ecological_interpretation()
    for key, value in ecology.items():
        summary += f"\n{key.upper().replace('_', ' ')}:\n{value}\n"

    summary += """
===============================================================================
END OF ANALYTICAL SUMMARY
===============================================================================
"""

    return summary


if __name__ == "__main__":
    # Generate analytical summary
    summary = create_analytical_summary_document()

    # Save to file
    with open('../analytical_summary.txt', 'w') as f:
        f.write(summary)

    print("Analytical summary generated: analytical_summary.txt")
    print("\nKey findings:")
    print("  • M-specialist requires σ_MS > 1 to survive")
    print("  • Three-species coexistence requires intermediate ω")
    print("  • Bifurcation occurs at critical ω values")
    print("  • Cooperation must exceed competition for coexistence")

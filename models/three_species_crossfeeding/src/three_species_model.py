"""
Three-Species Cross-Feeding Model

A mathematical model describing the dynamics of three strategists in a microbial community:
- S-specialist: Substrate specialists
- M-specialist: Metabolite specialists
- G-specialist: Generalists (can use both pathways)

The model integrates substrate utilization, metabolite cross-feeding, and synergistic
growth through a weighted combination of metabolic pathways.

Author: Jian Wang
Date: January 2026
"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import fsolve, root
from typing import Dict, Tuple, List, Optional
import warnings


class ThreeSpeciesModel:
    """
    Three-species cross-feeding model with substrate specialists, metabolite specialists,
    and generalists.

    Model equations:
    dN_S/dt = r_S * N_S * [1 + σ_{SM} * N_M/K_M + (1-ω) * σ_{SG} * N_G/K_G
                           - ω * α_{SG} * N_G/K_G - N_S/K_S]

    dN_M/dt = r_M * N_M * [-1 + σ_{MS} * N_S/K_S + ω * σ_{MG} * N_G/K_G
                           - (1-ω) * α_{MG} * N_G/K_G - N_M/K_M]

    dN_G/dt = r_G * N_G * [ω(1 - α_{GS}*N_S/K_S + σ_{GM}*N_M/K_M)
                           + (1-ω)(-1 - α_{GM}*N_M/K_M + σ_{GS}*N_S/K_S) - N_G/K_G]
    """

    def __init__(self, params: Optional[Dict] = None):
        """
        Initialize the three-species model with parameters.

        Parameters
        ----------
        params : dict, optional
            Model parameters. If None, uses default parameters.

        Parameter definitions:
        - r_S, r_M, r_G: Intrinsic growth rates (1/time)
        - K_S, K_M, K_G: Carrying capacities (cells/volume)
        - σ_{ij}: Synergistic interaction coefficients (dimensionless, > 0)
        - α_{ij}: Competition coefficients (dimensionless, > 0)
        - ω: Pathway weighting parameter (0 ≤ ω ≤ 1)
        """
        if params is None:
            self.params = self._default_parameters()
        else:
            self.params = params

    def _default_parameters(self) -> Dict:
        """
        Default parameter set representing moderate cross-feeding and competition.

        Returns
        -------
        dict
            Default parameter dictionary
        """
        return {
            # Growth rates
            'r_S': 1.0,    # S-specialist growth rate
            'r_M': 0.8,    # M-specialist growth rate (slightly slower)
            'r_G': 0.9,    # Generalist growth rate (intermediate)

            # Carrying capacities
            'K_S': 100.0,  # S-specialist carrying capacity
            'K_M': 100.0,  # M-specialist carrying capacity
            'K_G': 100.0,  # Generalist carrying capacity

            # Synergistic coefficients (cooperation)
            'sigma_SM': 0.5,  # S benefits from M's metabolites
            'sigma_MS': 0.6,  # M benefits from S's substrate breakdown
            'sigma_SG': 0.3,  # S benefits from G (substrate pathway)
            'sigma_GS': 0.4,  # G benefits from S (metabolite pathway)
            'sigma_MG': 0.3,  # M benefits from G (metabolite pathway)
            'sigma_GM': 0.4,  # G benefits from M (substrate pathway)

            # Competition coefficients
            'alpha_SG': 0.4,  # S competes with G (metabolite pathway)
            'alpha_MG': 0.4,  # M competes with G (substrate pathway)
            'alpha_GS': 0.3,  # G competes with S (substrate pathway)
            'alpha_GM': 0.3,  # G competes with M (metabolite pathway)

            # Pathway weighting
            'omega': 0.5,    # Equal weighting between pathways for generalist
        }

    def equations(self, t: float, N: np.ndarray) -> np.ndarray:
        """
        System of ODEs for the three-species model.

        Parameters
        ----------
        t : float
            Time point (not explicitly used in autonomous system)
        N : array-like, shape (3,)
            Population densities [N_S, N_M, N_G]

        Returns
        -------
        np.ndarray, shape (3,)
            Time derivatives [dN_S/dt, dN_M/dt, dN_G/dt]
        """
        N_S, N_M, N_G = N
        p = self.params

        # Prevent negative populations
        N_S = max(0, N_S)
        N_M = max(0, N_M)
        N_G = max(0, N_G)

        # S-specialist dynamics (Eq. 1)
        dN_S = p['r_S'] * N_S * (
            1
            + p['sigma_SM'] * N_M / p['K_M']
            + (1 - p['omega']) * p['sigma_SG'] * N_G / p['K_G']
            - p['omega'] * p['alpha_SG'] * N_G / p['K_G']
            - N_S / p['K_S']
        )

        # M-specialist dynamics (Eq. 2)
        dN_M = p['r_M'] * N_M * (
            -1
            + p['sigma_MS'] * N_S / p['K_S']
            + p['omega'] * p['sigma_MG'] * N_G / p['K_G']
            - (1 - p['omega']) * p['alpha_MG'] * N_G / p['K_G']
            - N_M / p['K_M']
        )

        # Generalist dynamics (Eq. 3)
        # Substrate pathway (weighted by ω)
        substrate_pathway = p['omega'] * (
            1
            - p['alpha_GS'] * N_S / p['K_S']
            + p['sigma_GM'] * N_M / p['K_M']
        )

        # Metabolite pathway (weighted by 1-ω)
        metabolite_pathway = (1 - p['omega']) * (
            -1
            - p['alpha_GM'] * N_M / p['K_M']
            + p['sigma_GS'] * N_S / p['K_S']
        )

        dN_G = p['r_G'] * N_G * (
            substrate_pathway + metabolite_pathway - N_G / p['K_G']
        )

        return np.array([dN_S, dN_M, dN_G])

    def simulate(self, N0: np.ndarray, t_span: Tuple[float, float],
                 t_eval: Optional[np.ndarray] = None, method: str = 'RK45') -> Dict:
        """
        Simulate the three-species model.

        Parameters
        ----------
        N0 : array-like, shape (3,)
            Initial conditions [N_S(0), N_M(0), N_G(0)]
        t_span : tuple
            Time span (t_start, t_end)
        t_eval : array-like, optional
            Time points at which to store solution
        method : str, default 'RK45'
            Integration method for solve_ivp

        Returns
        -------
        dict
            Solution dictionary with 't' (time) and 'N' (populations)
        """
        if t_eval is None:
            t_eval = np.linspace(t_span[0], t_span[1], 1000)

        sol = solve_ivp(
            self.equations,
            t_span,
            N0,
            method=method,
            t_eval=t_eval,
            dense_output=True
        )

        return {
            't': sol.t,
            'N_S': sol.y[0],
            'N_M': sol.y[1],
            'N_G': sol.y[2],
            'success': sol.success,
            'message': sol.message
        }

    def find_equilibria(self, n_attempts: int = 50,
                       search_range: Tuple[float, float] = (0, 200)) -> List[np.ndarray]:
        """
        Find equilibrium points of the system using multiple initial guesses.

        Parameters
        ----------
        n_attempts : int, default 50
            Number of random initial guesses to try
        search_range : tuple, default (0, 200)
            Range for random initial guesses

        Returns
        -------
        list of np.ndarray
            List of unique equilibrium points found
        """
        def equilibrium_equations(N):
            """Equilibrium condition: dN/dt = 0"""
            return self.equations(0, N)

        equilibria = []

        # Always check the trivial equilibrium (0, 0, 0)
        initial_guesses = [np.array([0.0, 0.0, 0.0])]

        # Add boundary equilibria
        K_S, K_M, K_G = self.params['K_S'], self.params['K_M'], self.params['K_G']
        initial_guesses.extend([
            np.array([K_S, 0.0, 0.0]),
            np.array([0.0, K_M, 0.0]),
            np.array([0.0, 0.0, K_G]),
            np.array([K_S, K_M, 0.0]),
            np.array([K_S, 0.0, K_G]),
            np.array([0.0, K_M, K_G]),
            np.array([K_S, K_M, K_G]),
        ])

        # Add random guesses
        for _ in range(n_attempts):
            N_guess = np.random.uniform(search_range[0], search_range[1], 3)
            initial_guesses.append(N_guess)

        for N_guess in initial_guesses:
            try:
                sol = fsolve(equilibrium_equations, N_guess, full_output=True)
                N_eq = sol[0]
                info = sol[1]

                # Check if solution converged and is valid
                if info['fvec'].dot(info['fvec']) < 1e-9:  # Residual is small
                    # Round to avoid numerical duplicates
                    N_eq_rounded = np.round(N_eq, decimals=6)

                    # Check if equilibrium is new and non-negative
                    if np.all(N_eq >= -1e-6):  # Allow small numerical errors
                        N_eq_rounded[N_eq_rounded < 0] = 0  # Clean up small negatives

                        is_duplicate = False
                        for existing_eq in equilibria:
                            if np.allclose(N_eq_rounded, existing_eq, atol=1e-5):
                                is_duplicate = True
                                break

                        if not is_duplicate:
                            equilibria.append(N_eq_rounded)
            except:
                continue

        return equilibria

    def jacobian(self, N: np.ndarray) -> np.ndarray:
        """
        Compute the Jacobian matrix at a given point.

        Parameters
        ----------
        N : array-like, shape (3,)
            Population densities [N_S, N_M, N_G]

        Returns
        -------
        np.ndarray, shape (3, 3)
            Jacobian matrix ∂f_i/∂N_j
        """
        N_S, N_M, N_G = N
        p = self.params

        # Extract parameters for readability
        r_S, r_M, r_G = p['r_S'], p['r_M'], p['r_G']
        K_S, K_M, K_G = p['K_S'], p['K_M'], p['K_G']
        omega = p['omega']

        # S-specialist equation components
        f_S = (1 + p['sigma_SM'] * N_M/K_M + (1-omega) * p['sigma_SG'] * N_G/K_G
               - omega * p['alpha_SG'] * N_G/K_G - N_S/K_S)

        # M-specialist equation components
        f_M = (-1 + p['sigma_MS'] * N_S/K_S + omega * p['sigma_MG'] * N_G/K_G
               - (1-omega) * p['alpha_MG'] * N_G/K_G - N_M/K_M)

        # Generalist equation components
        f_G = (omega * (1 - p['alpha_GS']*N_S/K_S + p['sigma_GM']*N_M/K_M)
               + (1-omega) * (-1 - p['alpha_GM']*N_M/K_M + p['sigma_GS']*N_S/K_S)
               - N_G/K_G)

        # Jacobian elements
        J = np.zeros((3, 3))

        # Row 1: ∂(dN_S/dt)/∂N_i
        J[0, 0] = r_S * (f_S - N_S/K_S)  # ∂/∂N_S
        J[0, 1] = r_S * N_S * p['sigma_SM'] / K_M  # ∂/∂N_M
        J[0, 2] = r_S * N_S * ((1-omega)*p['sigma_SG'] - omega*p['alpha_SG']) / K_G  # ∂/∂N_G

        # Row 2: ∂(dN_M/dt)/∂N_i
        J[1, 0] = r_M * N_M * p['sigma_MS'] / K_S  # ∂/∂N_S
        J[1, 1] = r_M * (f_M - N_M/K_M)  # ∂/∂N_M
        J[1, 2] = r_M * N_M * (omega*p['sigma_MG'] - (1-omega)*p['alpha_MG']) / K_G  # ∂/∂N_G

        # Row 3: ∂(dN_G/dt)/∂N_i
        J[2, 0] = r_G * N_G * (-omega*p['alpha_GS'] + (1-omega)*p['sigma_GS']) / K_S  # ∂/∂N_S
        J[2, 1] = r_G * N_G * (omega*p['sigma_GM'] - (1-omega)*p['alpha_GM']) / K_M  # ∂/∂N_M
        J[2, 2] = r_G * (f_G - N_G/K_G)  # ∂/∂N_G

        return J

    def stability_analysis(self, N_eq: np.ndarray) -> Dict:
        """
        Analyze stability of an equilibrium point.

        Parameters
        ----------
        N_eq : array-like, shape (3,)
            Equilibrium point to analyze

        Returns
        -------
        dict
            Dictionary containing:
            - 'eigenvalues': Complex eigenvalues of Jacobian
            - 'eigenvectors': Corresponding eigenvectors
            - 'stable': Boolean indicating stability
            - 'type': Classification of equilibrium point
        """
        J = self.jacobian(N_eq)
        eigenvalues, eigenvectors = np.linalg.eig(J)

        # Stability: all eigenvalues must have negative real parts
        stable = np.all(np.real(eigenvalues) < 0)

        # Classify equilibrium type
        real_parts = np.real(eigenvalues)
        imag_parts = np.imag(eigenvalues)

        if stable:
            if np.all(np.abs(imag_parts) < 1e-10):
                eq_type = "Stable node"
            else:
                eq_type = "Stable spiral"
        else:
            if np.any(real_parts > 0) and np.any(real_parts < 0):
                eq_type = "Saddle point"
            elif np.all(real_parts > 0):
                if np.all(np.abs(imag_parts) < 1e-10):
                    eq_type = "Unstable node"
                else:
                    eq_type = "Unstable spiral"
            else:
                eq_type = "Unknown"

        return {
            'equilibrium': N_eq,
            'jacobian': J,
            'eigenvalues': eigenvalues,
            'eigenvectors': eigenvectors,
            'stable': stable,
            'type': eq_type,
            'max_real_eigenvalue': np.max(real_parts)
        }

    def classify_equilibrium_ecology(self, N_eq: np.ndarray,
                                     threshold: float = 1e-3) -> str:
        """
        Classify equilibrium point ecologically.

        Parameters
        ----------
        N_eq : array-like, shape (3,)
            Equilibrium point [N_S, N_M, N_G]
        threshold : float, default 1e-3
            Threshold below which population is considered extinct

        Returns
        -------
        str
            Ecological classification
        """
        N_S, N_M, N_G = N_eq

        present = np.array([N_S, N_M, N_G]) > threshold

        if not np.any(present):
            return "Extinction"
        elif np.all(present):
            return "Three-species coexistence"
        elif np.sum(present) == 2:
            if present[0] and present[1]:
                return "S-M coexistence (no generalist)"
            elif present[0] and present[2]:
                return "S-G coexistence (no M-specialist)"
            else:  # M and G
                return "M-G coexistence (no S-specialist)"
        else:  # Only one species
            if present[0]:
                return "S-specialist dominance"
            elif present[1]:
                return "M-specialist dominance (unlikely - requires cross-feeding)"
            else:
                return "Generalist dominance"

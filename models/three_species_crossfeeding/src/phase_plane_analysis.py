"""
Phase Plane Analysis for Three-Species Model

Jeff Gore-style phase plane analysis including:
- 2D phase portraits (projections of 3D system)
- Nullcline calculations and visualization
- Vector field analysis
- Trajectory plotting
- Basin of attraction analysis

Author: Jian Wang
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D, proj3d
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Callable
from three_species_model import ThreeSpeciesModel


class Arrow3D(FancyArrowPatch):
    """Helper class for 3D arrows in matplotlib."""

    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)


class PhasePlaneAnalyzer:
    """
    Phase plane analysis tools for the three-species cross-feeding model.

    Since the system is 3-dimensional, we analyze 2D projections (phase planes)
    while treating the third variable as either constant or allowing it to vary.
    """

    def __init__(self, model: ThreeSpeciesModel):
        """
        Initialize phase plane analyzer.

        Parameters
        ----------
        model : ThreeSpeciesModel
            Instance of the three-species model
        """
        self.model = model
        self.setup_plotting_style()

    def setup_plotting_style(self):
        """Set up publication-quality plotting style."""
        sns.set_context("paper", font_scale=1.5)
        sns.set_style("whitegrid")
        plt.rcParams['figure.figsize'] = (10, 8)
        plt.rcParams['axes.labelsize'] = 14
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12
        plt.rcParams['legend.fontsize'] = 11

    def compute_nullclines_2D(self, species_pair: Tuple[int, int],
                             fixed_species_val: float,
                             grid_range: Tuple[float, float] = (0, 150),
                             n_points: int = 100) -> Dict[str, np.ndarray]:
        """
        Compute nullclines for a 2D projection of the 3D system.

        Parameters
        ----------
        species_pair : tuple of int
            Which two species to plot (0=S, 1=M, 2=G), e.g., (0, 1) for S-M plane
        fixed_species_val : float
            Value at which to fix the third species
        grid_range : tuple, default (0, 150)
            Range of values for grid
        n_points : int, default 100
            Number of grid points

        Returns
        -------
        dict
            Dictionary containing nullcline information
        """
        idx1, idx2 = species_pair
        fixed_idx = [i for i in range(3) if i not in species_pair][0]

        # Create grid
        x1 = np.linspace(grid_range[0], grid_range[1], n_points)
        x2 = np.linspace(grid_range[0], grid_range[1], n_points)
        X1, X2 = np.meshgrid(x1, x2)

        # Compute derivatives on grid
        dN1 = np.zeros_like(X1)
        dN2 = np.zeros_like(X2)

        for i in range(n_points):
            for j in range(n_points):
                N = np.zeros(3)
                N[idx1] = X1[i, j]
                N[idx2] = X2[i, j]
                N[fixed_idx] = fixed_species_val

                dN = self.model.equations(0, N)
                dN1[i, j] = dN[idx1]
                dN2[i, j] = dN[idx2]

        return {
            'X1': X1,
            'X2': X2,
            'dN1': dN1,
            'dN2': dN2,
            'x1': x1,
            'x2': x2,
            'species_pair': species_pair,
            'fixed_species': fixed_idx,
            'fixed_value': fixed_species_val
        }

    def plot_phase_portrait_2D(self, species_pair: Tuple[int, int],
                               fixed_species_val: float,
                               initial_conditions: Optional[List[np.ndarray]] = None,
                               t_max: float = 100,
                               grid_range: Tuple[float, float] = (0, 150),
                               n_grid: int = 20,
                               show_nullclines: bool = True,
                               show_equilibria: bool = True,
                               ax: Optional[plt.Axes] = None) -> plt.Figure:
        """
        Create a 2D phase portrait with nullclines, vector field, and trajectories.

        Parameters
        ----------
        species_pair : tuple of int
            Which two species to plot (0=S, 1=M, 2=G)
        fixed_species_val : float
            Value at which to fix the third species
        initial_conditions : list of array-like, optional
            List of initial conditions to simulate
        t_max : float, default 100
            Maximum simulation time
        grid_range : tuple, default (0, 150)
            Range for phase plane
        n_grid : int, default 20
            Number of grid points for vector field
        show_nullclines : bool, default True
            Whether to show nullclines
        show_equilibria : bool, default True
            Whether to mark equilibrium points
        ax : matplotlib.axes.Axes, optional
            Axes to plot on

        Returns
        -------
        matplotlib.figure.Figure
            Figure object
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 8))
        else:
            fig = ax.get_figure()

        idx1, idx2 = species_pair
        fixed_idx = [i for i in range(3) if i not in species_pair][0]

        species_names = ['S (substrate specialist)', 'M (metabolite specialist)', 'G (generalist)']
        species_abbrev = ['N_S', 'N_M', 'N_G']

        # Compute nullclines
        nullcline_data = self.compute_nullclines_2D(
            species_pair, fixed_species_val, grid_range, n_points=100
        )

        # Plot nullclines
        if show_nullclines:
            # Nullcline for species 1 (where dN1/dt = 0)
            ax.contour(nullcline_data['X1'], nullcline_data['X2'],
                      nullcline_data['dN1'], levels=[0],
                      colors='blue', linewidths=2, linestyles='--',
                      label=f'd{species_abbrev[idx1]}/dt = 0')

            # Nullcline for species 2 (where dN2/dt = 0)
            ax.contour(nullcline_data['X1'], nullcline_data['X2'],
                      nullcline_data['dN2'], levels=[0],
                      colors='red', linewidths=2, linestyles='--',
                      label=f'd{species_abbrev[idx2]}/dt = 0')

        # Vector field
        x1_vec = np.linspace(grid_range[0], grid_range[1], n_grid)
        x2_vec = np.linspace(grid_range[0], grid_range[1], n_grid)
        X1_vec, X2_vec = np.meshgrid(x1_vec, x2_vec)

        U = np.zeros_like(X1_vec)
        V = np.zeros_like(X2_vec)

        for i in range(n_grid):
            for j in range(n_grid):
                N = np.zeros(3)
                N[idx1] = X1_vec[i, j]
                N[idx2] = X2_vec[i, j]
                N[fixed_idx] = fixed_species_val

                dN = self.model.equations(0, N)
                U[i, j] = dN[idx1]
                V[i, j] = dN[idx2]

        # Normalize vectors for better visualization
        M = np.sqrt(U**2 + V**2)
        M[M == 0] = 1  # Avoid division by zero
        U_norm = U / M
        V_norm = V / M

        ax.quiver(X1_vec, X2_vec, U_norm, V_norm,
                 M, alpha=0.6, cmap='viridis', scale=25)

        # Plot trajectories
        if initial_conditions is not None:
            colors = plt.cm.tab10(np.linspace(0, 1, len(initial_conditions)))

            for ic_idx, N0 in enumerate(initial_conditions):
                # Ensure N0 has the correct dimension
                if len(N0) != 3:
                    N0_full = np.zeros(3)
                    N0_full[idx1] = N0[0]
                    N0_full[idx2] = N0[1]
                    N0_full[fixed_idx] = fixed_species_val
                else:
                    N0_full = N0

                sol = self.model.simulate(N0_full, (0, t_max))

                if sol['success']:
                    ax.plot(sol[f'N_{["S", "M", "G"][idx1]}'],
                           sol[f'N_{["S", "M", "G"][idx2]}'],
                           color=colors[ic_idx], linewidth=1.5, alpha=0.7)

                    # Mark initial condition
                    ax.plot(N0_full[idx1], N0_full[idx2], 'o',
                           color=colors[ic_idx], markersize=8,
                           markeredgecolor='black', markeredgewidth=1)

        # Find and mark equilibria
        if show_equilibria:
            equilibria = self.model.find_equilibria()

            for eq in equilibria:
                # Only plot if the fixed species is close to the specified value
                if np.abs(eq[fixed_idx] - fixed_species_val) < 10:
                    stability = self.model.stability_analysis(eq)

                    if stability['stable']:
                        marker = '*'
                        color = 'green'
                        size = 200
                    else:
                        marker = 'X'
                        color = 'red'
                        size = 150

                    ax.scatter(eq[idx1], eq[idx2], marker=marker,
                              s=size, c=color, edgecolors='black',
                              linewidths=2, zorder=10,
                              label='Stable' if stability['stable'] else 'Unstable')

        ax.set_xlabel(f'{species_names[idx1]} ({species_abbrev[idx1]})', fontsize=14)
        ax.set_ylabel(f'{species_names[idx2]} ({species_abbrev[idx2]})', fontsize=14)
        ax.set_xlim(grid_range)
        ax.set_ylim(grid_range)
        ax.set_title(f'Phase Portrait: {species_abbrev[idx1]}-{species_abbrev[idx2]} plane\n'
                    f'({species_abbrev[fixed_idx]} = {fixed_species_val:.1f}, ω = {self.model.params["omega"]:.2f})',
                    fontsize=16)

        # Remove duplicate labels in legend
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc='best', frameon=True)

        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        return fig

    def plot_3D_phase_space(self, initial_conditions: List[np.ndarray],
                           t_max: float = 100,
                           show_equilibria: bool = True,
                           elev: float = 20, azim: float = 45) -> plt.Figure:
        """
        Create a 3D phase space plot showing all three species.

        Parameters
        ----------
        initial_conditions : list of array-like
            List of initial conditions [N_S(0), N_M(0), N_G(0)]
        t_max : float, default 100
            Maximum simulation time
        show_equilibria : bool, default True
            Whether to mark equilibrium points
        elev : float, default 20
            Elevation angle for 3D view
        azim : float, default 45
            Azimuthal angle for 3D view

        Returns
        -------
        matplotlib.figure.Figure
            Figure object
        """
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        colors = plt.cm.tab10(np.linspace(0, 1, len(initial_conditions)))

        # Plot trajectories
        for ic_idx, N0 in enumerate(initial_conditions):
            sol = self.model.simulate(N0, (0, t_max))

            if sol['success']:
                ax.plot(sol['N_S'], sol['N_M'], sol['N_G'],
                       color=colors[ic_idx], linewidth=2, alpha=0.7,
                       label=f'IC {ic_idx+1}')

                # Mark initial condition
                ax.scatter(N0[0], N0[1], N0[2],
                          color=colors[ic_idx], s=100,
                          marker='o', edgecolors='black', linewidths=2)

                # Mark endpoint
                ax.scatter(sol['N_S'][-1], sol['N_M'][-1], sol['N_G'][-1],
                          color=colors[ic_idx], s=50,
                          marker='^', edgecolors='black', linewidths=1)

        # Find and mark equilibria
        if show_equilibria:
            equilibria = self.model.find_equilibria()

            for eq in equilibria:
                stability = self.model.stability_analysis(eq)

                if stability['stable']:
                    marker = '*'
                    color = 'green'
                    size = 300
                    label = 'Stable equilibrium'
                else:
                    marker = 'X'
                    color = 'red'
                    size = 200
                    label = 'Unstable equilibrium'

                ax.scatter(eq[0], eq[1], eq[2],
                          marker=marker, s=size, c=color,
                          edgecolors='black', linewidths=2, zorder=10)

        ax.set_xlabel('N_S (substrate specialist)', fontsize=12, labelpad=10)
        ax.set_ylabel('N_M (metabolite specialist)', fontsize=12, labelpad=10)
        ax.set_zlabel('N_G (generalist)', fontsize=12, labelpad=10)
        ax.set_title(f'3D Phase Space (ω = {self.model.params["omega"]:.2f})',
                    fontsize=16, pad=20)

        ax.view_init(elev=elev, azim=azim)
        ax.legend(loc='best', frameon=True)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        return fig

    def plot_timeseries(self, initial_conditions: List[np.ndarray],
                       t_max: float = 100) -> plt.Figure:
        """
        Plot population dynamics over time.

        Parameters
        ----------
        initial_conditions : list of array-like
            List of initial conditions
        t_max : float, default 100
            Maximum simulation time

        Returns
        -------
        matplotlib.figure.Figure
            Figure object
        """
        n_ics = len(initial_conditions)
        fig, axes = plt.subplots(n_ics, 1, figsize=(12, 4*n_ics))

        if n_ics == 1:
            axes = [axes]

        species_colors = {'N_S': 'blue', 'N_M': 'red', 'N_G': 'green'}
        species_labels = {
            'N_S': 'S (substrate specialist)',
            'N_M': 'M (metabolite specialist)',
            'N_G': 'G (generalist)'
        }

        for ic_idx, N0 in enumerate(initial_conditions):
            ax = axes[ic_idx]
            sol = self.model.simulate(N0, (0, t_max))

            if sol['success']:
                for species in ['N_S', 'N_M', 'N_G']:
                    ax.plot(sol['t'], sol[species],
                           color=species_colors[species],
                           linewidth=2, label=species_labels[species])

                ax.set_xlabel('Time', fontsize=12)
                ax.set_ylabel('Population density', fontsize=12)
                ax.set_title(f'Initial condition: N_S={N0[0]:.1f}, N_M={N0[1]:.1f}, N_G={N0[2]:.1f}',
                           fontsize=14)
                ax.legend(loc='best', frameon=True)
                ax.grid(True, alpha=0.3)

        plt.tight_layout()

        return fig

    def bifurcation_analysis_omega(self, omega_range: Tuple[float, float] = (0, 1),
                                   n_omega: int = 50,
                                   N0: Optional[np.ndarray] = None) -> Dict:
        """
        Analyze how system behavior changes with pathway weighting parameter ω.

        Parameters
        ----------
        omega_range : tuple, default (0, 1)
            Range of ω values to explore
        n_omega : int, default 50
            Number of ω values to test
        N0 : array-like, optional
            Initial condition for simulations

        Returns
        -------
        dict
            Bifurcation data
        """
        omega_values = np.linspace(omega_range[0], omega_range[1], n_omega)

        if N0 is None:
            N0 = np.array([50.0, 50.0, 50.0])

        results = {
            'omega': omega_values,
            'equilibria': [],
            'stable_equilibria': [],
            'final_states': {'N_S': [], 'N_M': [], 'N_G': []},
            'coexistence': []
        }

        original_omega = self.model.params['omega']

        for omega in omega_values:
            # Update ω
            self.model.params['omega'] = omega

            # Simulate to find final state
            sol = self.model.simulate(N0, (0, 500))
            final_state = np.array([sol['N_S'][-1], sol['N_M'][-1], sol['N_G'][-1]])

            results['final_states']['N_S'].append(final_state[0])
            results['final_states']['N_M'].append(final_state[1])
            results['final_states']['N_G'].append(final_state[2])

            # Find equilibria
            equilibria = self.model.find_equilibria(n_attempts=20)
            stable_eq = []

            for eq in equilibria:
                stability = self.model.stability_analysis(eq)
                if stability['stable']:
                    stable_eq.append(eq)

            results['equilibria'].append(equilibria)
            results['stable_equilibria'].append(stable_eq)

            # Check for coexistence (all three species > 1)
            coexistence = np.all(final_state > 1.0)
            results['coexistence'].append(coexistence)

        # Restore original ω
        self.model.params['omega'] = original_omega

        return results

    def plot_bifurcation_diagram(self, bifurcation_data: Dict) -> plt.Figure:
        """
        Plot bifurcation diagram showing how equilibria change with ω.

        Parameters
        ----------
        bifurcation_data : dict
            Output from bifurcation_analysis_omega()

        Returns
        -------
        matplotlib.figure.Figure
            Figure object
        """
        fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

        species_names = ['N_S', 'N_M', 'N_G']
        species_labels = ['S (substrate specialist)',
                         'M (metabolite specialist)',
                         'G (generalist)']
        colors = ['blue', 'red', 'green']

        omega = bifurcation_data['omega']

        for idx, (species, label, color) in enumerate(zip(species_names, species_labels, colors)):
            ax = axes[idx]

            # Plot final states from simulations
            ax.plot(omega, bifurcation_data['final_states'][species],
                   'o', color=color, markersize=4, alpha=0.7,
                   label='Simulation endpoint')

            # Plot stable equilibria
            for i, omega_val in enumerate(omega):
                stable_eq = bifurcation_data['stable_equilibria'][i]
                for eq in stable_eq:
                    species_idx = ['N_S', 'N_M', 'N_G'].index(species)
                    ax.plot(omega_val, eq[species_idx], '*',
                           color='darkgreen', markersize=8)

            ax.set_ylabel(f'{label}', fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.set_ylim(bottom=-5)

            if idx == 0:
                ax.set_title('Bifurcation Diagram: Effect of pathway weighting (ω)',
                           fontsize=14)

        axes[-1].set_xlabel('ω (pathway weighting parameter)', fontsize=12)
        axes[-1].set_xlim(omega[0], omega[-1])

        # Add shaded region for coexistence
        coexist_omega = omega[bifurcation_data['coexistence']]
        if len(coexist_omega) > 0:
            for ax in axes:
                for i in range(len(omega)-1):
                    if bifurcation_data['coexistence'][i]:
                        ax.axvspan(omega[i], omega[i+1],
                                 alpha=0.1, color='green', zorder=-1)

        plt.tight_layout()

        return fig

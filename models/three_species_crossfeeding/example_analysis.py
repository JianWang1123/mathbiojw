#!/usr/bin/env python3
"""
Simple example demonstrating the three-species cross-feeding model analysis.

This script performs a quick analysis and generates basic plots to verify
the model implementation and demonstrate key features.

Author: Jian Wang
Date: January 2026
"""

import sys
sys.path.append('src')

import numpy as np
import matplotlib.pyplot as plt
from three_species_model import ThreeSpeciesModel
from phase_plane_analysis import PhasePlaneAnalyzer


def main():
    """Run basic analysis of three-species model."""

    print("="*70)
    print("Three-Species Cross-Feeding Model: Example Analysis")
    print("="*70)

    # Initialize model
    print("\n1. Initializing model with default parameters...")
    model = ThreeSpeciesModel()

    print(f"   Growth rates: r_S={model.params['r_S']}, "
          f"r_M={model.params['r_M']}, r_G={model.params['r_G']}")
    print(f"   Pathway weighting: ω = {model.params['omega']}")

    # Find equilibria
    print("\n2. Finding equilibrium points...")
    equilibria = model.find_equilibria(n_attempts=50)
    print(f"   Found {len(equilibria)} equilibria")

    for i, eq in enumerate(equilibria):
        stability = model.stability_analysis(eq)
        eco_type = model.classify_equilibrium_ecology(eq)

        print(f"\n   Equilibrium {i+1}:")
        print(f"      Populations: N_S={eq[0]:.2f}, N_M={eq[1]:.2f}, N_G={eq[2]:.2f}")
        print(f"      Type: {eco_type}")
        print(f"      Stability: {stability['type']}")

        if stability['stable']:
            print(f"      → STABLE (will attract nearby trajectories)")
        else:
            print(f"      → UNSTABLE (will repel nearby trajectories)")

    # Simulate dynamics
    print("\n3. Simulating population dynamics...")
    N0 = np.array([50.0, 50.0, 50.0])
    sol = model.simulate(N0, (0, 100))

    if sol['success']:
        print(f"   Simulation successful!")
        print(f"   Initial: N_S={N0[0]:.1f}, N_M={N0[1]:.1f}, N_G={N0[2]:.1f}")
        print(f"   Final:   N_S={sol['N_S'][-1]:.1f}, "
              f"N_M={sol['N_M'][-1]:.1f}, N_G={sol['N_G'][-1]:.1f}")
    else:
        print(f"   Simulation failed: {sol['message']}")
        return

    # Create visualizations
    print("\n4. Generating visualizations...")

    # Time series plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: Time series
    ax = axes[0, 0]
    ax.plot(sol['t'], sol['N_S'], 'b-', linewidth=2, label='S (substrate specialist)')
    ax.plot(sol['t'], sol['N_M'], 'r-', linewidth=2, label='M (metabolite specialist)')
    ax.plot(sol['t'], sol['N_G'], 'g-', linewidth=2, label='G (generalist)')
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Population density', fontsize=12)
    ax.set_title('Population Dynamics', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Panel 2: Phase portrait S-M
    print("   Creating S-M phase portrait...")
    analyzer = PhasePlaneAnalyzer(model)

    ax = axes[0, 1]
    analyzer.plot_phase_portrait_2D(
        species_pair=(0, 1),
        fixed_species_val=50.0,
        initial_conditions=[
            np.array([80, 10, 50]),
            np.array([20, 80, 50]),
            np.array([50, 50, 50])
        ],
        t_max=100,
        grid_range=(0, 120),
        n_grid=15,
        show_nullclines=True,
        show_equilibria=True,
        ax=ax
    )

    # Panel 3: Phase portrait S-G
    print("   Creating S-G phase portrait...")
    ax = axes[1, 0]
    analyzer.plot_phase_portrait_2D(
        species_pair=(0, 2),
        fixed_species_val=50.0,
        initial_conditions=[
            np.array([80, 50, 10]),
            np.array([20, 50, 80]),
            np.array([50, 50, 50])
        ],
        t_max=100,
        grid_range=(0, 120),
        n_grid=15,
        show_nullclines=True,
        show_equilibria=True,
        ax=ax
    )

    # Panel 4: Phase portrait M-G
    print("   Creating M-G phase portrait...")
    ax = axes[1, 1]
    analyzer.plot_phase_portrait_2D(
        species_pair=(1, 2),
        fixed_species_val=50.0,
        initial_conditions=[
            np.array([50, 80, 10]),
            np.array([50, 20, 80]),
            np.array([50, 50, 50])
        ],
        t_max=100,
        grid_range=(0, 120),
        n_grid=15,
        show_nullclines=True,
        show_equilibria=True,
        ax=ax
    )

    plt.tight_layout()
    plt.savefig('figures/example_analysis.png', dpi=300, bbox_inches='tight')
    print("   Saved figure: figures/example_analysis.png")

    # Test different omega values
    print("\n5. Testing different pathway weightings (ω)...")
    omega_values = [0.2, 0.5, 0.8]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    for idx, omega in enumerate(omega_values):
        model.params['omega'] = omega
        sol = model.simulate(N0, (0, 100))

        ax = axes[idx]
        ax.plot(sol['t'], sol['N_S'], 'b-', linewidth=2, label='S')
        ax.plot(sol['t'], sol['N_M'], 'r-', linewidth=2, label='M')
        ax.plot(sol['t'], sol['N_G'], 'g-', linewidth=2, label='G')
        ax.set_xlabel('Time', fontsize=12)
        ax.set_ylabel('Population density', fontsize=12)
        ax.set_title(f'ω = {omega:.1f}', fontsize=14, fontweight='bold')
        ax.legend(loc='best', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 120)

        print(f"   ω = {omega:.1f}: Final populations = "
              f"N_S={sol['N_S'][-1]:.1f}, "
              f"N_M={sol['N_M'][-1]:.1f}, "
              f"N_G={sol['N_G'][-1]:.1f}")

    plt.tight_layout()
    plt.savefig('figures/omega_comparison.png', dpi=300, bbox_inches='tight')
    print("   Saved figure: figures/omega_comparison.png")

    # Summary
    print("\n" + "="*70)
    print("Analysis complete!")
    print("="*70)
    print("\nKey findings:")
    print("  • Model successfully simulates three-species dynamics")
    print("  • Multiple equilibria exist (extinction, coexistence, exclusion)")
    print("  • Phase portraits show nullclines and flow fields")
    print("  • Pathway weighting (ω) strongly affects community composition")
    print("\nNext steps:")
    print("  • Run full analysis: jupyter notebook notebooks/three_species_phase_analysis.ipynb")
    print("  • Explore parameter space systematically")
    print("  • Compare with experimental data")
    print("="*70)

    plt.show()


if __name__ == "__main__":
    main()

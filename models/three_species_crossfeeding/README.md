# Three-Species Cross-Feeding Model: Phase Plane Analysis

**Author:** Jian Wang
**Affiliation:** KU Leuven, Department of Chemical Engineering
**Date:** January 2026

---

## Overview

This repository contains a comprehensive mathematical and computational analysis of a three-species microbial community model featuring bidirectional cross-feeding, competition, and metabolic cooperation. The model describes the dynamics of:

1. **S-specialist** (Substrate specialists) - organisms that grow on primary substrate
2. **M-specialist** (Metabolite specialists) - obligate cross-feeders that depend on metabolic byproducts
3. **G-specialist** (Generalists) - flexible strategists that can utilize both pathways with tunable weighting

The analysis follows the rigorous phase plane methodology pioneered by Jeff Gore's laboratory at MIT, combining equilibrium analysis, stability theory, nullcline geometry, and bifurcation analysis to understand community assembly and coexistence.

---

## Model Equations

The dynamics are governed by the following system of ordinary differential equations:

### S-specialist (Substrate specialist)
```
dN_S/dt = r_S · N_S [1 + σ_SM · N_M/K_M + (1-ω) · σ_SG · N_G/K_G
                     - ω · α_SG · N_G/K_G - N_S/K_S]
```

**Interpretation:**
- Base growth rate: 1 (can grow independently on substrate)
- Benefits from M-specialist's presence (cooperation: σ_SM)
- Interaction with generalist depends on ω (cooperation when low ω, competition when high ω)
- Self-limitation at carrying capacity K_S

### M-specialist (Metabolite specialist)
```
dN_M/dt = r_M · N_M [-1 + σ_MS · N_S/K_S + ω · σ_MG · N_G/K_G
                     - (1-ω) · α_MG · N_G/K_G - N_M/K_M]
```

**Interpretation:**
- Base growth rate: -1 (obligate cross-feeder, cannot survive alone)
- Requires S-specialist for survival (cooperation: σ_MS)
- Interaction with generalist depends on ω
- Self-limitation at carrying capacity K_M

### G-specialist (Generalist)
```
dN_G/dt = r_G · N_G [ω(1 - α_GS · N_S/K_S + σ_GM · N_M/K_M)
                     + (1-ω)(-1 - α_GM · N_M/K_M + σ_GS · N_S/K_S)
                     - N_G/K_G]
```

**Interpretation:**
- **Substrate pathway** (weight ω): can grow independently, competes with S, benefits from M
- **Metabolite pathway** (weight 1-ω): obligate cross-feeder, competes with M, benefits from S
- ω is the **niche parameter** determining metabolic strategy

---

## Parameter Definitions

| Parameter | Description | Units | Typical Range |
|-----------|-------------|-------|---------------|
| r_S, r_M, r_G | Intrinsic growth rates | 1/time | 0.5 - 2.0 |
| K_S, K_M, K_G | Carrying capacities | cells/volume | 50 - 200 |
| σ_ij | Synergistic coefficients (cooperation) | dimensionless | 0.2 - 0.8 |
| α_ij | Competition coefficients | dimensionless | 0.2 - 0.6 |
| ω | Pathway weighting parameter | dimensionless | 0 - 1 |

**Key relationships:**
- σ_ij > 0: Species i benefits from species j (cross-feeding, facilitation)
- α_ij > 0: Species i is harmed by species j (competition, interference)
- ω = 0: Generalist uses only metabolite pathway (like M-specialist)
- ω = 1: Generalist uses only substrate pathway (like S-specialist)
- 0 < ω < 1: True generalist strategy

---

## Repository Structure

```
three_species_crossfeeding/
├── src/
│   ├── three_species_model.py          # Core model class with ODE system
│   └── phase_plane_analysis.py         # Phase plane tools and visualizations
├── notebooks/
│   └── three_species_phase_analysis.ipynb  # Comprehensive analysis notebook
├── figures/                             # Generated plots (created by notebook)
│   ├── timeseries_dynamics.png
│   ├── phase_plane_SM.png
│   ├── phase_plane_SG.png
│   ├── phase_plane_MG.png
│   ├── phase_space_3D.png
│   ├── bifurcation_omega.png
│   └── coexistence_map.png
├── data/                                # Data storage (for parameter sweeps)
├── requirements.txt                     # Python dependencies
└── README.md                            # This file
```

---

## Installation and Setup

### Prerequisites
- Python 3.8+
- Jupyter Notebook or JupyterLab

### Dependencies
```bash
pip install -r requirements.txt
```

Required packages:
- numpy
- scipy
- matplotlib
- seaborn
- jupyter

### Quick Start

1. Clone or navigate to this directory
2. Install dependencies: `pip install -r requirements.txt`
3. Launch Jupyter: `jupyter notebook notebooks/three_species_phase_analysis.ipynb`
4. Run all cells to reproduce the analysis

---

## Key Features

### 1. Model Implementation (`three_species_model.py`)

**ThreeSpeciesModel class:**
- `equations()`: ODE right-hand side for integration
- `simulate()`: Time-series integration using scipy.solve_ivp
- `find_equilibria()`: Numerical root-finding for all equilibrium points
- `jacobian()`: Analytical Jacobian matrix for stability analysis
- `stability_analysis()`: Eigenvalue decomposition and equilibrium classification
- `classify_equilibrium_ecology()`: Biological interpretation of equilibria

### 2. Phase Plane Analysis (`phase_plane_analysis.py`)

**PhasePlaneAnalyzer class:**
- `compute_nullclines_2D()`: Calculate nullclines for 2D projections
- `plot_phase_portrait_2D()`: Create phase portraits with vector fields
- `plot_3D_phase_space()`: Visualize trajectories in full 3D phase space
- `plot_timeseries()`: Population dynamics over time
- `bifurcation_analysis_omega()`: Explore parameter space
- `plot_bifurcation_diagram()`: Visualize bifurcations

### 3. Comprehensive Notebook (`three_species_phase_analysis.ipynb`)

Organized analysis following Gore lab style:
1. Model initialization and parameter exploration
2. Equilibrium finding and stability classification
3. Time series dynamics from diverse initial conditions
4. 2D phase portraits for all species pairs
5. 3D phase space visualization
6. Bifurcation analysis for pathway weighting (ω)
7. Sensitivity analysis for cooperation-competition balance
8. Biological interpretation and ecological insights

---

## Analytical Results

### Equilibrium Classification

The system exhibits multiple types of equilibria:

1. **Extinction equilibrium**: (0, 0, 0) - typically unstable
2. **Single-species equilibria**:
   - S-only: (K_S, 0, 0) - stable if generalist cannot invade
   - M-only: Not feasible (M is obligate cross-feeder)
   - G-only: (0, 0, K_G) - feasible only if ω > threshold
3. **Two-species equilibria**:
   - S-M coexistence: Cooperative cross-feeding
   - S-G coexistence: Depends on ω
   - M-G coexistence: Depends on ω
4. **Three-species coexistence**: Requires balanced parameters

### Coexistence Conditions

**Necessary conditions for three-species coexistence:**
1. Cooperation strength > Competition strength (σ > α)
2. Intermediate pathway weighting (ω ∈ [ω_min, ω_max])
3. M-specialist receives sufficient cross-feeding from S
4. Generalist's niche is sufficiently differentiated

**Stability requirements:**
- All eigenvalues of Jacobian must have negative real parts
- Local stability does not guarantee global stability
- Multiple stable equilibria can coexist (multistability)

---

## Biological Interpretations

### Ecological Insights

1. **Obligate vs. Facultative Cross-Feeding**:
   - M-specialist: Obligate cross-feeder (base growth = -1)
   - S-specialist: Facultative partner (can survive alone)
   - G-specialist: Metabolic flexibility enables coexistence

2. **Niche Differentiation**:
   - ω determines the generalist's metabolic niche
   - Intermediate ω creates sufficient differentiation
   - Extreme ω leads to competitive exclusion

3. **Cooperation-Competition Trade-off**:
   - Cross-feeding creates positive interactions
   - Resource competition creates negative interactions
   - Balance determines community assembly

4. **Priority Effects**:
   - Initial conditions can determine final community composition
   - Multistability enables alternative stable states
   - Relevant for community assembly and invasion resistance

### Experimental Connections

This model framework applies to:

- **Syntrophic partnerships** in anaerobic digestion (e.g., Syntrophomonas + Methanobacterium)
- **Cross-feeding in biofilms** (spatial structure + metabolic exchange)
- **Engineered consortia** for bioproduction (designed cooperation)
- **Natural soil communities** (carbon metabolism networks)
- **Gut microbiota** (primary degraders + secondary consumers + generalists)

---

## Example Usage

### Basic Simulation

```python
from three_species_model import ThreeSpeciesModel
import numpy as np

# Initialize model with default parameters
model = ThreeSpeciesModel()

# Define initial condition
N0 = np.array([50.0, 50.0, 50.0])  # Equal initial densities

# Simulate for 100 time units
sol = model.simulate(N0, (0, 100))

# Plot results
import matplotlib.pyplot as plt
plt.plot(sol['t'], sol['N_S'], label='S-specialist')
plt.plot(sol['t'], sol['N_M'], label='M-specialist')
plt.plot(sol['t'], sol['N_G'], label='Generalist')
plt.xlabel('Time')
plt.ylabel('Population density')
plt.legend()
plt.show()
```

### Equilibrium Analysis

```python
# Find all equilibria
equilibria = model.find_equilibria(n_attempts=100)

# Analyze each equilibrium
for eq in equilibria:
    stability = model.stability_analysis(eq)
    eco_type = model.classify_equilibrium_ecology(eq)

    print(f"Equilibrium: {eq}")
    print(f"Type: {eco_type}")
    print(f"Stable: {stability['stable']}")
    print(f"Eigenvalues: {stability['eigenvalues']}\n")
```

### Phase Portrait

```python
from phase_plane_analysis import PhasePlaneAnalyzer

# Create analyzer
analyzer = PhasePlaneAnalyzer(model)

# Plot S-M phase plane
fig = analyzer.plot_phase_portrait_2D(
    species_pair=(0, 1),  # S vs M
    fixed_species_val=50.0,  # G = 50
    initial_conditions=[
        np.array([80, 10, 50]),
        np.array([10, 80, 50])
    ],
    show_nullclines=True,
    show_equilibria=True
)
plt.show()
```

---

## Parameter Sensitivity

Key parameters that strongly influence coexistence:

1. **ω (pathway weighting)**: Most sensitive parameter
   - Controls generalist's metabolic strategy
   - Intermediate values promote coexistence
   - Bifurcation parameter

2. **σ_MS (M benefits from S)**: Critical for M survival
   - Must exceed threshold for M-specialist viability
   - Determines minimum S density needed

3. **Competition coefficients (α)**: Determine exclusion strength
   - High α → competitive exclusion
   - Balance with σ determines coexistence

4. **Growth rates (r_i)**: Set timescale
   - Relative values affect dominance
   - Fast growers can outcompete cooperators

---

## Future Extensions

### Theoretical
1. **Spatial structure**: Reaction-diffusion PDEs or agent-based models
2. **Stochastic dynamics**: Demographic noise, extinction risk
3. **Evolutionary dynamics**: Adaptive ω, mutation-selection balance
4. **Environmental fluctuations**: Variable substrate, temperature
5. **N-species generalization**: Metabolic network models

### Computational
1. **Parameter inference**: Fit to experimental data (Bayesian methods)
2. **Sensitivity analysis**: Sobol indices, variance decomposition
3. **Optimal control**: Engineering coexistence, maximizing productivity
4. **Machine learning**: Predict coexistence from parameters

### Experimental
1. **Validate with chemostat experiments**
2. **Measure interaction coefficients (σ, α)**
3. **Test bifurcation predictions**
4. **Engineer synthetic consortia with tunable ω**

---

## References

### Theoretical Foundations
- **Gore et al. (2009)** "Snowdrift game dynamics and facultative cheating in yeast" *Nature* 459:253-256
- **Momeni et al. (2013)** "Using artificial systems to explore the ecology and evolution of symbioses" *Cellular and Molecular Life Sciences* 70:1933-1948
- **Goldford et al. (2018)** "Emergent simplicity in microbial community assembly" *Science* 361:469-474

### Cross-Feeding Ecology
- **Estrela et al. (2022)** "Metabolic rules of microbial community assembly" *Nature Ecology & Evolution* 6:174-182
- **Pande et al. (2014)** "Fitness and stability of obligate cross-feeding interactions that emerge upon gene loss in bacteria" *ISME Journal* 8:953-962
- **Hoek et al. (2016)** "Resource availability modulates the cooperative and competitive nature of a microbial cross-feeding mutualism" *PLoS Biology* 14:e1002540

### Phase Plane Methods
- **Strogatz (2015)** *Nonlinear Dynamics and Chaos* (2nd ed.), Westview Press
- **Murray (2002)** *Mathematical Biology I: An Introduction* (3rd ed.), Springer
- **Bazykin (1998)** *Nonlinear Dynamics of Interacting Populations*, World Scientific

---

## Citation

If you use this model or code in your research, please cite:

```bibtex
@software{wang2026threespp,
  author = {Wang, Jian},
  title = {Three-Species Cross-Feeding Model: Phase Plane Analysis},
  year = {2026},
  url = {https://github.com/JianWang1123/mathbiojw},
  note = {Mathematical biology model for microbial community dynamics}
}
```

---

## Contact

**Jian Wang**
PhD Student, Mathematical Biology
KU Leuven, Department of Chemical Engineering
Email: [your-email@kuleuven.be]
Website: [https://JianWang1123.github.io/mathbiojw](https://JianWang1123.github.io/mathbiojw)

---

## License

This code is released under the MIT License. See LICENSE file for details.

---

## Acknowledgments

This work was inspired by the rigorous ecological modeling approaches developed in the laboratories of Jeff Gore (MIT), Alvaro Sanchez (Yale), and Wenying Shou (UCL). The mathematical framework builds on classic ecological theory while incorporating modern insights from microbial metabolic networks and synthetic ecology.

Special thanks to the mathematical biology and systems biology communities for developing the theoretical and computational tools that made this analysis possible.

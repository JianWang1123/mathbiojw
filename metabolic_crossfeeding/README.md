# Metabolic Cross-Feeding Model: Division of Labor

This repository contains a complete analytical and numerical analysis demonstrating that **metabolic cross-feeding favors the evolution of division of labor** in microbial communities.

## Key Result

In a two-species, two-amino-acid cross-feeding system:
- The **symmetric generalist equilibrium** (both species produce both amino acids equally) is evolutionarily **unstable**
- The **division of labor equilibrium** (each species specializes on one amino acid) is the **evolutionarily stable strategy (ESS)**

## Files

### Core Model
- `crossfeeding_model.py` - Main model implementation with:
  - `AnalyticalCrossFeedingModel`: Simplified model for analytical tractability
  - `FullChemostatModel`: Complete bioreactor model with explicit resource dynamics
  - Analytical derivations and game-theoretic analysis

### Visualization
- `visualizations.py` - Publication-quality figure generation
- `figures/` - Generated figures (PNG and PDF)

### Mathematical Derivation
- `division_of_labor_derivation.tex` - Complete LaTeX document with full analytical proofs

## Key Equations

### Fitness Function
```
W_A = (1 - f_A1 - f_A2) * g(P1, P2) - D
```
where:
- `f_A1, f_A2`: Investment fractions in amino acids F1, F2
- `g(P1, P2) = γ * P1^α * P2^α`: Growth rate
- `P1 = f_A1 + f_B1`, `P2 = f_A2 + f_B2`: Total amino acid production
- `D`: Dilution rate

### Selection Gradient
```
∂W_A/∂f_A1 = g * [α(1 - f_A1 - f_A2)/P1 - 1]
```

### Equilibria
- **Symmetric**: `f* = α / 2(α+1)` → For α=1: `f* = 1/4`
- **Division of Labor**: `f*_div = α / (α+1)` → For α=1: `f* = 1/2`

## Figures

1. **Model Schematic** - Bioreactor setup with cross-feeding interactions
2. **Fitness Landscape** - Direct costs vs indirect benefits of investment
3. **Phase Portrait** - Evolutionary trajectories toward division of labor
4. **Stability Analysis** - Eigenvalues and bifurcation diagram
5. **Game Theory** - Payoff matrix and Nash equilibria
6. **Eco-Evo Dynamics** - Full simulation of population and strategy evolution
7. **Summary** - Why division of labor evolves

## Usage

```bash
# Generate all figures
python visualizations.py

# Run analytical derivations
python crossfeeding_model.py
```

## Requirements

- numpy
- scipy
- matplotlib

## Citation

If you use this model, please cite:
```
Wang, J. (2024). Eco-evolutionary Dynamics of Metabolic Specialization:
From Symmetric Public Goods to Obligate Cross-Feeding Mutualism.
```

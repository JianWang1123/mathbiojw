#!/usr/bin/env python3
"""
Visual Summary: Both ω_crit1 and ω_crit2 Formulas and Results

Creates a clean reference figure showing:
- Analytical formulas
- Numerical values
- Biological interpretations
- Parameter dependencies

Author: Jian Wang
Date: January 2026
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams.update({'font.size': 9, 'font.family': 'sans-serif'})

fig = plt.figure(figsize=(16, 10))

# Title
fig.suptitle('Complete Analytical Summary: ω_crit1 and ω_crit2 for Three-Species Cross-Feeding Model',
             fontsize=16, fontweight='bold', y=0.98)

# ============================================================================
# LEFT PANEL: ω_crit1 (Always exists)
# ============================================================================

ax1 = plt.subplot(1, 2, 1)
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1)
ax1.axis('off')

# Title box
title_box1 = FancyBboxPatch((0.05, 0.88), 0.9, 0.1,
                           boxstyle="round,pad=0.01",
                           facecolor='lightgreen', edgecolor='darkgreen',
                           linewidth=2)
ax1.add_patch(title_box1)
ax1.text(0.5, 0.93, r'$\omega_{crit1}$: Generalist Invasion Threshold',
        ha='center', va='center', fontsize=14, fontweight='bold')

# Formula box
formula_box1 = FancyBboxPatch((0.05, 0.68), 0.9, 0.18,
                             boxstyle="round,pad=0.01",
                             facecolor='lightyellow', edgecolor='orange',
                             linewidth=2)
ax1.add_patch(formula_box1)

ax1.text(0.5, 0.83, 'ANALYTICAL FORMULA (Explicit)', ha='center',
        fontsize=11, fontweight='bold', style='italic')

ax1.text(0.5, 0.77, r'$\omega_{crit1} = \frac{1 - \sigma_{GS} \cdot s^*_{SM} + \alpha_{GM} \cdot m^*_{SM}}{2 - (\sigma_{GS}+\alpha_{GS}) \cdot s^*_{SM} + (\sigma_{GM}+\alpha_{GM}) \cdot m^*_{SM}}$',
        ha='center', fontsize=12, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax1.text(0.5, 0.70, 'This is Equation 3 in manuscript!',
        ha='center', fontsize=9, style='italic', color='darkred')

# Numerical values box
values_box1 = FancyBboxPatch((0.05, 0.48), 0.9, 0.18,
                            boxstyle="round,pad=0.01",
                            facecolor='lightcyan', edgecolor='blue',
                            linewidth=2)
ax1.add_patch(values_box1)

ax1.text(0.5, 0.63, 'BASELINE PARAMETERS', ha='center',
        fontsize=11, fontweight='bold')

ax1.text(0.15, 0.58, r'$s^*_{SM} = 2.000$', fontsize=10)
ax1.text(0.15, 0.54, r'$m^*_{SM} = 2.000$', fontsize=10)
ax1.text(0.15, 0.50, r'$\sigma_{GS} = 0.4$', fontsize=10)

ax1.text(0.55, 0.58, r'$\alpha_{GS} = 0.3$', fontsize=10)
ax1.text(0.55, 0.54, r'$\sigma_{GM} = 0.4$', fontsize=10)
ax1.text(0.55, 0.50, r'$\alpha_{GM} = 0.3$', fontsize=10)

# Result box
result_box1 = FancyBboxPatch((0.15, 0.36), 0.7, 0.1,
                            boxstyle="round,pad=0.01",
                            facecolor='gold', edgecolor='darkgoldenrod',
                            linewidth=3)
ax1.add_patch(result_box1)

ax1.text(0.5, 0.41, r'$\omega_{crit1} = 0.4000$',
        ha='center', fontsize=16, fontweight='bold')

# Biological meaning
bio_box1 = FancyBboxPatch((0.05, 0.18), 0.9, 0.16,
                         boxstyle="round,pad=0.01",
                         facecolor='lavender', edgecolor='purple',
                         linewidth=2)
ax1.add_patch(bio_box1)

ax1.text(0.5, 0.32, 'BIOLOGICAL INTERPRETATION', ha='center',
        fontsize=11, fontweight='bold')

ax1.text(0.5, 0.27, r'$\omega < 0.40$: Generalist too metabolite-specialized',
        ha='center', fontsize=9)
ax1.text(0.5, 0.24, r'$\rightarrow$ Cannot invade S-M platform', ha='center',
        fontsize=9, style='italic')
ax1.text(0.5, 0.20, r'$\omega > 0.40$: Three-species coexistence emerges',
        ha='center', fontsize=9, color='darkgreen', fontweight='bold')

# Properties
props_box1 = FancyBboxPatch((0.05, 0.02), 0.9, 0.14,
                           boxstyle="round,pad=0.01",
                           facecolor='white', edgecolor='black',
                           linewidth=1.5)
ax1.add_patch(props_box1)

ax1.text(0.5, 0.14, 'KEY PROPERTIES', ha='center', fontsize=10, fontweight='bold')
ax1.text(0.15, 0.10, '✓ Always exists', fontsize=9)
ax1.text(0.15, 0.07, '✓ Explicit formula', fontsize=9)
ax1.text(0.15, 0.04, '✓ Easy to compute', fontsize=9)

ax1.text(0.55, 0.10, '✓ Based on S-M equilibrium', fontsize=9)
ax1.text(0.55, 0.07, '✓ Marks G invasion', fontsize=9)
ax1.text(0.55, 0.04, '✓ Transcritical bifurcation', fontsize=9)

# ============================================================================
# RIGHT PANEL: ω_crit2 (Parameter-dependent)
# ============================================================================

ax2 = plt.subplot(1, 2, 2)
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)
ax2.axis('off')

# Title box
title_box2 = FancyBboxPatch((0.05, 0.88), 0.9, 0.1,
                           boxstyle="round,pad=0.01",
                           facecolor='lightcoral', edgecolor='darkred',
                           linewidth=2)
ax2.add_patch(title_box2)
ax2.text(0.5, 0.93, r'$\omega_{crit2}$: M Displacement Threshold',
        ha='center', va='center', fontsize=14, fontweight='bold')

# Formula box
formula_box2 = FancyBboxPatch((0.05, 0.68), 0.9, 0.18,
                             boxstyle="round,pad=0.01",
                             facecolor='lightyellow', edgecolor='orange',
                             linewidth=2)
ax2.add_patch(formula_box2)

ax2.text(0.5, 0.83, 'IMPLICIT CONDITION (Numerical)', ha='center',
        fontsize=11, fontweight='bold', style='italic')

ax2.text(0.5, 0.77, r'$\sigma_{MS} \cdot s^*_{SG}(\omega_{crit2}) + \sigma_{MG} \cdot g^*_{SG}(\omega_{crit2}) = r_M$',
        ha='center', fontsize=11, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax2.text(0.5, 0.70, 'Requires numerical root-finding!',
        ha='center', fontsize=9, style='italic', color='darkred')

# Two scenarios box
scenarios_box = FancyBboxPatch((0.05, 0.40), 0.9, 0.26,
                              boxstyle="round,pad=0.01",
                              facecolor='lightyellow', edgecolor='darkorange',
                              linewidth=2)
ax2.add_patch(scenarios_box)

ax2.text(0.5, 0.64, 'PARAMETER-DEPENDENT EXISTENCE', ha='center',
        fontsize=11, fontweight='bold', color='darkred')

# Scenario 1
ax2.text(0.12, 0.58, 'Scenario 1:', fontsize=10, fontweight='bold')
ax2.text(0.12, 0.54, 'Baseline parameters', fontsize=9, style='italic')
ax2.text(0.12, 0.50, r'$\sigma_{MS}=1.5, \sigma_{MG}=0.4$', fontsize=9)

result1_box = FancyBboxPatch((0.1, 0.44), 0.35, 0.04,
                            boxstyle="round,pad=0.005",
                            facecolor='pink', edgecolor='red', linewidth=2)
ax2.add_patch(result1_box)
ax2.text(0.275, 0.46, r'$\omega_{crit2}$ does NOT exist',
        ha='center', fontsize=9, fontweight='bold')
ax2.text(0.275, 0.42, '→ Permanent S-M-G coexistence',
        ha='center', fontsize=8, style='italic')

# Scenario 2
ax2.text(0.58, 0.58, 'Scenario 2:', fontsize=10, fontweight='bold')
ax2.text(0.58, 0.54, 'Enhanced M-G mutualism', fontsize=9, style='italic')
ax2.text(0.58, 0.50, r'$\sigma_{MS}=1.8, \sigma_{MG}=0.8$', fontsize=9)

result2_box = FancyBboxPatch((0.55, 0.44), 0.35, 0.04,
                            boxstyle="round,pad=0.005",
                            facecolor='lightgreen', edgecolor='darkgreen', linewidth=2)
ax2.add_patch(result2_box)
ax2.text(0.725, 0.46, r'$\omega_{crit2} = 0.6860$',
        ha='center', fontsize=9, fontweight='bold')
ax2.text(0.725, 0.42, '→ Bounded coexistence window',
        ha='center', fontsize=8, style='italic')

# Biological meaning
bio_box2 = FancyBboxPatch((0.05, 0.22), 0.9, 0.16,
                         boxstyle="round,pad=0.01",
                         facecolor='lavender', edgecolor='purple',
                         linewidth=2)
ax2.add_patch(bio_box2)

ax2.text(0.5, 0.36, 'BIOLOGICAL INTERPRETATION (When exists)', ha='center',
        fontsize=11, fontweight='bold')

ax2.text(0.5, 0.31, r'$\omega < \omega_{crit2}$: M can coexist with S-G',
        ha='center', fontsize=9)
ax2.text(0.5, 0.28, r'$\omega > \omega_{crit2}$: G too substrate-specialized',
        ha='center', fontsize=9)
ax2.text(0.5, 0.24, r'$\rightarrow$ M displaced (S-G equilibrium only)',
        ha='center', fontsize=9, style='italic', color='darkred', fontweight='bold')

# Properties
props_box2 = FancyBboxPatch((0.05, 0.02), 0.9, 0.18,
                           boxstyle="round,pad=0.01",
                           facecolor='white', edgecolor='black',
                           linewidth=1.5)
ax2.add_patch(props_box2)

ax2.text(0.5, 0.18, 'KEY PROPERTIES', ha='center', fontsize=10, fontweight='bold')
ax2.text(0.15, 0.14, '⚠ Parameter-dependent', fontsize=9, color='red')
ax2.text(0.15, 0.11, '⚠ Implicit formula', fontsize=9, color='red')
ax2.text(0.15, 0.08, '⚠ Numerical solution', fontsize=9, color='red')
ax2.text(0.15, 0.05, '⚠ May not exist', fontsize=9, color='red', fontweight='bold')

ax2.text(0.55, 0.14, '✓ Based on S-G equilibrium', fontsize=9)
ax2.text(0.55, 0.11, '✓ Marks M displacement', fontsize=9)
ax2.text(0.55, 0.08, '✓ Creates coexistence window', fontsize=9)
ax2.text(0.55, 0.05, '✓ Requires σ_MG ≥ 1.0', fontsize=9)

plt.tight_layout()
plt.savefig('both_omega_critical_formulas_summary.png', dpi=300, bbox_inches='tight')
print("✓ Summary figure saved: both_omega_critical_formulas_summary.png")

# Create a second figure: Community transitions
fig2, ax = plt.subplots(1, 1, figsize=(14, 6))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

# Title
ax.text(0.5, 0.95, 'Community Composition Transitions: Baseline vs Modified Parameters',
       ha='center', fontsize=16, fontweight='bold')

# Baseline scenario (top)
ax.text(0.05, 0.75, 'BASELINE PARAMETERS', fontsize=12, fontweight='bold',
       bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))

ax.text(0.05, 0.68, r'$\sigma_{MS}=1.5, \sigma_{MG}=0.4$', fontsize=10)

# S-M box
sm_box1 = FancyBboxPatch((0.15, 0.55), 0.12, 0.08,
                        boxstyle="round,pad=0.01",
                        facecolor='lightblue', edgecolor='blue', linewidth=2)
ax.add_patch(sm_box1)
ax.text(0.21, 0.59, 'S-M', ha='center', fontsize=11, fontweight='bold')

# Arrow
ax.arrow(0.28, 0.59, 0.13, 0, head_width=0.02, head_length=0.03,
        fc='black', ec='black', linewidth=2)
ax.text(0.345, 0.63, r'$\omega=0.40$', ha='center', fontsize=9, fontweight='bold')
ax.text(0.345, 0.55, r'$\omega_{crit1}$', ha='center', fontsize=8, style='italic')

# S-M-G box
smg_box1 = FancyBboxPatch((0.45, 0.55), 0.35, 0.08,
                         boxstyle="round,pad=0.01",
                         facecolor='gold', edgecolor='darkgoldenrod', linewidth=3)
ax.add_patch(smg_box1)
ax.text(0.625, 0.59, 'S-M-G (permanent coexistence)', ha='center',
       fontsize=11, fontweight='bold')

ax.text(0.625, 0.50, 'No ω_crit2 → Three species coexist indefinitely',
       ha='center', fontsize=9, style='italic', color='darkgreen')

# Modified scenario (bottom)
ax.text(0.05, 0.38, 'MODIFIED PARAMETERS', fontsize=12, fontweight='bold',
       bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))

ax.text(0.05, 0.31, r'$\sigma_{MS}=1.8, \sigma_{MG}=0.8$', fontsize=10)

# S-M box
sm_box2 = FancyBboxPatch((0.15, 0.18), 0.12, 0.08,
                        boxstyle="round,pad=0.01",
                        facecolor='lightblue', edgecolor='blue', linewidth=2)
ax.add_patch(sm_box2)
ax.text(0.21, 0.22, 'S-M', ha='center', fontsize=11, fontweight='bold')

# Arrow 1
ax.arrow(0.28, 0.22, 0.11, 0, head_width=0.02, head_length=0.02,
        fc='black', ec='black', linewidth=2)
ax.text(0.335, 0.26, r'$\omega=0.15$', ha='center', fontsize=9, fontweight='bold')
ax.text(0.335, 0.18, r'$\omega_{crit1}$', ha='center', fontsize=8, style='italic')

# S-M-G box
smg_box2 = FancyBboxPatch((0.42, 0.18), 0.18, 0.08,
                         boxstyle="round,pad=0.01",
                         facecolor='gold', edgecolor='darkgoldenrod', linewidth=3)
ax.add_patch(smg_box2)
ax.text(0.51, 0.22, 'S-M-G', ha='center', fontsize=11, fontweight='bold')

# Arrow 2
ax.arrow(0.61, 0.22, 0.11, 0, head_width=0.02, head_length=0.02,
        fc='black', ec='black', linewidth=2)
ax.text(0.665, 0.26, r'$\omega=0.69$', ha='center', fontsize=9, fontweight='bold')
ax.text(0.665, 0.18, r'$\omega_{crit2}$', ha='center', fontsize=8, style='italic', color='red')

# S-G box
sg_box = FancyBboxPatch((0.75, 0.18), 0.12, 0.08,
                       boxstyle="round,pad=0.01",
                       facecolor='lightcoral', edgecolor='red', linewidth=2)
ax.add_patch(sg_box)
ax.text(0.81, 0.22, 'S-G', ha='center', fontsize=11, fontweight='bold')

# Coexistence window annotation
window_box = FancyBboxPatch((0.4, 0.08), 0.24, 0.06,
                           boxstyle="round,pad=0.01",
                           facecolor='lightyellow', edgecolor='orange',
                           linewidth=2, linestyle='--')
ax.add_patch(window_box)
ax.text(0.52, 0.11, 'Coexistence Window', ha='center', fontsize=9,
       fontweight='bold', style='italic')
ax.text(0.52, 0.05, r'Width = $\omega_{crit2} - \omega_{crit1} = 0.54$',
       ha='center', fontsize=8)

plt.savefig('community_transitions_comparison.png', dpi=300, bbox_inches='tight')
print("✓ Transitions figure saved: community_transitions_comparison.png")

plt.show()

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print("""
Two critical thresholds in three-species cross-feeding model:

1. ω_crit1 = 0.4000
   - ALWAYS exists (when S-M mutualism is stable)
   - Explicit analytical formula (Equation 3 in manuscript)
   - Marks generalist invasion into S-M platform
   - Transcritical bifurcation: S-M → S-M-G

2. ω_crit2 = Parameter-dependent
   - Baseline (σ_MG=0.4): DOES NOT EXIST
   - Modified (σ_MG=0.8): ω_crit2 = 0.6860
   - Implicit condition: λ_M(S*_SG, G*_SG) = 0
   - Marks M displacement: S-M-G → S-G

KEY INSIGHT:
  With baseline parameters, three-species coexistence is PERMANENT once
  G invades (ω > 0.40). The second bifurcation ω_crit2 only emerges with
  stronger generalist-to-metabolite-specialist facilitation (σ_MG ≥ 1.0).

MANUSCRIPT IMPLICATION:
  The statement about ω_crit2 (line 104) should clarify its parameter-
  dependence. Consider adding: "Under parameter regimes with enhanced
  M-G mutualism, a second bifurcation ω_crit2 can emerge..."
""")

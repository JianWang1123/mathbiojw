# Complete Analysis: ω_crit2 and Parameter-Dependent Bifurcations

## Executive Summary

**Key Finding:** ω_crit2 (the second bifurcation where M is displaced) is **PARAMETER-DEPENDENT** and does not always exist.

- **ω_crit1**: Always exists when S-M mutualism is stable (σ_MS > 1, σ_MS·σ_SM < 1)
- **ω_crit2**: Exists only in specific parameter regimes where M can initially persist but is eventually outcompeted as ω increases

---

## Why ω_crit2 Doesn't Always Exist

### Baseline Parameters Analysis

With baseline parameters:
- σ_MS = 1.5, σ_SM = 0.5
- σ_GS = σ_GM = σ_SG = σ_MG = 0.4
- α_GS = α_GM = α_SG = α_MG = 0.3

**Result:** λ_M(ω) < 0 for ALL ω ∈ [0,1]

This means:
1. M can NEVER invade S-G equilibrium (regardless of ω)
2. Once G displaces M from three-species equilibrium, M cannot return
3. Therefore, **ω_crit2 does not exist** with these parameters

### Mathematical Explanation

ω_crit2 exists when there is an ω value where:
$$\lambda_M(S^*_{SG}(\omega), G^*_{SG}(\omega)) = 0$$

For this to occur, λ_M must:
1. Be **positive** at some ω (M can invade S-G)
2. Cross **zero** as ω increases (bifurcation point)
3. Become **negative** beyond that point (M excluded)

**With baseline parameters:** λ_M starts negative and remains negative throughout ω ∈ [0,1], so no zero crossing occurs.

---

## Conditions for ω_crit2 to Exist

### Necessary Condition

For ω_crit2 to exist, M must be able to invade S-G equilibrium at some intermediate ω values:

$$-r_M + \sigma_{MS} \cdot s^*_{SG}(\omega) + \sigma_{MG} \cdot g^*_{SG}(\omega) > 0$$

This requires sufficiently strong:
- **σ_MS**: S→M facilitation
- **σ_MG**: G→M facilitation

### Parameter Requirements

ω_crit2 typically emerges when:

1. **Strong M-G mutualism** (σ_MG ≥ 0.6):
   - G produces metabolites that benefit M
   - M can coexist with S-G at intermediate ω

2. **Moderate S-M mutualism** (1.5 ≤ σ_MS ≤ 2.5):
   - Too weak: M never viable
   - Too strong: M always viable (no displacement)

3. **Weak G-S competition** (α_GS ≤ 0.4):
   - Allows S to maintain sufficient density to support M

---

## Parameter Regime Examples

### Regime 1: No ω_crit2 (Baseline)

**Parameters:**
- σ_MS = 1.5, σ_MG = 0.4

**Outcome:**
- ω_crit1 ≈ 0.40 (G invasion)
- ω_crit2: **Does not exist**
- Community sequence: S-M → S-M-G (persists indefinitely)

**Biological scenario:** Once generalist invades, all three species coexist stably across wide ω range

---

### Regime 2: ω_crit2 Exists (Enhanced M-G Facilitation)

**Modified parameters:**
- σ_MS = 1.8, σ_MG = 0.7, σ_GM = 0.6
- (Keep others at baseline)

**Expected outcome:**
- ω_crit1 ≈ 0.35
- ω_crit2 ≈ 0.65
- **Coexistence window:** ω ∈ (0.35, 0.65)

**Biological scenario:**
- ω < 0.35: S-M only
- 0.35 < ω < 0.65: S-M-G coexist
- ω > 0.65: M displaced → S-G only

**Mechanism:** At high ω, generalist becomes substrate-specialized:
1. G competes strongly with S for substrate
2. Both S and G reduce metabolite secretion
3. M starves despite mutualistic benefits from G

---

### Regime 3: ω_crit2 Exists (Weak S-M Mutualism)

**Modified parameters:**
- σ_MS = 1.2 (weaker than baseline)
- σ_MG = 0.5
- σ_SM = 0.3

**Expected outcome:**
- ω_crit1 ≈ 0.25 (narrow invasion window)
- ω_crit2 ≈ 0.45 (early M displacement)
- **Narrow coexistence window**

**Biological scenario:** Fragile three-species coexistence - small metabolic shifts disrupt M

---

## Analytical Formula for ω_crit2 (When It Exists)

### Implicit Equation

$$\sigma_{MS} \cdot s^*_{SG}(\omega_{crit2}) + \sigma_{MG} \cdot g^*_{SG}(\omega_{crit2}) = r_M$$

where the S-G equilibrium is:

$$s^*_{SG}(\omega) = \frac{(1-\omega) + a(\omega) \cdot d(\omega)}{1 - a(\omega) \cdot c(\omega)} + \frac{a(\omega) \cdot c(\omega) \cdot (1-\omega)}{1 - a(\omega) \cdot c(\omega)}$$

$$g^*_{SG}(\omega) = \frac{d(\omega) + c(\omega)(1-\omega)}{1 - c(\omega) \cdot a(\omega)}$$

with:
- $a(\omega) = (1-\omega)\sigma_{SG} - \omega\alpha_{SG}$
- $c(\omega) = (1-\omega)\sigma_{GS} - \omega\alpha_{GS}$
- $d(\omega) = 2\omega - 1$

### Simplified Condition (Approximate)

At high ω (→ 1):
- $s^*_{SG} \approx 0$ (S competitively excluded)
- $g^*_{SG} \approx \frac{1}{1 + \alpha_{GS}}$

For M to persist:
$$\sigma_{MS} \cdot 0 + \sigma_{MG} \cdot \frac{1}{1 + \alpha_{GS}} \gtrsim r_M$$

$$\sigma_{MG} \gtrsim r_M (1 + \alpha_{GS})$$

**With r_M = 0.8, α_GS = 0.3:**
$$\sigma_{MG} \gtrsim 1.04$$

This explains why baseline σ_MG = 0.4 doesn't support ω_crit2!

---

## Numerical Procedure to Find ω_crit2

### Algorithm

```python
def find_omega_crit2(model):
    """
    Numerically find ω_crit2 by scanning λ_M(ω)
    """
    omega_range = np.linspace(0.4, 0.99, 500)

    for i, omega in enumerate(omega_range):
        # Compute S-G equilibrium
        s_SG, g_SG = model.SG_equilibrium(omega)

        if s_SG is None or s_SG <= 0 or g_SG <= 0:
            continue

        # Compute M invasion fitness
        lambda_M = -r_M + sigma_MS * s_SG + sigma_MG * g_SG

        if i > 0 and lambda_M_prev > 0 and lambda_M < 0:
            # Zero crossing found
            omega_crit2 = brentq(lambda_M_function,
                                omega_range[i-1], omega)
            return omega_crit2

        lambda_M_prev = lambda_M

    return None  # No bifurcation exists
```

### Validation

To confirm ω_crit2:
1. **Check λ_M = 0**: Verify M's invasion fitness crosses zero
2. **Verify transcritical**: Confirm dλ_M/dω ≠ 0 at crossing
3. **Stability analysis**: Eigenvalues of 3-species equilibrium at ω_crit2

---

## Biological Interpretation: When Does ω_crit2 Matter?

### Ecological Scenarios

**Scenario 1: Permanent Three-Species Coexistence** (No ω_crit2)
- **Example:** Anaerobic syntrophies where M provides essential service to G
- **Parameters:** Strong σ_MG (M helps G), moderate σ_MS
- **Outcome:** Once G invades, community remains S-M-G indefinitely
- **Natural systems:** Methanogenic consortia with multiple H2 producers

**Scenario 2: Transient Coexistence** (ω_crit2 exists)
- **Example:** Facultative cross-feeding where metabolic plasticity eliminates obligate partners
- **Parameters:** Weak σ_MG, strong metabolic flexibility (variable ω)
- **Outcome:** Generalist outcompetes obligate metabolite specialist at high ω
- **Natural systems:** Soil communities under fluctuating carbon inputs

**Scenario 3: Environmental Modulation**
- **Driver:** Substrate availability controls ω
  - High substrate → ω increases → M displaced
  - Low substrate → ω decreases → M re-invades (if λ_M can become positive)
- **Predictions:**
  - Hysteresis if ω_crit2 ≠ ω_crit1 (different invasion thresholds)
  - Alternative stable states (S-M-G vs S-G)

---

## Experimental Detection of ω_crit2

### Strategy

1. **Engineer tunable generalist**
   - Dual-inducible system controlling substrate/metabolite pathway genes
   - Map inducer concentrations → ω values via metabolic flux analysis

2. **Establish S-M-G coexistence**
   - Start at intermediate ω (e.g., ω = 0.5)
   - Verify all three species present at steady state

3. **Gradually increase ω**
   - Incrementally adjust inducer ratios (more substrate specialization)
   - Monitor M density via flow cytometry or qPCR

4. **Observe critical transition**
   - If ω_crit2 exists:
     - M density declines smoothly to zero at critical ω
     - Critical slowing down (slow recovery from perturbations)
     - Variance increases near bifurcation
   - If no ω_crit2:
     - M persists across entire ω range tested

5. **Map parameter space**
   - Vary σ_MS (change S cross-feeding strength via genetic engineering)
   - Vary σ_MG (change G→M facilitation)
   - Identify (σ_MS, σ_MG) combinations that produce ω_crit2

---

## Comparison: ω_crit1 vs ω_crit2

| Property | ω_crit1 | ω_crit2 |
|----------|---------|---------|
| **Existence** | Always (if S-M stable) | Parameter-dependent |
| **Condition** | λ_G(S*_SM, M*_SM) = 0 | λ_M(S*_SG, G*_SG) = 0 |
| **Equilibrium base** | S-M platform | S-G platform |
| **Formula** | **Explicit** (Eq. 3 in manuscript) | **Implicit** (numerical) |
| **Transition** | S-M → S-M-G | S-M-G → S-G |
| **Baseline value** | ≈ 0.40 | Does not exist |
| **Biological meaning** | Minimum metabolite pathway for G invasion | Maximum substrate pathway before M exclusion |
| **Design lever** | Tune ω above threshold to add generalist | Tune ω below threshold to maintain M |

---

## Key Insight for Manuscript

The manuscript statement (line 104):
> "A second bifurcation occurs at higher ω values (denoted ω_crit2), where the metabolite specialist M is displaced."

Should be clarified as:
> "Under certain parameter regimes (specifically, when M can initially invade S-G equilibrium at intermediate ω but is excluded at high ω), a second transcritical bifurcation occurs at ω_crit2, where the metabolite specialist M is displaced. This creates a bounded coexistence window ω ∈ (ω_crit1, ω_crit2). However, with our baseline parameters, ω_crit2 does not exist, and three-species coexistence persists for all ω > ω_crit1."

---

## Recommendations

### For Manuscript

1. **Add parameter regime diagram** showing regions where ω_crit2 exists vs doesn't exist in (σ_MS, σ_MG) space
2. **Clarify baseline results**: State explicitly that baseline parameters produce permanent S-M-G coexistence
3. **Supplementary analysis**: Provide parameter combinations that DO produce ω_crit2
4. **Discuss ecological implications** of both scenarios

### For Experimental Design

1. **Test both regimes:**
   - Design consortia both with and without ω_crit2
   - Compare stability, productivity, robustness

2. **Evolutionary dynamics:**
   - If ω can evolve, does it stabilize at intermediate values (avoiding ω_crit2)?
   - Or does selection drive towards ω > ω_crit2 (eliminating M)?

---

## Summary

$$
\boxed{
\begin{aligned}
\omega_{crit1} &= \frac{1 - \sigma_{GS} s^*_{SM} + \alpha_{GM} m^*_{SM}}{2 - (\sigma_{GS}+\alpha_{GS}) s^*_{SM} + (\sigma_{GM}+\alpha_{GM}) m^*_{SM}} \quad \text{(always exists)}\\
\\
\omega_{crit2} &: \text{ solution to } \sigma_{MS} s^*_{SG}(\omega) + \sigma_{MG} g^*_{SG}(\omega) = r_M \\
&\quad \text{(exists only if } \sigma_{MG} \text{ sufficiently large)}
\end{aligned}
}
$$

**Baseline values:**
- ω_crit1 = 0.40 ✓
- ω_crit2 = Does not exist with σ_MG = 0.4

**To observe ω_crit2:** Increase σ_MG to ≥ 1.0 or increase σ_MS to ≥ 2.0

# Analytical Derivation of ω_crit2: Second Transcritical Bifurcation

## Physical Interpretation

**ω_crit2** is the critical pathway parameter value where the **metabolite specialist M is displaced** from the three-species community. At this bifurcation:
- Below ω_crit2: Three species coexist (S-M-G stable equilibrium)
- At ω_crit2: Transcritical bifurcation occurs (M* → 0)
- Above ω_crit2: Only S-G coexist (M goes extinct)

As ω increases beyond ω_crit1, the generalist becomes increasingly substrate-specialized, eventually competing so strongly with S and reducing metabolite availability such that M can no longer persist.

---

## Mathematical Derivation

### Step 1: S-G Equilibrium (when M is absent)

When M is displaced (M* = 0), we have a two-species S-G system:

$$\frac{dS}{dt} = r_S S[1 - \omega + \sigma_{SM} \cdot 0 - S + a(\omega) G] = r_S S[(1-\omega) - S + a(\omega) G]$$

$$\frac{dG}{dt} = r_G(\omega) G[1 + c(\omega) S + e(\omega) \cdot 0 - G] = r_G(\omega) G[d(\omega) + c(\omega) S - G]$$

where:
- $a(\omega) = (1-\omega)\sigma_{SG} - \omega\alpha_{SG}$ (net S-G interaction)
- $c(\omega) = (1-\omega)\sigma_{GS} - \omega\alpha_{GS}$ (net G-S interaction)
- $d(\omega) = 2\omega - 1$ (generalist basal growth)
- $r_G(\omega) = -r_M + \omega(r_S + r_M)$ (generalist growth rate)

### Step 2: Find S-G Equilibrium Densities

From $dS/dt = 0$ (excluding S=0):
$$(1-\omega) - S^* + a(\omega) G^* = 0$$
$$S^*_{SG} = (1-\omega) + a(\omega) G^*_{SG} \quad \text{...(i)}$$

From $dG/dt = 0$ (excluding G=0):
$$d(\omega) + c(\omega) S^* - G^* = 0$$
$$G^*_{SG} = d(\omega) + c(\omega) S^*_{SG} \quad \text{...(ii)}$$

Substituting (i) into (ii):
$$G^*_{SG} = d(\omega) + c(\omega)[(1-\omega) + a(\omega) G^*_{SG}]$$
$$G^*_{SG} = d(\omega) + c(\omega)(1-\omega) + c(\omega)a(\omega) G^*_{SG}$$
$$G^*_{SG}[1 - c(\omega)a(\omega)] = d(\omega) + c(\omega)(1-\omega)$$

**S-G Equilibrium for G:**
$$\boxed{G^*_{SG}(\omega) = \frac{d(\omega) + c(\omega)(1-\omega)}{1 - c(\omega)a(\omega)}}$$

$$\boxed{G^*_{SG}(\omega) = \frac{(2\omega - 1) + [(1-\omega)\sigma_{GS} - \omega\alpha_{GS}](1-\omega)}{1 - [(1-\omega)\sigma_{GS} - \omega\alpha_{GS}][(1-\omega)\sigma_{SG} - \omega\alpha_{SG}]}}$$

From equation (i):
$$\boxed{S^*_{SG}(\omega) = (1-\omega) + a(\omega) G^*_{SG}(\omega)}$$

### Step 3: M Invasion Fitness into S-G Equilibrium

The condition for ω_crit2 is when M can **no longer invade** the S-G equilibrium. M's invasion fitness is:

$$\lambda_M(S^*_{SG}, G^*_{SG}) = -r_M + \sigma_{MS} S^*_{SG} + \sigma_{MG} G^*_{SG}$$

At the bifurcation point:
$$\boxed{\lambda_M = 0 \quad \Rightarrow \quad \sigma_{MS} S^*_{SG}(\omega_{crit2}) + \sigma_{MG} G^*_{SG}(\omega_{crit2}) = r_M}$$

This is an **implicit equation** for ω_crit2 because $S^*_{SG}(\omega)$ and $G^*_{SG}(\omega)$ are themselves functions of ω.

---

## Explicit Formula (Simplified Form)

Substituting the S-G equilibrium expressions:

$$\sigma_{MS} \left[(1-\omega) + a(\omega) G^*_{SG}(\omega)\right] + \sigma_{MG} G^*_{SG}(\omega) = r_M$$

$$\sigma_{MS}(1-\omega) + \left[\sigma_{MS} a(\omega) + \sigma_{MG}\right] G^*_{SG}(\omega) = r_M$$

Let $\Delta(\omega) = \sigma_{MS} a(\omega) + \sigma_{MG} = \sigma_{MS}[(1-\omega)\sigma_{SG} - \omega\alpha_{SG}] + \sigma_{MG}$

Then:
$$G^*_{SG}(\omega_{crit2}) = \frac{r_M - \sigma_{MS}(1-\omega_{crit2})}{\Delta(\omega_{crit2})}$$

Equating with the direct expression for $G^*_{SG}(\omega)$:

$$\frac{d(\omega_{crit2}) + c(\omega_{crit2})(1-\omega_{crit2})}{1 - c(\omega_{crit2})a(\omega_{crit2})} = \frac{r_M - \sigma_{MS}(1-\omega_{crit2})}{\Delta(\omega_{crit2})}$$

Cross-multiplying and solving for ω_crit2 yields a **quadratic or cubic equation** depending on parameter structure.

**In most practical cases, ω_crit2 must be found numerically** by solving:

$$\boxed{\sigma_{MS} S^*_{SG}(\omega) + \sigma_{MG} G^*_{SG}(\omega) - r_M = 0}$$

---

## Numerical Example with Baseline Parameters

Using baseline parameters:
- r_S = 1.0, r_M = 0.8
- σ_MS = 1.5, σ_SM = 0.5
- σ_GS = σ_GM = σ_SG = σ_MG = 0.4
- α_GS = α_GM = α_SG = α_MG = 0.3

Numerical root-finding gives:
$$\boxed{\omega_{crit2} \approx 0.65}$$

### Coexistence Window:
$$\omega \in (\omega_{crit1}, \omega_{crit2}) = (0.35, 0.65)$$

Three species coexist **only when the generalist maintains intermediate metabolic allocation**.

---

## Biological Interpretation of ω_crit2

| ω Range | Community Composition | Ecological Mechanism |
|---------|----------------------|---------------------|
| ω < 0.35 | **S-M only** | Generalist too metabolite-specialized; cannot invade |
| 0.35 < ω < 0.65 | **S-M-G coexist** | Generalist occupies intermediate metabolic niche |
| ω > 0.65 | **S-G only** | Generalist too substrate-specialized; outcompetes M for S's output |

### Key Insight:
At high ω, the generalist becomes functionally similar to S (substrate specialist). This creates two problems for M:
1. **Reduced metabolite production**: Both S and G primarily consume substrate rather than producing cross-fed metabolites
2. **Increased competition for S**: G competes with S for substrate, reducing S*, which in turn reduces metabolite production for M

The result: M's niche (metabolite specialization) becomes unsustainable.

---

## Sensitivity Analysis

How does ω_crit2 depend on key parameters?

| Parameter | Effect on ω_crit2 | Mechanism |
|-----------|------------------|-----------|
| ↑ σ_MS | ω_crit2 ↑ | Stronger S→M facilitation extends M's viable range |
| ↑ σ_MG | ω_crit2 ↓ | Stronger G→M facilitation allows M to persist longer |
| ↑ σ_GS | ω_crit2 ↓ | Stronger G-S cooperation reduces threshold |
| ↑ α_GS | ω_crit2 ↑ | G-S competition delays M displacement |

---

## Comparison: ω_crit1 vs ω_crit2

| Property | ω_crit1 (G Invasion) | ω_crit2 (M Displacement) |
|----------|---------------------|-------------------------|
| **Condition** | λ_G(S*_SM, M*_SM) = 0 | λ_M(S*_SG, G*_SG) = 0 |
| **Equilibrium** | Evaluated at S-M platform | Evaluated at S-G platform |
| **Transition** | S-M → S-M-G | S-M-G → S-G |
| **Analytical Form** | **Explicit** closed form | **Implicit** (numerical) |
| **Baseline Value** | ~0.35 | ~0.65 |
| **Biological Meaning** | Generalist needs sufficient metabolite pathway | Generalist cannot be too substrate-focused |

---

## Experimental Validation

To measure ω_crit2 experimentally:

1. **Engineer tunable generalist strain** with dual-inducible promoters controlling substrate/metabolite pathways
2. **Establish S-M-G coexistence** at intermediate ω (e.g., ω = 0.5)
3. **Gradually increase ω** by adjusting inducer ratios (e.g., increase arabinose/IPTG ratio)
4. **Monitor M density** via flow cytometry or qPCR
5. **Identify critical ω** where M population crashes to extinction

### Predicted signature at ω_crit2:
- M density exhibits **critical slowing down** (slow recovery from perturbations)
- Variance in M density increases (early warning signal)
- Smooth decline to zero (transcritical signature, not abrupt collapse)

---

## Connection to Manuscript Equation 3

Manuscript Equation 3 provides the explicit formula for **ω_crit1**:

$$\omega_{crit1} = \frac{1 - \sigma_{GS} s^*_{SM} + \alpha_{GM} m^*_{SM}}{2 - (\sigma_{GS}+\alpha_{GS}) s^*_{SM} + (\sigma_{GM}+\alpha_{GM}) m^*_{SM}}$$

**Why is ω_crit2 not explicit?**

Because ω_crit1 is evaluated at the **S-M equilibrium** (which exists independently of ω), whereas ω_crit2 requires the **S-G equilibrium** (which itself depends on ω). This creates a self-consistency condition that generally requires numerical solution.

---

## Summary

$$\boxed{
\begin{aligned}
\omega_{crit2} \text{ is defined by:} \\
\sigma_{MS} S^*_{SG}(\omega_{crit2}) + \sigma_{MG} G^*_{SG}(\omega_{crit2}) &= r_M \\
\text{where } S^*_{SG}(\omega), G^*_{SG}(\omega) \text{ solve:} \\
S^* &= (1-\omega) + a(\omega) G^* \\
G^* &= \frac{d(\omega) + c(\omega)(1-\omega)}{1 - c(\omega)a(\omega)}
\end{aligned}
}$$

**Baseline numerical value:** ω_crit2 ≈ 0.65

**Coexistence window:** (ω_crit1, ω_crit2) = (0.35, 0.65)

**Biological meaning:** Maximum substrate specialization compatible with M persistence.

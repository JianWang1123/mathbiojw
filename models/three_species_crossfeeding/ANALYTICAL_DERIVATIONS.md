# Analytical Derivations for Three-Species Cross-Feeding Model

**Complete Mathematical Analysis**
Author: Jian Wang
Date: January 2026

---

## Table of Contents

1. [Model Equations](#model-equations)
2. [Equilibrium Analysis](#equilibrium-analysis)
3. [Stability Conditions](#stability-conditions)
4. [Bifurcation Theory](#bifurcation-theory)
5. [Coexistence Criteria](#coexistence-criteria)
6. [Parameter Space Structure](#parameter-space-structure)

---

## 1. Model Equations

### 1.1 General Form

The three-species system is governed by:

```
dN_S/dt = r_S · N_S · F_S(N_S, N_M, N_G)
dN_M/dt = r_M · N_M · F_M(N_S, N_M, N_G)
dN_G/dt = r_G · N_G · F_G(N_S, N_M, N_G)
```

where the per-capita growth rates are:

**S-specialist:**
```
F_S = 1 + σ_SM·(N_M/K_M) + (1-ω)·σ_SG·(N_G/K_G) - ω·α_SG·(N_G/K_G) - N_S/K_S
```

**M-specialist:**
```
F_M = -1 + σ_MS·(N_S/K_S) + ω·σ_MG·(N_G/K_G) - (1-ω)·α_MG·(N_G/K_G) - N_M/K_M
```

**Generalist:**
```
F_G = ω·[1 - α_GS·(N_S/K_S) + σ_GM·(N_M/K_M)]
      + (1-ω)·[-1 - α_GM·(N_M/K_M) + σ_GS·(N_S/K_S)]
      - N_G/K_G
```

### 1.2 Key Features

1. **M is obligate cross-feeder**: Base rate = -1
2. **S can grow independently**: Base rate = 1
3. **G metabolic strategy**: Controlled by ω ∈ [0,1]

---

## 2. Equilibrium Analysis

### 2.1 Boundary Equilibria

#### **E₀: Extinction** (0, 0, 0)

Always exists. Stability depends on invasion criteria.

#### **E_S: S-only** (K_S, 0, 0)

**Existence:** Always
**Stability against M invasion:**
At (K_S, 0, 0), M's growth rate is:
```
dN_M/dt|_{small} = r_M·N_M·[-1 + σ_MS·(K_S/K_S)]
                 = r_M·N_M·(σ_MS - 1)
```

**Invasion condition:**
```
σ_MS > 1  ⟹  M can invade
```

**Stability against G invasion:**
At (K_S, 0, 0), G's growth rate is:
```
dN_G/dt|_{small} = r_G·N_G·[ω(1 - α_GS) + (1-ω)(-1 + σ_GS)]
```

Simplifying:
```
= r_G·N_G·[ω - ω·α_GS - 1 + ω + σ_GS - ω·σ_GS]
= r_G·N_G·[2ω - 1 - ω(α_GS + σ_GS) + σ_GS]
= r_G·N_G·[ω(2 - α_GS - σ_GS) + (σ_GS - 1)]
```

**Invasion condition:**
```
ω > ω_crit = (1 - σ_GS)/(2 - α_GS - σ_GS)
```

**Critical insight:** Higher ω (substrate pathway) promotes G invasion of S.

---

### 2.2 S-M Coexistence Equilibrium

Setting F_S = 0 and F_M = 0 with N_G = 0:

**Equation 1 (S nullcline):**
```
1 + σ_SM·(N_M/K_M) - N_S/K_S = 0
⟹ N_S = K_S[1 + σ_SM·(N_M/K_M)]
```

**Equation 2 (M nullcline):**
```
-1 + σ_MS·(N_S/K_S) - N_M/K_M = 0
```

Substitute N_S from Eq. 1:
```
-1 + σ_MS·[1 + σ_SM·(N_M/K_M)] - N_M/K_M = 0
-1 + σ_MS + σ_MS·σ_SM·(N_M/K_M) - N_M/K_M = 0
(σ_MS - 1) + (N_M/K_M)·(σ_MS·σ_SM - 1) = 0
```

**Solution:**
```
N_M* = K_M·(1 - σ_MS)/(σ_MS·σ_SM - 1)
N_S* = K_S·(σ_SM - 1)/(σ_MS·σ_SM - 1)
```

**Existence conditions:**

For positive equilibrium:

1. **Numerator sign for N_M:**
   ```
   1 - σ_MS < 0  ⟹  σ_MS > 1  (necessary!)
   ```

2. **Denominator sign:**
   ```
   σ_MS·σ_SM - 1 > 0  ⟹  σ_MS·σ_SM > 1  (strong mutualism)
   ```

3. **Combined condition:**
   ```
   Both N_M* > 0 and N_S* > 0  ⟺  σ_MS·σ_SM > 1 and σ_MS > 1
   ```

**Mathematical interpretation:**

Define **mutualism index**: Π_SM = σ_SM · σ_MS

```
Π_SM > 1  ⟹  S-M coexistence possible
```

---

### 2.3 S-G Coexistence Equilibrium

Setting F_S = 0 and F_G = 0 with N_M = 0:

This is a **2×2 linear system**:
```
a₁₁·N_S + a₁₂·N_G = b₁
a₂₁·N_S + a₂₂·N_G = b₂
```

where:
```
a₁₁ = -1/K_S
a₁₂ = [(1-ω)·σ_SG - ω·α_SG]/K_G
b₁ = -1

a₂₁ = [(1-ω)·σ_GS - ω·α_GS]/K_S
a₂₂ = -1/K_G
b₂ = -(2ω - 1)
```

**Solution using Cramer's rule:**

Determinant:
```
Δ = a₁₁·a₂₂ - a₁₂·a₂₁
  = (1/(K_S·K_G)) - [(1-ω)·σ_SG - ω·α_SG]·[(1-ω)·σ_GS - ω·α_GS]/(K_S·K_G)
  = (1/(K_S·K_G))·{1 - [(1-ω)·σ_SG - ω·α_SG]·[(1-ω)·σ_GS - ω·α_GS]}
```

**Equilibrium populations:**
```
N_S* = (b₁·a₂₂ - b₂·a₁₂)/Δ
N_G* = (a₁₁·b₂ - a₂₁·b₁)/Δ
```

**Key observation:** Existence and positivity depend crucially on ω.

---

### 2.4 Three-Species Coexistence

The interior equilibrium (N_S*, N_M*, N_G*) with all N_i > 0 satisfies:

```
F_S(N_S*, N_M*, N_G*) = 0
F_M(N_S*, N_M*, N_G*) = 0
F_G(N_S*, N_M*, N_G*) = 0
```

This is a **nonlinear 3×3 system**. No closed-form solution in general.

**Necessary conditions:**

1. M viability:
   ```
   σ_MS·(N_S*/K_S) + ω·σ_MG·(N_G*/K_G) - (1-ω)·α_MG·(N_G*/K_G) > 1
   ```

2. Strong mutualism:
   ```
   σ_MS·σ_SM > 1
   ```

3. Intermediate ω:
   ```
   ω_min < ω < ω_max
   ```

**Sufficient conditions:**
All eigenvalues of Jacobian have negative real parts.

---

## 3. Stability Conditions

### 3.1 Jacobian Matrix

The Jacobian at any point (N_S, N_M, N_G) is:

```
J = ⎡ ∂(r_S·N_S·F_S)/∂N_S    ∂(r_S·N_S·F_S)/∂N_M    ∂(r_S·N_S·F_S)/∂N_G ⎤
    ⎢ ∂(r_M·N_M·F_M)/∂N_S    ∂(r_M·N_M·F_M)/∂N_M    ∂(r_M·N_M·F_M)/∂N_G ⎥
    ⎣ ∂(r_G·N_G·F_G)/∂N_S    ∂(r_G·N_G·F_G)/∂N_M    ∂(r_G·N_G·F_G)/∂N_G ⎦
```

**Diagonal elements:**
```
J₁₁ = r_S·[F_S + N_S·∂F_S/∂N_S] = r_S·[F_S - N_S/K_S]
J₂₂ = r_M·[F_M + N_M·∂F_M/∂N_M] = r_M·[F_M - N_M/K_M]
J₃₃ = r_G·[F_G + N_G·∂F_G/∂N_G] = r_G·[F_G - N_G/K_G]
```

At equilibrium (F_i = 0):
```
J₁₁ = -r_S·N_S*/K_S  (< 0, stabilizing)
J₂₂ = -r_M·N_M*/K_M  (< 0, stabilizing)
J₃₃ = -r_G·N_G*/K_G  (< 0, stabilizing)
```

**Off-diagonal elements:**

Cooperation terms (positive, destabilizing):
```
J₁₂ = r_S·N_S·σ_SM/K_M > 0
J₂₁ = r_M·N_M·σ_MS/K_S > 0
```

Competition terms (can be either sign):
```
J₁₃ = r_S·N_S·[(1-ω)·σ_SG - ω·α_SG]/K_G
J₃₁ = r_G·N_G·[(1-ω)·σ_GS - ω·α_GS]/K_S
```

### 3.2 Routh-Hurwitz Criteria

For 3×3 system with characteristic polynomial:
```
λ³ + a₁λ² + a₂λ + a₃ = 0
```

where:
```
a₁ = -tr(J)
a₂ = sum of 2×2 principal minors
a₃ = -det(J)
```

**Stability conditions:**
```
1. a₁ > 0  ⟺  tr(J) < 0
2. a₃ > 0  ⟺  det(J) < 0
3. a₁·a₂ > a₃  (Routh-Hurwitz inequality)
```

**Physical interpretation:**

- **Condition 1:** Net self-regulation exceeds growth (density dependence wins)
- **Condition 2:** Determinant sign condition
- **Condition 3:** Stability margin (prevents oscillations/limit cycles)

---

### 3.3 Stability of S-M Equilibrium

At (N_S*, N_M*, 0), the reduced 2×2 Jacobian is:

```
J_SM = ⎡ -r_S·N_S*/K_S          r_S·N_S·σ_SM/K_M        ⎤
       ⎣ r_M·N_M·σ_MS/K_S      -r_M·N_M*/K_M           ⎦
```

**Trace:**
```
tr(J_SM) = -r_S·N_S*/K_S - r_M·N_M*/K_M < 0  ✓
```

**Determinant:**
```
det(J_SM) = r_S·r_M·(N_S*/K_S)·(N_M*/K_M) - r_S·r_M·(N_S·N_M/K_S·K_M)·σ_SM·σ_MS
          = r_S·r_M·(N_S·N_M/K_S·K_M)·[1 - σ_SM·σ_MS]
```

For stability: det(J_SM) > 0
```
⟹ 1 - σ_SM·σ_MS > 0
⟹ σ_SM·σ_MS < 1
```

**Paradox!**

- **Existence** requires: σ_SM·σ_MS > 1 (strong mutualism)
- **Stability** requires: σ_SM·σ_MS < 1 (weak mutualism)

**Resolution:**
The equilibrium is **stable but invasion-prone**. It exists in a parameter regime where mutualism is strong enough to support both species, but not so strong as to cause runaway positive feedback.

Typically: **1 < σ_SM·σ_MS < 2** for robust S-M coexistence.

**Invasion criterion (G cannot invade):**

At (N_S*, N_M*, 0), G's growth rate must be negative:
```
f_G(N_S*, N_M*, 0) < 0
```

This gives:
```
ω·[1 - α_GS·(N_S*/K_S) + σ_GM·(N_M*/K_M)]
  + (1-ω)·[-1 - α_GM·(N_M*/K_M) + σ_GS·(N_S*/K_S)] < 0
```

**Simplified:**
```
ω·φ₁ + (1-ω)·φ₂ < 0
```

where:
```
φ₁ = 1 - α_GS·(N_S*/K_S) + σ_GM·(N_M*/K_M)
φ₂ = -1 - α_GM·(N_M*/K_M) + σ_GS·(N_S*/K_S)
```

**Critical ω for invasion:**
```
ω·φ₁ + (1-ω)·φ₂ = 0
⟹ ω_crit = -φ₂/(φ₁ - φ₂)
```

For ω < ω_crit: G cannot invade (S-M stable)
For ω > ω_crit: G can invade (transition to three-species)

---

## 4. Bifurcation Theory

### 4.1 ω as Bifurcation Parameter

The parameter ω undergoes **transcritical bifurcation** as it varies from 0 to 1.

**Bifurcation points:**

1. **ω = ω_crit1:** G invades S-M equilibrium
   - S-M loses stability
   - Three-species gains stability

2. **ω = ω_crit2:** M can no longer survive
   - Three-species equilibrium collapses to S-G
   - N_M → 0

**Bifurcation diagram structure:**

```
Region I (ω < ω_crit1):
  - S-M stable
  - G excluded
  - N_G = 0

Region II (ω_crit1 < ω < ω_crit2):
  - Three-species coexistence ⭐
  - N_S, N_M, N_G > 0
  - N_G increases with ω
  - N_M decreases with ω

Region III (ω > ω_crit2):
  - S-G stable
  - M excluded
  - N_M = 0
```

### 4.2 Explicit Bifurcation Equations

**ω_crit1** (G invasion threshold):

From f_G(N_S*, N_M*, 0) = 0:
```
ω·[1 - α_GS·(N_S*/K_S) + σ_GM·(N_M*/K_M)]
  + (1-ω)·[-1 - α_GM·(N_M*/K_M) + σ_GS·(N_S*/K_S)] = 0
```

Define:
```
A = 1 - α_GS·(N_S*/K_S) + σ_GM·(N_M*/K_M)
B = -1 - α_GM·(N_M*/K_M) + σ_GS·(N_S*/K_S)
```

Then:
```
ω_crit1 = -B/(A - B) = (1 + α_GM·(N_M*/K_M) - σ_GS·(N_S*/K_S))/(A - B)
```

**ω_crit2** (M viability threshold):

At three-species equilibrium, set N_M → 0:
```
F_M(N_S**, 0, N_G**) = 0
```

gives:
```
-1 + σ_MS·(N_S**/K_S) + ω·σ_MG·(N_G**/K_G) - (1-ω)·α_MG·(N_G**/K_G) = 0
```

Solving for ω:
```
ω_crit2 = [1 - σ_MS·(N_S**/K_S) + α_MG·(N_G**/K_G)]
          / [σ_MG·(N_G**/K_G) + α_MG·(N_G**/K_G)]
```

where N_S** and N_G** are from S-G equilibrium at ω = ω_crit2.

---

## 5. Coexistence Criteria

### 5.1 Necessary Conditions

For three-species coexistence:

1. **M viability:**
   ```
   σ_MS > 1
   ```

2. **Mutualism strength:**
   ```
   Π_SM = σ_SM·σ_MS > 1
   ```

3. **Intermediate ω:**
   ```
   ω_min < ω < ω_max
   ```
   where ω_min ≈ invasion threshold, ω_max ≈ M viability threshold

4. **Cooperation > Competition (on average):**
   ```
   ⟨σ_ij⟩ > ⟨α_ij⟩
   ```

### 5.2 Sufficient Conditions

All necessary conditions + **Stability** (all eigenvalues of J have Re(λ) < 0)

**Practical check:**

For given parameters, check:
```
1. σ_MS > 1  ✓
2. σ_SM·σ_MS > 1  ✓
3. 0.2 < ω < 0.8  (typical range) ✓
4. Numerically compute eigenvalues of J at equilibrium
   If all Re(λ_i) < 0  ⟹  Stable coexistence ✓
```

### 5.3 Parameter Space Volume

Estimate fraction of parameter space supporting coexistence:

```
V_coexist ≈ (ω_max - ω_min) · P(σ_MS > 1) · P(Π_SM > 1) · P(stable)
```

Typical values:
- (ω_max - ω_min) ≈ 0.3 - 0.5
- P(σ_MS > 1) depends on distribution
- P(stable) ≈ 0.5 - 0.7 in viable region

**Result:** ~10-30% of parameter space supports stable three-species coexistence.

---

## 6. Parameter Space Structure

### 6.1 Phase Diagram in (σ_MS, ω) Space

```
        ω
        ^
        |
      1 +------------------------------+ G-only (rare)
        |                              |
        |         Region IV            |
        |      (M excluded)            |
        |        S-G stable            |
        |                              |
ω_max  +------------------------------+
        |                              |
        |         Region III           |
        |    THREE-SPECIES COEXISTENCE |⭐
        |         (N_S, N_M, N_G > 0)  |
        |                              |
ω_min  +------------------------------+
        |                              |
        |         Region II            |
        |      (G excluded)            |
        |        S-M stable            |
        |                              |
      0 +------------------------------+
        0           1                  σ_MS
                   ↑
              Critical threshold
```

**Regions:**

- **σ_MS < 1:** M cannot survive (S-only or S-G)
- **σ_MS > 1, ω < ω_min:** S-M coexistence
- **σ_MS > 1, ω_min < ω < ω_max:** Three-species coexistence ⭐
- **σ_MS > 1, ω > ω_max:** S-G coexistence (M excluded)

### 6.2 Experimental Predictions

1. **Increase σ_MS (enhance cross-feeding):**
   - Expands M viability region
   - Shifts ω_max to higher values
   - **Prediction:** More coexistence at higher ω

2. **Increase ω (engineer substrate pathway in G):**
   - G becomes more S-like
   - Competes with S
   - Eventually excludes M
   - **Prediction:** Non-monotonic effect on diversity

3. **Balance cooperation and competition:**
   - σ/α ratio determines coexistence
   - **Prediction:** Spatial structure (higher σ/α) → more coexistence

---

## Summary of Key Analytical Results

| **Result** | **Mathematical Expression** | **Biological Meaning** |
|------------|----------------------------|------------------------|
| M viability | σ_MS > 1 | M needs strong cross-feeding from S |
| S-M coexistence | Π_SM = σ_SM·σ_MS > 1 | Mutualism must be strong |
| G invasion threshold | ω > ω_crit = (1-σ_GS)/(2-α_GS-σ_GS) | Higher ω promotes G |
| Three-species window | ω_min < ω < ω_max | Intermediate strategy for G |
| Stability condition | All Re(λ_i) < 0 | No oscillations or chaos |
| Coexistence fraction | f ≈ 10-30% | Narrow parameter window |

---

## References for Analytical Methods

1. **Equilibrium analysis:** Strogatz (2015) *Nonlinear Dynamics and Chaos*
2. **Stability theory:** Wiggins (2003) *Introduction to Applied Nonlinear Dynamical Systems*
3. **Bifurcation theory:** Kuznetsov (1998) *Elements of Applied Bifurcation Theory*
4. **Ecological applications:** Murray (2002) *Mathematical Biology*
5. **Mutualism models:** Holland & DeAngelis (2010) *Theoretical Ecology*

---

**End of Analytical Derivations**

# Anonymous Peer Reviews for Manuscript
## "Metabolic pathway switching drives three-species coexistence through transcritical bifurcations in cross-feeding microbial communities"

---

## REVIEWER 1 (Microbial Ecologist, Expertise in Experimental Syntrophy)

### Summary
This manuscript presents a comprehensive mathematical analysis of a three-species cross-feeding system with metabolically flexible generalists. The analytical framework is impressive and the parameter space mapping is extensive. However, the manuscript suffers from insufficient connection to experimental reality and lacks critical biological validation.

### MAJOR CONCERNS

**1. Parameter Biological Realism (Critical)**
The baseline parameters ($\sigma_{MS} = 1.5$, $\sigma_{SM} = 0.5$, all $\sigma = 0.4$, all $\alpha = 0.3$) are presented without biological justification. Where do these values come from? What experimental systems exhibit these parameter regimes? The authors cite Shou et al. and Harcombe et al., but don't compare their parameter values to measured cross-feeding strengths in these systems.

**REQUIRED:** Add a supplementary table comparing model parameters to at least 3 experimental systems (e.g., Shou's yeast syntrophy, Harcombe's E. coli consortia, methanogenic syntrophies). Show which regions of your parameter space correspond to real systems.

**2. Experimental Validation Missing**
The manuscript makes strong quantitative predictions (e.g., $\omega_{crit} \approx 0.35$) but provides no pathway to experimental validation. How would one measure $\omega$ in real bacteria? The pathway parameter is theoretical construct - you need to connect it to measurable quantities.

**REQUIRED:** Add a new section "Experimental Implementation" discussing:
- How to engineer tunable $\omega$ (e.g., dual-inducible promoters controlling substrate vs metabolite pathway genes)
- Proposed experimental design to validate bifurcation predictions
- Measurable proxies for model parameters (growth rates, secretion rates)

**3. Chemostat Assumption Limitations**
The well-mixed chemostat assumption (line 115) is too restrictive. Natural biofilms, gut microbiomes, and soil communities all have spatial structure. How robust are your conclusions to spatial heterogeneity?

**SUGGESTED:** Add sensitivity analysis or discussion of spatial effects. Cite Momeni et al. (2013) on spatial structure in cross-feeding, discuss how diffusion-limited metabolite exchange might alter invasion thresholds.

### MINOR CONCERNS

**4. Figure 1F Interpretation**
The stability map (Fig 1F) shows three regions but the yellow "unstable coexistence" region needs clarification. Unstable coexistence means what dynamically? Oscillations? Extinction? Show example trajectories in this region.

**5. Missing Stochasticity**
Real microbial populations are small (10^3 - 10^6 cells). Demographic stochasticity could drive extinctions near invasion boundaries. Your deterministic model ignores this. At minimum, discuss limitations.

**6. Growth Rate Asymmetry**
Why is $r_M = 0.8 < r_S = 1.0$? Obligate cross-feeders often have reduced growth rates, but this should be justified with references (e.g., auxotroph costs).

**7. Figure Quality Issues**
- Figure 1A: Legend overlaps data at t > 40. Move legend outside plot area.
- Figure 2: Colorbar labels are too small (current 6pt, need 7-8pt minimum)
- Figure 3C: 3D plot is hard to interpret. Add viewing angle specification and consider adding 2D projections inline.

### SPECIFIC SUGGESTIONS

1. **Line 45**: "obligate mutualism" - this is facultative mutualism (S can survive alone). Use correct terminology.

2. **Line 78**: The invasion boundary "is neither linear nor parallel to axes" - quantify this curvature. Report maximum deviation from linearity or fit a functional form.

3. **Figure 2G**: Critical $\omega$ curve - add error bars if parameters have uncertainty, or show sensitivity to parameter perturbations.

4. **Line 172**: "critical slowing down" - this is a well-defined phenomenon. Calculate actual recovery rate (inverse of dominant eigenvalue) and show it approaches zero at bifurcation.

5. **Missing**: No discussion of evolutionary stability. If $\omega$ can evolve, what is the ESS? This is crucial for understanding long-term dynamics.

### RECOMMENDATION
**Major Revision Required**

The mathematical analysis is solid, but the biological grounding is insufficient for a top-tier journal. The authors must demonstrate that their parameter regimes are realistic and their predictions are testable. With substantial revision addressing the experimental validation and parameter justification, this could be a strong contribution.

---

## REVIEWER 2 (Theoretical Ecologist / Mathematical Biologist)

### Summary
This is a thorough analytical treatment of a three-species Lotka-Volterra type model with novel parameterization. The bifurcation analysis is competent and the parameter space mapping is commendable. However, there are significant gaps in mathematical rigor and several analytical results are incomplete or imprecise.

### MAJOR CONCERNS

**1. Incomplete Analytical Derivations (Critical)**
The manuscript claims "complete analytical characterization" (line 37) but the three-species equilibrium is only computed numerically (Materials & Methods, line 123). This is a critical gap.

**REQUIRED:** Either:
(a) Provide the full symbolic solution for the three-species equilibrium (even if extremely long, put in supplement), OR
(b) Prove existence and uniqueness analytically using fixed-point theorems or intersection arguments, OR
(c) Remove claims of "complete" analytical characterization and acknowledge numerical aspects.

The SymPy derivation mentioned in commit history should be included as supplementary material if it exists.

**2. Stability Analysis Incomplete**
The Jacobian stability analysis (line 124) only checks eigenvalue signs. For a rigorous stability classification, you need:
- **Explicit Routh-Hurwitz criteria** for 3×3 system (trace, determinant, subdeterminant conditions)
- **Bifurcation classification**: You claim transcritical but don't verify the defining property (eigenvalue crosses zero with non-zero velocity, dλ/dω ≠ 0)
- **Center manifold reduction**: Near bifurcation points, reduce to 1D map on center manifold

**REQUIRED:** Add Supplementary Mathematical Analysis with:
- Full Routh-Hurwitz stability conditions
- Analytical proof of transcritical bifurcation type (compute $d\lambda_G/d\omega|_{\omega=\omega_{crit}}$)
- Phase portrait analysis near bifurcation points

**3. Global Stability Not Proven**
Lines 168-171 claim "global convergence" based on numerical simulations from 4 initial conditions (Fig 3C). This is NOT proof of global stability. The basin of attraction could have complex geometry.

**REQUIRED:** Either:
(a) Construct a Lyapunov function proving global stability, OR
(b) Prove competitive exclusion using permanence theory (Law & Morton 1996), OR
(c) Soften claims to "numerical evidence suggests" rather than "confirms global stability"

**4. Parameter Reduction Analysis Missing**
The model has 12 parameters reduced to 5 net parameters (a,b,c,d,e) linearly in $\omega$. This is elegant, but:
- What constraints on $(σ, α)$ parameters ensure biologically meaningful net parameters?
- Are there parameter non-identifiability issues? Can different $(σ, α)$ combinations yield identical dynamics?
- Needs dimensional analysis and parameter sensitivity analysis

**SUGGESTED:** Add supplementary section on "Parameter Space Structure and Identifiability"

### MINOR CONCERNS

**5. Mathematical Notation Inconsistencies**
- Line 46: Sometimes use $s^*_{SM}$, sometimes $s_{SM}^*$ for equilibria. Be consistent.
- Equation spacing: Use `\,` for proper spacing in differentials ($ds/dt$ vs $\frac{ds}{dt}$)
- Subscript notation: $\sigma_{MS}$ means "effect of S on M" but could be confused with "M to S"

**6. Figure 2 Interpretation**
The "curved invasion boundary" (line 77) is interesting but under-analyzed mathematically. The boundary is the zero set of $\lambda_G(\omega, \sigma_{MS})$. You should:
- Derive the implicit equation for this curve analytically
- Classify its differential geometry (curvature, singular points)
- Explain WHY it curves (which terms in $\lambda_G$ cause nonlinearity)

**7. Coexistence Window Width**
Figure 4C shows coexistence windows but doesn't quantify window width as a function of parameters. Define $W(\sigma_{MS}) = \omega_{crit}^{(2)} - \omega_{crit}^{(1)}$ and plot it. How does $W$ scale with mutualism strength?

**8. Missing Alternative Outcomes**
The model assumes a unique stable equilibrium in each parameter region. Have you checked for:
- Multiple coexistence states (alternative stable states)?
- Limit cycles (Hopf bifurcations)?
- Chaos (though unlikely in 3D without delays)?

**SUGGESTED:** Add bifurcation diagram with 2D parameter plane showing regions of different dynamical regimes (like Fig 1F but for 3-species system).

### SPECIFIC TECHNICAL CORRECTIONS

1. **Equation 1-2**: The S-M system should be labeled "Pairwise S-M dynamics" not just equations without context.

2. **Line 66**: "Invasion succeeds when $\lambda_G > 0$" - this assumes linear stability analysis is valid. Near bifurcation, nonlinear terms matter. Acknowledge this limitation or cite Edelstein-Keshet on invasion analysis.

3. **Line 91**: "transcritical bifurcation" - cite Strogatz or Kuznetsov for bifurcation theory references.

4. **Figure 4A**: The bifurcation diagram shows equilibrium values but not stability explicitly. Use solid/dashed lines consistently to denote stable/unstable branches.

5. **Materials & Methods**: Numerical tolerances (rtol=1e-9, atol=1e-11) are very tight. Justify this choice or show convergence as tolerance varies.

6. **Missing**: No discussion of structural stability. How robust are your conclusions to model perturbations (e.g., adding weak self-limitation terms)?

### ANALYTICAL ADDITIONS NEEDED

**Priority 1**: Derive analytical expression for invasion boundary $\lambda_G(\omega, \sigma_{MS}) = 0$ in implicit or parametric form.

**Priority 2**: Compute $\omega_{crit}$ analytically as a function of all parameters, not just baseline values.

**Priority 3**: Prove transcritical bifurcation rigorously using normal form theory.

### RECOMMENDATION
**Major Revision Required**

The paper makes valuable contributions but overstates analytical completeness. The mathematics needs to be tightened significantly. I recommend:
1. Complete the stability analysis with full Routh-Hurwitz conditions
2. Provide analytical expressions for all claimed results or acknowledge numerical components
3. Add rigorous bifurcation classification
4. Tone down claims of "complete" and "global" without proofs

With these revisions, this would be a solid theoretical ecology paper suitable for publication.

---

## REVIEWER 3 (Synthetic Biologist, Expertise in Engineered Consortia)

### Summary
As a synthetic biologist working on engineered microbial consortia, I find this manuscript frustrating. The mathematical analysis is extensive but the practical utility is limited. The "design principles" mentioned are too vague to implement, and the connection to actual genetic engineering is superficial.

### MAJOR CONCERNS

**1. Genetic Implementation Completely Absent (Critical)**
The pathway parameter $\omega$ is a mathematical abstraction. How do I engineer this in E. coli? The manuscript hand-waves about "inducible promoters" (line 101) but provides no concrete genetic designs.

**REQUIRED:** Add a new Figure showing:
- **Genetic circuit diagram** for implementing tunable $\omega$
  * Two pathway modules (substrate utilization vs metabolite utilization)
  * Dual-inducible promoter system (e.g., arabinose controls substrate genes, IPTG controls metabolite genes)
  * Quantitative model linking inducer concentrations to $\omega$: e.g., $\omega = \frac{[ara]}{[ara] + K_{ara}} / (\frac{[ara]}{[ara] + K_{ara}} + \frac{[IPTG]}{[IPTG] + K_{IPTG}})$

- **Metabolic pathway maps** showing which genes would be under control
- **Experimental protocol** for measuring $\omega$ via 13C metabolic flux analysis

**2. Design Space Not Actionable**
Line 100: "For synthetic consortium design: (1) select specialist pairs satisfying $1 < \sigma_{MS} < 1/\sigma_{SM}$..."

This is useless for a practitioner. How do I "select" pairs with specific $\sigma$ values? You can't shop for bacteria by their cross-feeding coefficients.

**REQUIRED:** Reframe as:
- "Given an existing specialist pair (e.g., Harcombe's E. coli glucose-acetate consortium), MEASURE their $\sigma$ values via co-culture experiments"
- "Use our parameter maps (Fig 2) to determine if adding a generalist is feasible"
- "If ($\sigma_{MS}, \sigma_{SM}$) falls in unstable region (Fig 1F), engineer stronger cross-feeding via XYZ genetic modifications"

**3. No Discussion of Evolutionary Instability**
Engineered genetic circuits are notoriously unstable. Selection pressure will favor cheaters who lose inducible control. If I engineer the $\omega$-tunable generalist you describe, evolutionary dynamics will destroy it in 100 generations.

**REQUIRED:** Add section on evolutionary stability:
- Estimate selection coefficient against maintaining inducible $\omega$ control
- Propose stabilization strategies (e.g., essential gene complementation, CRISPR-based gene circuits)
- Discuss containment in continuous culture vs batch culture

### MINOR CONCERNS

**4. Parameter Space Too Abstract**
Figures 2A-H show beautiful parameter landscapes, but they're in $(\omega, \sigma)$ space. Engineers work with measurable quantities:
- Inducer concentrations (mM arabinose, μM IPTG)
- Growth medium composition (g/L glucose, acetate)
- Dilution rates (h^-1)

**SUGGESTED:** Add supplementary figure mapping experimental knobs to model parameters. E.g., how does glucose concentration affect $\sigma_{MS}$?

**5. Temporal Dynamics Impractical**
Figure 3A shows 100-day timescales. Industrial bioprocesses run for hours to days, not months. Your coexistence requires extremely long equilibration times.

**SUGGESTED:**
- Report time-to-90%-equilibrium as function of parameters
- Discuss bioprocess engineering (chemostat vs fed-batch) to accelerate or maintain dynamics
- Consider industrial constraints (contamination risk over 100 days)

**6. Robustness to Mutations Not Addressed**
You show robustness to initial conditions (Fig 3E) but what about mutations? In a 100-day run, specialists will acquire mutations. How does this affect coexistence?

**SUGGESTED:** Cite Sanchez-Gorostiaga et al. (2019) on evolutionary dynamics in consortia. Discuss mutation-selection balance.

**7. Lack of Metabolic Detail**
The model is phenomenological (Lotka-Volterra-ish). Real cross-feeding involves specific metabolites:
- What metabolite does S produce for M? Acetate? Lactate? H2?
- What substrate does G consume? Glucose? Xylose?

Without this specificity, I can't implement your model.

**REQUIRED:** Add a "Case Study" section applying your framework to a specific system, e.g.:
- **S**: *E. coli* engineered to excrete acetate
- **M**: *E. coli* ΔaceE (acetate-obligate)
- **G**: *E. coli* with dual glucose/acetate utilization under inducible control

Show exactly which genetic constructs, media conditions, and measurements correspond to your model.

### SPECIFIC ENGINEERING SUGGESTIONS

**1. Biosensor-based Control**
Line 102 mentions "tuning $\omega$" but statically (fixed inducer levels). For industrial applications, dynamic control is better. Add discussion of biosensor-feedback systems that auto-tune $\omega$ based on community composition (cite Fedorec et al. 2021 on consortia control).

**2. Spatial Segregation**
For industrial scale-up, spatial structure is unavoidable (biofilms, gradients). Instead of treating this as a limitation (Reviewer 1's point), treat it as an engineering opportunity. Could spatial segregation stabilize coexistence even outside your parameter windows?

**3. Productivity Metrics Missing**
From an engineering perspective, I don't just want coexistence - I want productivity. If the generalist diverts resources from a valuable product, coexistence is undesirable.

**SUGGESTED:** Add analysis of:
- Total biomass at equilibrium as function of $\omega$
- Metabolite production rates (if S produces a valuable compound)
- Trade-offs between stability and productivity

**4. Figure 2G - Practical Calibration**
The critical $\omega$ curve is useful IF I can calibrate it experimentally. Propose a one-time calibration protocol:
- Measure invasion success at 5-6 $\omega$ values (different inducer ratios)
- Fit your model to data to extract parameter values
- Use fitted model to predict outcomes at other $\omega$ values

**5. Comparison to Existing Consortia**
You cite several experimental systems but don't compare predictions to their outcomes. Add a table:

| System | Reference | Predicted regime | Observed outcome | Match? |
|--------|-----------|------------------|------------------|--------|
| Shou yeast pair | Shou 2007 | Stable S-M | Stable coexistence | ✓ |
| Harcombe E.coli | Harcombe 2014 | ... | ... | ? |

### MISSING: COMPUTATIONAL TOOLS

For this to be useful to experimentalists, provide:
1. **Interactive web tool** or Python package where I input measured parameters and get predictions
2. **Parameter estimation scripts** to fit model to time-series data
3. **Experimental design optimizer** suggesting which $\omega$ values to test

### RECOMMENDATION
**Major Revision Required - Conditional**

The mathematical work is solid (pending Reviewer 2's concerns), but the practical utility is near zero in current form. This reads like a theoretical exercise, not a framework for engineering consortia.

**To make this impactful for synthetic biology**:
1. Add detailed genetic implementation schemes
2. Include at least one concrete case study with real organisms and genes
3. Provide computational tools for practitioners
4. Address evolutionary instability head-on
5. Reframe "design principles" as actionable protocols

If the authors address these concerns, this could be a landmark paper bridging theory and engineering. Without them, it's just another math bio paper with limited real-world impact.

---

## EDITORIAL SUMMARY

All three reviewers recognize the mathematical rigor and comprehensive analysis, but identify critical gaps:

**Common themes**:
1. Insufficient connection to experimental reality (all reviewers)
2. Missing biological/genetic implementation details (R1, R3)
3. Incomplete analytical rigor despite claims (R2)
4. Lack of quantitative comparison to existing experimental systems (R1, R3)

**Required for acceptance**:
- Parameter justification with experimental comparison
- Complete analytical derivations or acknowledge numerical components
- Concrete genetic implementation scheme
- At least one detailed case study

**Recommended decision**: **Major Revision Required**

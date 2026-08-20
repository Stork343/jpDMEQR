# Pilot-gate theory decision: variance, bandwidth and Dantzig direction

**Status:** authoritative theory-side decision for the failed P01--P04 pilot gate  
**Date:** 2026-08-20  
**Applies before:** any new pilot used to satisfy the freeze gate, any `B=500/1000` main run, and confirmatory empirical inference  
**Current final-run authorisation:** **NO**

This decision responds to the three groups of questions raised after the diagnostic pilot. It is prospective: no constant below is selected to improve empirical coverage. The changes are required because the present implementation is internally inconsistent with the unsmoothed profile target/influence representation and because the current bandwidth lies on a non-vanishing Bahadur-remainder boundary.

---

## 1. Decision on the sandwich variance

### 1.1 Primary variance meat must use the **unsmoothed quantile score**

The current method specification has an internal mismatch:

- the inferential target is the **unsmoothed** regularised profile target `beta_star`;
- the intended influence representation in Section 7 already uses

\[
 g_i^\star=-m_i^{-1}X_i^{\mathsf T}\psi_\tau(r_i^\star),
\]

with the ordinary check-score `psi_tau`;
- but the current Section 8 / code estimates the meat from the smoothed score `psi_{tau,h}`.

For inference on the unsmoothed target, the primary cluster meat is henceforth

\[
\widehat g_i^{(0)}
=-\frac{1}{m_i}X_i^{\mathsf T}
\psi_\tau(\widehat r_i),
\qquad
\widehat{\bar g}^{(0)}=n^{-1}\sum_i\widehat g_i^{(0)},
\]

and, for coordinate `k`,

\[
\widehat\sigma_{0,k}^{2}
=\frac1n\sum_{i=1}^{n}
\left[
\widehat\omega_k^{\mathsf T}
\{\widehat g_i^{(0)}-\widehat{\bar g}^{(0)}\}
\right]^2,
\]

\[
\widehat{\operatorname{se}}(\widetilde\beta_k)
=\widehat\sigma_{0,k}/\sqrt n.
\]

Use `psi_tau(u)=tau-I(u<0)`. The value at exactly zero is immaterial under the continuous-error assumptions; code must nevertheless use a deterministic convention.

The **bread remains the corrected smoothed effective Hessian**. Smoothing is used to obtain differentiable optimisation and a stable local-curvature estimator; it does not change the first-order score defining the asymptotic variance of the unsmoothed target.

This is the clustered/profile analogue of the variance used in convolution-smoothed quantile-regression inference: the smoothed loss supplies the estimator/bread, while the first-order limiting variance is based on the ordinary quantile score. It is a theory correction, not a coverage calibration.

### 1.2 Smoothed-meat sandwich is retained only as a diagnostic

Define

\[
\widehat g_i^{(h)}
=-m_i^{-1}X_i^{\mathsf T}\psi_{\tau,h}(\widehat r_i).
\]

The current sandwich based on `g_i^(h)` is retained as `sandwich_smoothed_meat` for diagnostics and sensitivity, but it is no longer the primary Wald variance for `beta_star`.

If the scientific target were instead the fixed-`h` smoothed population minimiser `beta_h_star`, the first-order finite-`h` variance would be

\[
H_h^{-1}\,\operatorname{Var}\{g_{i,h}(\beta_h^\star)\}\,H_h^{-1}.
\]

There is **no universal scalar additive correction** that converts this into the variance for the unsmoothed target in the clustered/profile model. Therefore do not multiply the old SE by an empirical factor and do not add a fitted coverage correction.

### 1.3 Exact first-order smoothing deficit for the iid scalar score

For a symmetric kernel and a smooth scalar error density, the finite-`h` score second moment has the expansion

\[
E\{\psi_{\tau,h}(\varepsilon)^2\}
=\tau(1-\tau)-c_K h f_\varepsilon(0)+O(h^2),
\]

where

\[
c_K=2\int_0^1 F_K(-v)\{1-F_K(-v)\}\,dv.
\]

For the Epanechnikov kernel used by the reference implementation,

\[
c_K=\frac{9}{35}.
\]

This identity is a useful unit-test/diagnostic. It shows why a smoothed meat is systematically smaller at finite bandwidth, but it also shows why a single iid correction such as `+(9/35) h f(0)` is **not** an adequate clustered-profile variance formula. The unsmoothed cluster score automatically preserves unequal cluster sizes, within-cluster dependence, and the fitted profile residual structure.

### 1.4 Cluster-size imbalance does not require an additional scalar correction

The estimand uses equal cluster weighting and within-cluster averaging. Therefore the score contribution is exactly `m_i^{-1} X_i' psi`, and the meat is the empirical covariance of those cluster contributions. Random unequal `m_i` is already represented in the distribution of `g_i`.

No extra factor involving `mean(m_i)`, `E(1/m_i)`, or an effective row count is to be added.

The optional `n/(n-1)` correction remains a labelled degrees-of-freedom sensitivity only. It cannot explain, and is not intended to repair, a 40--50% SE deficit.

---

## 2. Decision on the nonlinear/Bahadur remainder and bandwidth

### 2.1 `h = n^{-1/3}` is no longer admissible as the primary inferential bandwidth

The current rule

\[
h=c_h n^{-1/3}
\]

is on the boundary at which

\[
nh^3=1
\]

for every `n`. Consequently, the proof-relevant nonlinear remainder associated with estimating the local curvature is not forced to vanish as `n` increases. This is consistent with the diagnostic observation that moving from `n=80` to `n=120` did not repair the first-order approximation.

The high-dimensional debiased SQR literature likewise requires the bandwidth to lie strictly between the `n^{-1/3}` and `n^{-1/4}` scales for the Bahadur remainder to vanish; an explicit necessary condition is `n h^3 -> infinity`.

### 2.2 New primary bandwidth

Use

\[
\boxed{h=c_h n^{-3/10}},
\qquad
c_h\in\{0.75,1,1.25\},
\]

with `c_h=1` as the primary inferential setting.

This exponent simultaneously gives

\[
\sqrt n\,h^2=n^{-1/10}\to0
\]

for smoothing bias and

\[
nh^3=n^{1/10}\to\infty
\]

for the nonlinear/Bahadur remainder.

All method specifications, registry helper functions, POP-H construction, tests and manuscript formulas must be synchronised to this rule before the next pilot.

### 2.3 Do not use the heuristic `O(1/(nh))` calculation as a variance correction

A Taylor remainder for the profile score is tensor-valued and its stochastic order depends on local kernel-derivative concentration, estimation error and high-dimensional sparsity. The quantity `1/(nh)` is not a general finite-sample variance term that can simply be added to the sandwich.

For the next pilot, the nonlinear effect is assessed through the **measured Bahadur decomposition**, not through an analytic SE inflation factor:

\[
R_{B,k}
=\sqrt n(\widetilde\beta_k-\beta_k^\star)
-\frac1{\sqrt n}\sum_i\xi_{ik}^{\mathrm{oracle}},
\]

using the population/Oracle direction where available. Record the mean, SD, median absolute value and 90th percentile of `R_B,k`, scaled by the first-order SD.

### 2.4 Pilot thresholds are not relaxed in advance

The current `[0.80,1.20]` baseline `mean(SE)/empirical SD` diagnostic remains in force after the theory corrections. The broad `B=200` coverage diagnostic also remains in force.

If the corrected unsmoothed meat, corrected bandwidth and corrected precision direction still fail those thresholds, the pilot fails and the source of the remaining discrepancy must be decomposed. Do not widen the thresholds to make the existing `n=80--120` cells pass.

---

## 3. Decision on the Dantzig tolerance

### 3.1 The multiplier `0.5` is not a theoretical lower bound

The old grid

```text
c_mu in {0.5, 1, 2, 4}
```

was an initial numerical grid, not a theorem-imposed lower limit. The observed P01 value `mu ~= 0.294` is too coarse to diagnose the inverse row accurately.

Replace the frozen candidate grid by

```text
c_mu in {0.10, 0.25, 0.50, 1, 2, 4}
```

with

\[
\mu(c)=c\left\{\sqrt{\log p/(nh)}+h^2\right\},
\qquad h=c_hn^{-3/10}.
\]

No candidate may be chosen by true direction error, coefficient truth, empirical bias or coverage.

### 3.2 Feasibility alone is not an adequate precision-direction gate

The constraint

\[
\|\widehat H\widehat\omega_k-e_k\|_\infty\le\mu
\]

only certifies Dantzig feasibility. It does **not** certify that `omega_hat` is close enough to the population inverse row for debiasing.

The next implementation must additionally record:

- `||omega_hat||_1` and `||omega_hat||_2`;
- adjacent-grid direction stability;
- POP-H relative `l1/l2` direction error in simulations where the asset exists;
- cosine similarity to the POP-H direction;
- the inverse-defect remainder bound defined below.

### 3.3 Exact deterministic debiasing-defect bound

Let

\[
\Delta=\widehat\beta-\beta^\star,
\qquad
\delta_k=\|\widehat H\widehat\omega_k-e_k\|_\infty.
\]

The empirical inverse-defect term satisfies exactly

\[
\left|
(e_k-\widehat H\widehat\omega_k)^{\mathsf T}\Delta
\right|
\le
\delta_k\|\Delta\|_1
\le
\mu\|\Delta\|_1.
\]

After root-`n` scaling,

\[
\sqrt n\left|
(e_k-\widehat H\widehat\omega_k)^{\mathsf T}\Delta
\right|
\le
\sqrt n\,\mu\|\Delta\|_1.
\]

With `||Delta||_1 = O_p{s sqrt(log p/n)}`, this contribution is

\[
O_p\{s\,\mu\sqrt{\log p}\}
=
O_p\left\{
\frac{s\log p}{\sqrt{nh}}+s h^2\sqrt{\log p}
\right\}.
\]

This is the relevant bound for judging the Dantzig grid. The fact that `mu` converges more slowly than `n^{-1/2}` is not itself a defect; what must vanish is the **product** with the l1 estimation error after root-`n` scaling.

If `omega_k^star` is feasible and the population inverse has bounded row-`l1` norm `M`, standard CLIME/Dantzig arguments give, up to constants,

\[
\|\widehat\omega_k-\omega_k^\star\|_\infty
\lesssim M\mu,
\]

and under row sparsity the `l1` error is of the corresponding sparse order. Thus a large feasible `mu` directly permits a large debiasing error.

### 3.4 Primary non-truth-based selection rule for `c_mu`

For the final method, do **not** select `c_mu` by empirical coverage or POP-H error. Use a deterministic cluster-level inverse-Hessian cross-validation rule over the fixed grid.

For each target coordinate and each candidate `c`:

1. fit the profile estimator once on the analysis sample;
2. partition independent clusters into four deterministic folds;
3. form training and validation averages of the per-cluster effective-Hessian contributions at the fitted value;
4. solve the Dantzig row on the training Hessian using `mu(c)`;
5. evaluate on the held-out Hessian the inverse quadratic loss

\[
L_{f,k}(c)
=\frac12\omega_{-f,k}(c)^{\mathsf T}
H_f\omega_{-f,k}(c)-e_k^{\mathsf T}\omega_{-f,k}(c);
\]

6. average over folds;
7. choose the **largest** feasible multiplier within one standard error of the minimum validation loss;
8. solve the row once more on the full Hessian using the selected multiplier.

All folds and tie rules are deterministic and stored. This rule uses neither simulation truth nor interval coverage. Because the candidate set is finite and every candidate has the same theoretical rate, this tuning changes constants, not the required order of `mu`.

For a temporary diagnostic implementation, the expanded grid may first be run with every candidate reported side-by-side, but no final pilot gate is issued until the non-truth-based selection rule above is implemented.

### 3.5 New simulation-only row-accuracy gate

Where a POP-H asset is available, add the following theory diagnostics:

\[
E_{2,k}=\frac{\|\widehat\omega_k-\omega_k^{\mathrm{pop}}\|_2}
{\|\omega_k^{\mathrm{pop}}\|_2},
\]

\[
E_{1,k}=\frac{\|\widehat\omega_k-\omega_k^{\mathrm{pop}}\|_1}
{1+\|\omega_k^{\mathrm{pop}}\|_1},
\]

and the normalized inverse-defect bound

\[
D_k=
\frac{\sqrt n\,\delta_k\|\widehat\beta-\beta^\star\|_1}
{\sigma_{0,k}^{\mathrm{pop}}}.
\]

These are **diagnostics/gates**, never tuning criteria. For P01, the next pilot is not accepted if the inverse-direction diagnostics are grossly non-first-order (in particular, if the median `D_k` exceeds `0.5` or its 90th percentile exceeds `1`). A failure triggers further precision-method review rather than a wider Wald interval.

---

## 4. Asset reuse and a newly identified POP-H bandwidth bug

### 4.1 If only `mu` changes

Mathematically, neither the unsmoothed profile target nor the population effective direction depends on the analysis Dantzig tolerance. Therefore a **pure `mu`-grid change** does not alter the numerical target/direction objects.

However, the current freeze implementation keys assets to the full registry checksum and current commit. Such assets will be marked stale even when their mathematical dependencies are unchanged. Do not bypass that protection silently. Either rebuild them or replace the full-registry identity by a documented dependency-specific asset hash before reuse.

### 4.2 Because this decision also changes `h`, population-direction assets must be rebuilt

The present decision changes the inferential bandwidth from `n^{-1/3}` to `n^{-3/10}`. Population effective directions depend on the analysis curvature and hence on `h`. Consequently:

- **profile-target assets:** target the unsmoothed parameter and can in principle be revalidated/reused if their DGP/target construction is unchanged;
- **POP-H / population-direction assets:** must be rebuilt under the new analysis bandwidth.

### 4.3 Fix `approximate_population_direction_v2` before rebuilding

The current function overwrites

```text
cfg$n_clusters <- n_population
```

and then computes

```text
h_analysis <- h_multiplier * cfg$n_clusters^(-1/3)
```

so the stored population direction is evaluated at a bandwidth determined by `n_population`, not by the scenario's analysis sample size. This is not the intended finite-`h` oracle direction.

Before rebuilding POP-H assets, preserve

```text
n_analysis <- config$n_clusters
```

and use

\[
h_{\mathrm{analysis}}
=c_h n_{\mathrm{analysis}}^{-3/10}.
\]

The large `n_population` controls Monte Carlo accuracy of `H_pop`; it must not change the bandwidth of the finite-sample scenario being approximated.

The profile-target approximation may continue to use its own small target-approximation bandwidth because its purpose is to approximate the unsmoothed target, but that target bandwidth must be explicitly independent of the analysis `c_mu` and should not be conflated with the POP-H bandwidth.

---

## 5. Effect on pilot gate and Module A

### 5.1 Pilot gate

After implementing all changes above, rerun P01--P04 from scratch. Do not reuse the previous pilot gate result.

Retain the existing baseline diagnostics:

```text
mean(SE) / empirical SD in [0.80, 1.20]
95% coverage diagnostic in [0.88, 0.98] at B=200, interpreted with MCSE
Dantzig feasibility >= 0.95
profile identity thresholds unchanged
```

Add the inverse-defect/row-accuracy diagnostics from Section 3.5 and separate reporting of `unsmoothed_meat_SE` versus `smoothed_meat_SE`.

### 5.2 Module A `n=100`

Do **not** delete or raise the `n=100` cells before the corrected pilot. They are valuable as finite-sample theory-scaling/stress cells.

If corrected inference is still poor at `n=100`, report it as a finite-sample boundary and base positive near-nominal claims on the sample sizes where the data support them. Never describe 0.55 SE/SD or 70% coverage as "close to nominal" and never remove the `n=100` cell solely because it performs poorly.

---

## 6. Required implementation sequence before another pilot

1. Synchronise `docs/METHOD_SPECIFICATION.md`:
   - unsmoothed cluster meat is primary;
   - smoothed meat is sensitivity;
   - `h=c_h n^{-3/10}`;
   - expanded `c_mu` grid and the inverse-Hessian CV rule.
2. Update profile inference code to compute/store both meats and use the unsmoothed meat for primary intervals.
3. Update all bandwidth helper functions from `n^{-1/3}` to `n^{-3/10}` for analysis/inference; keep target-approximation bandwidth separately named.
4. Fix the POP-H scenario-bandwidth bug described in Section 4.3.
5. Expand `c_mu` to `{0.10,0.25,0.50,1,2,4}` and implement the cluster-level inverse-Hessian CV/tie rule.
6. Add per-cluster effective-Hessian contributions needed by the CV rule.
7. Add the exact inverse-defect and POP-H row diagnostics.
8. Update registry/config checksums and invalidate/rebuild POP-H assets.
9. Rerun strict geometry because the implementation commit changed.
10. Rerun micro-preflight.
11. Rerun P01--P04 with `B=200`.
12. Issue a new pilot-gate decision only from the corrected pipeline.

No `B=500/1000` main run is authorised until this sequence passes.

---

## 7. Bottom-line adjudication

1. **Sandwich:** current smoothed-meat primary variance is not aligned with the unsmoothed target/influence representation. Replace it by the unsmoothed cluster-score meat. Do not use a fitted multiplicative SE factor.
2. **Second order:** `h=n^{-1/3}` sits on `nh^3=1` and is not an acceptable primary debiasing bandwidth. Replace by `h=n^{-3/10}` and measure the Bahadur remainder directly.
3. **Cluster sizes:** `m_i^{-1}` weighting is already correct for the equal-cluster target; no extra cluster-size scalar correction.
4. **Dantzig:** `c_mu=0.5` is not a lower bound. Expand the grid to `{0.10,0.25,0.50,1,2,4}`, use a non-truth-based inverse-Hessian CV rule, and add a remainder-relevant row-accuracy gate.
5. **Pilot thresholds:** keep them; fix the estimator/variance/precision approximation instead of moving the goalposts.
6. **Assets:** a pure `mu` change would not alter targets/directions mathematically, but this decision changes `h`, so POP-H directions must be rebuilt. Profile targets may be revalidated if their dependency hash is unchanged.
7. **A-module n=100:** retain as a stress/scaling cell and report the observed finite-sample boundary if it remains difficult.

References informing the bandwidth/variance decision include He, Pan, Tan and Zhou (2023, *Journal of Econometrics*, DOI 10.1016/j.jeconom.2021.07.010) and Yan, Wang and Zhang (2023, *JMLR* 24:245), whose debiased convolution-SQR Bahadur theory requires `n h^3 -> infinity` and whose Wald variance is based on the ordinary quantile-score limiting variance rather than a finite-`h` smoothed-score variance.
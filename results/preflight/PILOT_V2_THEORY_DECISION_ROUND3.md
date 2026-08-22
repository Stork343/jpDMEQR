# Pilot v2 theory decision — round 3

**Status:** authoritative theory-side adjudication after the failed P05 calibration gate  
**Date:** 2026-08-23  
**Evidence basis:** `results/preflight/PILOT_V2_P05_CALIBRATION_GATE.md` and the non-conflicting round-1/round-2 evidence  
**Supersedes:** the first-stage penalty rule in `build_tuning_for_analysis_v2`; no other mathematical object is changed unless stated below  
**Applies before:** any new P01--P06 pilot, any final `B=500/1000` registry, and confirmatory GSE65391 inference  
**Final-run authorisation:** **NO**

The P05 evidence identifies a decisive failure upstream of the precision and variance layers. With the frozen unweighted penalty

\[
\lambda_\beta=\sqrt{\log(p)/n},
\]

`max|profile score(0)| < lambda_beta`, so the zero vector satisfies the first-stage KKT condition. The median first-stage error consequently remains `||beta0||_1`, and a one-step correction is being asked to repair an order-one error. The resulting failure is not evidence against the exact profile score, corrected profile Hessian, unsmoothed cluster meat, or population target. It is evidence that the first-stage penalty was not calibrated to the scale of the **actual cluster-profile score**.

This decision fixes that issue prospectively. No constant below is chosen by coefficient truth, empirical coverage, or the observed P05 bias.

---

## 1. Q1 — first-stage penalty

## 1.1 Revoke the raw unweighted rule

The following rule is revoked for every proposed-method fit:

```text
lambda_beta = 1 * sqrt(log(p)/n)
penalty_factor_j = 1 for every penalised coordinate
```

The expression has the correct broad stochastic order only when the score coordinates have a common unit scale. The repository score is

\[
\widehat g_h(\beta)
=-\frac1n\sum_{i=1}^n
\frac1{m_i}X_i^{\mathsf T}
\psi_{\tau,h}\{r_i(\beta)\},
\]

so its coordinate scale depends on within-cluster averaging, covariate scale, quantile level, nuisance profiling and the working penalty. A constant multiplying `sqrt(log p/n)` cannot be transported from an independently normalised model without score loadings.

## 1.2 New primary rule: cluster self-normalised profile-score penalty

Let `P` be the set of penalised coordinates and `U` the set of always-included coordinates. Write `p_P=|P|`. Use the same analysis bandwidth, nuisance penalty and exact profile score as the final estimator.

For a provisional coefficient vector `b`, profile each cluster nuisance effect and define the **smoothed cluster score contribution**

\[
s_i(b)
=-\frac1{m_i}X_i^{\mathsf T}
\psi_{\tau,h}\{r_i(b)\}\in\mathbb R^p.
\]

For each `j in P`, define

\[
\bar s_j(b)=n^{-1}\sum_{i=1}^n s_{ij}(b),
\qquad
\widehat\ell_j(b)
=\left[
\frac1n\sum_{i=1}^n
\{s_{ij}(b)-\bar s_j(b)\}^2
\right]^{1/2}.
\]

The loading is the empirical standard deviation of one independent-cluster score contribution. It is **not** the standard deviation of the sample mean and therefore is divided by `sqrt(n)` only in the scalar penalty level below.

Freeze

\[
\alpha_{\lambda,n}
=\frac{0.10}{\log\{\max(n,3)\}},
\qquad
q_{\lambda,n}
=\Phi^{-1}\left(1-rac{\alpha_{\lambda,n}}{2p_P}\right),
\]

\[
\lambda_{0,n}
=1.10\,\frac{q_{\lambda,n}}{\sqrt n}.
\]

No additional `sqrt{tau(1-tau)}` is inserted because the loading is computed from the actual `psi_{tau,h}` cluster score and already contains its quantile-specific scale.

### Two-pass loading calibration within the first stage

1. **Null-profile loadings.** Set `b^(0)=0`, profile the nuisance effects at `b^(0)`, and compute `ell_j^(0)=ell_j(b^(0))`.
2. **Preliminary fit.** Fit the penalised profile estimator with coordinate penalty
   \[
   \lambda_j^{(0)}=\lambda_{0,n}\widehat\ell_j^{(0)},
   \qquad j\in P,
   \]
   and zero penalty on `U`. Denote the result by `b^(1)`.
3. **One loading update.** Recompute `ell_j^(1)=ell_j(b^(1))` and set
   \[
   \widehat\ell_j^{\mathrm{final}}
   =\max\{\widehat\ell_j^{(0)},\widehat\ell_j^{(1)}\}.
   \]
   The maximum prevents an in-sample fitted-residual squeeze from reducing the penalty loading.
4. **Final first-stage fit.** Warm-start at `b^(1)` and refit with
   \[
   \boxed{
   \lambda_j
   =\lambda_{0,n}\widehat\ell_j^{\mathrm{final}},
   \qquad j\in P.
   }
   \]

This is exactly two penalty-loading passes. Do not keep iterating until a favourable support is obtained.

### Implementation through the current API

The current optimiser accepts a scalar `lambda_beta` and a vector `penalty_factor`. Implement the final fit as

```text
lambda_beta = lambda_0_n
penalty_factor_j = ell_j_final * base_penalty_factor_j
```

where `base_penalty_factor_j=0` on unpenalised coordinates and `1` on ordinary penalised coordinates.

### Degenerate loading rule

Before fitting, outcome-blind zero-variance design columns must already have been removed. Let `ell_ref` be the median positive loading among `P`. If any retained penalised coordinate has

```text
ell_j < 1e-8 * ell_ref
```

or a non-finite loading, return a `lambda_loading_degenerate` failure. Do not silently replace it by an arbitrary floor. Such a failure indicates a design/preprocessing problem.

## 1.3 Status of the old multiplier grid

`lambda_beta_multipliers` is no longer a hidden selection grid around the raw rate. Its interpretation is changed to a **prespecified sensitivity multiplier** applied to the calibrated coordinate penalties:

\[
\lambda_j(a)=a\lambda_{0,n}\widehat\ell_j^{\mathrm{final}}.
\]

The primary proposed method always uses `a=1`. Other values may appear only in explicitly registered tuning-sensitivity cells; no value is selected by truth, bias or coverage. Primary P01--P06 and primary main-registry rows should carry only the primary multiplier `1` after the registry is versioned.

## 1.4 Why this rule is theory-aligned

The relevant Lasso condition is domination of the coordinatewise empirical score at the smoothed target, after standardising by the coordinate score scale. The new rule implements that principle for the exact independent-cluster profile score. It replaces the unjustified implicit loading `ell_j=1` by an estimated cluster loading.

The normal critical value is a fixed Bonferroni/self-normalised approximation. It is computationally deterministic and does not add a per-replication multiplier bootstrap to the large simulation. A future cluster multiplier critical value would require a separate prospective decision; it is not needed for this round.

## 1.5 First-stage KKT and convergence specification

Let the effective coordinate penalty be

\[
p_j=\lambda_{0,n}\widehat\ell_j^{\mathrm{final}}
\]

for `j in P`, and let `p_ref` be the median positive `p_j`. Define coordinate KKT residuals

\[
r_j=
\begin{cases}
|g_j+p_j\operatorname{sign}(\widehat\beta_j)|,
&j\in P,\ |\widehat\beta_j|>10^{-10},\\
\max(0,|g_j|-p_j),
&j\in P,\ |\widehat\beta_j|\le10^{-10},\\
|g_j|,&j\in U.
\end{cases}
\]

Define the scale-normalised residual

\[
R_{\mathrm{KKT}}
=\max\left\{
\max_{j\in P}\frac{r_j}{p_j},
\max_{j\in U}\frac{r_j}{p_{\mathrm{ref}}}
\right\}.
\]

A final first-stage fit is accepted only if all of the following hold:

```text
R_KKT <= 1e-3
max_nuisance_gradient <= 1e-7
last beta_change <= 1e-7 * max(1, ||beta_hat||_inf)
profile objective and every coefficient are finite
```

The reference solver budget becomes

```text
max_iter >= 2000
max_backtrack >= 50
beta_tol = 1e-7
```

with warm starts from the preliminary fit. Reaching the iteration limit without all acceptance conditions is a `penalised_fit` failure. The previous `0.25` raw-multiplier fit is not accepted merely because it produced a nonzero vector; it must be rerun under the new penalty and convergence contract.

## 1.6 First-stage theory statement

Maintain the sparse first-stage rate, but state it correctly. Let `beta_h^star` be the population minimiser of the smoothed profile risk at the analysis bandwidth. Under:

1. the standardised score-domination event
   \[
   \max_{j\in P}
   \frac{|g_{h,j}(\beta_h^\star)|}{\ell_j^\star}
   \le \frac{\lambda_{0,n}}{1+\eta}
   \]
   for some fixed `eta>0`;
2. loadings bounded above and below and uniformly comparable to their population counterparts;
3. profile restricted strong convexity on the Lasso cone;
4. sparsity and moment conditions giving `s log p / n -> 0`;
5. numerical KKT error `o_p(lambda_0_n)`;

we retain

\[
\|\widehat\beta-\beta_h^\star\|_1
=O_p\left\{s\sqrt{\frac{\log p}{n}}\right\},
\qquad
\|\widehat\beta-\beta_h^\star\|_2
=O_p\left\{\sqrt{\frac{s\log p}{n}}\right\}.
\]

For the unsmoothed target,

\[
\|\widehat\beta-\beta^\star\|_1
=O_p\left\{s\sqrt{\frac{\log p}{n}}
+\|\beta_h^\star-\beta^\star\|_1\right\},
\]

which is `O_p{s sqrt(log p/n)+s h^2}` under the stated target-approximation conditions.

The constant in an `O_p` rate is not one, so `s sqrt(log p/n)` is not a hard numerical upper bound for one simulation cell. However, an exactly zero estimator with non-shrinking order-one error across `n` is incompatible with the intended first-order regime and must be fixed before assessing debiasing.

---

## 2. Q2 — Dantzig tolerance after the first-stage failure

## 2.1 Do not change the anchor rate in round 3

The round-2 anchor remains

\[
\boxed{
\mu(c)=c\left\{
\sqrt{\frac{\log p}{nh}}+h^2
\right\},
\qquad h=c_hn^{-3/10},
}
\]

with

```text
c in {0.02,0.05,0.10,0.25,0.50,1,2,4}.
```

The observation that `sqrt(n) mu` does not vanish is not by itself a contradiction. The inverse-defect term is controlled by

\[
\sqrt n\,\mu\|\widehat\beta-\beta^\star\|_1,
\]

not by `sqrt(n) mu` alone. Under the restored first-stage rate, its leading order is

\[
O_p\left\{
\frac{s\log p}{\sqrt{nh}}
+s h^2\sqrt{\log p}
\right\},
\]

which can vanish under the manuscript's sparsity/dimension conditions. With `beta_hat=0`, the required first-stage factor is order one and no admissible practical precision row can repair the expansion.

The leading `sqrt{log p/(n h)}` scale is also the natural max-norm fluctuation scale of a kernel Hessian. More refined high-dimensional SQR bounds contain additional plug-in and nonlinear terms; the present evidence does not show that replacing the anchor by `n^{-alpha}` would be theoretically correct.

## 2.2 Keep the round-2 selection criterion conditionally

Retain the round-2 non-truth selection rule:

- `n<200`: deterministic two-fold cluster CV;
- `n>=200`: deterministic four-fold cluster CV;
- primary validation loss
  \[
  \|H_{\mathrm{val}}\widehat\omega_{\mathrm{train},k}-e_k\|_\infty;
  \]
- choose the smallest mean held-out defect;
- numerical ties choose the smaller `mu`;
- re-solve on the full Hessian.

The P05 result obtained with an empty first stage is not evidence against this loss. It shows that inverse-Hessian validation is intentionally insensitive to coefficient truth and cannot diagnose a bad first-stage estimator.

## 2.3 Mandatory lambda--mu sequencing

The precision stage may begin only after the final score-calibrated first stage passes the KKT/nuisance acceptance rules. Then:

1. recompute all profiled nuisance effects at the accepted final `beta_hat`;
2. recompute the complete effective Hessian at that `beta_hat`;
3. construct CV fold Hessians from those final per-cluster Hessian contributions;
4. run the frozen `mu` selector;
5. compute the one-step estimator and both sandwich meats.

A Hessian or `mu` selected at the preliminary fit, at `beta=0`, or from a stale first-stage object is invalid.

## 2.4 Conditional review trigger

After the lambda-corrected mechanism run, reconsider the `mu` anchor only if the following remain poor jointly:

- the actual normalized inverse-defect term `A_k`;
- POP-H relative row errors and cosine similarity;
- the total normalized Bahadur remainder;
- scaling from P05 to P06.

`D_k` alone remains a conservative Hölder diagnostic and cannot trigger a rate change. No `mu` change is authorised in round 3.

## 2.5 Assets

Changing the practical first-stage penalty does not change:

- the unsmoothed regularised profile target;
- the population effective direction at fixed DGP, `tau`, `Lambda` and analysis `h`.

Profile-target and POP-H assets may therefore be reused when their dependency hashes exclude analysis `lambda` and match all mathematical dependencies. The practical pilot outputs, Hessians, selected precision rows and intervals must all be regenerated.

---

## 3. Q3 — oracle coverage and the variance layer

## 3.1 Interpretation of the P05 oracle result

The TRUE-SUPPORT result at `n=200` has negligible bias, `SE/SD=0.80--0.95`, and coverage `0.87--0.92`. With `B=200`, the binomial Monte Carlo standard error is about `0.019--0.024` over that coverage range. A coordinate estimate of `0.87` is only `0.01` below the diagnostic lower boundary and does not identify a missing deterministic variance factor.

This is acceptable finite-sample behaviour for reporting as a bounded-cluster boundary/calibration result, provided P06 and the final scaling cells show improvement or stability. It is not evidence for exact nominal calibration at `n=200`.

## 3.2 Do not reopen the primary variance formula

The primary variance remains:

```text
smoothed corrected profile bread
+ unsmoothed fitted cluster-score meat
```

No scalar multiplier, coverage-fitted adjustment, or new main variance formula is authorised. The delete-one-cluster jackknife, KC/CR2-style and MD/CR3-style quantities remain prespecified diagnostics/sensitivities. They may be reported in the supplement but are not selected coordinate by coordinate or DGP by DGP.

A high-order correction can be promoted only after a separate prospective review establishes a profile-specific formula and demonstrates stable improvement across the registered error, random-effect and cluster-size settings. P05 alone is insufficient.

## 3.3 Gate interpretation

The bands

```text
SE/SD in [0.80, 1.20]
coverage in [0.88, 0.98]
```

remain unchanged for the **practical proposed-method P05 calibration gate** after the lambda repair. They are freeze/debug criteria, not theorem statements and not publication promises.

TRUE-SUPPORT and POP-H methods are mechanism diagnostics and do not independently fail the freeze because one oracle coordinate is `0.01` below the coverage band. Their results must remain visible and be used to interpret the practical method.

P06 remains the scaling confirmation. The manuscript may say that the `n=200` oracle shows mild finite-sample undercoverage and that calibration improves with the number of independent clusters only if the P06/main results support that statement.

---

## 4. Required versioning and rerun

The lambda rule changes every practical proposed-method fit. Therefore all existing practical P01--P06 outputs are superseded for freeze purposes.

Required actions:

1. create a new versioned pilot registry, e.g. `config/simulation_pilot_v3.csv`;
2. assign new run IDs with a `pilot_v3_` prefix;
3. preserve all old raw results and reports as historical evidence;
4. implement and unit-test the score-loading calculation, two-pass fit and normalized KKT;
5. run a deterministic lambda mechanism test verifying that the raw unit-loading rule and new loaded rule differ as intended;
6. run P05 with `B=50` first, retaining the full lambda-loading/KKT/first-stage diagnostics;
7. if the first-stage mechanism is healthy, rerun P01--P06 from scratch at the frozen replication counts;
8. rerun the practical P05 hard gate and P06 scaling assessment;
9. only then resume freeze-preflight.

The old P06 may finish for historical diagnosis, but it cannot satisfy the round-3 gate because it uses the superseded lambda rule.

## 5. Required new outputs

For every proposed-method replication save:

```text
lambda_rule
lambda_alpha
lambda_normal_quantile
lambda_safety_constant
lambda_base
lambda_loading_pass0_min/median/max
lambda_loading_pass1_min/median/max
lambda_coordinate_min/median/max
zero_profile_score_max
zero_kkt_ratio_max
preliminary_kkt_normalized
final_kkt_normalized
final_kkt_absolute
first_stage_iterations
first_stage_beta_change
first_stage_nonzero_count
first_stage_l1_error        # simulations only
first_stage_l2_error        # simulations only
```

The first-stage errors are post-fit simulation diagnostics and never tune the penalty.

## 6. Final authorisation state

- Exact profile score/Hessian/target: unchanged.
- Primary unsmoothed sandwich meat: unchanged.
- Analysis bandwidth: unchanged at `c_h n^{-3/10}`.
- First-stage lambda: **changed to cluster self-normalised score loading**.
- Dantzig anchor/grid/defect-CV rule: unchanged from round 2.
- Oracle variance formula: unchanged.
- Final simulation and confirmatory application inference: still blocked.

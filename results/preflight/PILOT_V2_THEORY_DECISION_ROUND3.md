# Pilot v2 theory decision — round 3

**Status:** authoritative theory-side adjudication after the failed P05 calibration gate  
**Date:** 2026-08-23  
**Evidence:** `results/preflight/PILOT_V2_P05_CALIBRATION_GATE.md`  
**Executable specification:** `docs/METHOD_SPECIFICATION_ROUND3_AMENDMENT.md`  
**Final-run authorisation:** **NO**

P05 identifies a failure upstream of the precision and variance layers. Under the frozen rule

\[
\lambda_\beta=\sqrt{\log(p)/n},
\]

with unit penalty factors, the zero vector satisfies the first-stage KKT condition because `max|profile score(0)| < lambda_beta`. The median first-stage error therefore remains order one, and the one-step correction is being asked to repair an invalid starting estimator. This does not invalidate the exact profile score, corrected profile Hessian, unsmoothed cluster meat, population target, or POP-H definition.

No rule below is selected by coefficient truth, support truth, empirical bias, or coverage.

---

## Q1. First-stage lambda rule

### Decision Q1(a): revoke the raw unit-loaded penalty

The raw rule `sqrt(log p/n)` implicitly assumes unit coordinate score scale. That assumption is false for the repository's exact score

\[
\widehat g_h(\beta)
=-\frac1n\sum_{i=1}^n
m_i^{-1}X_i^{\mathsf T}
\psi_{\tau,h}\{r_i(\beta)\}.
\]

The penalty must be calibrated to independent-cluster score variation.

### New primary rule

Let `P` be the penalised coordinates, `U` the unpenalised coordinates, and `p_P=|P|`. At a provisional vector `b`, define the cluster score contribution

\[
s_i(b)=-m_i^{-1}X_i^{\mathsf T}
\psi_{\tau,h}\{r_i(b)\}.
\]

For `j in P`, define

\[
\bar s_j(b)=n^{-1}\sum_i s_{ij}(b),
\qquad
\widehat\ell_j(b)
=\left[n^{-1}\sum_i
\{s_{ij}(b)-\bar s_j(b)\}^2\right]^{1/2}.
\]

Freeze

\[
\alpha_{\lambda,n}=0.10/\log\{\max(n,3)\},
\]

\[
q_{\lambda,n}
=\Phi^{-1}\left(
1-\frac{\alpha_{\lambda,n}}{2p_P}
\right),
\qquad
\lambda_{0,n}=1.10\,q_{\lambda,n}/\sqrt n.
\]

No extra `sqrt{tau(1-tau)}` is inserted because the loading is calculated from the actual quantile-specific profile score.

Use exactly two loading passes:

1. set `b^(0)=0`, profile nuisance effects, and compute `ell_j^(0)`;
2. fit a preliminary estimator with coordinate penalty `lambda_0,n * ell_j^(0)`;
3. compute `ell_j^(1)` at the preliminary fit and set
   \[
   \ell_j^{\mathrm{final}}=\max\{\ell_j^{(0)},\ell_j^{(1)}\};
   \]
4. warm-start and refit with
   \[
   \boxed{p_j=\lambda_{0,n}\ell_j^{\mathrm{final}}}
   \]
   on `P`, with zero penalty on `U`.

Through the current optimiser, use

```text
lambda_beta = lambda_0_n
penalty_factor_j = ell_j_final * base_penalty_factor_j
```

where `base_penalty_factor_j` is zero for unpenalised coordinates and one for ordinary penalised coordinates.

The two-pass construction is fixed. Do not keep updating until a preferred support appears.

### Degenerate loadings

After outcome-blind zero-variance filtering, let `ell_ref` be the median positive loading. A non-finite loading or `ell_j < 1e-8 * ell_ref` causes an explicit `lambda_loading_degenerate` failure. Do not silently floor the loading.

### Status of `lambda_beta_multipliers`

The field is reinterpreted as an explicit sensitivity multiplier around the calibrated coordinate penalty:

\[
p_j(a)=a\lambda_{0,n}\ell_j^{\mathrm{final}}.
\]

The primary method uses `a=1`. Other values appear only in separately registered tuning-sensitivity cells. There is no lambda CV in the primary method and no unrecorded selection from the old grid.

### Decision Q1(b): numerical acceptance

Let `p_j` be the final coordinate penalty. The coordinate KKT residual is

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

Let `p_ref` be the median positive `p_j` and define

\[
R_{\mathrm{KKT}}
=\max\left\{
\max_{j\in P}r_j/p_j,
\max_{j\in U}r_j/p_{\mathrm{ref}}
\right\}.
\]

A final first-stage fit is accepted only if

```text
R_KKT <= 1e-3
max_nuisance_gradient <= 1e-7
last beta_change <= 1e-7 * max(1, ||beta_hat||_inf)
all coefficients and objective values are finite
```

Use at least

```text
max_iter = 2000
max_backtrack = 50
beta_tol = 1e-7
```

with a warm start from the preliminary fit. Reaching the iteration limit without satisfying the contract is a `penalised_fit` failure. The old non-converged `0.25`-multiplier fit is not acceptable evidence for the revised estimator.

### Decision Q1(c): first-stage rate

Maintain the theoretical rate, stated relative to the smoothed population profile target `beta_h_star`. Under standardised score domination, loadings uniformly comparable to their population scales, profile restricted strong convexity, sparsity/moment conditions, and numerical KKT error `o_p(lambda_0,n)`, retain

\[
\|\widehat\beta-\beta_h^\star\|_1
=O_p\{s\sqrt{\log p/n}\},
\]

\[
\|\widehat\beta-\beta_h^\star\|_2
=O_p\{\sqrt{s\log p/n}\}.
\]

The score-domination event is

\[
\max_{j\in P}
\frac{|g_{h,j}(\beta_h^\star)|}{\ell_j^\star}
\leq\frac{\lambda_{0,n}}{1+\eta}
\]

for fixed `eta>0`. For the unsmoothed target,

\[
\|\widehat\beta-\beta^\star\|_1
=O_p\{s\sqrt{\log p/n}
+\|\beta_h^\star-\beta^\star\|_1\}.
\]

Under the target approximation assumption, the second term is of order `s h^2`.

The expression `s sqrt(log p/n)` is an order, not a unit-constant finite-sample upper bound. An exactly zero estimator with non-shrinking order-one error across `n`, however, is incompatible with the intended first-order regime.

---

## Q2. Mu anchor and lambda--mu interaction

### Decision Q2(a): do not change the mu rate in round 3

Keep

\[
\mu(c)=c\left\{
\sqrt{\log p/(nh)}+h^2
\right\},
\qquad
h=c_hn^{-3/10},
\]

and

```text
c in {0.02,0.05,0.10,0.25,0.50,1,2,4}.
```

The fact that `sqrt(n) mu` does not vanish is not the relevant requirement. The inverse-defect contribution is controlled by

\[
\sqrt n\,\mu\|\widehat\beta-\beta^\star\|_1.
\]

After restoring the first-stage rate, its leading order is

\[
O_p\left\{
\frac{s\log p}{\sqrt{nh}}
+s h^2\sqrt{\log p}
\right\},
\]

which can vanish under the manuscript's sparsity and dimension conditions. With an empty first stage, the multiplying estimation error is order one, so the current P05 result cannot identify a defective mu rate.

The leading `sqrt{log p/(n h)}` term is also the natural max-norm fluctuation scale of a kernel Hessian. More complete debiased-SQR bounds contain additional plug-in and nonlinear terms; the current evidence does not justify replacing the anchor by an arbitrary `n^{-alpha}` rule.

### Decision Q2(b): keep the round-2 non-truth selector

Retain:

- two-fold cluster CV for `n<200`;
- four-fold cluster CV for `n>=200`;
- held-out loss `||H_val omega_train-e_k||_inf`;
- minimum mean held-out defect;
- smaller `mu` under a numerical tie;
- a full-Hessian resolve at the selected value.

The loss is intended to assess the inverse equation and therefore is not supposed to diagnose whether the first-stage coefficient estimator is empty.

### Decision Q2(c): mandatory sequencing

The precision stage may run only after the final score-calibrated first stage passes its KKT and nuisance criteria. Then recompute:

1. every nuisance profile at the accepted final `beta_hat`;
2. every per-cluster Hessian contribution;
3. the full effective Hessian and fold Hessians;
4. the mu path and selected precision row;
5. the one-step estimator and both sandwich meats.

A Hessian inherited from `beta=0`, the preliminary fit, or a stale object is invalid.

Reopen the mu rate only if the lambda-corrected diagnostics remain poor jointly in actual `A_k`, POP-H row errors, total Bahadur remainder, and P05-to-P06 scaling. `D_k` alone remains a conservative bound and is not a hard gate.

Changing practical lambda does not change the profile target or POP-H direction. Those assets may be reused under matching dependency hashes; every practical fit, Hessian, selected row, interval, and summary must be regenerated.

---

## Q3. Oracle coverage and variance gate

### Decision Q3(a): acceptable as mild finite-sample undercoverage

At `B=200`, coverage estimates in the range `0.87--0.92` have binomial Monte Carlo standard errors of roughly `0.019--0.024`. With negligible oracle bias and `SE/SD=0.80--0.95`, two coordinates at `0.87` do not identify a missing deterministic variance factor. The result is acceptable for the paper as a transparent finite-sample boundary/calibration result, conditional on P06 and the main scaling sequence.

### Decision Q3(b): no new production correction

The primary variance remains

```text
corrected smoothed profile bread
+ unsmoothed fitted cluster-score meat.
```

No scalar multiplier, coverage-fitted correction, or alternative main sandwich is authorised. Delete-one-cluster jackknife, KC/CR2-style, and MD/CR3-style quantities remain diagnostics/sensitivities. A correction may be promoted only after a separate prospective derivation and stable evidence across the registered error, random-effect, and cluster-size settings.

### Decision Q3(c): gate interpretation

The practical proposed-method P05 gate remains

```text
SE/SD in [0.80,1.20]
coverage in [0.88,0.98].
```

The bands are freeze/debug criteria, not theorem statements and not publication promises. TRUE-SUPPORT and POP-H are mechanism diagnostics and do not independently fail the freeze because an oracle coordinate is `0.01` below the band. P06 remains the scaling confirmation.

The manuscript may state that the `n=200` oracle exhibits mild undercoverage and that calibration improves with the number of independent clusters only if P06 and the final scaling results support that wording.

---

## Required versioning and execution

The lambda change invalidates every practical P01--P06 output for freeze purposes.

1. Create a new versioned pilot registry, such as `config/simulation_pilot_v3.csv`.
2. Use new `pilot_v3_*` run IDs.
3. Preserve old raw outputs and reports as historical evidence.
4. Implement and test the score loadings, two-pass fit, penalty-factor mapping, and normalised KKT.
5. Run P05 with `B=50` as a lambda-mechanism check first.
6. If the first stage is non-degenerate and passes the numerical contract, rerun P01--P06 from scratch at the frozen replication counts.
7. Apply the unchanged practical P05 gate and P06 scaling assessment.
8. Resume strict geometry, micro-preflight, and freeze-preflight only after the new pilot decision.

The old P06 may finish for historical diagnosis but cannot satisfy the round-3 freeze gate.

For every proposed-method replication save the lambda rule, critical value, safety constant, loading ranges, coordinate penalty ranges, zero-profile score/KKT ratios, preliminary/final normalised KKT, absolute KKT, iteration count, final coefficient change, support size, and simulation-only first-stage `l1/l2` errors.

---

## Authorisation state

- Exact profile score/Hessian/target: unchanged.
- Primary unsmoothed sandwich meat: unchanged.
- Analysis bandwidth: unchanged at `c_h n^{-3/10}`.
- First-stage lambda: changed to cluster self-normalised profile-score loadings.
- Dantzig anchor/grid/selector: unchanged from round 2.
- Oracle variance formula: unchanged.
- Final simulation and confirmatory application inference: still blocked.

# Method specification — round-five pilot amendment

**Status:** authoritative amendment pending incorporation into `docs/METHOD_SPECIFICATION.md`  
**Date:** 2026-08-24  
**Source decision:** `results/preflight/PILOT_V4_THEORY_DECISION_ROUND5.md`  
**Scope:** post-L1 profile refit, inferential starting estimator, post-refit numerical acceptance, and round-five mechanism diagnostics

This amendment does **not** alter the exact profile score, corrected effective Hessian, unsmoothed regularised profile target, round-three score-loaded L1 selection rule, round-four `h_est`, round-two Dantzig program, round-four `h_inf`, or the primary unsmoothed sandwich meat.

It changes the coefficient vector at which the inferential reprofile/one-step stage begins.

## 1. Distinguish the selection estimator from the inferential starting estimator

The round-four penalised estimator remains

\[
\widehat\beta^{L1}
\]

and is used for high-dimensional selection/screening and its audit fields.

It is no longer passed directly to the inferential reprofile.

Let

```text
P = penalised coordinates
U = always-included unpenalised coordinates
T = prespecified inference target coordinates from the registry
```

and define

\[
\widehat S
=\{j\in P:|\widehat\beta^{L1}_j|>10^{-10}\},
\]

\[
\boxed{
S_R=\widehat S\cup U\cup T.
}
\]

`T` is fixed before simulation/data analysis. It is never chosen from coefficient truth or p-values.

## 2. Refit-dimension gate

Let `d_R=|S_R|`. The post-refit is attempted only when

\[
\boxed{
d_R\le
\left\lfloor
\frac{n}{\log\{\max(n,3)\}}
\right\rfloor.
}
\]

If the condition fails, return

```text
failure_stage = post_refit_dimension
status = failed
```

and retain the failure in the primary denominator. There is no silent fallback to the penalised starting estimator.

This gate is outcome-independent given the selected set and guarantees that the unpenalised refit dimension remains `o(n)` in the intended regime.

## 3. Post-L1 profile refit at the inference bandwidth

The refit uses

\[
h_{\rm inf}=c_h n^{-3/10}
\]

with the already frozen primary `c_h=1`.

For a vector `b in R^{d_R}`, define the reduced profiled criterion

\[
\widehat Q^{R}_{h_{\rm inf}}(b)
=\frac1n\sum_{i=1}^n
\min_{\gamma_i\in\mathbb R^q}
\left[
\frac1{m_i}\sum_{j=1}^{m_i}
\rho_{\tau,h_{\rm inf}}
\{Y_{ij}-X_{ij,S_R}^{\mathsf T}b-Z_{ij}^{\mathsf T}\gamma_i\}
+\frac12\gamma_i^{\mathsf T}\Lambda\gamma_i
\right].
\]

The post-refit is

\[
\boxed{
\bar\beta_{S_R}
=\arg\min_{b\in\mathbb R^{d_R}}
\widehat Q^R_{h_{\rm inf}}(b),
}
\]

with **zero L1 penalty**. Set

\[
\bar\beta_j=0,\qquad j\notin S_R.
\]

Warm-start `b` from the corresponding coordinates of `beta^L1`, but the final solution is the unpenalised reduced-model optimum.

The nuisance ridge `Lambda` is unchanged and remains part of every profiled cluster criterion.

## 4. Post-refit numerical acceptance

At the final reduced-model solution require:

```text
max absolute reduced-model profile gradient <= 1e-7
max nuisance gradient <= 1e-7
last beta change <= 1e-8 * max(1, ||beta_refit||_inf)
all coefficients/objectives finite
```

and, for the symmetric reduced effective Hessian,

```text
minimum eigenvalue > 1e-8
condition number < 1e10
```

Reference solver budget:

```text
max_iter >= 1000
max_backtrack >= 50
```

Failure to satisfy any condition gives

```text
failure_stage = post_refit
status = failed
```

with no silent fallback.

## 5. Full inferential reprofile after post-refit

After `beta_refit` passes:

1. embed `beta_refit` in the full `p`-vector with zeros outside `S_R`;
2. discard all reduced-model nuisance/Hessian objects for the main precision calculation;
3. using the full design matrix `X`, reprofile every nuisance effect at `beta_refit` with `h_inf`;
4. recompute the exact full profile score;
5. recompute the full effective Hessian and all deterministic CV-fold Hessians;
6. run the unchanged Dantzig selector;
7. form the unchanged one-step coordinate estimator;
8. compute the unchanged unsmoothed primary meat and smoothed diagnostic meat.

Thus the primary practical estimator is

\[
\widetilde\beta_k
=\bar\beta_k
-\widehat\omega_k^{\mathsf T}
\widehat g_{h_{\rm inf}}(\bar\beta).
\]

The first-stage L1 vector remains in the result object as `beta_l1`; the inferential start is stored separately as `beta_refit`.

## 6. Theoretical interpretation

Let `S_star` denote the active coordinates of the relevant sparse approximation and assume, for the strengthened post-refit statement,

\[
S_\star\subseteq S_R,
\qquad d_R=O(s).
\]

Under local reduced-model profile curvature, regular density/design conditions, bounded cluster sizes and the same smoothing-target conditions used elsewhere, the post-refit has schematic low-dimensional rate

\[
\|\bar\beta-\beta^\star\|_2
=O_p\left\{
\sqrt{\frac{d_R}{n}}+\sqrt{d_R}\,h_{\rm inf}^2
\right\},
\]

and

\[
\boxed{
\|\bar\beta-\beta^\star\|_1
=O_p\left\{
\frac{d_R}{\sqrt n}+d_Rh_{\rm inf}^2
\right\}.
}
\]

If `S_star` is not contained in `S_R`, add an omitted-signal approximation term such as

\[
\|\beta^\star_{S_\star\setminus S_R}\|_1
\]

under the corresponding reduced-model expansion. The exact constant/order for misspecified selected models must be stated as an assumption/lemma in the manuscript rather than silently omitted.

The post-refit therefore reduces shrinkage bias when the selected model has sure inclusion and controlled size. It does not make model-selection mistakes disappear.

This is the profile analogue of the post-L1 quantile refit studied by Belloni and Chernozhukov in high-dimensional sparse quantile regression.

## 7. Consequence for the inverse-defect term

The algebraic precision remainder is

\[
R_{\rm inv,k}
=(e_k-\widehat H\widehat\omega_k)^{\mathsf T}
(\bar\beta-\beta^\star).
\]

With the frozen Dantzig residual scale, under sure inclusion and `d_R=O(s)`,

\[
\sqrt n|R_{\rm inv,k}|
\le
\sqrt n\,\mu
\|\bar\beta-\beta^\star\|_1
\]

has schematic order

\[
O_p\left\{
 d_R\mu
+\sqrt n\,d_R\mu h_{\rm inf}^2
\right\}.
\]

This is the round-five mechanism target. Do not change the `mu` rate before testing whether the large penalised-start `Delta` was the dominant finite-sample factor.

## 8. Dantzig rule remains unchanged

Keep

\[
\mu(c)=c\left\{
\sqrt{\log p/(nh_{\rm inf})}+h_{\rm inf}^2
\right\},
\]

with

```text
c_mu in {0.02,0.05,0.10,0.25,0.50,1,2,4}
```

and the existing two-/four-fold held-out inverse-defect CV selector.

A smaller asymptotic rate is not authorised in round five. The leading kernel-Hessian concentration term remains `sqrt(log p/(n h_inf))`; existing iid debiased-SQR theory uses a sufficient tolerance containing this term plus additional plug-in terms rather than a faster rate.

## 9. POSTREFIT-EXACT-H diagnostic

For simulations only, add a data-selected reduced-model exact-inverse diagnostic without truth.

At `beta_refit`, compute the reduced effective Hessian

\[
H_R=\widehat H_{\rm eff,S_R,S_R}(\bar\beta;h_{\rm inf}).
\]

When the numerical gates in Section 4 hold, define for target `k in S_R`

\[
\omega^{R,\rm exact}_{k,S_R}=H_R^{-1}e_{k,S_R},
\qquad
\omega^{R,\rm exact}_{k,S_R^c}=0.
\]

Form the corresponding one-step coordinate using the **full-p score** evaluated at `beta_refit`.

Label this method/diagnostic

```text
POSTREFIT-EXACT-H
```

and never use its output to select `mu`, lambda, support, or a manuscript primary result.

It separates residual global-Dantzig regularisation from support/refit error.

## 10. Required round-five outputs

For every practical proposed-method replication store:

```text
beta_l1
beta_refit
selected_support_size
refit_set_size
refit_contains_targets
post_refit_status
post_refit_iterations
post_refit_gradient_max
post_refit_nuisance_gradient_max
post_refit_beta_change
post_refit_hessian_min_eigenvalue
post_refit_hessian_condition_number
l1_error_lasso              # simulations only
l1_error_refit              # simulations only
l2_error_lasso              # simulations only
l2_error_refit              # simulations only
omitted_active_count        # simulations only
A_k_lasso_start             # simulations only if retained
A_k_refit_start             # simulations only
bahadur_refit_scaled
```

Also retain the full `mu` candidate path for the `B=50` mechanism probe, but none of the truth-based fields may alter the selected candidate.

## 11. Round-five mechanism probe and gate

Before any full pilot, run

```text
P05-v5: n=200,p=500,s=5,tau=0.5,q=1,B=50
P06-v5: n=400,p=500,s=5,tau=0.5,q=1,B=50
```

with nested clusters/common random numbers where feasible.

The `B=50` probe is diagnostic only.

Compare:

- L1-start versus post-refit error;
- support inclusion and refit dimension;
- practical Dantzig versus POP-H versus POSTREFIT-EXACT-H versus TRUE-SUPPORT;
- actual `A_k`, `D_k`, row errors and Bahadur remainder;
- bias/coverage with explicit MCSE.

Only after review of this report may the full P01--P06 or P07 run begin.

## 12. Assets and invalidation

- profile-target assets: reusable under matching mathematical dependency hashes;
- POP-H assets: reusable because `h_inf`, target and DGP are unchanged, provided the dependency hash excludes the practical starting-estimator rule;
- round-four practical results: historical mechanism evidence only;
- all post-refits, full inferential reprofiles, Dantzig rows, intervals and summaries: regenerate under round five.

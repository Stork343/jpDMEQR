# Method specification — round-four pilot amendment

**Status:** authoritative amendment pending incorporation into `docs/METHOD_SPECIFICATION.md`  
**Date:** 2026-08-24  
**Source decision:** `results/preflight/PILOT_V3_THEORY_DECISION_ROUND4.md`  
**Scope:** first-stage curvature bandwidth, lambda/bandwidth sequencing, RSC diagnostics, and pilot scaling

This amendment changes the **first-stage estimation bandwidth only**. It does not alter the exact profile score, corrected effective Hessian, unsmoothed profile target, round-3 score-loaded lambda, one-step formula, round-2 Dantzig program, inference bandwidth, or primary unsmoothed sandwich meat.

## 1. Two bandwidths are now distinct objects

### 1.1 First-stage estimation bandwidth

Let `p_P` be the number of penalised coordinates in the fitted design. Define

\[
\boxed{
h_{\mathrm{est}}
=c_{\mathrm{est}}
\left\{\frac{\log(p_P\vee2)}{n}\right\}^{1/4}}
\]

with primary

\[
\boxed{c_{\mathrm{est}}=1}.
\]

This bandwidth is used only for:

- round-3 score-loading pass 0;
- preliminary penalised profile fit;
- round-3 score-loading pass 1;
- final penalised first-stage profile fit;
- first-stage KKT and curvature diagnostics.

The primary procedure does not select `c_est` from the data. Values `0.75` and `1.25` may appear only in separately registered bandwidth-sensitivity cells.

### 1.2 Inference bandwidth

Keep

\[
\boxed{
h_{\mathrm{inf}}=c_hn^{-3/10}}
\]

with primary `c_h=1` and the already frozen inferential bandwidth sensitivity settings.

This bandwidth is used for:

- the profile nuisance refit at the accepted `beta_hat`;
- exact inferential profile score;
- effective Hessian and fold Hessians;
- Dantzig direction selection;
- one-step correction.

The primary sandwich meat remains unsmoothed and therefore has no smoothing bandwidth.

## 2. Round-3 lambda is retained and evaluated at h_est

At a provisional coefficient `b`, compute

\[
s_i^{\mathrm{est}}(b)
=-m_i^{-1}X_i^{\mathsf T}
\psi_{\tau,h_{\mathrm{est}}}\{r_i^{\mathrm{est}}(b)\}.
\]

The round-3 loading is

\[
\widehat\ell_j(b)
=\left[n^{-1}\sum_i
\{s_{ij}^{\mathrm{est}}(b)-\bar s_j^{\mathrm{est}}(b)\}^2
\right]^{1/2}.
\]

Keep exactly

\[
\alpha_{\lambda,n}=0.10/\log\{\max(n,3)\},
\]

\[
q_{\lambda,n}
=\Phi^{-1}\left(1-
\frac{\alpha_{\lambda,n}}{2p_P}
\right),
\qquad
\lambda_{0,n}=1.10q_{\lambda,n}/\sqrt n,
\]

and the two-pass loading rule

\[
\ell_j^{\rm final}=\max\{\ell_j^{(0)},\ell_j^{(1)}\}.
\]

The final first-stage coordinate penalty is

\[
\boxed{p_j=\lambda_{0,n}\ell_j^{\rm final}}
\]

for penalised coordinates. The primary sensitivity multiplier remains `a=1`.

Do not lower `a`, perform lambda CV, or use support/coverage information in the primary procedure.

## 3. First-stage numerical acceptance is unchanged

The round-3 contract remains:

```text
weighted normalized KKT <= 1e-3
max nuisance gradient <= 1e-7
last beta change <= 1e-7 * max(1, ||beta_hat||_inf)
max_iter >= 2000
max_backtrack >= 50
beta_tol = 1e-7
all coefficients/objective values finite
```

Failure to satisfy the contract is a `penalised_fit` failure.

## 4. First-stage theoretical rate under h_est

The relevant local-RSC sufficient bandwidth order for high-dimensional smoothed quantile estimation has the form

\[
h_{\mathrm{est}}
\gtrsim
\max\left\{
C_1\sqrt{s\log p/n},
C_2s\log p/n
\right\}
\]

up to density/design/kernel constants. The inference-oriented rule `n^{-3/10}` is below this order in the current P05/P06 finite-sample cells.

The new bandwidth satisfies

\[
h_{\mathrm{est}}^2
=\sqrt{\log(p_P)/n}.
\]

Under the usual sparse regime

\[
s\lesssim\sqrt{n/\log p_P},
\]

`h_est` is compatible with the first-stage RSC lower order up to constants while tending to zero.

With the round-3 score-domination/loading assumptions, profile RSC, and numerical KKT error `o_p(lambda_0,n)`, write the first-stage bound relative to the unsmoothed target as

\[
\|\widehat\beta-\beta^\star\|_1
=O_p\left\{
 s\sqrt{\frac{\log p}{n}}
+s h_{\mathrm{est}}^2
\right\}.
\]

Because `h_est^2=sqrt(log p/n)`, this remains

\[
\boxed{
\|\widehat\beta-\beta^\star\|_1
=O_p\{s\sqrt{\log p/n}\}
}
\]

up to constants. The l2 statement is updated analogously with its smoothing-target term.

Do not infer a finite-sample sample-size threshold from this order without measuring the profile-specific curvature constants.

## 5. Hard separation between estimation and inference objects

After the final first-stage fit passes:

1. keep only `beta_hat` and first-stage audit fields;
2. discard first-stage nuisance profiles, score, Hessian, fold Hessians and curvature matrices for inferential use;
3. reprofile nuisance effects at the same `beta_hat` using `h_inf`;
4. recompute exact score/Hessian/fold Hessians at `h_inf`;
5. run the unchanged round-2 Dantzig selector;
6. form the unchanged one-step estimator;
7. compute the unchanged unsmoothed primary meat and smoothed diagnostic meat.

Code must expose both fields explicitly:

```text
h_est
h_inf
```

and tests must fail if an `h_est` Hessian is passed into the inferential precision/one-step pipeline.

## 6. First-stage RSC diagnostics

For simulations, save the following post-fit, truth-permitted diagnostics. They may not tune `h_est`, lambda, or mu:

```text
h_est
h_inf
rsc_order_lower = sqrt(s*log(p)/n)
h_est_over_rsc_order
lambda_min_Hest_SS_population
lambda_min_Hest_SS_target_sample
lambda_min_Hest_SS_fitted_sample
condition_Hest_SS_population
condition_Hest_SS_target_sample
condition_Hest_SS_fitted_sample
cone_curvature_proxy_min
```

`SS` denotes the true active block, used only in simulations. The cone proxy uses a fixed deterministic bank of sparse/cone directions generated from the replication seed before outcomes are inspected.

The population first-stage curvature may be estimated using an independent large Monte Carlo sample. This is a separate diagnostic object and is not a POP-H replacement.

## 7. Dantzig and variance rules are unchanged

Keep

\[
\mu(c)=c\{\sqrt{\log p/(nh_{\mathrm{inf}})}+h_{\mathrm{inf}}^2\},
\]

with the existing grid and two-/four-fold held-out inverse-defect selector.

Run precision estimation only after recomputing the inferential Hessian at `h_inf`.

Keep the primary variance as corrected smoothed profile bread plus unsmoothed fitted cluster-score meat. No empirical multiplier or new production small-sample correction is authorised.

## 8. Pilot v4 mechanism probe

Create a new versioned pilot registry. Before any full pilot run, execute:

```text
P05-v4: n=200, p=500, s=5, tau=0.50, q=1, B=50
P06-v4: n=400, p=500, s=5, tau=0.50, q=1, B=50
```

with identical DGPs and the round-4 method.

The probe is diagnostic. Compare first-stage error, RSC diagnostics, active-coordinate `A_k`, POP-H row error and total Bahadur remainder across `n`. Do not issue a formal coverage gate from `B=50`.

Where feasible, use common random numbers/nested clusters across P05/P06 to sharpen the scaling comparison, without changing either marginal DGP.

## 9. Pilot gate and conditional P07

The full practical P05 calibration bands remain

```text
SE/SD in [0.80,1.20]
coverage in [0.88,0.98].
```

They are freeze/debug gates, not theorem guarantees.

If the round-4 P05/P06 probe exhibits clear improvement in first-stage/RSC/Bahadur diagnostics but the subsequent full P05 gate remains below the bands, authorise the following scaling diagnostic:

```text
P07: n=800, p=500, s=5, tau=0.50, q=1,
     same baseline Gaussian random-intercept DGP,
     B=200.
```

P07 does not automatically replace P05 as the gate. A later theory/project decision is required before narrowing the manuscript's practical inference regime to larger `n`.

## 10. Assets and invalidation

- Profile-target assets: reusable under matching dependency hashes.
- POP-H assets: reusable if their dependency hash depends on `h_inf` and not on first-stage-only `h_est`.
- All round-3 practical P01--P06 fits: historical only; invalid for freeze.
- All practical first-stage fits, inferential reprofiles, Hessians, directions, intervals and summaries: rerun under pilot v4.

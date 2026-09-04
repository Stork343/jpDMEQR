# Method specification — round-six pilot amendment

**Status:** authoritative amendment pending incorporation into `docs/METHOD_SPECIFICATION.md`  
**Date:** 2026-08-25  
**Source decision:** `results/preflight/PILOT_V5_THEORY_DECISION_ROUND6.md`  
**Scope:** variance-mechanism decomposition, conditional pilot gate interpretation, and GSE65391 inferential scope

This amendment does **not** change the point estimator, exact profile score/Hessian, target, round-three lambda, round-four bandwidths, round-five post-L1 refit, Dantzig program, or the current primary asymptotic sandwich formula. It adds a required diagnostic decomposition before the production variance is declared final.

## 1. Current primary variance remains provisional

For coordinate `k`, the current primary root-`n` variance scale is

\[
\widehat V_k^{\rm fit,fit}
=\frac1n\sum_{i=1}^n
\left[
\widehat\omega_k^{\mathsf T}
\{\widehat g_i^{(0)}-\overline g^{(0)}\}
\right]^2,
\]

with standard error

\[
\widehat{\operatorname{se}}(\widetilde\beta_k)
=\sqrt{\widehat V_k^{\rm fit,fit}/n}.
\]

No scalar or additive finite-sample correction is part of the primary method in round six.

## 2. Required T=L+Q+R decomposition

In simulations with a population direction asset, define

\[
T_k=\sqrt n(\widetilde\beta_k-\beta_k^\star).
\]

Let

\[
G^\star=n^{-1/2}\sum_i g_i^\star
\]

be the target-level ordinary quantile cluster-score sum under the already frozen target/nuisance convention. Let `omega_pop,k` be the population effective direction at `h_inf`.

Define

\[
L_k=-\omega_k^{\mathrm{pop}\,\mathsf T}G^\star,
\]

\[
Q_k=-(\widehat\omega_k-\omega_k^{\mathrm{pop}})^\mathsf T G^\star,
\]

and

\[
R_k=T_k-L_k-Q_k.
\]

Across replications verify

\[
\operatorname{Var}(T)
=\operatorname{Var}(L)+\operatorname{Var}(Q)+\operatorname{Var}(R)
+2\operatorname{Cov}(L,Q)
+2\operatorname{Cov}(L,R)
+2\operatorname{Cov}(Q,R).
\]

This decomposition separates the population influence term, finite-sample direction-estimation contribution and all remaining higher-order terms.

## 3. Required variance ladder

For each replication compute:

\[
V^{\rm fit,fit}
=\frac1n\sum_i
[\widehat\omega^\mathsf T(g_i^{\rm fit}-\bar g^{\rm fit})]^2,
\]

\[
V^{\rm fit,target}
=\frac1n\sum_i
[\widehat\omega^\mathsf T(g_i^\star-\bar g^\star)]^2,
\]

\[
V^{\rm pop,fit}
=\frac1n\sum_i
[\omega^{\rm pop\,\mathsf T}(g_i^{\rm fit}-\bar g^{\rm fit})]^2,
\]

and

\[
V^{\rm pop,target}
=\frac1n\sum_i
[\omega^{\rm pop\,\mathsf T}(g_i^\star-\bar g^\star)]^2.
\]

All are root-`n` variance scales. The standard error of a coefficient divides the square root by `sqrt(n)`.

The diagnostic ratios are

\[
M_{\rm fit}
=E_{MC}(V^{\rm fit,fit})/\operatorname{Var}_{MC}(L+Q),
\]

\[
M_{\rm pop,fit}
=E_{MC}(V^{\rm pop,fit})/\operatorname{Var}_{MC}(L),
\]

\[
M_{\rm pop,target}
=E_{MC}(V^{\rm pop,target})/\operatorname{Var}_{MC}(L),
\]

and

\[
M_T=\operatorname{Var}_{MC}(T)/\operatorname{Var}_{MC}(L+Q).
\]

Use 5,000 paired bootstrap resamples over the replication index with seed `20260825` to obtain percentile 95% Monte Carlo intervals for these ratios.

## 4. Mechanism classification

- **V1: higher-order remainder.** `M_fit` is statistically compatible with 1, but `M_T>1` and the `R`/covariance terms explain the excess variance. Retain the current asymptotic sandwich and report small/moderate-`n` undercoverage as a finite-sample remainder phenomenon.
- **V2: direction plug-in deficit.** `M_fit` is below 1 while the population-direction meat ratios are compatible with 1. Do not use a scalar correction; a direction-aware variance procedure requires another method amendment.
- **V3: fitted-score squeeze.** `M_pop,fit` is below 1 while `M_pop,target` is compatible with 1. A profile-specific leverage/replication correction must be derived before a full pilot.
- **V4: mixed/inconclusive.** Increase P06-v5 diagnostic replication count from 50 to 100 and repeat the same ladder; do not resolve ambiguity with the full pilot.

The classification uses Monte Carlo uncertainty of the mechanism ratios, not interval coverage.

## 5. Conditional pilot-gate semantics

Before the variance mechanism is classified, P01--P06 and P07-B=200 remain blocked.

If V1 is supported and the production variance remains unchanged:

1. P07-v5 `B=200` is run first;
2. its calibration metrics are treated as finite-sample performance estimates, not as tuning inputs;
3. if no new implementation anomaly appears, the full versioned P01--P06 `B=200` pilot is authorised under the same estimator and variance.

Under V1, the coverage/SE-SD bands remain visible calibration diagnostics but are no longer binary implementation-correctness requirements for intentionally small-`n` stress cells. Under V2--V4 they remain blockers because the variance mechanism is not closed.

## 6. Application scope at approximately 129 independent subjects

### 6.1 High-dimensional layer

Gene-level high-dimensional Wald intervals and p-values are not confirmatory at this sample size under the current evidence.

Permitted outputs:

- subject-held-out prediction/pinball loss;
- selection stability;
- exploratory coefficient paths/rankings;
- exploratory intervals clearly labelled as non-confirmatory.

Prohibited outputs:

- confirmatory gene-level Wald claims;
- high-dimensional FDR statements based on the current pointwise intervals;
- claims that simulation has established nominal calibration at `n≈129`.

### 6.2 Prespecified low-dimensional layer

For clinical covariates and externally frozen immune modules, use a direct unpenalised profile fit with fixed-effect dimension

\[
d_{app}\le\lfloor n/\log\{\max(n,3)\}\rfloor.
\]

Do not use the high-dimensional Dantzig direction in this layer. Report:

- ordinary unsmoothed cluster-sandwich intervals;
- full delete-one-subject jackknife intervals;
- both estimates/intervals side by side.

Claim strength remains conditional on the P-app low-dimensional calibration cells below.

## 7. Application-matched simulation cells

Freeze the eligible GSE65391 subject visit-count histogram/multiset before fitting any application outcome model. Only aggregate cluster-size information is used.

### 7.1 PAPP-HD

```text
n=129
p=500
s=5
signal=0.75
rho_x=0.5
tau=0.50
q=1 random intercept
cluster sizes=frozen empirical GSE65391 visit-count design
Gaussian quantile-centred errors
B=200
```

Methods: PROFILE-DQR, POSTREFIT-EXACT-H, TRUE-SUPPORT; POP-H is optional if a valid asset is built prospectively. This cell documents the high-dimensional boundary and is not required to pass the practical coverage band.

### 7.2 PAPP-LD25 / PAPP-LD50 / PAPP-LD75

```text
n=129
d=20 total fixed-effect coordinates
s=5
tau in {0.25,0.50,0.75}
q=1 random intercept
cluster sizes=same frozen empirical design
AR1 rho_x=0.5
Gaussian quantile-centred errors
zero L1 penalty
B=200 per tau
```

For each cell compare the ordinary cluster sandwich and full delete-one-cluster jackknife. Both are prespecified and reported; neither is chosen retrospectively by coverage.

## 8. Provenance and freeze manifest

The final pilot manifest must contain a `mechanism_evidence` section with run ID, commit, config hash, replication count, report checksum and `formal_gate=false` for:

- P05/P06-v4 B=50;
- P05/P06-v5 B=50;
- P07-v5 B=50;
- round-six variance ladder;
- any P06-v5 extension used to resolve V4.

Formal pilot/freeze execution occurs only after the production variance is frozen.

## 9. Final execution order

1. round-six variance ladder;
2. if required, variance/direction method amendment and dedicated mechanism probe;
3. freeze final method code/config;
4. unit/package tests;
5. strict geometry on final commit;
6. micro-preflight on final registry;
7. authorised B=200 pilot/scaling cells;
8. final pilot decision manifest including all failures;
9. freeze-preflight and matching commit/config checksum;
10. final B=500/1000 registry.

No main run may start earlier.

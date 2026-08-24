# Pilot v4 theory decision — round 5

**Status:** authoritative theory-side adjudication after the paired P05-v4/P06-v4 mechanism probe  
**Date:** 2026-08-24  
**Evidence:** `results/preflight/PILOT_V4_PAIRED_PROBE.md` at commit `991d52616c9fec6b55d3e518b7b4221d01efe072`  
**Executable specification:** `docs/METHOD_SPECIFICATION_ROUND5_AMENDMENT.md`  
**Final-run authorisation:** **NO**

The paired probe establishes two facts simultaneously:

1. the round-four dual-bandwidth first stage is working as intended: at `n=400` the selected support is exact in all 50 replications, the RSC diagnostics are healthy and improving, and the penalised first-stage error decreases with `n`;
2. exact support does **not** imply a small initial estimation error: the median penalised `l1_error` at `n=400` is still about `2.27`. Consequently the inverse-defect term
   \[
   (e_k-\widehat H\widehat\omega_k)^\mathsf T
   (\widehat\beta-\beta^\star)
   \]
   can remain large even when variable selection is perfect.

The report's phrase “mu-slack dominates” is therefore only partly correct. The observed term is the **product** of an approximate inverse defect and a still heavily shrunken penalised starting estimator. POP-H reduces practical row-estimation error but still leaves active-coordinate `A_k` around `2--3.2`, while TRUE-SUPPORT is near calibrated. This pattern does not support shrinking the theoretical `mu` rate in isolation.

No rule below is chosen from empirical coverage.

---

## 1. Q1 — is the mu reopen condition met, and should the mu rate change?

### Decision 1.1 — the review condition is met, but a smaller asymptotic mu rate is **not authorised**

Round four required another precision review if `A_k`, POP-H row errors, total Bahadur remainder and P05-to-P06 scaling remained poor after the first-stage repair. That review condition is met.

However, the resulting decision is **not** to replace

\[
\mu(c)=c\left\{\sqrt{\frac{\log p}{n h_{\rm inf}}}+h_{\rm inf}^2\right\}
\]

by an arbitrary faster `n^{-alpha}` rate.

The leading term `sqrt{log p/(n h)}` is the max-norm fluctuation scale of a kernel Hessian. In the iid debiased SQR theory of Yan, Wang and Zhang (JMLR 2023, 24:245), the sufficient CLIME/Dantzig tolerance bound contains this term **plus additional positive plug-in/nonlinear terms** involving `s`, `h`, `log p` and `log(p vee n)`. That theory therefore provides no basis for declaring a smaller order than `sqrt{log p/(n h)}` to be uniformly feasible for the population inverse. Their numerical implementation tunes the constant by cross-validation rather than replacing the theoretical order.

The current cluster-profile problem is more complicated, not less. A smaller rate would need a new concentration proof for the exact profiled Hessian. The present simulation cannot supply that proof.

### Decision 1.2 — reduce the other factor in the inverse-defect product: add a post-L1 profile refit

The primary finite-sample modification is a **post-penalisation profile refit**. This is a bias-reduction step, not a coverage correction. Post-L1 quantile refitting has a direct precedent in high-dimensional quantile regression: Belloni and Chernozhukov analyse ordinary quantile regression on the model selected by the first-step L1 procedure and show that the post-selection estimator retains near-oracle rates while reducing shrinkage bias under suitable conditions.

The exact profile version is specified in `docs/METHOD_SPECIFICATION_ROUND5_AMENDMENT.md`.

The penalised estimator remains the screening/selection estimator. The unpenalised post-refit becomes the starting estimator for the one-step inference.

### Decision 1.3 — Dantzig anchor, candidate grid and defect-CV rule remain frozen for the round-five probe

Keep:

```text
h_inf = c_h*n^(-3/10), primary c_h=1
c_mu in {0.02,0.05,0.10,0.25,0.50,1,2,4}
n<200: deterministic 2-fold cluster CV
n>=200: deterministic 4-fold cluster CV
validation loss = ||H_val omega_train - e_k||_inf
select minimum mean defect; numerical tie -> smaller mu
```

After the post-refit, every nuisance profile and every full/fold Hessian used by the precision stage must be recomputed at the post-refit coefficient vector.

### Decision 1.4 — why post-refit addresses the measured rate bottleneck

Let `S_R` be the selected post-refit set, with dimension `d_R`, and suppose for interpretation that the true active set is contained in `S_R` and `d_R=O(s)`. A low-dimensional unpenalised profile refit at `h_inf` has the schematic rate

\[
\|\bar\beta-\beta^\star\|_1
=O_p\left\{\frac{d_R}{\sqrt n}+d_R h_{\rm inf}^2\right\}
\]

under the corresponding local profile curvature and smooth-target conditions. If relevant active variables are omitted, add the approximation term

\[
\|\beta^\star_{S^\star\setminus S_R}\|_1.
\]

Thus, under sure inclusion, the normalized inverse-defect contribution becomes

\[
\sqrt n\,\mu\|\bar\beta-\beta^\star\|_1
=O_p\left\{
 d_R\mu
 +\sqrt n\,d_R\mu h_{\rm inf}^2
\right\},
\]

rather than the penalised-start order

\[
O_p\{s\mu\sqrt{\log p}\}.
\]

This removes a `sqrt(log p)` shrinkage factor from the leading algebraic bound when model inclusion is successful. It directly targets the mechanism seen in P06, where support is exact but the penalised coefficient error remains large.

This is the quantity that must be tested before any change to the `mu` rate.

---

## 2. Q2 — should the full P01--P06 be skipped and P07 run now?

### Decision — skip the full pilot for now, but do **not** run P07-B=200 under the superseded starting estimator

The full P01--P06 run under the current penalised-start one-step method is not informative enough to justify its cost. It is therefore **not authorised**.

P07 with `B=200` under the same current method is also not authorised. It would mainly trace the already identified product `mu * ||Delta||_1` without testing the newly identified shrinkage component.

Instead run a paired round-five mechanism probe after implementing the post-refit:

```text
P05-v5: n=200, p=500, s=5, tau=0.50, q=1, B=50
P06-v5: n=400, p=500, s=5, tau=0.50, q=1, B=50
```

Use the same nested-cluster/common-random-number construction as round four where feasible.

P07 is conditionally authorised **only after** the P05/P06-v5 probe if:

- the post-refit substantially reduces `l1_error` and `A_k`;
- P06 moves materially toward POP-H/TRUE-SUPPORT calibration;
- the P05-to-P06 direction is coherent;
- but P06 still does not enter the practical gate.

In that case, first run a **P07-v5 B=50 scaling probe**, not `B=200`. A later decision will determine whether a full P07 is worth the cost.

---

## 3. Q3 — exact support but incomplete one-step correction

### Decision — the observed x00002 bias is explained by the current first-order inverse-defect remainder; it is not evidence of a new leading term

Exact support at `n=400` says that the locations of the nonzero coefficients are recovered. It does not say that the penalised coefficient vector is close enough to the target. The paired probe reports median penalised `l1_error ~= 2.27`, so the factor

\[
\Delta=\widehat\beta-\beta^\star
\]

is still large.

For a one-step estimator,

\[
\widetilde\beta_k-\beta_k^\star
= -\omega_k^\mathsf T g(\beta^\star)
+ (e_k-\widehat H\omega_k)^\mathsf T\Delta
+ R_{\rm nonlinear,k},
\]

up to sign convention. Therefore a large `A_k` at exact support is fully compatible with a non-negligible **remainder** produced by shrinkage times approximate inverse defect. It is not part of the desired first-order influence term, and it is not asymptotically acceptable while `A_k` is many standard deviations.

The post-refit probe is the discriminating experiment:

- if post-refit makes `Delta` small and practical `A_k`/Bahadur remainder collapse toward POP-H/TRUE-SUPPORT, the diagnosis is closed without changing `mu`;
- if `Delta` becomes small but practical `A_k` remains large, reopen the precision method/tolerance;
- if both practical and POP-H remain poor after post-refit, reopen `h_inf`/population-direction theory rather than merely shrinking the Dantzig constant.

### Required simulation-only exact-support diagnostic

For the round-five probe add a reduced-model exact-inverse diagnostic on the **data-selected refit set**:

1. form the post-refit set without using truth;
2. compute the reduced inferential Hessian at the post-refit estimate;
3. invert that reduced Hessian when numerically nonsingular;
4. extend the requested row by zero outside the refit set;
5. report its one-step bias/coverage as `POSTREFIT-EXACT-H`.

This is a mechanism diagnostic, not a primary method and not a tuning criterion. It distinguishes residual Dantzig regularisation from support/refit error without using the true support.

---

## 4. Primary post-L1 profile refit rule

The precise implementation is in the round-five method amendment. The essential definition is:

1. run the unchanged round-four penalised first stage at `h_est`;
2. let
   \[
   \widehat S=\{j\in P:|\widehat\beta^{L1}_j|>10^{-10}\};
   \]
3. let `T` be the target coordinates fixed by the registry before the data are generated, and `U` the always-included unpenalised coordinates;
4. define
   \[
   S_R=\widehat S\cup T\cup U;
   \]
5. require
   \[
   |S_R|\le\left\lfloor\frac{n}{\log\{\max(n,3)\}}\right\rfloor;
   \]
6. at **`h_inf`**, minimise the profiled smooth quantile criterion over coefficients in `S_R` with **zero L1 penalty**, setting all other fixed-effect coefficients to zero;
7. embed the result in `R^p`; this `beta_refit` becomes the starting estimator for the full-p score/Hessian/Dantzig/one-step stage.

No target truth, support truth, coefficient truth or coverage enters `S_R`.

If the refit dimension, curvature or numerical convergence gate fails, the practical method fails that replication. There is no silent fallback to the penalised start.

---

## 5. Pilot-gate status

The practical bands remain unchanged:

```text
SE/SD in [0.80,1.20]
coverage in [0.88,0.98]
```

They remain freeze/debug criteria rather than theorem promises.

No full pilot gate is issued from `B=50`.

After the P05/P06-v5 mechanism probe:

- if P06 becomes near calibrated while P05 remains poor because active variables are omitted from `S_R`, stop for a scope decision on the `n=200` boundary and GSE65391 high-dimensional inference;
- if both improve coherently, run the versioned full P01--P06 pilot;
- if neither improves, stop for another theory review before any P07 or main simulation.

---

## 6. Invalidation and assets

The post-refit changes the practical inference starting estimator. Therefore all round-four practical P01--P06 results are historical mechanism evidence only.

Profile-target and POP-H assets remain mathematically unchanged and may be reused under dependency hashes that correctly exclude the practical starting-estimator rule.

All practical post-refits, inferential nuisance profiles, Hessians, precision rows, one-step estimates, variances and summaries must be regenerated.

---

## 7. Authorisation state

- DGP: unchanged.
- Round-three score-loaded L1 selection rule: unchanged.
- `h_est`: unchanged.
- `h_inf`: unchanged.
- New step: **post-L1 unpenalised profile refit at `h_inf` on the data-selected set plus prespecified target/unpenalised coordinates**.
- Dantzig anchor/grid/CV selector: unchanged pending the post-refit probe.
- Primary variance: unchanged.
- Full P01--P06: blocked.
- P07-B=200: blocked.
- P05/P06-v5 B=50 paired mechanism probe: authorised after implementation/tests.
- Final simulation and confirmatory high-dimensional application inference: blocked.

# Pilot v3 theory decision — round 4

**Status:** authoritative theory-side adjudication after `PILOT_V3_P05_MECHANISM_CHECK.md`  
**Date:** 2026-08-24  
**Evidence:** `results/preflight/PILOT_V3_P05_MECHANISM_CHECK.md`  
**Executable specification:** `docs/METHOD_SPECIFICATION_ROUND4_AMENDMENT.md`  
**Final-run authorisation:** **NO**

Round 3 successfully removed the empty-fit pathology: the score-loaded first stage is nonzero and satisfies its numerical KKT contract. The remaining P05 failure is concentrated on active coordinates. The new evidence does **not** justify lowering the lambda multiplier, changing the AR(1) DGP, or changing the inferential bandwidth. It reveals that the same bandwidth is currently being asked to satisfy two different requirements:

1. enough local observations / empirical curvature for a high-dimensional penalised first stage;
2. sufficiently small smoothing bias and nonlinear remainder for the final debiased inference.

Those requirements need not use the same bandwidth.

---

## 1. Q1 — first-stage RSC / curvature regime

### 1.1 Population identification is not rejected by the frozen P05 DGP

Keep the P05/P06 DGP unchanged:

```text
rho_x = 0.5 AR(1)
signal = 0.75
m_i in {3,...,8}
q = 1 random intercept nuisance
Gaussian baseline error
```

An AR(1) design with `rho_x=0.5` is positive definite and is not a near-singular design. Under the frozen baseline, the fixed-effect covariates are centred and independent of the cluster effect/random-intercept direction. Profiling a low-dimensional random intercept changes the effective curvature but does not by itself imply loss of population slope identification.

The B=50 mechanism check is therefore **not** evidence that the DGP should be weakened by reducing `rho_x`, increasing the signal, enlarging `m_i`, or removing the nuisance effect. No R1--R10 DGP change is authorised in round 4.

### 1.2 The current single-bandwidth RSC guarantee is not certifiable at n=200 or n=400

The relevant high-dimensional smoothed-quantile theory requires the first-stage empirical loss to satisfy local restricted strong convexity. For an l1-penalised smoothed QR estimator, the known sufficient bandwidth window has the schematic lower bound

\[
h_{\rm est}
\gtrsim
\max\left\{
C_1\sqrt{\frac{s\log p}{n}},
C_2\frac{s\log p}{n}
\right\},
\]

with constants depending on the feature scale, density near the target quantile and kernel. This lower bound is a **curvature/concentration** requirement: when `h` is too small, too few residuals contribute materially to the empirical kernel Hessian, even though the population Hessian remains of constant order.

For the P05 orders, ignoring unknown constants,

\[
\sqrt{\frac{s\log p}{n}}
=\sqrt{\frac{5\log(500)}{200}}
\approx0.394,
\]

while the current inference bandwidth is

\[
h_{\rm inf}=200^{-3/10}\approx0.204.
\]

At P06,

\[
\sqrt{\frac{5\log(500)}{400}}
\approx0.279,
\qquad
h_{\rm inf}=400^{-3/10}\approx0.166.
\]

Thus neither P05 nor P06 lies in the known order-level sufficient RSC window if the **inference bandwidth is also used for the first-stage lasso**. This does not prove empirical RSC fails in every replication; it means the retained first-stage rate cannot currently be justified from the standard local-RSC argument at those sample sizes.

With fixed `s=5,p=500`, solving only the order equality

\[
n^{-3/10}=\sqrt{s\log(p)/n}
\]

gives

\[
n\asymp\{s\log p\}^{5/2}\approx5.4\times10^3,
\]

before unknown constants are considered. Therefore it is not theoretically reasonable to expect the **single-bandwidth** rule to move from an uncertified RSC regime at `n=200` to a clearly certified one merely by increasing to `n=400`.

This number is not a claimed finite-sample threshold. The exact onset depends on profile-specific constants that are not identified by the existing B=50 run. It is an order calculation explaining why P06 under the old single-bandwidth construction is not the correct diagnostic.

### 1.3 Required RSC diagnostics

The next mechanism run must record, for simulations only:

- minimum eigenvalue of the population first-stage Hessian on the true active set;
- minimum eigenvalue of the sample first-stage Hessian on the true active set at the frozen target;
- the same quantity at the final first-stage estimate;
- condition numbers of those active-set Hessians;
- a fixed random-cone curvature audit using deterministic sparse/cone directions, reported as a diagnostic only;
- the ratio `h_est / sqrt(s log(p)/n)` using truth only as a post-fit simulation diagnostic.

True support/sparsity may not enter bandwidth or lambda selection. These diagnostics distinguish population identification from finite-sample RSC concentration.

---

## 2. Q2 — lambda sensitivity multiplier a

### Decision: keep a=1 as the primary rule in round 4

Do **not** reduce the primary multiplier to `0.5` or `0.25`, and do not add CV/BIC selection of `a` in this round.

The round-3 penalty is already calibrated to the independent-cluster profile-score scale. The P05 coordinate penalty around `0.069` is not abnormally large relative to a high-dimensional maximum-score penalty; the failure was obtained while using an inference-oriented bandwidth that is below the known first-stage RSC sufficient order.

Reducing `a` before repairing the curvature regime would confound two mechanisms and would also weaken the score-domination event used in the first-stage oracle inequality.

Therefore retain

\[
p_j
=\lambda_{0,n}\widehat\ell_j^{\rm final}
\]

with primary sensitivity multiplier

\[
\boxed{a=1}.
\]

The existing registered alternative multipliers remain sensitivity-only settings. They are not used to select the primary method.

Reopen `a` only after the dual-bandwidth mechanism check if **all** of the following occur jointly:

1. the first-stage RSC diagnostics are healthy;
2. the final first-stage numerical KKT contract passes;
3. active-coordinate l1/l2 error still fails to shrink from P05 to P06;
4. the actual inverse-defect term `A_k` remains large after the improved first stage.

No automatic `a` change is authorised by this trigger; it requires another theory decision.

---

## 3. Q3 — dual bandwidth

### 3.1 Correct interpretation of smoothing

The empirical Hessian contains `K_h(r)=K(r/h)/h`. Its **expectation** is a density-weighted covariance and remains of constant order as `h -> 0`; it is not attenuated by a factor `1/h`. The finite-sample difficulty is that a smaller bandwidth concentrates curvature on fewer residuals, increasing the stochastic error of the empirical Hessian and weakening local RSC certification.

Thus the current evidence does not motivate a smaller first-stage bandwidth. It motivates a **larger first-stage bandwidth** while retaining the inference bandwidth required by the Bahadur analysis.

### 3.2 Frozen dual-bandwidth rule

Introduce two named bandwidths.

#### First-stage estimation bandwidth

Let `p_P` be the number of penalised coordinates actually entering the fit. Define

\[
\boxed{
h_{\rm est}
=c_{\rm est}\left\{\frac{\log(p_P\vee2)}{n}\right\}^{1/4}}
\]

with primary

\[
\boxed{c_{\rm est}=1}.
\]

The values `0.75` and `1.25` may be registered later as explicit bandwidth-sensitivity cells, but they are not selected from the data in the primary pilot.

This rule is non-oracle: it does not use the true sparsity `s`.

Its square is

\[
h_{\rm est}^2
=\sqrt{\frac{\log p_P}{n}},
\]

so the first-stage smoothing-target error is of the same order as the coordinate stochastic penalty scale. Under the usual sparse regime

\[
s\lesssim\sqrt{n/\log p_P},
\]

the order `h_est=(log p/n)^{1/4}` is compatible with the local-RSC lower scale `sqrt{s log p/n}` up to constants, while still tending to zero.

For P05,

\[
h_{\rm est}\approx0.420,
\]

which is of the same order as `sqrt{5 log(500)/200}=0.394`. For P06,

\[
h_{\rm est}\approx0.353,
\]

compared with the RSC scale `0.279`.

#### Inference bandwidth

Keep the existing rule unchanged:

\[
\boxed{
h_{\rm inf}=c_h n^{-3/10}},
\qquad c_h=1
\]

for the primary pilot, with the already frozen sensitivity multipliers for the bandwidth module.

### 3.3 Exact estimation/inference sequence

1. Compute round-3 score loadings using `h_est`.
2. Run both round-3 loading passes and the final penalised first-stage fit using `h_est`.
3. Apply the round-3 weighted KKT/nuisance convergence contract to that fit.
4. After the final first stage passes, **discard all first-stage profile nuisance/Hessian objects for inferential purposes**.
5. At the accepted `beta_hat`, reprofile every cluster nuisance effect using `h_inf`.
6. Recompute the exact smoothed profile score and effective Hessian using `h_inf`.
7. Construct fold Hessians using `h_inf`, run the unchanged round-2 Dantzig selector, and form the one-step estimator.
8. Use the ordinary unsmoothed fitted cluster score for the primary sandwich meat, as already frozen.

No `h_est` Hessian, nuisance fit, Dantzig row, or sandwich object may leak into the final inferential layer.

### 3.4 First-stage rate under dual bandwidth

The first-stage rate is now stated as

\[
\|\widehat\beta-\beta^\star\|_1
=O_p\left\{
 s\sqrt{\frac{\log p}{n}}
+s h_{\rm est}^2
\right\}.
\]

Since

\[
h_{\rm est}^2=\sqrt{\log p/n},
\]

this remains

\[
\boxed{
\|\widehat\beta-\beta^\star\|_1
=O_p\left\{s\sqrt{\frac{\log p}{n}}\right\}
}
\]

up to constants under the target-approximation conditions. Thus the larger estimation bandwidth does not alter the retained first-stage order.

The debiasing expansion and inference target continue to use `h_inf`; no theorem should identify the two bandwidths after this amendment.

### 3.5 Asset implications

The unsmoothed profile target does not depend on either analysis bandwidth. Existing profile-target assets remain reusable under their dependency hash.

POP-H is an inferential-direction asset and continues to use

\[
h_{\rm inf}=c_hn^{-3/10}.
\]

It does **not** need rebuilding merely because `h_est` changes, provided its dependency hash correctly excludes first-stage-only bandwidth fields and all inferential dependencies are unchanged.

A new optional first-stage curvature diagnostic asset may be created, but it is not substituted for POP-H.

---

## 4. Q4 — P05/P06 B=50 scaling probe

### Decision: authorised after implementation, but rerun P05 as well

Because `h_est` changes the practical estimator, the existing round-3 P05 B=50 output is historical mechanism evidence only. A P06-only probe would compare two different estimators and is therefore invalid.

After implementing and testing the dual-bandwidth rule, authorise a paired mechanism probe:

```text
P05-v4: n=200, p=500, s=5, B=50
P06-v4: n=400, p=500, s=5, B=50
```

using identical frozen DGP definitions, `a=1`, the round-3 lambda loading rule, `h_est=(log p_P/n)^(1/4)`, and `h_inf=n^(-3/10)`.

Where the simulation infrastructure permits, use the same replicate-level base seed and make P05 the first 200 clusters of the corresponding P06 DGP draw. This common-random-number coupling is optional but preferred because it sharpens the scaling comparison without changing either marginal DGP.

The B=50 run is a **mechanism/scaling diagnostic**, not a formal coverage gate. Record at least:

- first-stage l1/l2 errors;
- active/null coordinate bias;
- active-selection count/TPR as diagnostic only;
- active-set Hessian/RSC diagnostics listed in Section 1.3;
- `A_k`, `D_k`, POP-H row errors and total Bahadur remainder;
- lambda loadings and KKT diagnostics;
- practical and TRUE-SUPPORT SE/SD and coverage, with MCSE clearly labelled as B=50 exploratory precision.

Do not launch the full P01--P06 pilot until the paired B=50 report is reviewed.

---

## 5. Q5 — interpretation if n=200 still fails after round 4

### 5.1 P05 remains a hard debugging gate for the current manuscript scope

Do not demote the practical P05 calibration gate in advance. Keep

```text
SE/SD in [0.80,1.20]
coverage in [0.88,0.98]
```

for the full practical P05 pilot after the dual-bandwidth mechanism has been validated.

The reason is substantive: the planned empirical application has only on the order of 10^2 independent subjects. If the proposed high-dimensional inference is materially invalid even at `n=200` under the benign baseline DGP, the current manuscript cannot simultaneously claim practical small-cluster inference and use a roughly `n=129` application for confirmatory high-dimensional Wald inference.

These bands remain debugging/freeze criteria, not theorem statements and not publication promises.

### 5.2 Conditional P07 authorisation

If the round-4 P05/P06 B=50 probe shows all of the following:

1. first-stage l1/l2 error decreases materially with `n`;
2. active-set/RSC diagnostics improve or remain healthy;
3. `A_k` and total Bahadur remainder decrease;
4. P06 is materially closer to calibration than P05;

but the full P05 gate subsequently remains below the frozen calibration bands, then a new scaling cell is authorised:

```text
P07: n=800, p=500, s=5, tau=0.50, q=1,
     identical baseline Gaussian random-intercept DGP,
     B=200.
```

P07 is a **scaling diagnostic**, not an automatic replacement for P05. It may show where first-order behaviour becomes practically adequate.

If P05 fails while P06/P07 pass, a separate theory/project decision is required before changing the manuscript scope to an `n>=400` or `n>=800` inference regime. Such a scope change would also restrict the GSE65391 analysis: high-dimensional confirmatory Wald inference at approximately 129 subjects would no longer be supported, although prediction/stability and low-dimensional confirmatory analyses could remain.

### 5.3 Reporting language

Until the new evidence exists, do not write that inference is 'near nominal' at `n=200` or that it 'approaches nominal' with `n`. If the final scaling results support monotone improvement, acceptable language is:

> The smallest-cluster regime exhibits visible finite-sample distortion, while calibration improves as the number of independent clusters increases.

The exact boundary sample size must be result-derived from the frozen scaling experiment, not declared from the order calculation in Section 1.2.

---

## 6. Versioning and invalidation

The new `h_est` rule changes every practical proposed-method first stage. Therefore:

- all round-3 practical P01--P06 results are invalid for freeze purposes;
- they remain immutable historical mechanism evidence;
- create a new versioned registry, e.g. `config/simulation_pilot_v4.csv`;
- use new `pilot_v4_*` run IDs;
- do not overwrite round-3 raw directories;
- target assets remain reusable under dependency hashes;
- POP-H assets remain reusable if their inferential dependency hash is unchanged;
- every practical first-stage fit, inferential reprofile, Hessian, selected precision row, interval and summary must be regenerated.

---

## 7. Authorisation state

- DGP: unchanged.
- Exact score/Hessian/target: unchanged.
- Round-3 cluster score-loaded lambda: unchanged.
- Primary lambda multiplier: `a=1`, unchanged.
- First-stage bandwidth: **changed** to `h_est=(log(p_P)/n)^(1/4)`.
- Inference bandwidth: unchanged at `h_inf=n^(-3/10)`.
- Dantzig anchor/grid/selector: unchanged.
- Primary unsmoothed sandwich: unchanged.
- Paired P05/P06 B=50 mechanism probe: authorised after implementation/tests.
- Full P01--P06 pilot: blocked pending review of that probe.
- Final simulation and confirmatory application inference: blocked.

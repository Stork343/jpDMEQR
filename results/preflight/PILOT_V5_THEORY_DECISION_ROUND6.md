# Pilot v5 theory decision — round 6

**Status:** authoritative theory-side adjudication after the P05/P06/P07-v5 mechanism sequence  
**Date:** 2026-08-25  
**Evidence:** `results/preflight/PILOT_V5_PAIRED_PROBE.md`, including the P07-v5 addendum at commit `ca74a77cef863959b0b207b7db560966b7900305`  
**Executable specification:** `docs/METHOD_SPECIFICATION_ROUND6_AMENDMENT.md`  
**Final-run authorisation:** **NO**

Round five closed the first-stage shrinkage and inverse-defect mechanisms. At `n=800`, the post-refit is accurate, active-coordinate bias is essentially zero, and the actual inverse-defect diagnostic is small relative to earlier rounds. The remaining practical discrepancy is variance calibration: fitted-direction Wald standard errors remain below the repeated-sampling SD while POP-H/TRUE-SUPPORT controls are much closer to calibration.

The remaining discrepancy cannot yet be labelled a universal “meat bias.” The fitted practical direction is estimated from the same clusters and is random; the usual plug-in meat treats the realised direction as fixed. A second-order contribution from direction estimation and its covariance with the score can therefore appear in repeated-sampling variance while being absent from the within-replication plug-in covariance. Conversely, a nonlinear remainder of the complete one-step estimator can inflate `Var(T)` even when the fitted-direction meat correctly estimates the variance of its own first-order term. These mechanisms require different remedies.

No variance multiplier or correction below is selected from coverage.

---

## 1. Q1 — variance layer

### Decision 1.1 — the n=800 discrepancy is not acceptable as “calibrated,” but it is not enough to authorise a correction

The observed `SE/SD` range `0.72--0.87` at `n=800` corresponds to a material variance deficit. It may be reported later as finite-sample behaviour, but it is not sufficient evidence for claiming well-calibrated practical high-dimensional Wald inference at that sample size.

The `B=50` coverage values `0.84--0.94` have substantial Monte Carlo uncertainty and therefore are not used to estimate a correction factor or to change the nominal gate.

### Decision 1.2 — no analytic higher-order production correction is authorised in round 6

Do **not** add any of the following to the primary variance at this stage:

- a scalar inflation chosen from `SE/SD`;
- an additive `mu`-slack variance term;
- a generic Kauermann--Carroll/CR2 multiplier;
- a generic Mancl--DeRouen/CR3 multiplier;
- a degrees-of-freedom factor based on selected-model size;
- a correction based on the realised Dantzig norm alone.

The `mu` inverse-defect term is a remainder/bias term, not an additive variance component. KC/MD/CR-type corrections are designed for low-dimensional estimating equations with a specified leverage operator; they do not automatically account for a high-dimensional, data-dependent Dantzig direction selected from the same sample. Promoting one without a profile/direction derivation would be a new unproved method.

### Decision 1.3 — mandatory T=L+Q+R variance ladder before any full pilot

Run the variance mechanism diagnostic on the existing **P06-v5 `B=50` replications**, where the POP-H asset is already available. P05-v5 is a secondary replication of the same diagnostic. P07 does not need a new POP-H asset for this round.

For replication `r` and target coordinate `k`, let

\[
T_{rk}=\sqrt n(\widetilde\beta_{rk}-\beta_k^\star).
\]

Let `omega_pop,k` be the population direction at the scenario's `h_inf`. Let

\[
g_{ir}^{\star}
=-m_i^{-1}X_i^{\mathsf T}\psi_\tau(r_{ir}^{\star})
\]

be the ordinary quantile cluster score at the frozen target and target-level profile nuisance map, using the same target convention already used by the POP-H/Bahadur diagnostics. Define

\[
G_{rk}^{\star}=n^{-1/2}\sum_{i=1}^{n}g_{ir}^{\star}.
\]

The population first-order component is

\[
L_{rk}=-\omega_{k}^{\mathrm{pop}\,\mathsf T}G_r^{\star}.
\]

The finite-sample **direction-estimation component** is

\[
Q_{rk}
=-(\widehat\omega_{rk}-\omega_k^{\mathrm{pop}})^{\mathsf T}G_r^{\star}.
\]

Define the residual component

\[
R_{rk}=T_{rk}-L_{rk}-Q_{rk}.
\]

Across the `B=50` replications report

\[
\operatorname{Var}(T),\quad
\operatorname{Var}(L),\quad
\operatorname{Var}(Q),\quad
\operatorname{Var}(R),
\]

and every covariance term required to verify

\[
\operatorname{Var}(T)
=\operatorname{Var}(L)+\operatorname{Var}(Q)+\operatorname{Var}(R)
+2\operatorname{Cov}(L,Q)
+2\operatorname{Cov}(L,R)
+2\operatorname{Cov}(Q,R).
\]

The numerical identity error must be below `1e-10` after using the same stored replication values.

### Decision 1.4 — four within-replication meat levels

For each replication also compute the root-`n` variance scale

\[
\widehat V^{\rm fit,fit}_{rk}
=\frac1n\sum_i
\left[\widehat\omega_{rk}^{\mathsf T}
\{g_{ir}^{\rm fit}-\bar g_r^{\rm fit}\}\right]^2,
\]

where `g_i^fit` is the ordinary quantile score at the fitted post-refit residuals. This is the current primary variance scale.

Compute additionally

\[
\widehat V^{\rm fit,target}_{rk}
=\frac1n\sum_i
\left[\widehat\omega_{rk}^{\mathsf T}
\{g_{ir}^{\star}-\bar g_r^{\star}\}\right]^2,
\]

\[
\widehat V^{\rm pop,fit}_{rk}
=\frac1n\sum_i
\left[\omega_k^{\mathrm{pop}\,\mathsf T}
\{g_{ir}^{\rm fit}-\bar g_r^{\rm fit}\}\right]^2,
\]

and

\[
\widehat V^{\rm pop,target}_{rk}
=\frac1n\sum_i
\left[\omega_k^{\mathrm{pop}\,\mathsf T}
\{g_{ir}^{\star}-\bar g_r^{\star}\}\right]^2.
\]

These isolate, respectively, the fitted direction, fitted residual-score plug-in, and population-direction effects.

Compare the Monte Carlo means of these quantities with

\[
\operatorname{Var}(L),\qquad
\operatorname{Var}(L+Q),\qquad
\operatorname{Var}(T).
\]

All comparisons are on the root-`n` variance scale; do not mix these quantities with the reported standard error, which divides by `sqrt(n)`.

### Decision 1.5 — Monte Carlo uncertainty for mechanism ratios

For each coordinate compute at least the following ratios:

\[
M_{\rm fit}
=\frac{E_{MC}(\widehat V^{\rm fit,fit})}
{\operatorname{Var}_{MC}(L+Q)},
\]

\[
M_{\rm pop,fit}
=\frac{E_{MC}(\widehat V^{\rm pop,fit})}
{\operatorname{Var}_{MC}(L)},
\]

\[
M_{\rm pop,target}
=\frac{E_{MC}(\widehat V^{\rm pop,target})}
{\operatorname{Var}_{MC}(L)},
\]

and

\[
M_T=\frac{\operatorname{Var}_{MC}(T)}
{\operatorname{Var}_{MC}(L+Q)}.
\]

Estimate Monte Carlo uncertainty using **5,000 paired nonparametric bootstrap resamples of the replication index**, seed `20260825`. Report percentile 95% intervals. The bootstrap is a diagnostic for Monte Carlo uncertainty only; it does not alter any model interval.

### Decision 1.6 — interpretation branches

The next action is determined by the decomposition, not by coverage.

**V1 — higher-order estimator remainder:** if the 95% interval for `M_fit` includes 1 while `M_T` is clearly above 1 and the `R`/covariance terms explain the excess variance, retain the present asymptotic sandwich. The practical undercoverage is then a finite-sample remainder problem, not a missing meat factor.

**V2 — fitted-direction plug-in deficit:** if the upper 95% endpoint of `M_fit` is below 1 while `M_pop,fit` and/or `M_pop,target` are compatible with 1, the data-dependent precision direction contributes repeated-sampling variability not captured by the current within-replication plug-in meat. Stop before a full pilot. The next candidate remedy must be a direction-aware procedure (for example a prospectively derived cross-fitted direction variance or a replication/jackknife variance), not a scalar factor.

**V3 — fitted-score squeeze:** if `M_pop,fit` is below 1 but `M_pop,target` is compatible with 1, the fitted residual/profile-score plug-in is the principal variance deficit. Stop before a full pilot; a profile-specific leverage/leave-cluster correction must be derived and tested prospectively.

**V4 — mixed/inconclusive:** if more than one mechanism is material or the intervals cannot distinguish the branches, add **50 new P06-v5 diagnostic replications only** (for a total `B=100`) using the existing frozen method and rerun the same ladder. Do not run the full pilot to resolve an inconclusive variance mechanism.

No branch automatically authorises an empirical correction.

---

## 2. Q2 — full pilot and P07-B=200

### Decision 2.1 — neither is authorised before the variance ladder

Do not run the full P01--P06 `B=200` pilot and do not run P07 `B=200` until Section 1 is classified.

### Decision 2.2 — conditional authorisation after branch V1 only

If branch **V1** is supported and no production variance formula changes, then:

1. run **P07-v5 at B=200** first using the unchanged estimator/variance, extending the existing nested/common-random-number sequence where feasible;
2. report the `B=200` coverage, `SE/SD`, bias, `A_k`, Bahadur remainder and MCSE without changing thresholds;
3. if the run introduces no new numerical/mechanistic anomaly, the full versioned P01--P06 `B=200` pilot is authorised automatically under the same frozen method.

In branch V1, the calibration bands remain reported diagnostics but cease to be a binary **implementation-correctness** blocker for the small-`n` stress cells. The paper must report the finite-sample undercoverage and the `n` scaling rather than asserting nominal calibration at those cells.

If branch V2, V3 or V4 occurs, P07-B=200 and the full pilot remain blocked until the variance/direction mechanism receives another prospective method decision.

### Decision 2.3 — why the full pilot is deferred

P01--P04 are still required for the final manuscript and freeze record because they test other quantiles/error structures. Their cost is justified only after the variance rule that will be used to analyse them is frozen. Running them before Section 1 is resolved risks another complete rerun.

---

## 3. Q3 — GSE65391 inferential scope and P-app

### Decision 3.1 — no confirmatory high-dimensional Wald claim at approximately 129 subjects

The current simulation evidence does **not** support a statement that the practical high-dimensional Wald procedure is well calibrated at an independent-cluster sample size near 129. Therefore the GSE65391 application must not use gene-level PROFILE-DQR Wald intervals/p-values/FDR as confirmatory evidence.

This restriction is prospective and follows the simulation evidence; it does not prevent high-dimensional estimation, prediction or exploratory ranking.

### Decision 3.2 — authorised application structure

The application is split into two layers.

**High-dimensional exploratory layer**

- PROFILE-DQR selection and post-refit may be used;
- subject-level held-out pinball loss/prediction is primary;
- report selection stability, coefficient trajectories across quantiles and exploratory gene rankings;
- any gene-level Wald interval is labelled exploratory/diagnostic only and is not used for confirmatory significance or FDR control.

**Low-dimensional prespecified layer**

- use only clinical covariates and externally frozen immune modules specified before outcome association fitting;
- fit the unpenalised low-dimensional profile model directly; do not invoke the high-dimensional Dantzig direction for this layer;
- require total fixed-effect dimension
  \[
  d_{app}\le\left\lfloor n/\log\{\max(n,3)\}\right\rfloor;
  \]
- report the ordinary unsmoothed cluster-sandwich interval and a delete-one-subject jackknife interval side by side;
- no module is declared confirmatory solely because one of the two intervals excludes zero until the P-app low-dimensional calibration cell below is reviewed.

The low-dimensional layer is therefore prespecified inferential work, but its final claim strength is conditional on the application-matched simulation audit.

### Decision 3.3 — authorise application-matched simulation cells

Freeze the audited GSE65391 eligible-subject visit-count **histogram/multiset only** as the cluster-size design; no subject identifiers or outcomes enter the simulation.

Authorise the following cells before confirmatory application interpretation:

**PAPP-HD**

```text
n = 129
p = 500
s = 5
signal = 0.75
rho_x = 0.5
q = 1 random intercept
tau = 0.50
cluster sizes = frozen empirical eligible-visit-count multiset/histogram
Gaussian quantile-centred error
B = 200
```

Run `PROFILE-DQR`, `POSTREFIT-EXACT-H` and `TRUE-SUPPORT` (plus POP-H only if an asset is already economical to build). This cell documents the high-dimensional boundary; it is not required to pass the practical coverage band.

**PAPP-LD25 / PAPP-LD50 / PAPP-LD75**

```text
n = 129
d = 20 total fixed-effect coordinates
s = 5
tau in {0.25,0.50,0.75}
q = 1 random intercept
same frozen empirical cluster-size design
same baseline Gaussian/AR1 design
no L1 penalty; direct low-dimensional profile fit
B = 200 per tau
```

For each cell compare the ordinary cluster sandwich with the full delete-one-cluster jackknife of the low-dimensional estimator. These cells determine how strongly the prespecified module layer can be interpreted. No variance rule is selected by which one has better coverage after the fact; both are reported.

The P-app cells are application-scope diagnostics and do not replace the main simulation modules.

---

## 4. Q4 — freeze path and records

### Decision 4.1 — mechanism probes are formal provenance, not formal calibration gates

The following existing runs must be included in the eventual `pilot_gate.json` (or equivalent manifest) under a `mechanism_evidence` section:

- P05/P06-v4 `B=50`;
- P05/P06-v5 `B=50`;
- P07-v5 `B=50`;
- the round-six variance ladder and any P06 diagnostic extension.

For each store run ID, commit, registry/config hash, `B`, role, report path and checksum. Mark `formal_gate=false` for all mechanism probes.

### Decision 4.2 — refresh timing

Use the following order:

1. complete the round-six variance ladder;
2. if it triggers a method/variance change, implement that change and rerun its dedicated mechanism probe;
3. freeze the final method code/config version;
4. rerun unit tests and strict geometry on that final commit;
5. rerun micro-preflight on the final versioned registry;
6. execute the authorised `B=200` pilot/scaling cells;
7. write the final pilot decision/manifest, including failed calibration cells rather than suppressing them;
8. run freeze-preflight and create the matching commit/config checksum manifest;
9. only then start the final `B=500/1000` registry.

If round six branch V1 leaves the mathematical estimator unchanged, strict geometry need not be rerun merely to execute the variance ladder; it must still be refreshed on the final commit immediately before the formal pilot/freeze sequence.

### Decision 4.3 — final main-run conditions

The final registry remains blocked until all of the following are satisfied:

- the round-six variance mechanism is classified and the production variance rule is frozen;
- all required benchmark adapters/fidelity evidence pass;
- all required profile-target/POP-H assets pass dependency checks;
- source/unit/package checks pass;
- strict geometry passes on the final commit;
- micro-preflight passes on the final versioned registry;
- the authorised formal pilot/scaling runs are complete and fully accounted for;
- the pilot manifest records both successes and finite-sample failures honestly;
- the manuscript/application scope is synchronised with the observed inferential regime;
- freeze-preflight emits a matching commit/config checksum manifest.

A small-`n` coverage failure is not silently removed from the registry. Whether it is a method blocker or a reported finite-sample limitation is determined by the mechanism classification above, not by widening the coverage band.

---

## 5. Authorisation state

- Exact profile score/Hessian/target: unchanged.
- Round-3 score-loaded lambda: unchanged.
- Round-4 dual bandwidth: unchanged.
- Round-5 post-L1 profile refit: unchanged.
- Dantzig anchor/grid/selector: unchanged.
- Primary asymptotic variance: provisionally unchanged pending the round-six ladder.
- Full P01--P06 pilot: blocked pending the ladder.
- P07-B=200: blocked pending the ladder; conditionally authorised only under V1.
- GSE65391 high-dimensional confirmatory Wald inference: not authorised.
- GSE65391 exploratory high-dimensional analysis: authorised after ordinary empirical gates.
- PAPP application-matched diagnostics: authorised.
- Final `B=500/1000` simulation: blocked.

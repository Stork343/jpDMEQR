# Pilot v2 theory decision — round 2

**Status:** authoritative theory-side adjudication after the corrected P01 decomposition  
**Date:** 2026-08-21  
**Supersedes:** the Dantzig selection/gating portions of `results/preflight/PILOT_GATE_THEORY_DECISION.md`; all non-conflicting decisions in that document remain in force  
**Applies before:** any new pilot freeze decision, any `B=500/1000` main run, and confirmatory GSE65391 inference  
**Final-run authorisation:** **NO**

This decision responds to the second-round questions raised after implementing the first theory correction: unsmoothed score meat, `h=c_h n^{-3/10}`, expanded Dantzig grid, dependency-hash asset reuse and corrected POP-H analysis bandwidth. The reported P01 evidence is: the oracle/true-support estimator has negligible bias but `mean(SE)/empirical SD` remains approximately `0.72--0.75`; the prior one-SE/largest-in-band Dantzig rule selects coarse rows (typically `c_mu>=0.5`) and produces normalized inverse-defect bounds around `10--11`; `c_mu=0.10` remains feasible. These facts must be treated as diagnostic evidence, not as values to which a correction is calibrated.

---

## 1. A1: what the remaining oracle variance shortfall does and does not establish

### Decision A1.1 — do not attribute the remaining 27% to a scalar missing factor

The remaining shortfall after replacing the smoothed meat by the unsmoothed meat is **not sufficient to identify one mechanism from the ratio alone**. In particular, it does not justify multiplying standard errors by `1/0.73`, adding a fixed second-order variance term, or introducing a cluster-size factor.

The three proposed explanations have different status:

1. **Fitted-score / fitted-curvature small-sample bias is plausible and must be measured.** Small-sample sandwich estimators are well known to be downward biased because estimating equations and fitted residuals are evaluated after parameter fitting; leverage corrections such as Kauermann--Carroll and Mancl--DeRouen were developed for precisely this phenomenon. In this project, however, the relevant estimating equation is a profiled cluster M-estimating equation, so an off-the-shelf GEE multiplier cannot be assumed valid without a profile-specific derivation.
2. **A genuine nonlinear/Bahadur remainder is also plausible.** Even with an exact/population precision row, the empirical estimator contains randomness from the fitted profile curvature, nuisance profiling and evaluation away from the population target. A non-negligible remainder raises the empirical SD of the estimator but is absent from a first-order sandwich.
3. **The baseline DGP is not the leading theoretical explanation.** P01 has independent clusters, non-informative bounded cluster size and a prespecified within-feature AR(1) design. Unequal `m_i` is already part of the distribution of the equal-cluster-weight score `m_i^{-1} X_i' psi`. No additional factor in `mean(m_i)` or `E(1/m_i)` is authorised. A balanced-`m` sensitivity is nevertheless useful to verify that imbalance is not driving the finite-sample gap.

### Decision A1.2 — run an oracle variance ladder before adding any finite-sample correction

For each target coordinate in the low-dimensional TRUE-SUPPORT P01 diagnostic, construct the following objects using the same Monte Carlo replications.

Let

\[
T_k=\sqrt n(\widetilde\beta_k-\beta_k^\star).
\]

Let `omega_pop,k` be the population effective direction at the scenario's **analysis** bandwidth, and let `g_i^star` denote the corresponding cluster score evaluated at the frozen target and the target-level profiled nuisance map. Define the first-order oracle term

\[
L_k=-\frac{1}{\sqrt n}\sum_{i=1}^n
\omega_{k}^{\mathrm{pop}\,\mathsf T}g_i^\star,
\qquad
R_k=T_k-L_k.
\]

The diagnostic must report

\[
\operatorname{Var}(T_k),\quad
\operatorname{Var}(L_k),\quad
\operatorname{Var}(R_k),\quad
2\operatorname{Cov}(L_k,R_k),
\]

and verify numerically

\[
\operatorname{Var}(T_k)
=\operatorname{Var}(L_k)+\operatorname{Var}(R_k)+2\operatorname{Cov}(L_k,R_k).
\]

In addition, report four sandwich levels:

- population/oracle first-order variance from the POP-H asset;
- sample score covariance evaluated at the frozen target with the population direction;
- sample score covariance evaluated at fitted profile residuals with the population direction;
- fitted-score covariance with the sample exact inverse row in TRUE-SUPPORT.

Interpretation is frozen as follows:

- a large gap between target-score and fitted-score meat is evidence of plug-in/residual-score squeeze;
- agreement of oracle first-order variance with `Var(L_k)` but a large gap between `Var(L_k)` and `Var(T_k)` is evidence of a nonlinear/higher-order contribution;
- both gaps may coexist;
- only after this ladder is available may a finite-sample covariance correction be promoted to the main method.

### Decision A1.3 — add a leave-one-cluster-out diagnostic, not a coverage-fitted multiplier

For TRUE-SUPPORT only, add a delete-one-cluster jackknife variance diagnostic by refitting the low-dimensional profile estimator after deleting each cluster. This is a mechanism diagnostic because it estimates the variability of the complete finite-sample estimator rather than rescaling the sandwich to match coverage.

Also implement low-dimensional profile-score analogues of the Kauermann--Carroll / CR2-type and Mancl--DeRouen / CR3-type leverage corrections as **diagnostics**, not yet as the primary high-dimensional variance estimator. For a low-dimensional profile estimating equation with cluster Hessian contribution `H_i` and `A=sum_i H_i`, the relevant leverage operator is of the form

\[
\mathcal H_i=H_iA^{-1}.
\]

A KC/CR2-style score correction uses an inverse square-root of `I-mathcal H_i`; an MD/CR3-style correction uses an inverse. The precise symmetric implementation and numerical stabilisation must be documented and tested. Because cluster sizes are small, KC/CR2 is the preferred first diagnostic and MD/CR3 is a sensitivity check; neither is selected by coverage.

No high-dimensional production correction is authorised until these diagnostics show which component is responsible for the residual shortfall.

---

## 2. A2: status of the `[0.80,1.20]` SE/SD band

### Decision A2.1 — the band is a diagnostic target, not a theorem

The interval `[0.80,1.20]` is not an asymptotic theorem and is not guaranteed at `n=80` or `n=120`. The corrected evidence now shows that even the low-dimensional oracle/TRUE-SUPPORT estimator can lie below the band. Therefore failure of this ratio at `n=80` cannot, by itself, be interpreted as a defect in the high-dimensional Dantzig step.

### Decision A2.2 — do not relax the band and do not manufacture a multiplier

Keep `[0.80,1.20]` unchanged wherever it is used as an asymptotic-calibration diagnostic. Do not replace it by a wider interval and do not choose a finite-sample multiplier to force the ratio into it.

### Decision A2.3 — separate the stress pilot from the first-order calibration gate

P01--P04 remain mandatory **stress/diagnostic** cells. Their convergence, score/Hessian identities, failure accounting and variance-decomposition outputs remain hard requirements, but P01 at `n=80` no longer serves as the sole hard asymptotic-calibration gate.

Before the final freeze, add two prospectively fixed calibration cells aligned with Module A:

```text
P05: n=200, p=500, s=5, tau=0.50, q=1, baseline Gaussian random intercept
P06: n=400, p=500, s=5, tau=0.50, q=1, same DGP
```

Use the same estimator, tuning rules and target definitions as the main study. The `[0.80,1.20]` SE/SD and `[0.88,0.98]` broad coverage diagnostics are evaluated at P05 as the hard first-order calibration gate; P06 verifies the expected scaling direction. If P05 fails after the precision and variance mechanisms are repaired, the main run remains blocked. P01--P04 are still reported and may show substantial finite-sample distortion.

This change does **not** lower a threshold after observing results; it separates an intentionally difficult small-`n` stress test from the regime used to authorise asymptotic claims.

---

## 3. B1: revoke the one-SE/largest-in-band Dantzig rule

### Decision B1.1 — the previous tie rule was unsuitable for a debiasing direction

The previous rule, 'choose the largest feasible multiplier within one standard error of the minimum inverse-quadratic CV loss', is revoked. It borrowed a parsimony convention from regularisation model selection, but here increasing `mu` directly enlarges the allowed inverse defect

\[
\|\widehat H\widehat\omega_k-e_k\|_\infty.
\]

When validation loss is flat, the old rule systematically selects the coarsest admissible row, exactly the wrong tie direction for debiasing.

### Decision B1.2 — new non-truth selection rule

Use deterministic **two-fold cluster cross-validation** for `n<200` and four-fold cluster cross-validation for `n>=200`. Two folds at `n=80` leave approximately 40 validation clusters rather than 20 and reduce the variance of the row-validation criterion.

For each candidate `c_mu`, solve the row on the training-cluster Hessian. The primary validation criterion is the held-out inverse defect

\[
L_{f,k}^{\mathrm{defect}}(c)
=\|H_f\widehat\omega_{-f,k}(c)-e_k\|_\infty.
\]

The inverse-quadratic loss

\[
L_{f,k}^{\mathrm{quad}}(c)
=\frac12\widehat\omega_{-f,k}(c)^\mathsf T
H_f\widehat\omega_{-f,k}(c)-e_k^\mathsf T\widehat\omega_{-f,k}(c)
\]

is retained as a secondary diagnostic.

Select the candidate with the **smallest mean held-out inverse defect** across folds. If several candidates are tied within a fixed numerical tolerance, choose the **smallest `mu`**. Do not use a one-SE enlargement. Re-solve the selected row on the full sample.

This rule targets the population equation `H omega=e_k` directly and uses neither simulation truth nor empirical coverage.

---

## 4. B2: interpretation of the large `D_k`

### Decision B2.1 — `D_k=10--11` under the old selection rule does not prove an intrinsic small-n impossibility

The quantity

\[
D_k=\frac{\sqrt n\,\delta_k\|\widehat\beta-\beta^\star\|_1}
{\sigma_{0,k}^{\mathrm{pop}}}
\]

is a **Hölder upper bound**, not the actual inverse-defect remainder. It can be very conservative because it replaces the signed inner product by the product of an infinity norm and an `l1` norm. A median value around 10 is evidence that the current certificate is uselessly loose/coarse, but it is not proof that the actual debiasing remainder is ten standard deviations.

### Decision B2.2 — demote the old `D_k<0.5 / Q90<1` rule from hard gate to review diagnostic

The numerical thresholds `median(D_k)<0.5` and `Q90(D_k)<1` were not theorem-derived. They are no longer hard freeze gates.

In simulations, compute the **actual** normalized inverse-defect term

\[
A_k=
\frac{\sqrt n\,
|(e_k-\widehat H\widehat\omega_k)^\mathsf T
(\widehat\beta-\beta^\star)|}
{\sigma_{0,k}^{\mathrm{pop}}},
\]

as well as `D_k`, POP-H relative row errors, cosine similarity and the total normalized Bahadur remainder. `A_k` and `D_k` use truth only as post-fit diagnostics; neither may select `mu`.

A large `D_k` with small `A_k` is interpreted as looseness of the norm bound. A large `A_k` confirms a genuine precision-induced first-order failure.

The new precision review trigger is based on the joint pattern of `A_k`, POP-H row error and total Bahadur remainder, not `D_k` alone.

---

## 5. B3: expand the candidate grid downward

### Decision

For the second diagnostic pilot use

```text
c_mu in {0.02, 0.05, 0.10, 0.25, 0.50, 1, 2, 4}
```

with

\[
\mu(c)=c\left\{\sqrt{\log p/(nh)}+h^2\right\},
\qquad h=c_h n^{-3/10}.
\]

The expansion is prospective and is justified by two facts: `c=0.10` remains feasible in P01, and the previous rule systematically selected substantially larger values despite poor inverse accuracy. Adding smaller fixed constants does not alter the theoretical order of `mu`; it changes only the constant.

If `0.02` or `0.05` is infeasible for a particular row, record infeasibility and continue. Do not enlarge the grid below `0.02` after viewing coverage. Any further expansion requires another committed theory decision based on row-equation diagnostics.

---

## 6. C1: profile-target asset reuse

### Decision — confirmed

The profile-target assets approximate the **unsmoothed** regularised profile parameter. They are not finite-`h` POP-H direction assets and do not need to use the analysis bandwidth `n_analysis^{-3/10}`.

The target-construction routine may continue to use its own vanishing numerical approximation bandwidth (including the existing `n_population^{-1/3}`-scale rule), because its purpose is to make smoothing error negligible while solving the unsmoothed population target. The `n h^3 -> infinity` condition pertains to the debiased finite-sample Bahadur expansion, not to a deterministic high-accuracy population-target approximation.

Reuse is authorised when the dependency hash confirms that the DGP, `tau`, `Lambda`, nuisance specification, target submodel and target-construction code are unchanged. Changes only to the analysis `mu`, variance estimator or analysis bandwidth do not mathematically change `beta_star`.

---

## 7. C2: POP-H definition

### Decision — confirmed

The intended population direction for a finite-sample scenario is

\[
\omega_{k}^{\mathrm{pop}}(h_{\mathrm{analysis}})
=H_{\mathrm{pop}}
\{\beta^\star,h_{\mathrm{analysis}}\}^{-1}e_k,
\]

with

\[
h_{\mathrm{analysis}}=c_h n_{\mathrm{analysis}}^{-3/10}.
\]

The 100,000 population clusters are used only to estimate this expectation accurately; they do not determine `h_analysis`. The current asset fields `n_analysis`, `h_analysis`, `sigma0_pop` and `Sigma0_population` are the intended definition and should remain top-level auditable fields.

The population first-order variance for the unsmoothed-target diagnostic uses the same population direction together with the population covariance of the ordinary quantile cluster score.

---

## 8. D1/D2/D3: what to do after the precision failure and what scaling to expect

### D1 — order of response

Do **not** declare `n=80,p=200` intrinsically below the first-order regime until the defective `mu` selector is replaced. The required order is:

1. replace the largest-in-band rule by the two-/four-fold held-out inverse-defect rule above;
2. expand the candidate grid to `0.02--4`;
3. run a short P01 mechanism diagnostic with all candidates retained;
4. inspect actual `A_k`, POP-H row errors and the variance ladder;
5. only then classify the `n=80` cell as a finite-sample boundary if first-order diagnostics remain poor even at the accurately selected row.

### D2 — expected scaling is slow at these n values

Using the deterministic bound

\[
\sqrt n\,\delta_k\|\Delta\|_1
\lesssim \sqrt n\,\mu\,s\sqrt{\log p/n},
\]

and `h=n^{-3/10}`, the dominant fixed-`p,s` component of the bound behaves approximately as

\[
n^{-7/20}=n^{-0.35}
\]

(up to logarithmic and model-dependent constants). Therefore increasing `n` from 80 to 120 is expected to reduce this bound by only about 13%; 80 to 200 by about 27%; and 80 to 400 by about 43%. A value around 10 cannot plausibly fall below 1 through the 80-to-120 sample-size increase alone. This is another reason to repair row selection before interpreting the pilot as a sample-size impossibility.

The high-dimensional debiased SQR literature is also substantially more comfortable at larger sample sizes: its published simulations use `n=500`, while its Bahadur theorem requires a bandwidth range consistent with `n h^3 -> infinity`. Small-`n` behaviour therefore must be treated as a finite-sample stress problem, not assumed to inherit the asymptotic approximation automatically.

### D3 — four-fold CV is too noisy for P01

At `n=80`, four folds leave only about 20 validation clusters per fold. The observed flat validation loss is consistent with this criterion being noise-dominated. The frozen rule is changed to two folds for `n<200`. Do not use more folds to solve this problem.

A repeated two-fold CV may be run once as a sensitivity diagnostic, but the final method uses one deterministic split rule and does not average over many CV randomisations in every Monte Carlo replication.

---

## 9. E1/E2/E3: rerun cost, solver and Module A reporting

### E1 — staged rerun, not immediate full B=200

Before rerunning all P01--P04 at `B=200`, run a **P01 mechanism check with B=50** under the new `mu` rule. It must retain the full candidate path and include the variance ladder, TRUE-SUPPORT delete-one-cluster jackknife diagnostic, POP-H row diagnostics and both sandwich meats. This B=50 run is a debugging/mechanism study and cannot satisfy the final pilot gate.

If the mechanism check behaves coherently, rerun P01--P04 from scratch at `B=200`, add P05/P06, and then issue the formal pilot decision.

### E2 — HiGHS is authorised as a numerical solver, conditional on parity tests

Changing the LP solver does not change the statistical estimator when the same Dantzig linear program is solved to the same tolerance. Therefore HiGHS is authorised for production use if it passes a frozen solver-parity suite against CLARABEL on representative rows spanning low/high `p` and small/large `mu`.

Require at least 20 representative rows and verify:

- both solvers return a feasible row or agree on infeasibility;
- primal infinity residual is at most `mu+1e-6`;
- relative difference in the `l1` objective is at most `1e-5` (with a small absolute tolerance near zero);
- resulting one-step coordinates agree within a committed numerical tolerance;
- repeated runs are deterministic.

If these pass, HiGHS may be the default LP solver; CLARABEL remains a reference audit solver. No further theory approval is required for a solver-only change.

The two-fold rule at small `n` and elimination of the one-SE search also reduce the number of row solves. POP-H and exact TRUE-SUPPORT directions do not require CV row selection.

### E3 — retain `n=100` as a boundary cell

Keep the Module A `n=100` rows. The manuscript wording must be result-contingent. The prespecified language is:

> `n=100` is included as a finite-cluster boundary cell; calibration is evaluated through the full `n=100,200,400` scaling sequence.

Only if the results support it may the paper subsequently state that inference 'approaches nominal calibration as the number of clusters increases.' Do not call the `n=100` result 'near nominal' unless its reported interval/SE diagnostics actually justify that phrase.

---

## 10. F1: implication for GSE65391 with about 129 subjects

### Decision — the application remains viable, but the inferential layer must be tiered

A subject count around 129 is in the same finite-cluster range in which the pilot has exposed first-order calibration risk. This does **not** invalidate the application, because the confirmatory layer is deliberately low-dimensional and the exploratory high-dimensional layer can be separated from confirmatory inference.

The empirical protocol is amended conceptually as follows:

1. Prediction and selection/stability analyses may proceed once the computational method passes its correctness gates.
2. Confirmatory low-dimensional module/covariate inference requires an empirical-design small-sample variance audit. Report the ordinary unsmoothed-meat sandwich together with a prespecified delete-one-subject jackknife and, if the low-dimensional profile KC/CR2 diagnostic is validated in simulation, the corresponding leverage-corrected interval.
3. No correction is selected because it gives a preferred p-value or significance pattern. The primary finite-sample rule must be fixed from the simulation diagnostic before C3 associations are inspected.
4. Exploratory gene-level debiased intervals at `n≈129` are permitted only if the selected precision rows show adequate held-out inverse-defect stability and the simulation scaling evidence supports first-order behaviour in a comparable regime. Otherwise the gene-level analysis is limited to held-out prediction, selection stability and explicitly exploratory rankings; no confirmatory gene-level Wald/FDR claim is made.

Thus the GSE65391 application is not abandoned, but confirmatory inference remains blocked until this second-round pilot diagnosis is resolved.

---

## 11. Frozen decisions from round 1 that remain unchanged

The following remain authoritative:

- primary target: unsmoothed regularised profile parameter;
- exact profile score uses original `X`, not `X_tilde`;
- Schur-complement Hessian includes the ridge-curvature identity term;
- primary inferential bandwidth `h=c_h n^{-3/10}`;
- primary asymptotic sandwich meat is the ordinary quantile cluster score; smoothed meat is diagnostic;
- no scalar `m_i` correction and no coverage-fitted SE multiplier;
- POP-H direction uses the finite-sample scenario's analysis bandwidth;
- profile-target assets may be reused by dependency hash when mathematical dependencies are unchanged;
- final `B=500/1000` simulation and confirmatory empirical inference remain blocked until a fresh pilot/freeze manifest passes.

---

## 12. Immediate implementation order

1. Commit/retain the full P01 decomposition evidence file used for this round-two review. If `results/preflight/PILOT_V2_P01_DECOMPOSITION.md` is not in the branch, add it before the next gate.
2. Replace the old one-SE/largest-in-band `mu` rule with the held-out inverse-defect selector: two folds for `n<200`, four folds otherwise, minimum mean defect, numerical ties to the smaller `mu`.
3. Expand `c_mu` to `{0.02,0.05,0.10,0.25,0.50,1,2,4}`.
4. Demote `D_k` from hard gate to conservative diagnostic; add actual normalized inverse-defect `A_k` and the variance-ladder decomposition.
5. Add low-dimensional TRUE-SUPPORT delete-one-cluster jackknife and KC/CR2 + MD/CR3 diagnostics; do not promote a correction to the main method yet.
6. Add P05 (`n=200,p=500,s=5`) and P06 (`n=400,p=500,s=5`) baseline calibration cells. Keep P01--P04 as stress diagnostics.
7. Add and pass the CLARABEL--HiGHS solver parity test before switching the production LP solver.
8. Run the P01 B=50 mechanism check first. Do not issue a freeze decision from it.
9. If coherent, rerun P01--P06 at the frozen pilot replication counts and evaluate the formal gate; the hard SE/SD and broad coverage calibration gate is assessed at P05, with P06 as scaling confirmation.
10. Only after that gate passes may the freeze manifest be regenerated and the main registry launched.

No step above may use coverage, p-values or truth to choose a tuning constant or finite-sample covariance correction.
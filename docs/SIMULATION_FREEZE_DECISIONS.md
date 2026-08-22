# Simulation freeze decisions

**Project:** `jpDMEQR` SJS reconstruction  
**Status:** authoritative pre-freeze decision record  
**Applies to:** `config/simulation_main.csv`, `config/simulation_pilot.csv`, all `R/*_v2.R` simulation/benchmark code, simulation scripts, and the numerical sections of the manuscript  
**Final-scale simulation authorised:** **No.** Final-scale execution remains blocked until every gate in Section 15 passes.

## 1. Purpose and authority

This document resolves items R1--R10 in `docs/SIMULATION_DESIGN_REVIEW.md`. It fixes the target definitions, data-generating mechanisms, benchmark obligations, population-target standard, screening comparisons, numerical checks, and final-run sequence before any `B=500/1000` experiment begins.

When this document conflicts with an earlier simulation note, use the following order:

1. `docs/METHOD_SPECIFICATION.md` for mathematical definitions;
2. this document for frozen simulation decisions;
3. `docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md` and `config/benchmark_requirements.csv` for comparator obligations;
4. `docs/RESULTS_CONTRACT.md` for output schemas;
5. `docs/SIMULATION_PROTOCOL.md` for the broader experimental rationale;
6. the frozen registry files;
7. manuscript prose;
8. legacy code.

A registry promise must either be implemented and reported or removed by a committed protocol amendment before the registry checksum is frozen. It may not be silently omitted after results are seen.

## 2. Target notation

For a fixed positive-definite working nuisance penalty matrix `Lambda`, define

\[
q_i(\beta,\gamma;D_i)
=\frac{1}{m_i}\sum_{j=1}^{m_i}
\rho_\tau(Y_{ij}-X_{ij}^{\mathsf T}\beta-Z_{ij}^{\mathsf T}\gamma)
+\frac12\gamma^{\mathsf T}\Lambda\gamma,
\]

\[
\gamma_i^\dagger(\beta;D_i)
\in\arg\min_\gamma q_i(\beta,\gamma;D_i),
\qquad
Q_\Lambda(\beta)
=E\left[q_i\{\beta,\gamma_i^\dagger(\beta;D_i);D_i\}\right],
\]

and

\[
\beta^\star(\Lambda)
\in\arg\min_\beta Q_\Lambda(\beta).
\]

The registry field `target_mode` has only two admissible values:

- `structural`: evaluation uses the generated slope vector `beta0`; the DGP must satisfy the structural-target conditions in Section 3;
- `profile_mc`: evaluation uses a separately generated high-accuracy Monte Carlo approximation to `beta_star(Lambda)`.

No simulation may evaluate coverage against `beta0` merely because the DGP contains a variable named `beta0`.

## 3. R1: nuisance-ridge sensitivity and the E04 target

### Decision

Keep E03--E06 on the **structural slope target** under the frozen orthogonal baseline DGP. In particular, E04 with `lambda_gamma=4` remains `target_mode=structural`. Delete any note claiming that the slope target automatically changes merely because `Lambda` is large.

### Reason

In the E-module baseline:

- each row of `X_i` is centred;
- `X_i` is independent of `(m_i,Z_i,b_i,epsilon_i)`;
- cluster size is non-informative;
- no fixed intercept is included in the penalised high-dimensional block;
- the fitted nuisance design contains the generated random-effect direction.

At `beta=beta0`, the profiled nuisance solution depends on `(Z_i,b_i,epsilon_i,m_i,Lambda)` but not on `X_i`. The population score is therefore

\[
\nabla Q_\Lambda(\beta^0)
=-E\left[
\frac{1}{m_i}X_i^{\mathsf T}
\psi_\tau\{Z_i(b_i-\gamma_i^\dagger)+\varepsilon_i\}
\right]=0
\]

for any fixed positive-definite `Lambda`. Convexity and identification then give the structural slope target. A large `Lambda` can materially alter nuisance shrinkage, effective curvature, interval length, and finite-sample bias without changing the population slope target in this orthogonal design.

### Required prospective audit

Before final E-module execution, run a low-dimensional population audit for

```text
Lambda multiplier in {0.25, 0.5, 1, 2, 4}
```

using at least 100,000 clusters and four independent repeats. Save `beta_star_mc-beta0`, coordinatewise Monte Carlo standard errors, KKT residuals, and between-repeat differences. The structural label passes only when any observed displacement is compatible with target-MC and numerical error and does not persist as population size increases.

If the DGP is changed so that `X` is correlated with random effects, cluster size is informative, the nuisance direction is omitted, or a fixed intercept is targeted, the structural argument no longer applies and `profile_mc` is required.

## 4. R2: weak- and strong-signal scenarios

### Decision

Add four registered signal-strength scenarios around the common baseline:

```text
n = 200
p = 1000
s = 8
q = 1
tau in {0.25, 0.50}
signal in {0.40, 1.10}
error = Gaussian
m_i in {3,...,8}
B = 500
```

The default signal remains `0.75` elsewhere.

### Required reporting

For each signal level report estimation error, TPR, FDP, selected size, active-coordinate power, null type-I error, bias, empirical SD, mean SE, coverage, interval length, convergence, and Dantzig feasibility. Weak-signal losses must not be hidden by conditioning on successful selection.

## 5. R3: within-cluster dependence and cluster-size imbalance

### 5.1 Gaussian-copula dependence

Retain the dependence sensitivity experiment. For cluster `i`, generate latent Gaussian variables

\[
(U_{i1},\ldots,U_{im_i})^{\mathsf T}
\sim N(0,R_i),
\qquad
(R_i)_{jk}=0.4^{|j-k|}.
\]

Set `V_ij=Phi(U_ij)` and transform through the selected marginal CDF:

\[
\varepsilon_{ij}=F^{-1}(V_{ij})-F^{-1}(\tau).
\]

Required registered cases:

1. Gaussian marginal, `tau=0.50`;
2. rescaled `t3` marginal, `tau=0.25`.

Clusters remain independent. The errors remain marginally quantile-centred and independent of `X`, so the fixed slope target remains structural. Report cluster-score skewness/kurtosis, interval coverage, and the difference between cluster-sandwich and naive observation-level uncertainty.

### 5.2 Non-informative cluster-size imbalance

Retain a separate imbalance setting. Draw

\[
m_i\in\{2,\ldots,12\},
\qquad
P(m_i=m)\propto \pi(1-\pi)^{m-2},
\quad \pi=0.35,
\]

independently of all covariates, random effects, and outcomes. This is a truncated shifted-geometric design and remains `target_mode=structural`.

Every replication records the realised minimum, maximum, mean, standard deviation, median, quartiles, and number of subjects at each cluster size.

## 6. R4: misspecified and boundary DGPs

The following scenarios are retained and use `target_mode=profile_mc`.

### 6.1 Omitted random slope

Generate random intercept and random time slope, fit only a random intercept, and approximate the resulting profile target. Use `p=300` for the final target-MC scenario so the active coordinates plus designated null controls can be optimised accurately. This lower `p` is a planned target-approximation exception, not part of the core high-dimensional scaling claim.

### 6.2 Random scale driven by the random intercept

Generate

\[
\varepsilon_{ij}
=\exp\{0.35\,\operatorname{clip}(b_{i0}/\sigma_{b0},-3,3)\}
\{U_{ij}-F_U^{-1}(\tau)\}.
\]

The fitted location-only nuisance model is misspecified. Evaluate coverage against the profile-MC target and report `beta_star-beta0` separately.

### 6.3 Correlation between `X` and random effects

Generate a baseline Gaussian design and replace the first predictor by

\[
X_{ij,1}
=\sqrt{1-\rho_{Xb}^2}\,U_{ij,1}
+\rho_{Xb}\,b_{i0}/\sigma_{b0},
\qquad \rho_{Xb}=0.40.
\]

Re-standardise the population predictor to unit variance. This scenario is `profile_mc` because the structural-score orthogonality is broken.

### 6.4 Informative cluster size

Generate an independent cluster latent variable `u_i~N(0,1)` and define

\[
\operatorname{logit}(p_i)
=-0.4+0.8\{b_{i0}/\sigma_{b0}+u_i\}/\sqrt2,
\qquad
m_i=2+\operatorname{Binomial}(10,p_i).
\]

Make `X_{ij,1}` contain the component `0.5 u_i` before standardisation. This creates dependence among cluster size, covariates, and cluster heterogeneity. The target is the profile minimiser under the **observed-cluster distribution**, with equal cluster weighting in the population criterion.

### 6.5 Nonlinear random-effect contribution

Generate

\[
Y_{ij}
=X_{ij}^{\mathsf T}\beta^0+Z_{ij}^{\mathsf T}b_i
+0.5\left\{(b_{i0}/\sigma_{b0})^2-1\right\}t_{ij}
+\varepsilon_{ij},
\]

while fitting only the linear working nuisance design. Use `profile_mc` and report target displacement.

### 6.6 Working-penalty sensitivity

Scalar changes to a fixed positive-definite `Lambda` are tuning/curvature sensitivity under the orthogonal baseline and are not automatically labelled model misspecification. A separate non-diagonal `Lambda` check is allowed for `q=2`, but its target label follows the structural conditions in Section 3 rather than the magnitude of the penalty alone.

## 7. R5: benchmark set and comparison boundaries

### 7.1 Required practical methods

The final registry must contain, where their task and dimension are appropriate:

- `PROFILE-DQR`;
- `POOLED-QR-LASSO`;
- `SQR-DEBIASED-IID`;
- `LQMM`;
- `BIAS-ADJ-LQMM`;
- `DOUBLE-PEN-QLMM`;
- `QGEE-SCAD`;
- `QIF-SEE`.

### 7.2 Required oracle/mechanism methods

- `PROFILE-DQR-TRUE-SUPPORT`;
- `PROFILE-DQR-TRUE-NUISANCE`;
- `PROFILE-DQR-POP-H`;
- `PROFILE-DQR-SPLIT`.

### 7.3 Excluded method

`BAYES-MIXED-LASSO` is excluded from final quantile-inference claims. It may appear only in a separately labelled mean-prediction supplement and is never a coverage comparator.

### 7.4 Target and metric restrictions

- Coverage is compared only for methods with a documented interval procedure and a commensurate target.
- `QGEE-SCAD` and `QIF-SEE` are primarily estimation/selection comparators when their marginal target differs from the profile target.
- Low-dimensional methods receive either a common independently screened set in the practical table or the true support in an explicitly oracle table.
- A missing or failed adapter may not be replaced by a method with a similar name.

The exact acceptance requirements are in `docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md` and `config/benchmark_requirements.csv`.

## 8. R6: population profile-target standard

For every `profile_mc` scenario, the final target asset must use

```text
n_population >= 100000 independent clusters
independent_repeats >= 4
```

The four repeats use deterministic, non-overlapping seeds. Save:

- the estimate from every repeat;
- the average target;
- coordinatewise between-repeat MCSE;
- maximum pairwise difference;
- convergence, KKT, nuisance-gradient, and identity diagnostics;
- DGP/config snapshot;
- implementation commit;
- file and configuration checksums.

Target uncertainty is propagated through a sensitivity audit by reevaluating coverage at the estimated target plus and minus `1.96 * target_mc_se` coordinatewise. A target asset is invalid if any repeat fails, if the target design does not cover all fitted columns, or if between-repeat variation is not small relative to simulation standard errors.

Development runs may use fewer clusters/repeats only when labelled `development`; such assets cannot satisfy the final gate.

## 9. R7: Module A scaling grid

The primary rate grid is

```text
p in {500, 1000, 2000}
n in {100, 200, 400}
s = 5
tau = 0.50
q = 1
```

A secondary sparsity grid uses

```text
p = 1000
n in {100, 200, 400}
s = 10
```

The registry must contain complete `n` sequences for every rate slope that will be reported. The manuscript may state `s in {5,8,10}` only as a cross-module range; Module A itself uses the two sequences above.

The `p=2000` high-cost cells may use `B=500`; other core rate cells use `B=1000`. This exception must be visible in the registry and result table.

## 10. R8: Module F screening comparisons

Module F includes five explicitly labelled routes:

1. no screening, where computationally feasible;
2. outcome-blind, cluster-weighted variance/MAD filtering;
3. CA-IQR-SIS as a heuristic;
4. independent subject-split quantile-score screening;
5. oracle support.

Rules:

- variance/MAD filtering uses no outcome information;
- split screening uses disjoint screening and inference clusters;
- designated inference coordinates may be forced into the inferential design only when this is declared before simulation; sure-inclusion is reported separately;
- same-sample CA-IQR-SIS followed by ordinary Wald intervals is exploratory, not confirmatory;
- no-screening failures at ultra-high dimension remain in the denominator and are reported as computational boundary results.

## 11. R9: numerical and tuning supplements

### 11.1 Strict profile geometry

The strict geometry audit covers at least 20 deterministic seeds and crosses:

```text
tau in {0.25, 0.50, 0.75}
q in {1, 2}
h multiplier in {0.75, 1.00, 1.25}
balanced and unbalanced cluster sizes
diagonal and non-diagonal positive-definite Lambda
finite-difference steps around 1e-5 with sensitivity checks
```

`CVXR` and an approved solver are mandatory. Dantzig validation may not be skipped with a warning.

Default thresholds remain those in `docs/METHOD_SPECIFICATION.md`; a relaxation requires a committed numerical note and a new preflight run.

### 11.2 Module G

In addition to `n`, `p`, `q`, and selected dimension, vary:

- number of requested target coordinates;
- number of worker processes.

Report total runtime, runtime per target coordinate, peak memory, speedup, parallel efficiency, ordering determinism, and failure rates.

### 11.3 Module E

Use the full frozen nuisance-penalty grid:

```text
Lambda multiplier in {0.25, 0.5, 1, 2, 4}
```

and retain the frozen grids for `h`, `lambda_beta`, and Dantzig tolerance. No tuning value is chosen by empirical coverage.

## 12. R10: finite-sample theory audit and replication counts

### 12.1 Replication counts

- deterministic geometry grid: at least 20 seeds;
- preflight micro-registry: small deterministic debugging run;
- pilot P01--P04: `B=200` each;
- core inference/rate scenarios: `B=1000`;
- secondary robustness and misspecification: `B=500`;
- computation module: at least 20 timed repeats, normally 50;
- population targets: at least 100,000 clusters and four repeats.

### 12.2 Required diagnostics

Where oracle quantities are computationally available, record:

1. smoothed-target minus unsmoothed-target error;
2. target MCSE;
3. penalised estimation error;
4. effective-Hessian max-norm and operator-norm errors;
5. smallest eigenvalue and condition number;
6. Schur/ridge identity discrepancy;
7. requested and selected Dantzig tolerance;
8. Dantzig residual, solver status, and row norms;
9. precision-direction error;
10. one-step correction magnitude;
11. scaled Bahadur remainder;
12. cluster-score skewness, kurtosis, and QQ diagnostics;
13. empirical SD, mean/median SE, and their ratio;
14. 90% and 95% coverage, interval length, null type-I error, and active power;
15. realised cluster-size summaries;
16. runtime, memory, iteration counts, KKT residual, nuisance gradients, and failure stage;
17. paired method differences and paired MCSE;
18. log--log rate slopes with uncertainty.

Failed fits, infeasible precision rows, and warnings remain in the primary denominator. Successful-fit-only summaries are secondary.

## 13. Registry exceptions that must be documented

- Misspecified `profile_mc` cells may use `p=300` to make population optimisation reliable.
- The larger-cluster cell may use `n=150`; its expected row count is comparable to higher-`n` small-cluster cells.
- Low-dimensional LQMM-type methods use a common independently screened set or oracle support and are never described as full-`p` competitors.
- Dense-precision scenarios are assumption-violation stress tests, not expected-validity settings.

## 14. Prohibited changes after freeze

After the final manifest is created, do not:

- alter scenario rows, method lists, target labels, seeds, tuning grids, or replication counts;
- overwrite target assets or raw result directories;
- remove a comparator because it performs well, poorly, or slowly;
- recalibrate a correction or standard error using coverage;
- change a `profile_mc` target after inspecting method performance;
- switch a failed empirical or simulation outcome to a more favourable alternative without a protocol amendment.

Any required post-freeze change creates a new versioned registry and run identifier. The previous run remains immutable.

## 15. Mandatory execution sequence and final-run gate

Execute in this order:

1. apply the complete ordinary-file freeze patch;
2. confirm this document, `docs/EMPIRICAL_FREEZE_DECISIONS.md`, `docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md`, `config/benchmark_requirements.csv`, and the complete `AGENTS.md` are present;
3. implement and test all required benchmark adapters;
4. parse all R files and run package/unit tests;
5. validate the registry and manuscript/config consistency;
6. run strict profile geometry on the current commit;
7. build and validate every required population target asset;
8. run the micro-preflight registry;
9. run P01--P04 and evaluate the frozen pilot gate;
10. run the freeze-preflight script and create a commit/config checksum manifest;
11. start `B=500/1000` final execution only when the manifest matches the current commit and registry checksum.

The main runner must refuse to start when any required comparator is unimplemented, any strict geometry check is stale or failed, any `profile_mc` asset is missing or below final accuracy, the pilot gate has not passed, or the freeze manifest does not match the current commit.

## 16. Current authorisation

The simulation design is now specified, but **final-scale computation is not yet authorised**. Adapter implementation, strict numerical validation, target construction, micro-preflight, and pilot gates remain prospective requirements. No manuscript performance claim may be generated from a development run.
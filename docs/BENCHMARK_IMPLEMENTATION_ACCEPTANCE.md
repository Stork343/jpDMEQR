# Benchmark implementation and acceptance contract

**Status:** authoritative comparator implementation contract  
**Applies to:** every method listed in `config/simulation_main.csv`, `config/simulation_pilot.csv`, `config/simulation_preflight.csv`, or an empirical comparison table  
**Final-run policy:** a required method is not available merely because an adapter function with its name exists. It must pass the evidence and fidelity gates in this document.

## 1. Purpose

The simulation study compares procedures with different targets, dimensions, and inferential capabilities. This document prevents three recurrent errors:

1. silently replacing a published method with an easier approximation;
2. reporting coverage for a method that has no matching interval procedure or no commensurate target;
3. declaring an adapter complete without deterministic tests, convergence diagnostics, or a reference reproduction.

`config/benchmark_requirements.csv` is the machine-readable companion. When the CSV and this document disagree, this document governs and the CSV must be corrected before the registry is frozen.

## 2. General adapter API

Every adapter has the interface

```r
fit_benchmark_<id>_v2(
  train,
  tau,
  target_coords,
  tuning,
  seed,
  context = list(),
  control = list()
)
```

and returns a named list containing at least

```text
method_id
status                   ok / warning / failed / not_implemented
failure_stage
failure_message
beta_hat                 named vector on the supplied design
beta_tilde               optional named vector
se                        optional named vector
ci_lower, ci_upper        optional named vectors
selected                  feature identifiers
prediction_function       or a documented serialisable prediction object
runtime_sec
converged
kkt_residual              when defined
warning_messages
implementation_version
reference_identifier
adapter_fidelity_status
```

Additional method-specific objects are retained in `fit_object` or a serialisable diagnostic object.

### 2.1 Status rules

- `ok`: all method-specific stopping and output checks pass;
- `warning`: a finite result exists but a declared warning criterion fails;
- `failed`: the method was attempted and failed at a named stage;
- `not_implemented`: no faithful adapter exists.

An adapter must not catch an error, insert zeros/NA-free placeholders, and return `ok`. A method must not call another adapter and merely change `method_id` unless the equivalence is mathematically exact and documented as a tested limiting case.

## 3. Fairness rules shared by all methods

1. Use the same generated data and common random numbers within a replication.
2. All tuning uses training clusters only.
3. Subject/cluster splits, not row splits, are used for clustered methods.
4. A low-dimensional method receives either:
   - the same independently selected feature set as other practical low-dimensional methods; or
   - the true support in an explicitly labelled oracle table.
5. Feature names and coefficient order must be preserved.
6. Prediction for a new cluster may not use generated latent random effects. Conditional follow-up prediction is a separate task.
7. Oracle information is passed only through the `context` object and is never available to practical methods.
8. Failure and warning counts remain in the primary denominator.
9. No adapter tuning uses empirical coverage, bias, true support, test outcomes, or post hoc correction constants.
10. A comparator may be excluded only by a committed pre-freeze amendment, not after its performance is observed.

## 4. Evidence required before a method passes

Every required adapter has all of the following.

### 4.1 Source note

Create

```text
docs/benchmark_notes/<method_id>.md
```

containing:

- full reference and DOI/paper identifier;
- exact equations/algorithm sections implemented;
- source or reference-code location and version, when available;
- tuning rules and working-correlation/random-effect choices;
- target interpretation;
- known deviations from the source;
- metrics for which comparison is permitted;
- dependency and license information.

### 4.2 Deterministic unit test

Create

```text
tests/testthat/test-benchmark-<method_id>.R
```

covering deterministic output, coefficient-name alignment, finite predictions, declared convergence, and failure classification.

### 4.3 Limiting-case test

At least one analytically interpretable or simpler limiting case must pass. Examples include negligible random-effect variance, one observation per cluster, true support supplied, or zero penalty.

### 4.4 Reference-fidelity check

For a published method requiring reimplementation, satisfy at least one of:

1. reproduce a numerical example/table/figure from the paper within stated Monte Carlo or numerical tolerance;
2. match maintained reference code on the same deterministic input;
3. obtain two independent implementations that agree within a frozen tolerance.

A vague similarity in qualitative performance is not a fidelity check.

### 4.5 Output-schema test

The adapter must pass `docs/RESULTS_CONTRACT.md`, including method-level and coordinate-level rows, solver/runtime fields, warning/failure fields, target label, and implementation version.

### 4.6 Evidence manifest

Write a machine-readable result under

```text
results/preflight/benchmarks/<method_id>/acceptance.json
```

or an equivalent CSV/RDS manifest containing:

```text
method_id
commit_sha
reference_identifier
unit_test_pass
limiting_case_pass
fidelity_check_pass
schema_check_pass
deterministic_seed_pass
allowed_metrics
created_utc
notes
```

The final preflight reads this evidence; it does not infer completion from function names.

## 5. Method-specific requirements

## 5.1 `PROFILE-DQR`

### Role

Proposed practical method; regularised profile target; high-dimensional coordinate inference.

### Mathematical invariants

The adapter must use the exact objects in `docs/METHOD_SPECIFICATION.md`:

\[
\widehat g_h(\beta)
=-\frac1n\sum_i\frac1{m_i}X_i^{\mathsf T}\psi_{\tau,h}\{r_i(\beta)\},
\]

\[
\widehat H_{\mathrm{eff}}
=\frac1n\sum_i
\left[X_i^{\mathsf T}W_iX_i
-X_i^{\mathsf T}W_iZ_i
(Z_i^{\mathsf T}W_iZ_i+\Lambda)^{-1}
Z_i^{\mathsf T}W_iX_i\right],
\]

and the residualised identity must include `A_i' Lambda A_i`.

### Acceptance

- strict gradient/Hessian/Schur tests pass on the current commit;
- nuisance gradients and profile KKT residuals meet frozen tolerances;
- every requested Dantzig row reports the requested/selected `mu`, residual, solver, status, and row norms;
- one-step correction and cluster sandwich use the same cluster score;
- deterministic reference and optimised implementations agree within committed tolerances;
- no truth-based calibration constant exists in code.

Allowed metrics: estimation, selection, prediction, coordinate coverage, theory diagnostics, computation.

## 5.2 `POOLED-QR-LASSO`

### Role

Classical high-dimensional negative control that ignores cluster heterogeneity.

### Requirements

- fit an l1-penalised quantile regression to pooled rows;
- use training-cluster tuning even though the objective is pooled;
- return named coefficients and new-cluster fixed-part predictions;
- do not reuse the proposed profile nuisance or cluster sandwich while calling the result pooled QR-Lasso.

Coverage is not a primary metric unless a separately implemented naive iid interval is explicitly labelled. Estimation, selection, and prediction are allowed.

Limiting check: when cluster effects are absent and rows are independent, results should agree with the underlying quantile-Lasso implementation.

## 5.3 `SQR-DEBIASED-IID`

### Role

Modern debiased convolution-smoothed high-dimensional quantile regression under an iid working assumption; used to isolate the cost of ignoring clustering.

### Non-acceptance of the current shortcut

Creating a zero nuisance column, assigning each row its own cluster, calling `PROFILE-DQR`, and changing the method label is **not by itself a fidelity implementation**. The final adapter must expose the iid smoothed-QR objective, iid score/Hessian, precision-direction estimator, one-step update, and iid variance formula directly and document their correspondence to the cited method.

### Required fidelity checks

1. Formula-level mapping to Yan, Wang and Zhang (2023), JMLR 24(245), paper `22-1217`.
2. Observation-level sample size `N` is used for the iid estimator and variance; no cluster sandwich is used.
3. No nuisance ridge or profiled nuisance optimisation appears in the iid objective.
4. Bandwidth, l1 penalty, and precision tolerance follow a frozen rule compatible with the source method.
5. In the limiting design with one observation per independent cluster and no nuisance column, the iid implementation agrees with the corresponding reduced proposed implementation up to documented numerical tolerance.
6. On genuinely clustered data, the adapter intentionally ignores clustering and is labelled `naive_iid` in target/variance metadata.
7. Reproduce a source simulation cell or maintained reference-code output before final acceptance.

Allowed metrics: estimation, selection, prediction, and naive-iid coordinate coverage, always labelled as such.

## 5.4 `LQMM`

### Role

Classical direct linear quantile mixed-model comparator.

### Requirements

- use `lqmm` with the same random-intercept/intercept-slope working design as declared for the scenario;
- fit only a common independently screened set or true support because the method is not a full-`p` high-dimensional comparator;
- save quadrature/control settings, convergence messages, coefficient names, and runtime;
- prediction for new subjects uses fixed effects only unless conditional prediction is separately declared.

Allowed metrics: low-dimensional estimation, prediction, and fixed-effect intervals where the package returns stable inferential output. Practical and oracle-support uses must be separate.

Limiting check: negligible random-effect variance approaches ordinary low-dimensional quantile regression.

## 5.5 `BIAS-ADJ-LQMM`

### Role

Many-small-cluster bias-adjusted comparator based on Battagliola et al. (2022), DOI `10.1016/j.ecosta.2021.07.003`.

### Fidelity requirements

- implement the estimator, bias adjustment, and resampling/Monte Carlo scheme exactly as specified in the paper or official code;
- do not replace it with a generic cluster bootstrap around `lqmm`;
- preserve the paper's conditioning, random-effect specification, number of resamples, and interval construction;
- tuning/resampling seeds are deterministic and independent of simulation truth;
- record the unadjusted estimate, estimated bias, adjusted estimate, standard error/interval, number of successful resamples, and failure reasons;
- use a common low-dimensional screened/oracle set.

### Acceptance

- reproduce at least one published simulation/example result within its Monte Carlo uncertainty, or match official/reference code;
- demonstrate reduced small-cluster bias relative to unadjusted LQMM in a frozen diagnostic cell without using that cell to tune the adjustment;
- return failure rather than an adjusted estimate when the required resampling success fraction is not met.

Allowed metrics: low-dimensional estimation, bias, prediction where defined, and coverage only when the faithful interval procedure is implemented.

## 5.6 `DOUBLE-PEN-QLMM`

### Role

Direct penalised mixed-quantile comparator based on Li, Liu and Luo (2020), DOI `10.1007/s11424-020-9065-4`.

### Fidelity requirements

- implement the paper's stated objective, fixed-effect penalty, random-effect penalty/selection mechanism, optimisation sequence, and tuning criterion;
- distinguish random-effect shrinkage/selection in the source method from the proposed fixed positive-definite ridge profile penalty;
- do not call the proposed estimator with different penalty constants and relabel it;
- report fixed and random selected sets, objective value, convergence, tuning values, and runtime;
- use only dimensions for which the source method can be faithfully executed, with any prescreening explicitly common across methods.

### Acceptance

- zero/random-effect limiting cases behave as specified by the source objective;
- reproduce a published numerical cell or match reference code;
- coefficient and random-effect dimensions are verified on deterministic examples;
- no inferential metric is reported unless the source method supplies and the adapter faithfully implements it.

Allowed metrics: estimation, selection, prediction; coverage only if a source-supported interval is implemented.

## 5.7 `QGEE-SCAD`

### Role

Recent ultra-high-dimensional longitudinal quantile estimating-equation comparator based on Zu, Lian, Green and Yu (2023), DOI `10.1080/01621459.2022.2128806`.

### Fidelity requirements

- implement the paper's quantile estimating equations, working-correlation treatment, SCAD penalty/algorithm, screening stage if required, and tuning rule;
- preserve the marginal/estimating-equation target; do not reinterpret it as the profile target;
- all correlation and tuning estimation uses training clusters only;
- report selected set, coefficient estimate, working-correlation choice, convergence, and any source-supported standard error;
- low-dimensional or screened dimensions must follow the source method and be declared.

### Acceptance

- reproduce a published simulation cell or match reference code;
- independence working correlation reduces to the corresponding quantile estimating-equation special case;
- cluster permutation leaves results invariant up to numerical tolerance;
- inference, if implemented, is labelled with the method's own target.

Allowed metrics: estimation, selection, prediction, and source-supported inference on its own target. Same-target coverage against `beta_star` is permitted only in DGPs where target equivalence is justified in advance.

## 5.8 `QIF-SEE`

### Role

Recent SJS longitudinal quantile variable-selection comparator based on Bhattacharya, Bhuiyan and Chatla (2026), DOI `10.1111/sjos.70077`.

### Fidelity requirements

- implement the QIF basis/extended score construction, smooth-threshold estimating equations, variable-selection update, and tuning/selection rule specified by the paper;
- record basis matrices, selected set, tuning values, convergence, and objective/estimating-equation residual;
- preserve its longitudinal marginal target and do not present it as a profile-debiased estimator;
- use training clusters only for tuning and any preliminary estimator.

### Acceptance

- reproduce a published numerical cell or match author/reference code;
- independence-basis limiting case behaves as documented;
- deterministic input/seed gives identical selected set and coefficients;
- cluster-order permutation does not change results.

Allowed metrics: estimation, selection, prediction. Coordinate coverage is excluded unless a separately documented source-supported inferential procedure is implemented and its target is commensurate.

## 5.9 `PROFILE-DQR-TRUE-SUPPORT`

### Role

Oracle selection diagnostic.

### Requirements

- restrict the design to the true active set plus prespecified active/null inference coordinates;
- do not penalise the true-support coefficients unless the diagnostic is explicitly named differently;
- use the same nuisance profiling, score, Hessian, Dantzig, and sandwich rules as `PROFILE-DQR`;
- map coefficients back to the full named design without turning absent features into false selections.

Allowed metrics: estimation, prediction, coordinate coverage, and the oracle gap relative to the practical method. It is never a practical comparator.

## 5.10 `PROFILE-DQR-TRUE-NUISANCE`

### Role

Oracle diagnostic for the cost of estimating cluster-specific nuisance effects.

### Scope restriction

Use only in correctly specified structural scenarios where generated random effects correspond to the fitted nuisance design. Do not use this oracle as a `profile_mc` comparator under omitted/nonlinear nuisance misspecification.

### Requirements

- subtract the generated `Z_i b_i` from the response or equivalently use generated nuisance effects in residuals;
- do not optimise a nuisance profile or apply nuisance ridge in this oracle fit;
- use the fixed-effect smoothed score/Hessian after the true nuisance is removed;
- use cluster-level score aggregation and the declared precision/one-step procedure;
- label the target as structural-known-nuisance.

Acceptance: matches low-dimensional quantile regression on `Y-Zb` in a deterministic case, is invariant to cluster relabelling, and never reads generated nuisance values in practical adapters.

## 5.11 `PROFILE-DQR-POP-H`

### Role

Oracle diagnostic for the cost of estimating the effective precision direction.

### Requirements

- retain the practical penalised profile estimate and empirical score;
- replace only the estimated Dantzig direction by a scenario/coordinate-specific high-accuracy population direction;
- construct the population effective Hessian at the frozen target and working `Lambda` using an independent population sample;
- save the population Hessian, inverse/direction, condition number, MC error/repeat information, target coordinates, DGP/config checksum, and commit;
- never use the same Monte Carlo sample for target construction and replication evaluation.

Acceptance: `H_pop * omega_pop` is numerically equal to the canonical vector on the stored target subspace; independent population repeats agree within frozen tolerance; the adapter fails if the asset does not match the current scenario/config/commit.

Allowed metrics: coordinate bias, SE, coverage, and oracle precision gap.

## 5.12 `PROFILE-DQR-SPLIT`

### Role

Independent subject-split screening followed by profile estimation/inference.

### Requirements

- split clusters, not rows, with deterministic saved roles;
- screen only on the screening clusters;
- fit and infer only on disjoint inference clusters;
- retain a prespecified target coordinate in the inferential design only under a declared forced-inclusion rule; report whether it would have been selected without forcing;
- map the selected design and targets consistently to test data;
- never reuse inference-cluster outcomes in screening/tuning.

Acceptance: zero overlap between role sets, deterministic split, no full-data preprocessing leakage, and agreement with direct `PROFILE-DQR` on the same selected/inference subset.

## 5.13 `ORACLE-LQMM`

Optional oracle classical comparator. It uses true support, is separated from practical LQMM, and follows the LQMM acceptance rules. It is not a substitute for `BIAS-ADJ-LQMM` or `DOUBLE-PEN-QLMM`.

## 5.14 `BAYES-MIXED-LASSO`

Excluded from final quantile benchmark requirements. A Gaussian mixed mean model may be retained only as a separately labelled supplemental mean-prediction benchmark. It has no permission for quantile coverage, quantile target, or primary superiority claims.

## 6. Metric-permission matrix

| Method | Estimation | Selection | Prediction | Coordinate coverage |
|---|---:|---:|---:|---:|
| PROFILE-DQR | yes | yes | yes | yes |
| POOLED-QR-LASSO | yes | yes | yes | no by default |
| SQR-DEBIASED-IID | yes | yes | yes | yes, labelled naive iid |
| LQMM | yes | no/limited | yes | source/package supported only |
| BIAS-ADJ-LQMM | yes | no | optional | faithful procedure only |
| DOUBLE-PEN-QLMM | yes | yes | yes | source supported only |
| QGEE-SCAD | yes | yes | yes | own target/source procedure only |
| QIF-SEE | yes | yes | yes | no by default |
| PROFILE-DQR-TRUE-SUPPORT | yes | oracle | yes | yes |
| PROFILE-DQR-TRUE-NUISANCE | yes | oracle | yes | yes, structural scenarios only |
| PROFILE-DQR-POP-H | diagnostic | no | no | yes |
| PROFILE-DQR-SPLIT | yes | yes | yes | yes, independent split |
| ORACLE-LQMM | yes | oracle | yes | package supported only |
| BAYES-MIXED-LASSO | supplemental | supplemental | mean prediction only | no |

## 7. Final benchmark gate

Before final-scale simulation, the preflight must:

1. parse `config/benchmark_requirements.csv`;
2. identify every method listed in the frozen registry;
3. verify that every `required_if_listed` method has an adapter and acceptance manifest on the current commit;
4. run each deterministic acceptance test;
5. confirm the adapter's allowed metrics match the requested output;
6. confirm reference identifiers and implementation versions are present;
7. stop on `not_implemented`, missing fidelity evidence, target mismatch, or schema failure.

A development/pilot run may explicitly allow missing methods for debugging the proposed estimator, but the output must list every omission and may not be used for final comparative claims.

## 8. Remaining implementation order

After the ordinary-file freeze patch is applied, implement in this order:

1. `SQR-DEBIASED-IID` formula-level fidelity test and any required adapter correction;
2. `BIAS-ADJ-LQMM`;
3. `DOUBLE-PEN-QLMM`;
4. `QGEE-SCAD`;
5. `QIF-SEE`.

The patch-provided `PROFILE-DQR-TRUE-NUISANCE`, `PROFILE-DQR-POP-H`, and `PROFILE-DQR-SPLIT` implementations still require their own acceptance tests and manifests before the final gate passes.
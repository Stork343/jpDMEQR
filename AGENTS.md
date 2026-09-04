# AGENTS.md

## 1. Purpose

This file governs Codex, autonomous coding agents, and human contributors working in this repository. The project is a statistical-methodology manuscript and reproducibility package. Mathematical consistency, target clarity, honest failure reporting, and prospective experimental governance take precedence over speed, cosmetic refactoring, or obtaining favourable numerical results.

This is a repository-level instruction file. Read it before modifying code, configurations, manuscript sources, tests, or data pipelines.

## 2. Current project state

Target journal: **Scandinavian Journal of Statistics**.

Primary methodological identity:

```text
many independent clusters
+ bounded cluster sizes
+ high-dimensional fixed effects
+ low-dimensional cluster-specific nuisance effects
+ exact ridge-profile geometry
+ row-wise debiased coordinate inference
+ cluster-level sandwich uncertainty
```

Current execution status:

- the method specification is authoritative;
- simulation decisions R1--R10 are frozen in ordinary text files;
- final benchmark adapters are not all accepted;
- final population target assets have not all passed their gate;
- final `B=500/1000` simulation is **not authorised**;
- GSE65391 download/build/audit may proceed;
- confirmatory/full empirical fitting is **not authorised** until the empirical gates pass.

Do not describe development outputs as final results.

## 3. Mandatory reading order

Before starting any task, read the files relevant to it in this order.

### 3.1 All method/simulation work

1. `AGENTS.md`
2. `docs/METHOD_SPECIFICATION.md`
3. `docs/SIMULATION_FREEZE_DECISIONS.md`
4. `docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md`
5. `config/benchmark_requirements.csv`
6. `docs/RESULTS_CONTRACT.md`
7. `docs/SIMULATION_PROTOCOL.md`
8. `config/simulation_pilot.csv`
9. `config/simulation_main.csv`
10. current `R/*_v2.R` and `scripts/simulation/*`
11. manuscript simulation/theory source

### 3.2 All empirical work

1. `AGENTS.md`
2. `docs/METHOD_SPECIFICATION.md`
3. `docs/EMPIRICAL_FREEZE_DECISIONS.md`
4. `docs/EMPIRICAL_PROTOCOL_GSE65391.md`
5. `config/application.yml`
6. `config/module_signatures.csv` when created
7. `docs/RESULTS_CONTRACT.md`
8. current `R/*application*_v2.R` and `scripts/application/*`
9. manuscript application source

### 3.3 Authority order when files disagree

1. `docs/METHOD_SPECIFICATION.md` for mathematical objects;
2. `docs/SIMULATION_FREEZE_DECISIONS.md` or `docs/EMPIRICAL_FREEZE_DECISIONS.md` for frozen design decisions;
3. `docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md` and `config/benchmark_requirements.csv` for comparator obligations;
4. `docs/RESULTS_CONTRACT.md` for output schemas;
5. the corresponding protocol;
6. frozen configuration files;
7. current v2 reference code;
8. manuscript prose;
9. legacy package code and `legacy-pre-sjs-reconstruction`.

Never copy an inferential formula from legacy code without checking it against the method specification.

## 4. File-delivery and patch integrity rules

- Commit ordinary UTF-8 source files directly to the branch.
- Do not deliver required research documents or code as split base64 fragments, segmented archives, opaque blobs, or pasted compressed payloads.
- Do not attempt to reconstruct or apply an archive when any segment is missing or its checksum is unavailable.
- A patch is complete only when all intended files are visible in the repository tree, parseable, and covered by a commit/PR diff.
- Large binary data and generated result objects remain outside Git; source/config/document files remain ordinary files.
- If a partial archive already exists, mark it unusable or remove it after the equivalent ordinary files are committed.

## 5. Non-negotiable mathematical definitions

We observe independent clusters

\[
D_i=\{(Y_{ij},X_{ij},Z_{ij}):j=1,\ldots,m_i\},
\qquad i=1,\ldots,n,
\]

with bounded `m_i`, high-dimensional fixed-effect design `X`, fixed low-dimensional nuisance design `Z`, and a fixed positive-definite working penalty matrix `Lambda`.

All primary stochastic rates, Monte Carlo sample sizes, and standard errors use `n`, the number of independent clusters. Do not treat the total row count as the number of independent observations.

### 5.1 Complete cluster criterion

\[
q_{i,h}(\beta,\gamma)
=\frac1{m_i}\sum_{j=1}^{m_i}
\rho_{\tau,h}(Y_{ij}-X_{ij}^{\mathsf T}\beta-Z_{ij}^{\mathsf T}\gamma)
+\frac12\gamma^{\mathsf T}\Lambda\gamma.
\]

The nuisance penalty is part of the criterion being profiled. It may not be removed when deriving the score/Hessian.

### 5.2 Exact profile score

For the nuisance minimiser `gamma_hat_i(beta)`, residual vector `r_i(beta)`, and smoothed score `psi_{tau,h}`:

\[
\widehat g_h(\beta)
=-\frac1n\sum_{i=1}^n
\frac1{m_i}X_i^{\mathsf T}\psi_{\tau,h}\{r_i(\beta)\}.
\]

The score uses the original fixed-effect design `X_i`.

**Forbidden legacy score:**

```text
-n^{-1} sum_i m_i^{-1} X_tilde_i' psi_i
```

is not the score of the complete ridge-profiled criterion. Do not use it in KKT expansions, one-step updates, influence functions, or sandwich variances.

### 5.3 Exact effective Hessian

Let

\[
W_i=\frac1{m_i}\operatorname{diag}\{\phi_{\tau,h}(r_{ij})\},
\qquad
B_i=Z_i^{\mathsf T}W_iZ_i+\Lambda,
\]

\[
A_i=B_i^{-1}Z_i^{\mathsf T}W_iX_i.
\]

Then

\[
\widehat H_{\mathrm{eff}}
=\frac1n\sum_i
\left[
X_i^{\mathsf T}W_iX_i
-X_i^{\mathsf T}W_iZ_iB_i^{-1}Z_i^{\mathsf T}W_iX_i
\right].
\]

With `X_tilde_i=X_i-Z_i A_i`, the equivalent identity is

\[
\widehat H_{\mathrm{eff}}
=\frac1n\sum_i
\left[
\widetilde X_i^{\mathsf T}W_i\widetilde X_i
+A_i^{\mathsf T}\Lambda A_i
\right].
\]

The term `A_i' Lambda A_i` is mandatory.

### 5.4 Row-wise precision and one-step update

For target coordinate `k`:

\[
\widehat\omega_k
\in\arg\min_\omega\|\omega\|_1
\quad\text{subject to}\quad
\|\widehat H_{\mathrm{eff}}\omega-e_k\|_\infty\le\mu.
\]

\[
\widetilde\beta_k
=\widehat\beta_k
-\widehat\omega_k^{\mathsf T}
\widehat g_h(\widehat\beta).
\]

The smallest feasible value on the frozen deterministic `mu` grid is used. Record requested and selected `mu`, residual, solver/status, and row norms.

### 5.5 Cluster sandwich

The fitted cluster score is

\[
\widehat g_i
=-\frac1{m_i}X_i^{\mathsf T}
\psi_{\tau,h}(\widehat r_i).
\]

The coordinate variance is built from centred cluster-level projected scores and divided by `n` for the finite-sample standard error. The score in the estimator, influence representation, and sandwich must be the same object.

## 6. Target rules

The primary inferential target is the unsmoothed regularised profile target

\[
\beta^\star(\Lambda)
\in\arg\min_\beta
E\left[
\min_\gamma
\left\{
\frac1{m_i}\sum_j\rho_\tau(Y_{ij}-X_{ij}^{\mathsf T}\beta-Z_{ij}^{\mathsf T}\gamma)
+\frac12\gamma^{\mathsf T}\Lambda\gamma
\right\}
\right].
\]

A simulation uses `target_mode=structural` only when equality with the generated slope is justified prospectively. The orthogonal baseline conditions are:

- centred `X` independent of `(m,Z,b,epsilon)`;
- non-informative cluster size;
- no targeted fixed intercept in the high-dimensional block;
- generated nuisance directions included in the fitted nuisance design;
- quantile-centred errors.

Under these conditions, a fixed positive-definite `Lambda` may change curvature and finite-sample nuisance shrinkage without changing the population slope target. See `docs/SIMULATION_FREEZE_DECISIONS.md` for the E-module decision and required audit.

Use `target_mode=profile_mc` under omitted/nonlinear nuisance structure, `X`--random-effect correlation, informative cluster size, random-scale misspecification, or any other setting where structural equality is not established.

Do not silently substitute `beta0` for `beta_star`.

## 7. Prohibited statistical and computational practices

- No correction multiplier chosen to make empirical coverage nominal.
- No post hoc standard-error inflation/deflation.
- No tuning of `h`, `lambda_beta`, `Lambda/lambda_gamma`, `mu`, screen dimension, or adapter options using truth, empirical bias, empirical coverage, or test outcomes.
- No silent omission of failed/nonconverged replications or infeasible Dantzig rows.
- No encoding failed estimates as zeros.
- No changing a comparator after observing its performance.
- No row-level split for clustered data.
- No full-data outcome-guided feature filtering before validation/sample splitting.
- No same-sample heuristic screening followed by confirmatory pointwise Wald/FDR claims.
- No manual insertion or transcription of numerical results into LaTeX.
- No overwriting immutable raw run directories.
- No replacing a published comparator with a convenient approximation while retaining its name.
- No presenting a marginal, conditional-mean, or oracle target as the proposed profile target.

## 8. Simulation design governance

The exact R1--R10 decisions are in `docs/SIMULATION_FREEZE_DECISIONS.md`.

Mandatory elements include:

- complete Module A `n` sequences for `p=500,1000,2000` and the secondary `s=10` sequence;
- weak/strong signals `0.40/1.10`;
- eight quantile-centred error families and the missing quantile cells;
- Gaussian-copula within-cluster dependence;
- truncated-geometric non-informative cluster-size imbalance;
- `X`--random-effect correlation, informative size, nonlinear nuisance, omitted slope, and random-scale profile-target scenarios;
- complete `Lambda` grid `0.25,0.5,1,2,4`;
- no-screening, outcome-blind variance filtering, CA-IQR-SIS, independent split screening, and oracle support;
- target-coordinate and worker-count computation scaling;
- common random numbers and paired method comparisons.

### 8.1 Replication counts

- pilot P01--P04: `B=200`;
- core inference/rate cells: `B=1000`;
- secondary robustness/misspecification: `B=500`;
- computation: at least 20 timed repeats, normally 50;
- final population targets: at least 100,000 clusters and four independent repeats.

Every reported Monte Carlo mean/proportion includes MCSE; coverage/type-I error also include a binomial interval.

### 8.2 Population target assets

A final `profile_mc` asset must contain scenario/config/commit identity, four estimates, average target, coordinatewise target MCSE, between-repeat variation, KKT/convergence/identity diagnostics, and checksums. A development asset is never accepted as final.

### 8.3 Freeze rule

Once a final checksum manifest exists, do not change registry rows, seeds, methods, tuning grids, target labels, or replication counts. A necessary change creates a new versioned registry and run ID.

## 9. Benchmark adapter governance

Read `docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md` and `config/benchmark_requirements.csv` before modifying any comparator.

### 9.1 No-name-only implementations

A method is not implemented because an adapter with its ID exists. It must have:

- a source note under `docs/benchmark_notes/`;
- deterministic unit tests;
- a limiting-case test;
- a paper/reference-code fidelity check when reimplemented;
- output-schema validation;
- an acceptance manifest on the current commit.

### 9.2 Required remaining work after the ordinary-file patch

1. `SQR-DEBIASED-IID` formula-level fidelity implementation/test;
2. `BIAS-ADJ-LQMM`;
3. `DOUBLE-PEN-QLMM`;
4. `QGEE-SCAD`;
5. `QIF-SEE`.

Patch-provided `PROFILE-DQR-TRUE-NUISANCE`, `PROFILE-DQR-POP-H`, and `PROFILE-DQR-SPLIT` still require tests/manifests.

### 9.3 Metric restrictions

- QIF/SCAD longitudinal methods are estimation/selection comparators unless their own interval procedure is faithfully implemented.
- Coverage comparisons require a declared target and interval formula.
- LQMM-type methods use common screened/oracle low-dimensional sets.
- `BAYES-MIXED-LASSO` is excluded from final quantile claims.
- Oracle methods must be visibly labelled and excluded from practical superiority claims.

## 10. Required validation sequence

Do not start a pilot or main simulation merely because code parses.

### Stage 0: source integrity

1. confirm all required ordinary files exist;
2. ensure no required patch depends on missing archive segments;
3. parse all R files/scripts/tests;
4. run package/unit checks.

### Stage 1: registry and contract

1. validate required columns and unique IDs;
2. validate scenario promises against the freeze decision;
3. validate method IDs against `benchmark_requirements.csv`;
4. validate output schemas.

### Stage 2: strict geometry

Run the strict geometry audit on the current commit. It covers at least 20 deterministic seeds, all target quantiles, `q=1/2`, `h` multipliers `0.75/1/1.25`, balanced/unbalanced clusters, non-diagonal `Lambda`, finite-difference sensitivity, and mandatory Dantzig validation.

Default acceptance thresholds are those in `docs/METHOD_SPECIFICATION.md`. No Dantzig check may be skipped because CVXR/solver is absent.

### Stage 3: benchmark gate

Every method required by a selected registry row must pass its method-specific acceptance evidence. `not_implemented` blocks final execution.

### Stage 4: target assets

Build and validate all `profile_mc` assets at final scale.

### Stage 5: micro-preflight

Run the small registry containing every new DGP/screening/target branch. Check deterministic seeds, data dimensions, target loading, adapter dispatch, and result schema.

### Stage 6: pilot

Run P01--P04 with `B=200`. Evaluate convergence, identity error, Dantzig feasibility, bias trend, SE/SD ratio, coverage with MCSE, Bahadur trend, and failure classes. Pilot thresholds are debugging gates, not result calibration targets.

### Stage 7: final freeze

Create a manifest containing current commit, registry checksum, benchmark evidence, strict geometry evidence, target checksums, pilot decision, package/session/hardware details, and worker count.

### Stage 8: final run

The main runner starts only when the manifest matches the current commit and registry checksum. Development bypass flags are prohibited for manuscript results.

## 11. Coding conventions

- Primary language: R.
- Use explicit dimension validation and informative errors.
- Use deterministic seeds derived from `experiment_id`, replicate, method, and role; never from process scheduling.
- Preserve cluster identifiers explicitly; do not assume row order.
- Parallel output is sorted deterministically by experiment, replicate, method, and coordinate.
- Save convergence codes, KKT residuals, iterations, runtime, solver status, warnings, and failure stage.
- Retain unrounded raw values. Display values below `1e-4` as zero only in presentation layers.
- Final tables normally use three decimals; numerical identity diagnostics use enough precision to audit tolerances.
- Use `ggplot2` for final figures.
- Add new dependencies to `scripts/00_install_dependencies.R`, `DESCRIPTION` as appropriate, and a dependency note.
- Do not add network downloads to unit tests.
- Tests must be deterministic and small enough for CI.
- Performance optimisation comes after correctness; retain a reference implementation for equivalence checks.

## 12. Failure taxonomy and outputs

Every attempted method/replication has a result row, including failure/warning cases. Use explicit stages such as:

```text
data_generation
screening
tuning
nuisance_fit
penalised_fit
profile_geometry
precision_solver
variance_estimation
prediction
benchmark_not_implemented
dependency
schema
unknown
```

Raw outputs include:

- method-level replication table;
- coordinate-level inference table;
- failure/warning log;
- screening-role/selected-feature record;
- exact config snapshot;
- implementation commit;
- session/package/hardware manifest;
- target-asset identity;
- checksum manifest.

See `docs/RESULTS_CONTRACT.md`.

## 13. Empirical data rules

Read `docs/EMPIRICAL_FREEZE_DECISIONS.md` before application work.

### 13.1 Allowed now

- download official GSE65391 files;
- build expression/phenotype objects;
- reconcile subject/visit identifiers;
- audit C3 and cohort gates;
- audit batches, replicates, annotation, missingness, and cluster sizes;
- generate non-sensitive manifests and summary reports.

### 13.2 Not allowed until gates pass

- full/confirmatory expression--C3 fitting;
- choosing an alternative outcome after C3 gate failure;
- selecting modules after viewing C3 associations;
- full-data feature filtering/standardisation for cross-validation;
- same-subject screening and inference;
- empirical manuscript claims.

### 13.3 Storage

Do not commit raw GEO files, Series Matrix files, large RDS objects, or subject-level clinical data.

Use:

```text
data/raw/GSE65391/
data/interim/GSE65391/
data/derived/GSE65391/
```

These are gitignored. Commit only code, empty directory READMEs, checksums/manifests without sensitive row-level content, and aggregate audits.

### 13.4 Frozen outcome/model boundary

- primary outcome: `log(C3)` subject to the C3 gate;
- primary quantiles: `0.25,0.50,0.75`;
- sensitivity quantile: `0.10` if local-curvature gate passes;
- primary random effect: subject intercept;
- sensitivity: intercept plus centred-time slope;
- Model A: demographic/time/treatment/batch adjustment;
- Model B: Model A plus SLEDAI;
- processed Series Matrix primary; raw route sensitivity;
- confirmatory immune-module file must be nonempty and checksum-frozen before fitting.

## 14. Git and collaboration rules

- Work on a named branch; do not force-push shared branches unless explicitly authorised.
- Make coherent commits with messages describing the statistical/computational change.
- Do not mix mathematical-target changes with formatting refactors.
- Update tests and relevant governance docs in the same PR as a behavioural change.
- Never merge a PR with failed required CI or an unresolved target/fidelity issue.
- Record deviations from source papers in method notes and PR descriptions.
- Do not mark an Issue complete until acceptance evidence exists.
- Generated large results are stored outside Git or as approved release artifacts; their manifest/checksum is committed.

## 15. Manuscript rules

- The manuscript is generated from source; do not hard-code simulation or application results.
- Preserve explicit distinction among structural, profile, marginal, iid, mean, and oracle targets.
- Do not claim uniform dominance.
- Report assumption-violation scenarios and failures honestly.
- Main text contains the minimum interpretable results; detailed grids/diagnostics go to Supporting Information.
- Any change in method formula must first update `docs/METHOD_SPECIFICATION.md`, receive theory review, and then propagate to code/tests/manuscript.
- Screening is a computational supplement, not the principal theoretical contribution.
- Same-sample exploratory gene findings are not confirmatory/FDR-controlled.

## 16. Definition of completion

### 16.1 Adapter completion

An adapter is complete only when its source note, code, deterministic tests, limiting-case test, fidelity evidence, schema test, and acceptance manifest pass on the current commit.

### 16.2 Pilot completion

A pilot is complete only when raw results, failures, coordinate metrics, MCSE summaries, theory diagnostics, and a written gate decision exist. A completed pilot does not authorise post hoc tuning.

### 16.3 Final simulation completion

Final simulation is complete only when the frozen manifest matches all outputs, every registered replication is accounted for, summaries are reproducible from raw files, paired uncertainty is reported, and the manuscript tables/figures regenerate without manual edits.

### 16.4 Empirical completion

The application is complete only when the cohort/outcome gate, leakage audit, fold freeze, module freeze, method gate, repeated validation, confirmatory layer, exploratory split layer, sensitivity analyses, and provenance outputs all pass.

## 17. Immediate execution sequence

The next agent must perform the following sequence exactly:

1. verify the four governance files and this complete `AGENTS.md` are ordinary repository files;
2. ignore/remove the incomplete segmented patch archive;
3. apply the complete ordinary-file code/config patch;
4. implement the remaining adapters in the order specified in `docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md`;
5. add source notes, fidelity tests, and acceptance manifests;
6. run source/registry/schema checks;
7. run strict geometry;
8. build final-scale `profile_mc` targets;
9. run micro-preflight;
10. run P01--P04 and issue a pilot gate decision;
11. run freeze preflight;
12. start final simulation only after a matching manifest exists.

GSE65391 download/build/audit may run in parallel. Confirmatory/full empirical fitting waits for the empirical gate.

## 18. Escalation rule

Stop and request theory/project review when a task would require changing:

- the profile score or Hessian;
- the target definition or target_mode of a frozen scenario;
- the registered method set;
- the DGP equations;
- the pilot/final acceptance thresholds;
- the primary empirical outcome, quantiles, cohort, modules, or preprocessing route;
- a published comparator's algorithm beyond a documented faithful implementation.

Do not resolve these changes by improvisation after results are available.
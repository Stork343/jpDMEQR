# AGENTS.md

## Purpose

This file governs work by Codex and other coding agents in this repository. The project is a statistical-methodology manuscript, so mathematical consistency and reproducibility take precedence over speed or cosmetic refactoring.

## Authority order

When files disagree, use the following precedence:

1. `docs/METHOD_SPECIFICATION.md`
2. `docs/RESULTS_CONTRACT.md`
3. `docs/SIMULATION_PROTOCOL.md`
4. `docs/EMPIRICAL_PROTOCOL_GSE65391.md`
5. current files under `R/*_v2.R`
6. SJS manuscript source
7. legacy package code and the `legacy-pre-sjs-reconstruction` branch

Never copy an inferential formula from legacy code without checking it against the method specification.

## Non-negotiable mathematical rules

1. The sample criterion is averaged over independent clusters. Within each cluster, the loss is averaged over that cluster's observations.
2. The nuisance penalty matrix `Lambda` is positive definite and is part of the criterion being profiled.
3. The exact profile score is based on the original fixed-effect design:

   ```text
   g_hat(beta) = -n^{-1} sum_i m_i^{-1} X_i' psi_h(r_i(beta)).
   ```

4. The effective Hessian is the Schur complement of the nuisance block.
5. In residualised-design form, the Hessian must include `A_i' Lambda A_i`.
6. Do not use `X_tilde' psi` as the score of the complete ridge-profiled criterion.
7. The debiased update, influence function and sandwich variance must all use the same score definition.
8. The inferential target is the unsmoothed regularised profile target unless a simulation scenario explicitly declares another target.
9. In misspecified scenarios, do not evaluate coverage against `beta0` unless equality with the profile target has been established.

## Prohibited practices

- Do not introduce constants chosen to make empirical coverage close to nominal.
- Do not multiply the debiasing correction or standard errors by post hoc factors.
- Do not choose `h`, `lambda_beta`, `lambda_gamma` or `mu` using the true parameter, empirical bias, empirical coverage or test-set outcomes.
- Do not silently omit failed replications, nonconverged fits or infeasible Dantzig programs.
- Do not use observation-level random splits for clustered data.
- Do not perform outcome-guided filtering on the full empirical dataset before cross-validation or sample splitting.
- Do not report same-sample screening followed by ordinary pointwise Wald p-values as confirmatory FDR-controlled inference.
- Do not hard-code numerical results into the manuscript.

## Required validation sequence

Before changing the full simulation pipeline:

1. Run `scripts/simulation/00_validate_profile_geometry.R`.
2. Confirm analytic profile gradient against central finite differences.
3. Confirm analytic profile Hessian against finite differences of the profile score.
4. Confirm the Schur-complement identity, including the ridge-curvature term.
5. Confirm `max(abs(H omega_k - e_k)) <= mu + tolerance` for each Dantzig row.
6. Add or update tests under `tests/testthat/`.

No pilot or main simulation should run when these tests fail.

## Coding conventions

- Primary language: R.
- Functions must validate dimensions and fail with informative messages.
- Use deterministic seeds. Derive replication seeds from `experiment_id` and `replicate`; never depend on process scheduling.
- Use cluster identifiers explicitly; never assume rows are already sorted.
- Record convergence code, KKT residual, iterations, runtime and failure stage.
- Values with absolute magnitude below `1e-4` may be displayed as zero, but retain unrounded values in raw outputs.
- Tables use three decimal places unless a diagnostic requires more precision.
- Use `ggplot2` for final figures.
- Parallel execution must return results in deterministic experiment/replicate order.
- New dependencies must be added to `scripts/00_install_dependencies.R` and documented.

## Data rules

- Do not commit GEO raw files, Series Matrix files, large RDS objects or subject-level clinical data.
- Store downloads under `data/raw/GSE65391/` and derived objects under `data/derived/GSE65391/`; both are gitignored.
- Write a manifest containing source URL, file size, download time and MD5 checksum.
- Treat healthy technical replicates as QC material, not longitudinal SLE outcome observations.
- All folds are formed by subject and saved before model fitting.
- Batch normalization/correction must be outcome-blind and must not leak validation-subject information into training.

## Deliverable contract

A completed experiment must produce:

- one raw row per method, experiment and replication;
- a separate coordinate-level inference table;
- a failure log including successful fits with warnings;
- a session-information file;
- the exact configuration row used;
- Monte Carlo standard errors in summary tables;
- paired method-difference summaries when methods use common random numbers.

See `docs/RESULTS_CONTRACT.md` for schemas.

## Immediate priorities

1. Optimize and independently test `R/profile_v2.R` without changing its mathematical target.
2. Implement benchmark adapters marked `implementation_required` in `docs/BENCHMARK_MATRIX.md`.
3. Run the frozen pilot registry and inspect every acceptance metric.
4. Download and audit GSE65391; record the actual eligible C3 cohort.
5. Freeze preprocessing and experiment configurations before running final simulations.

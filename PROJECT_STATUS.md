# Project status

Last structural update: 2026-08-17

## Target journal and paper identity

The working target is the *Scandinavian Journal of Statistics*. The paper is positioned as a method for high-dimensional profile inference under many independent, small clusters and cluster-specific incidental nuisance effects.

## Completed

- SJS/Wiley NJDv5 main manuscript and Supporting Information working drafts.
- Separation of the unsmoothed inferential target from the computational smoothing bandwidth.
- Corrected profile score for the complete ridge-regularised criterion.
- Corrected Schur-complement Hessian and ridge-curvature identity.
- Coordinate-wise Dantzig/CLIME debiasing specification.
- Cluster-level sandwich variance specification.
- Non-empty illustrative bandwidth/sparsity regime in the working theory.
- Legacy implementation preserved on `legacy-pre-sjs-reconstruction`.
- Reference v2 R implementation and deterministic geometry tests.
- Frozen pilot and main simulation registries.
- GSE65391 download, metadata parsing and outcome-audit pipeline.
- Output schemas, failure policy and Codex instructions.

## Open scientific gates

### Gate G1: primitive theory audit

The manuscript currently states high-level profile regularity and remainder conditions. These must be derived from primitive design, density, smoothing, sparsity and nuisance-map conditions. In particular, the plug-in effective-Hessian rate and nonlinear profile remainder need a complete empirical-process proof.

### Gate G2: reference implementation

The v2 implementation must pass gradient, Hessian, ridge-identity and Dantzig feasibility tests over a grid of quantiles, bandwidths, random-effect dimensions and unbalanced cluster sizes. The legacy inferential functions cannot be used as substitutes.

### Gate G3: pilot simulation

The pilot must meet the prespecified rules in `docs/SIMULATION_PROTOCOL.md` without any coverage calibration. If it fails, diagnose the estimator, target approximation, tuning or variance formula. Do not tune directly to coverage.

### Gate G4: benchmark completeness

The final study requires classical, direct recent, marginal-longitudinal and oracle comparators. Adapters that are not yet executable are explicitly marked in `docs/BENCHMARK_MATRIX.md` and must be implemented or their omission justified before the final run.

### Gate G5: GSE65391 audit

The primary `log(C3)` outcome is provisional. The application proceeds only if the frozen audit thresholds are met. Actual sample sizes, missingness, visit distribution, batch structure and probe annotation must be recorded before the model is fitted.

### Gate G6: final numerical run

Final tables require at least 1000 replications for core inference scenarios and 500 for secondary robustness scenarios, with Monte Carlo uncertainty reported. All scripts and configuration files must be frozen before the run.

## Immediate work queue

1. Run `Rscript scripts/simulation/00_validate_profile_geometry.R`.
2. Profile and optimize the reference implementation while preserving numerical equivalence.
3. Implement/verify the pooled SQR, LQMM, double-penalised QLMM, longitudinal QGEE-SCAD and oracle adapters.
4. Run pilot configurations P01--P04 with `B=200`.
5. Build and independently repeat the `profile_mc` target approximations before misspecified scenarios are run.
6. Download the processed GSE65391 Series Matrix and generate the data audit report.
7. Decide whether C3 passes the outcome gate and freeze the empirical cohort definition.
8. Update Sections 5--6 of the manuscript only after the results files have passed validation.

## Branch policy

- `main`: current SJS reconstruction and active code.
- `legacy-pre-sjs-reconstruction`: exact state before the reconstruction push.
- Feature branches: use `codex/<task-name>` and open a pull request for nontrivial changes.

# jpDMEQR: Debiased Profile Inference for High-Dimensional Mixed-Effects Quantile Regression

This repository is the working research repository for the manuscript

> **Debiased Profile Inference for High-Dimensional Mixed-Effects Quantile Regression with Many Small Clusters**

currently being reconstructed for the *Scandinavian Journal of Statistics* (SJS).

## Current status

The repository contains a complete SJS/Wiley working manuscript, a corrected mathematical specification of the profile score and Schur-complement Hessian, a preregistered simulation protocol, a reference R implementation, and a reproducible data-ingestion/audit pipeline for GEO series **GSE65391**.

The project is **not submission-ready**. Numerical results in the working manuscript remain intentionally blank until the corrected implementation passes the geometry tests, the pilot simulations pass the prespecified acceptance rules, and the GSE65391 outcome audit is complete. No coverage-calibration constants or truth-dependent tuning are permitted.

The previous implementation and manuscript logic are preserved on branch:

```text
legacy-pre-sjs-reconstruction
```

Do not use that branch for new inferential results. It is retained only for auditability and regression testing.

## Start here

Codex and other coding agents should read these files in order:

1. [`AGENTS.md`](AGENTS.md)
2. [`PROJECT_STATUS.md`](PROJECT_STATUS.md)
3. [`docs/METHOD_SPECIFICATION.md`](docs/METHOD_SPECIFICATION.md)
4. [`docs/SIMULATION_PROTOCOL.md`](docs/SIMULATION_PROTOCOL.md)
5. [`docs/EMPIRICAL_PROTOCOL_GSE65391.md`](docs/EMPIRICAL_PROTOCOL_GSE65391.md)
6. [`docs/RESULTS_CONTRACT.md`](docs/RESULTS_CONTRACT.md)

The mathematical specification in `docs/METHOD_SPECIFICATION.md` is authoritative when the legacy package code and the new manuscript differ.

## Repository layout

```text
AGENTS.md                          Codex/agent operating rules
PROJECT_STATUS.md                  completed work, open gates, next tasks
config/                            frozen experiment and application settings
docs/                              method, simulation, benchmark, data and output contracts
manuscript/                        SJS main paper, supplement and support scripts
R/profile_v2.R                     corrected reference estimator and inference
R/simulation_v2.R                  data-generating mechanisms and target utilities
R/benchmark_adapters_v2.R          competitor/oracle adapter interface
R/metrics_v2.R                     finite-sample and Monte Carlo diagnostics
R/application_v2.R                 GEO metadata, QC and subject-split helpers
scripts/simulation/                geometry, target, pilot, main and summary runners
scripts/application/               download, build, audit, fit and summarize GSE65391
data/README.md                     data provenance and non-commit rules
results/README.md                  output schema and directory rules
```

## Minimal setup

R 4.4 or newer is recommended. Install the declared dependencies with:

```bash
Rscript scripts/00_install_dependencies.R
```

Run the non-negotiable profile-geometry checks before any Monte Carlo work:

```bash
Rscript scripts/simulation/00_validate_profile_geometry.R
```

Run the pilot registry:

```bash
Rscript scripts/simulation/01_run_pilot.R \
  --config=config/simulation_pilot.csv \
  --jobs=4
```

The pilot is a gate, not a source of final paper tables. Main experiments must not start until the acceptance rules in `docs/SIMULATION_PROTOCOL.md` are met. Before any registry row with `target_mode=profile_mc` is run, build its target files:

```bash
Rscript scripts/simulation/04_build_profile_targets.R \
  --config=config/simulation_main.csv \
  --n-population=10000 --repeats=2
```

The command above is a development run. Frozen paper targets require at least 100,000 population clusters and four independent repeats after the target routine has been optimized.

## GSE65391 application

The official GEO accession page is:

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65391
```

Processed Series Matrix data and optional raw/non-normalized files can be downloaded by:

```bash
Rscript scripts/application/00_download_GSE65391.R --processed=true --raw=false
Rscript scripts/application/01_build_GSE65391.R
Rscript scripts/application/02_audit_GSE65391.R
Rscript scripts/application/03_fit_GSE65391.R --mode=smoke --fold=1 --repetition=1
Rscript scripts/application/04_summarise_GSE65391.R
```

The primary outcome is provisionally `log(C3)` at the same visit. It is used only if the preregistered missingness, sample-size and variation gates pass. The script must stop rather than silently substitute another outcome.

## Manuscript

The SJS manuscript source and Supporting Information are under `manuscript/`. The source follows the Wiley NJDv5 layout. Font files are intentionally not stored in this repository. Retrieve the non-font support files with `manuscript/fetch_wiley_template_files.sh`, obtain the current `Fonts/` directory from Wiley locally, and then run `manuscript/compile.sh`.

## Reproducibility principles

- Independent sampling unit: subject/cluster, not observation.
- Primary asymptotic scale: number of independent clusters `n`.
- Primary target: the unsmoothed regularised profile parameter `beta_star`.
- The profile score uses the original fixed-effect design `X`, not the residualised design.
- The effective Hessian includes the ridge-curvature term.
- All cross-validation and sample splitting are performed at the subject level.
- No empirical coverage tuning, post hoc standard-error inflation, or truth-dependent method selection.
- Every failure is retained and classified; no failed replication is silently dropped.
- Raw GEO data and large derived matrices are not committed.

## Historical code

The original package functions (`JP_DME_QR`, `Debias_Inference`, `CA_IQR_SIS`) remain in the repository so earlier work can be reproduced. They are **legacy interfaces** until they are reconciled with `R/profile_v2.R` and the new test suite. New experiments must call the `*_v2` functions or source the v2 files explicitly.

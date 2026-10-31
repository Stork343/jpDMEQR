# SoyNAM application

The study estimates canopy--yield profile associations after year-specific genomic adjustment and evaluates the additional prediction information in canopy coverage. Its design is recorded in `PROTOCOL_zh.md` and `config/application_soynam_20261031.yml`.

Run all commands from the repository root. Set `TMPDIR` to the absolute path of `results/scratch/runtime`; set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `VECLIB_MAXIMUM_THREADS` to one. The source revisions and file objects are listed in `results/preflight/soynam_20261031/source_manifest.csv`. `fetch_sources.R` retrieves those exact objects from the official CRAN mirrors. `SOYNAM_GIT_TRANSPORT=ssh` selects GitHub's SSH endpoint on port 443 when that transport is required.

The reproduction order is:

1. `Rscript research/soynam_application_20261031/fetch_sources.R`
2. `Rscript research/soynam_application_20261031/audit_source.R`
3. `Rscript research/soynam_application_20261031/prepare.R`
4. `Rscript research/soynam_application_20261031/check_pair_qp.R`
5. `Rscript research/soynam_application_20261031/check_pooled_compiled.R`
6. `Rscript research/soynam_application_20261031/run_inference_queue.R`
7. `Rscript research/soynam_application_20261031/run_prediction_queue.R`
8. `Rscript research/soynam_application_20261031/summarise_inference.R`
9. `Rscript research/soynam_application_20261031/summarise_prediction.R`
10. `Rscript research/soynam_application_20261031/descriptive.R`
11. `Rscript research/soynam_application_20261031/export_manuscript.R`
12. `Rscript research/soynam_application_20261031/plot_results.R`
13. `Rscript research/soynam_application_20261031/check_delivery.R`

The preparation and individual run directories are immutable. Completed runs are reused by the queue. Any existing incomplete case requires inspection before resuming. Prediction queue arguments may identify an already running case and its process identifier; the queue verifies that process before monitoring it.

The inference output contains all 21 settings and quantiles, with eight coefficient systems per fit. The prediction output contains three methods at each of three quantiles. Complete lines, sample roles and internal validation predictions are retained. All reported intervals and losses are computed from these saved objects. The manuscript export requires the complete set of results.

The independent numerical checks use actual field observations. `check_pair_qp.R` compares the structured solver with Clarabel. `check_basis_precision.R` verifies original-coordinate gradients on saved selected models. `check_pooled_compiled.R` checks the pooled LP against quantreg and the R/HiGHS path in `pooled_lp.R`, including deterministic repetition and coefficient mapping. Details are recorded in `docs/benchmark_notes/SOYNAM-POOLED-LASSO.md`.

Each successful case is copied to the configured independent backup directory. The Git repository contains the code and aggregate outputs. Full genotype matrices, fitted objects and line-level predictions are stored under the ignored data and result directories.

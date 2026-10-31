# Application dependencies

The SoyNAM application uses R 4.3.2 with `data.table`, `yaml`, `jsonlite`, `Rcpp`, `RcppEigen`, `Matrix`, `glmnet`, `mgcv`, `highs`, `filelock`, `processx` and `ps`. The independent convex checks also use `clarabel` and `quantreg`. Dependency declarations are included in `DESCRIPTION` and `scripts/00_install_dependencies.R`.

The fitted spatial calibration uses `mgcv` 1.9-1. Pooled prediction uses HiGHS 1.15.1 through its official C++ API. The native dependency is installed by the repository's dependency script. Every inference request retains its R session and implementation commit; prediction requests retain the sample roles and implementation commit.

Run scripts from the repository root with `TMPDIR` set to `results/scratch/runtime` and numerical thread variables set to one. The raw data and full result objects remain in ignored directories and are copied to the independent project backup. Ordinary source files and aggregate tables are versioned in Git.

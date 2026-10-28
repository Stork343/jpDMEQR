# Repository and package checks

`scripts/research/check_repository_package_20261022.R` runs the complete repository test inventory and an `rcmdcheck` package build. Each invocation uses a new run identifier under `results/scratch/repository_package_checks_20261022`. The `unit`, `check` and `all` phases are available.

```sh
TMPDIR=results/scratch Rscript scripts/research/check_repository_package_20261022.R . new_run all
```

A final invocation supplies the exact current commit as its fourth argument. The entry rejects a changed commit or uncommitted source files when this argument is present. The request records the actual commit and source state for development runs.

## Unit processes

Every `test-*.R` file runs in its own R process with the installed local `jpDMEQR` 0.4.0 package. Both registered native entry points must be loaded. The worker uses the package library under `results/scratch/profile_pair_package_20261022/library` and retains every test result, source error, assertion table, session record and console log.

The initial process configuration is `compiled`, `batch`, `Eigen_single_thread`, `aggregate`, `newton_armijo`, `cluster_rms`, `dual_warm`, one CV worker, Hessian refresh 20 and shared first-stage fitting. Tests of a reference implementation declare their own local configuration. Dedicated parallel tests explicitly release the inherited CV-worker environment variable before setting their tested worker counts.

`unit_file_summary.csv` includes every scheduled file. Counts come from the actual `testthat` results, including warnings and skips. A process that produces no result is recorded with unavailable assertion counts and an unsuccessful process status. Numerical failures retain their original conditions. The empty-fold OSQP case also saves its actual data, matrices, solver outputs and strict residual decisions.

`assertion_pass` requires no failed assertions, errors or skipped tests. `strict_pass` additionally requires no warnings. `unit_warning_details.csv` retains the source warning text, test and line. The unit certificate reads `assertion_pass` together with these warning details; the package check requires zero ERROR and zero WARNING.

## Scratch paths

The entry maintains `/Users/houjian/.codex/jpdmeqr-scratch-20261022` as a symbolic link to this project's ignored `results/scratch` directory. Both paths are checked and recorded. Child processes use this absolute path without spaces so that changes to the current working directory preserve the physical scratch location. UTF-8 locale settings accompany the complete inherited process environment.

## Package source and check

`prepare_package_source_20261022.R` stages the package using R's `tools::inRbuildignore` implementation and the standard R build exclusions. It preserves every repository test file. The staging directory and the build destination are separate children of the ignored scratch directory. The staged file inventory is retained.

The package worker runs installation, examples, tests, vignette code, vignette rebuilding and the PDF manual through `rcmdcheck`. Its output includes the complete check object, stdout, stderr and a manifest with every ERROR, WARNING and NOTE. Optional suggested packages that are unavailable remain recorded as NOTE messages. The validated C++14 compiler declaration is retained.

The installed paired API test runs actual clustered data through the exported fitter, paired components and inference functions. Its correction and centred cluster sandwich are checked independently from the returned cluster scores.

## Source warning records

`check_native_warning_records_20261022.R` runs the actual BIAS-ADJ-LQMM DGP used by its source checks. The training data use seed 202, with the runner's deterministic subject roles and fitting seed. A direct source call uses the same fitted subjects and controls. The check compares coefficients and intervals exactly, and verifies each warning's message, call, class, sequence and recorded stage against the conditions emitted during fitting. It saves the replication and coordinate CSV files with the complete comparison object.

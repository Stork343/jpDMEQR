# Joint profile manuscript reproduction

Run scripts from the `jpDMEQR` repository root with `TMPDIR=results/scratch`. This directory supplies the fixed-design application interface, prediction learner, paired simulation ablation and manuscript builders. The general high-dimensional estimator is loaded by `research/score_validity_20261026/load_joint.R`; its local package library is prepared by the adjacent `setup.R`.

## Numerical definitions

`fixed_inference.R` fits the selected complete profile at both bandwidths and returns their paired coefficients and matched cluster covariance. `prediction.R` selects an unsmoothed joint search model using five inner subject folds. Random-intercept search uses the primal-dual solver. `conic_search.R` expresses random-slope search as a sparse quadratic programme with the same check-loss objective and prediction ridge. Its dual certificate is recomputed after solving.

Clarabel is already a declared project dependency. The additional scripts use `dplyr`, `ggplot2`, `statmod`, `processx` and the existing preprocessing dependencies. R scripts fail at a failed numerical certificate and preserve a failure record. Dispatchers execute independent processes and retain each exit code.

## Independent simulation

`run_ablation.R WORKER WORKERS` executes the registered 24-setting design on independent seed base 1800000000 with the design adjustment removed. Four disjoint worker indices were used. Every unit is stored under `results/scratch/joint_profile_manuscript_20261026/ablation/`. Full-method records are under `results/scratch/score_validity_20261026/design_independent/`.

The ablation retains the outcome and profile selections and computes the exact Gaussian selected-target displacement from each saved support. The manuscript builder requires all 1440 matched dataset keys before aggregation. Coordinatewise binomial intervals and replication-level Monte Carlo errors are both retained.

## Application execution

`prepare_modules.R` prepares the reference module designs and all subject deletion designs. `launch_modules.R` dispatches their 40 model, quantile and sensitivity settings. `launch_jobs.R processing_modules` dispatches the twelve alternative-processing settings. Their immutable outputs are under `results/raw/application/joint_profile_20261026/` and `results/raw/application/joint_profile_processing_20261026/`.

`run_prediction_fold.R PROCESSING REPETITION FOLD` calculates an outer role. `launch_jobs.R prediction` schedules the frozen roles, except the separately executed reference role 1 1. Passing `slope_only` as the fourth argument calculates the six random-slope cases for an outer role in `reference_conic_slope/`; `launch_jobs.R slope_prediction` dispatches all 25 such roles. The final manifest uses all 750 random-intercept fits and the 150 specified conic random-slope fits. Earlier numerical attempts remain in their original directories and are not duplicated in the final prediction table.

`launch_jobs.R exploratory` runs all independent gene roles. Each role selects exactly three additions to its clinical basis on one subject group and estimates their selected-model profile on the other. There are 360 completed model, quantile, processing and role combinations, stored under `results/raw/application_exploratory/joint_profile_20261026/`.

Every application uses the passed cohort audit in `data/derived/GSE65391/audit_20261022_bootstrap_v1/`. Probe choices, gene filtering, module construction and scaling use their specified training subjects. The two processing sensitivities reuse their saved outcome-blind transformations with identical subject roles.

## Verification

```bash
TMPDIR=results/scratch Rscript research/joint_profile_manuscript_20261026/validate_fixed.R
TMPDIR=results/scratch Rscript research/joint_profile_manuscript_20261026/validate_conic.R
TMPDIR=results/scratch Rscript research/score_validity_20261026/check_reproduction.R
```

The fixed-design checks compare 18 coordinate estimates and standard errors with independently requested selected-coordinate fits. The conic checks compare 18 genuine quadratic and primal-dual solutions through their objectives and strong-convexity gap bounds. The reproduction check repeats two saved full high-dimensional fits on their original data seeds.

`population_roots.R` computes deterministic population roots and checks the integral by independent Gauss–Legendre quadrature. `build_inference_assets.R` generates the simulation, module and gene tables. Simultaneous module intervals are derived from saved coefficients and standard errors with family size 30, or 10 for the lower-decile sensitivity. `build_prediction_assets.R` requires all 900 final fits and complete five-repetition subject vectors. It preserves comparator convergence states and constructs conditional subject-resampling ranges.

The editable manuscript and generated fragments are under `manuscript/`. `check_manuscript.py` performs source and PDF object checks without visual inspection. The numerical source tables remain under `results/tables/joint_profile_manuscript_20261026/`.

## Remaining study calculations

The separate 78-setting full simulation and the complete published-method prediction comparison are pending. Their frozen registries and retained records remain available for those calculations. The current article reports the completed independent evaluation, its paired ablation and the completed new-method application results; comparator tables state their actual availability.

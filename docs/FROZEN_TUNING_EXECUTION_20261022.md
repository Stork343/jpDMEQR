# Frozen penalty tuning

`scripts/simulation/_run_registry_helpers_part06.R` supplies the penalty selection used by the simulation runner. The registry's method and scenario definitions determine every candidate before outcome evaluation. The outer prediction dataset supplies the reported prediction metrics. Penalty cross-validation uses the training subjects recorded in the replication object.

Cross-validation operates on the supplied feature design after the registered screening step. F04 retains its declared full-training-subject CA-IQR-SIS heuristic, uses independently generated test clusters for performance evaluation and labels its Wald output exploratory. Independent screening uses its existing separate subject roles.

## Paired profile method

Ordinary paired fits retain the configured multiplier `a=1`. E07 evaluates the eight distinct values in `lambda_beta_multipliers`. Each candidate receives five deterministic subject folds. A fold calls `fit_profile_pair_v2` through the existing data wrapper, including both loading passes, the final L1 fit, the prescribed refit and the two inferential profiles. It predicts with the fitted fixed-effect coefficient vector. The heldout loss averages pinball loss within each subject and then gives equal weight to every validation subject. Precision estimation follows penalty selection on the full training dataset.

The true-support adapter fits the union of its prescribed support and target coordinates with zero L1 penalty. Its E07 diagnostic record retains every configured multiplier and fold. Candidate ties select the first value on the increasing grid.

## Comparator definitions

POOLED uses every distinct registry multiplier with the existing anchor `sqrt(log(max(p,2))/n_training_clusters)`. Its objective remains the pooled check-loss sum plus the raw coordinate L1 penalty. SQR uses its source Gaussian convolution initial estimator within each fold. Both its bandwidth and its penalty are recalculated from that fold's observation count `N_training`. Its precision calculation occurs on the selected full training fit.

DOUBLE-PEN-QLMM receives `config$lambda_beta_multipliers` through the runner control. Its random-effect multiplier grid and source SIC or GACV block selection remain defined by the adapter. The replication diagnostic retains the actual source grids and each iteration's candidate criteria when that computation succeeds.

## Diagnostics and failures

A candidate is eligible only when all five folds succeed with finite validation loss. Its criterion is the total heldout subject loss divided by the number of training subjects. Ineligible candidates retain their fold statuses and numerical diagnostics. An empty eligible set terminates the method at stage `tuning`.

`method_diagnostics[[method_id]]$tuning` records the configured and distinct grids, selected multiplier, fold seed, subject roles, fold losses, candidate eligibility and training coefficient summaries. The final training record includes the method seed, all training subject identifiers, selected multiplier, status and fitted numerical summaries. The total method runtime includes penalty selection.

The replication table records the selected source penalty for SQR and DOUBLE-PEN-QLMM and the source bandwidth for SQR. The final tuning record also retains the original raw anchor.

## Reproduction

Run the actual DGP checks from the repository directory with a new output identifier:

```sh
TMPDIR=results/scratch Rscript scripts/research/audit_frozen_tuning_20261022.R . results/scratch/profile_pair_tuning_20261022/new_run
```

The output contains the assertion table, complete tuning records, session information and a machine-readable manifest. The checks cover all E07 multipliers, both comparator grids, repeatable selection, training-only fits under changed heldout outcomes, the source bandwidth and penalty formulas, and genuine iteration-budget failures.

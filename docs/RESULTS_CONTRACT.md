# Results and logging contract

## 1. Run directories

Every run is immutable and uses a timestamped/versioned root:

```text
results/raw/simulation/<run_id>/
results/raw/application/<run_id>/
```

A run directory contains:

```text
config_snapshot.csv           frozen registry rows used by the run
config_sha256.txt
implementation_commit.txt
implementation_branch.txt
sessionInfo.txt
hardware.json
solver_configuration.json
run_request.json              requested experiment ids, seed, replication cap, final flag
replication_metrics.csv       one row per experiment x replicate x method
coordinate_metrics.csv        one row per experiment x replicate x method x target coordinate
theory_diagnostics.csv
screening_records.csv         screening role, selected features, split overlap counts
runner_failures.csv           present when at least one task failed to run
replication_objects/          one RDS object per experiment x replicate
RUNNING.lock                  pid, start time, commit, resumed object count
INCOMPLETE.txt                present only while a shard is unfinished
completed_manifest.json       written once every requested task finished
```

The runner writes one object per finished `experiment_id x replicate` into `replication_objects/` before continuing, so a shard interrupted by a reboot keeps the replications it already completed. Resuming with `--resume` accepts an existing run directory when its recorded registry checksum, implementation commit, seed, replication cap, final flag, benchmark policy and experiment set match the current request, its lock belongs to a process that is no longer running, and the stored objects carry the current record-layer columns; `completed_manifest.json` records how many objects were reused. A completed run directory is never resumed and never overwritten. A corrected run receives a new `run_id` and a note explaining the change.

A directory that merges shards is derived and reproducible: it carries `merge_manifest.json` with the shard identifiers and row counts, rewrites `run_request.json` with the union of the shard experiment sets, and rejects duplicated rows across shards.

采用 `JPDMEQR_RESULT_STORAGE=streaming` 时，`raw_objects.rds` 保存
`replication_file_index_v1` 索引，每次重复的对象位于
`replication_objects/`。使用 `read_registry_replication_v2()` 读取真实对象。
移动整个运行目录后，传入当前 `raw_objects.rds` 路径或明确的
`base_directory`。索引保留失败对象及其对应的重复编号。

数值库发出的 warning 保存在每种方法的 `native_warning_records`，
内容包括出现次序、完整消息和函数调用。控制台同时保留这些 warning。
拟合返回可用结果时，记录 `status=warning`；拟合失败时，保存实际
condition 与已有计算记录。

## 2. Simulation replication schema

`method_version` and `score_definition` identify the equations used by each
method. The paired method records the stochastic and bandwidth contributions
to its initial penalty as `lambda_stochastic` and `lambda_mean` before applying
the common coordinate loadings and the registered sensitivity multiplier.
`lambda_mean_multiplier` records the prospectively chosen coefficient of the
bandwidth contribution. Each replication RDS also contains `method_diagnostics`,
including all precision cross-validation candidates, excluded candidates,
selected directions, and any partial fit retained after a fitting failure.

`peak_memory_mb` is the largest sampled process RSS during a method evaluation.
`process_peak_wset_mb` is available on Windows. `process_peak_rss_mb` records
the process-lifetime maximum RSS reported by `ps_memory_full_info()` when the
platform supplies `maxrss`. Unsupported platform metrics are recorded as `NA`.

One row per `experiment_id x replicate x method`, with these columns in this order:

```text
run_id
experiment_id
module
replicate
seed
method_id
method_version
score_definition
status
failure_stage
failure_message
warning_messages
reference_identifier
adapter_fidelity_status
n_clusters
n_clusters_original
total_rows
p
p_original
s
tau
q
m_rule
m_min
m_mean
m_median
m_max
error_dist
error_dependence
design_type
random_effect_dist
x_b_corr
informative_size
nonlinear_re_strength
screen_method
screen_dimension
screen_sure_inclusion
forced_target_count
h
lambda_beta
lambda_0_n
lambda_stochastic
lambda_mean
lambda_mean_multiplier
lambda_gamma
selected_size
selected_names
tpr
fdp
fpr
exact_support
l1_error
l2_error
max_active_error
pinball_loss
pinball_loss_intercept
pinball_skill_intercept
pinball_skill_pooled
cluster_loss_q10
cluster_loss_q50
cluster_loss_q90
quantile_calibration
objective
runtime_sec
first_stage_sec
refit_sec
inference_sec
component_timing_available
peak_memory_mb
process_peak_wset_mb
process_peak_rss_mb
outer_iterations
nuisance_iterations
dantzig_feasible_fraction
warning_count
converged
kkt_residual
profile_identity_error
max_nuisance_gradient
lambda_rule
lambda_alpha
lambda_normal_quantile
lambda_safety_constant
lambda_base
lambda_loading_pass0_min
lambda_loading_pass0_median
lambda_loading_pass0_max
lambda_loading_pass1_min
lambda_loading_pass1_median
lambda_loading_pass1_max
lambda_coordinate_min
lambda_coordinate_median
lambda_coordinate_max
zero_profile_score_max
zero_kkt_ratio_max
preliminary_kkt_normalized
final_kkt_normalized
final_kkt_absolute
first_stage_iterations
first_stage_beta_change
first_stage_nonzero_count
selected_support_size
refit_set_size
refit_contains_targets
post_refit_status
post_refit_iterations
post_refit_gradient_max
post_refit_nuisance_gradient_max
post_refit_beta_change
post_refit_hessian_min_eigenvalue
post_refit_hessian_condition_number
l1_error_lasso
l1_error_refit
l2_error_lasso
l2_error_refit
omitted_active_count
target_type
target_mc_se
implementation_commit
requested_workers
h_est
h_inf
rsc_order_lower
h_est_over_rsc_order
lambda_min_Hest_SS_population
lambda_min_Hest_SS_target_sample
lambda_min_Hest_SS_fitted_sample
condition_Hest_SS_population
condition_Hest_SS_target_sample
condition_Hest_SS_fitted_sample
cone_curvature_proxy_min
```

`module`, `error_dist`, `design_type` and the remaining scenario descriptors repeat the registry row so that a raw table can be read without the configuration file. `m_min`, `m_mean`, `m_median` and `m_max` summarise the realised cluster sizes of the analysis sample. `reference_identifier` and `adapter_fidelity_status` are written by the benchmark adapters and carry their acceptance note and fidelity status. The `lambda_*` and `post_refit_*` fields are the round-3 to round-5 acceptance diagnostics of the proposed route and stay missing for adapters that do not compute them. `l1_error_lasso`, `l1_error_refit`, `l2_error_lasso`, `l2_error_refit` and `omitted_active_count` separate selection error from refit error. `requested_workers` is the worker count the registry row asked for, which enters the computation summaries. `h_est`, `h_inf`, `rsc_order_lower` and the four `lambda_min_Hest_SS_*`/`condition_Hest_SS_*` fields record the two bandwidths, the rate order and the restricted-eigenvalue diagnostics.

`pinball_loss_intercept` is the held-out check loss of the training-sample `tau`-quantile of the response. `pinball_skill_pooled` and `pinball_skill_intercept` equal `1 - pinball_loss / baseline_loss` with the baseline paired inside the same experiment and replicate. `cluster_loss_q10`, `cluster_loss_q50` and `cluster_loss_q90` are quantiles of the per-cluster held-out check loss. `first_stage_sec`, `refit_sec` and `inference_sec` partition `runtime_sec`; `component_timing_available` states whether all three were recorded. `peak_memory_mb` is the largest resident set size sampled at the start of the fit, at the end of the fit and at the end of the replication; `process_peak_wset_mb` is the process working-set high-water mark at the end of the replication. `nuisance_iterations` totals the Newton iteration count over the cluster nuisance solves in the first stage and the post-refit, and `outer_iterations` totals the outer iterations of both stages.

`n_clusters` and `p` describe the analysis sample that the fitted method actually used, while `n_clusters_original` and `p_original` carry the registry values of the scenario. The two differ in screening cells, where an independent split keeps half of the clusters for fitting and reduces the fixed-effect dimension to `screen_dimension`; manuscript tables report the registry values and treat the analysis sample sizes as diagnostics.

`dantzig_feasible_fraction` is the share of target coordinates whose row-wise precision program was feasible. It is recorded only for methods that solve that program; for a comparator whose interval route does not, the field is missing rather than zero, and the coordinate table carries `feasible = FALSE` with a missing `dantzig_mu` as a placeholder of that route.

Missing values require a reason in `failure_stage` or an explicit `not_applicable` field. Do not encode failed metrics as zero.

## 3. Coordinate-level schema

One row per requested coordinate, with these columns in this order:

```text
run_id
experiment_id
module
replicate
seed
method_id
method_version
score_definition
coordinate
coordinate_type
target_type
target_value
target_mc_se
beta_hat
beta_tilde
estimated_se
estimated_se_smoothed
ci_level
ci_lower
ci_upper
covered
ci_level_secondary
ci_lower_secondary
ci_upper_secondary
covered_secondary
interval_length_secondary
covered_target_minus_1.96mcse
covered_target_plus_1.96mcse
wald_z
p_value
interval_length
feasible
dantzig_mu
mu_source
dantzig_residual
omega_l1
omega_l2
adjacent_stability
mu_cv_min_defect
mu_cv_defect
mu_cv_quad
mu_cv_folds_used
E1_pop
E2_pop
cosine_pop
D_k
A_k
A_k_lasso_start
bahadur_scaled
bahadur_refit_scaled
ladder_L_k
ladder_se1
ladder_se2
ladder_se3
ladder_se4
fs_se_cr0
fs_se_cr3
fs_se_kc
fs_leverage_max
status
implementation_commit
```

`coordinate_type` is `active`, `null` or `always_included` and uses the generated support of the scenario. `estimated_se` is the reported standard error and `estimated_se_smoothed` the same quantity at the smoothed working criterion, kept for the sensitivity diagnostic. `covered_target_minus_1.96mcse` and `covered_target_plus_1.96mcse` repeat the coverage comparison at the lower and upper edge of the target Monte Carlo interval. `feasible`, `dantzig_mu`, `dantzig_residual`, `omega_l1` and `omega_l2` describe the row-wise precision program; `mu_source` records how the tolerance was selected, and the `mu_cv_*` fields record the cross-validation diagnostics behind that choice. `adjacent_stability` records whether the neighbour on the `mu` grid gave the same interval decision. `E1_pop`, `E2_pop`, `cosine_pop`, `D_k`, `A_k` and `A_k_lasso_start` compare the fitted direction with the population direction. `bahadur_scaled`, `bahadur_refit_scaled` and the `ladder_se*`/`ladder_L_k` fields carry the standardisation ladder, and `fs_se_cr0`, `fs_se_cr3`, `fs_se_kc` and `fs_leverage_max` carry the first-stage sensitivity variants.

The interval columns without a suffix are at level `ci_level` = 0.95; `ci_level_secondary` is 0.90 and `covered_secondary` records the same comparison at that level. `covered`, `covered_secondary` and `p_value` remain missing when the target is unavailable, method/target equality is not justified, or the method returns no interval. A comparator that produces point estimates and selections without an interval procedure therefore carries missing standard error, interval bounds, coverage indicators and p-value, and is excluded from coverage denominators rather than counted as non-covering.

已经登记区间计算的 `PROFILE-DQR` 系列、`SQR-DEBIASED-IID` 和
`BIAS-ADJ-LQMM` 即使拟合失败，也为每个指定坐标保存一行缺失区间。
无条件覆盖率以全部登记重复为分母，失败区间计为未覆盖；成功拟合条件下的
覆盖率另行保存。当前 `LQMM` 和 `ORACLE-LQMM` adapter 返回点估计，
其比较指标按照相应接口记录。

## 4. Theory-diagnostic schema

One row per `experiment_id x replicate x method`, with these columns in this order:

```text
run_id
experiment_id
replicate
seed
method_id
method_version
score_definition
target_type
target_mc_se
hessian_max_error
hessian_operator_error
hessian_min_eigenvalue
hessian_condition_number
beta_h_target_distance
profile_identity_error
max_nuisance_gradient
kkt_residual
max_dantzig_residual
max_omega_l1
precision_row_error
one_step_correction_max
max_E1_pop
max_E2_pop
max_D_k
median_D_k
q90_D_k
max_A_k
median_A_k
max_bahadur_scaled
scaled_bahadur_remainder
score_skewness
score_excess_kurtosis
n_clusters
p
```

`hessian_max_error` and `hessian_operator_error` compare the fitted effective Hessian with its population version, `precision_row_error` and `one_step_correction_max` compare the fitted precision row and the debiasing step with their population versions, and `beta_h_target_distance` is the distance between the smoothed working target and the unsmoothed profile target. `max_dantzig_residual` and `max_omega_l1` summarise the row-wise precision program, and the `max_*`/`median_*` quantities over `E1_pop`, `E2_pop`, `D_k`, `A_k` and `bahadur_scaled` are the coordinate-wise extremes of those diagnostics. Fields stay missing where the oracle quantity is unavailable, which happens for adapters that do not target the profile estimand.

## 5. Screening-record schema

One row per experiment and replicate, describing how the analysis sample was formed:

```text
run_id
experiment_id
replicate
seed
screen_method
screen_engine
screen_cluster_count
fit_cluster_count
overlap_count
requested_dimension
realised_dimension
sure_inclusion
forced_target_count
selected_names
forced_target_names
```

The three cluster counts describe the design of the screening step. When `screen_method` is `none`, no screen set exists: `screen_cluster_count` is zero, `fit_cluster_count` is the analysis cluster count, and `overlap_count` is zero. Under independent split screening, `screen_cluster_count` and `fit_cluster_count` are the disjoint screen and fit groups, and `overlap_count` must be zero; the runner refuses a split whose two groups intersect. Under the same-sample screens (`variance_filter` and `ca_iqr_sis`), both counts equal the analysis cluster count and `overlap_count` equals it as well, which records that the screen and the fit used the same clusters; these cells are reported as a screening supplement and not as confirmatory pointwise claims. `sure_inclusion` records whether every generated active feature entered the selected set, which is a diagnostic rather than a screening guarantee.

## 6. Failure classes

Use one of:

```text
data_generation
screening
penalised_fit
nuisance_fit
kkt_not_met
hessian_nonfinite
precision_infeasible
precision_solver
variance_nonpositive
prediction
benchmark_not_implemented
resource_limit
unknown
```

A fit with a warning but usable output has `status=warning` and remains in summaries. A failed fit remains in the denominator for the failure-rate table. Coverage summaries use a clearly reported denominator and include both unconditional and successful-fit-only versions; the unconditional version is primary.

## 7. Monte Carlo summaries

Every reported mean includes:

- number scheduled;
- number successful;
- mean/median;
- standard deviation/IQR;
- Monte Carlo standard error;
- paired difference and paired MCSE where relevant.

Coverage/type-I/power includes a binomial interval. A difference of a few percentage points is not described as meaningful when it is within Monte Carlo noise.

The simulation summaries cover the following families, each with its own denominator:

- bias, empirical standard deviation, mean and median estimated standard error, their ratio, interval length at 0.95 and 0.90, coverage at both levels, type-I error on null coordinates, and power on active coordinates;
- studentised statistic quantiles with a normal reference, and null p-value uniformity reported as binned counts together with a Kolmogorov--Smirnov statistic and its Monte Carlo standard error;
- selection summaries: true-positive rate, false-discovery proportion, false-positive rate, exact support recovery, selected size, screening sure-inclusion, and the across-replication selection frequency of each active coordinate;
- prediction summaries: held-out check loss, skill against the pooled quantile-regression baseline and against the intercept-only quantile predictor, quantile calibration, and the cluster-level loss quantiles;
- theory summaries: smoothing-target distance, effective-Hessian maximum error, Schur identity discrepancy, precision-row error, Dantzig residual, scaled Bahadur remainder, cluster-score skewness and excess kurtosis, and log-log slopes of the principal error terms against `n` with their Monte Carlo standard errors;
- computation summaries for `module=G_computation` only: median, interquartile range and maximum of `runtime_sec`, `first_stage_sec`, `refit_sec`, `inference_sec` and `peak_memory_mb`, together with iteration counts and failure rates. These rows are never combined with statistical performance tables.

The summary directory `results/tables/<run_id>/` contains:

```text
method_summary.csv
coverage_summary.csv
coverage_by_scenario.csv
inference_summary.csv
null_pvalue_uniformity.csv
studentised_quantiles.csv
selection_stability.csv
theory_summary.csv
rate_slopes.csv
computation_summary.csv
paired_pinball_differences.csv
paired_method_differences.csv
registry_accounting.csv
```

`coverage_by_scenario.csv` aggregates the target coordinates of each scenario and reports the replicate-level coverage with its Monte Carlo standard error, so that correlated coordinates inside one replication are not treated as independent units. `registry_accounting.csv` records the scheduled, recorded, successful, warning, failed and not-implemented counts of every registered `experiment_id x method_id` cell; a final run must show zero missing replications. `rate_slopes.csv` is written only for cells that contain at least three distinct cluster counts, and `computation_summary.csv` only for the computation module; their absence states that the corresponding design is not part of the run.

## 8. Application outputs

The GSE65391 pipeline must create:

```text
download_manifest.csv
metadata_dictionary.csv
cohort_flow.csv
outcome_audit.csv
visit_distribution.csv
missingness_by_variable.csv
batch_qc.csv
subject_folds.csv
preprocessing_manifest.yml
fold_level_predictions.csv
fold_level_metrics.csv
confirmatory_module_inference.csv
exploratory_split_inference.csv
stability_selection.csv
fit_diagnostics.csv
application_failure_log.csv
```

No manuscript table may precede the cohort-flow and outcome-audit files.

## 9. Manuscript generation

Final LaTeX tables/figures must be generated by scripts from validated summary files. At minimum, the build script checks:

- no placeholder `pending` remains in a final numerical table;
- table denominators match the run manifest;
- every reported row is bound to a registry `experiment_id`, and no row label exists without a matching registry row;
- the replication count printed in the manuscript equals the frozen registry count;
- all plotted confidence intervals use the declared level;
- method labels match `docs/BENCHMARK_MATRIX.md`;
- no legacy calibration constants appear in code or output;
- all displayed values have an underlying unrounded raw value.

`scripts/simulation/10_build_manuscript_assets.R` writes the table environments into `manuscript/generated/` and the figures into `manuscript/generated/figures/`. The table files are `tab_simulation_main.tex`, `tab_simulation_methods.tex`, `tab_simulation_inference.tex`, `tab_simulation_theory.tex`, `tab_simulation_rates.tex`, `tab_simulation_computation.tex`, `tab_simulation_failures.tex` and `tab_simulation_selection.tex`; the figure files are `fig_coverage_clusters.pdf`, `fig_se_sd_ratio.pdf`, `fig_rate_slopes.pdf`, `fig_type1_power.pdf`, `fig_null_pvalue_uniformity.pdf`, `fig_studentised_quantiles.pdf` and `fig_computation_time.pdf`. A non-finite summary entry is printed as `--` together with the caption's statement of what the entry marks, which is distinct from an unfilled placeholder.

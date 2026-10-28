source("research/submission_completion_20261027/load.R")
source("research/submission_completion_20261027/quantile_highs.R")
reference_procedure <- fit_benchmark_frozen_tuning_v2
fit_benchmark_frozen_tuning_v2 <- function(method_id, train, tau, target_coords, tuning, seed,
                                          context = list(), control = list()) {
  control$pooled_optimizer <- "dual_lp"
  reference_procedure(method_id, train, tau, target_coords, tuning, seed, context, control)
}
output <- "results/scratch/submission_completion_20261027/pooled_dual_check"
if (dir.exists(output)) stop("The complete pooled dual comparison already exists.")
dir.create(output, recursive = TRUE)
registry <- submission_registry("main")
config <- config_row_to_list_v2(registry[registry$experiment_id == "A01", ])
config$methods <- "POOLED-QR-LASSO"
candidate <- run_one_replication_v2(".", config, 1L, base_seed = 202610272L,
  final = TRUE, allow_missing_benchmarks = FALSE, run_id = "submission_pooled_dual_check")
saveRDS(candidate, file.path(output, "result.rds"))
reference <- readRDS("results/scratch/submission_completion_20261027/comparator_preflight/POOLED-QR-LASSO.rds")
first <- reference$method_diagnostics[["POOLED-QR-LASSO"]]
second <- candidate$method_diagnostics[["POOLED-QR-LASSO"]]
difference <- max(abs(c(first$beta_hat, first$prediction_intercept) - c(second$beta_hat, second$prediction_intercept)))
if (candidate$metrics$status != reference$metrics$status || difference > 2e-5)
  stop("The complete pooled quantile procedure changed across LP formulations.")
record <- data.frame(maximum_coefficient_difference = difference,
  reference_seconds = reference$metrics$runtime_sec, candidate_seconds = candidate$metrics$runtime_sec)
write.csv(record, file.path(output, "checks.csv"), row.names = FALSE)
print(record)

source("research/submission_completion_20261027/load_comparators.R")
source("research/submission_completion_20261027/dispatch_state.R")
output <- "results/scratch/submission_completion_20261027/pooled_procedure_check"
if (dir.exists(output)) stop("The complete pooled comparison already exists.")
dir.create(output, recursive = TRUE)
registry <- submission_registry("main")
config <- config_row_to_list_v2(registry[registry$experiment_id == "A01", ])
config$methods <- "POOLED-QR-LASSO"
reference_environment <- new.env(parent = .GlobalEnv)
reference_environment$fit_benchmark_frozen_tuning_v2 <- submission_comparator_reference
reference_runner <- run_one_replication_v2
environment(reference_runner) <- reference_environment
reference <- reference_runner(".", config, 1L, base_seed = 202610272L,
  final = TRUE, allow_missing_benchmarks = FALSE, run_id = "submission_pooled_reference_check")
submission_atomic_rds(reference, file.path(output, "reference.rds"))
candidate <- run_one_replication_v2(".", config, 1L, base_seed = 202610272L,
  final = TRUE, allow_missing_benchmarks = FALSE, run_id = "submission_pooled_numerical_check")
saveRDS(candidate, file.path(output, "result.rds"))
first <- reference$method_diagnostics[["POOLED-QR-LASSO"]]
second <- candidate$method_diagnostics[["POOLED-QR-LASSO"]]
difference <- max(abs(first$beta_hat - second$beta_hat))
if (candidate$metrics$status != reference$metrics$status || difference > 2e-5)
  stop("The complete pooled quantile procedure changed after column generation.")
record <- data.frame(maximum_coefficient_difference = difference,
  reference_seconds = reference$metrics$runtime_sec, candidate_seconds = candidate$metrics$runtime_sec,
  reference_solver = "quantreg", candidate_solver = "HiGHS_dual_LP")
write.csv(record, file.path(output, "checks.csv"), row.names = FALSE)
print(record)

source("research/submission_completion_20261027/load_comparators.R")
output <- "results/scratch/submission_completion_20261027/qgee_procedure_check"
if (dir.exists(output)) stop("The complete QGEE comparison already exists.")
dir.create(output, recursive = TRUE)
reference <- readRDS("results/scratch/submission_completion_20261027/comparator_preflight/QGEE-SCAD.rds")
config <- reference$config
config$methods <- "QGEE-SCAD"
candidate <- run_one_replication_v2(".", config, 1L, base_seed = 202610272L,
  final = TRUE, allow_missing_benchmarks = FALSE, run_id = "submission_qgee_procedure_check")
saveRDS(candidate, file.path(output, "result.rds"))
first <- reference$method_diagnostics[["QGEE-SCAD"]]
second <- candidate$method_diagnostics[["QGEE-SCAD"]]
difference <- max(abs(c(first$beta_hat, first$prediction_intercept) - c(second$beta_hat, second$prediction_intercept)))
if (candidate$metrics$status != reference$metrics$status || difference > 2e-5)
  stop("The complete QGEE procedure changed after numerical optimisation.")
record <- data.frame(maximum_coefficient_difference = difference,
  reference_seconds = reference$metrics$runtime_sec, candidate_seconds = candidate$metrics$runtime_sec,
  status = candidate$metrics$status, converged = candidate$metrics$converged)
write.csv(record, file.path(output, "checks.csv"), row.names = FALSE)
print(record)

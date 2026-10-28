source("research/submission_completion_20261027/load_comparators.R")
source("research/submission_completion_20261027/dispatch_state.R")
output <- "results/scratch/submission_completion_20261027/sqr_procedure_check"
if (dir.exists(output)) stop("The complete SQR comparison already exists.")
dir.create(output, recursive = TRUE)
registry <- submission_registry("main")
config <- config_row_to_list_v2(registry[registry$experiment_id == "A01", ])
config$methods <- "SQR-DEBIASED-IID"
reference_environment <- new.env(parent = .GlobalEnv)
reference_environment$fit_benchmark_frozen_tuning_v2 <- submission_comparator_reference
reference_runner <- run_one_replication_v2
environment(reference_runner) <- reference_environment
reference <- reference_runner(".", config, 1L, base_seed = 202610272L,
  final = TRUE, allow_missing_benchmarks = FALSE, run_id = "submission_sqr_reference_check")
submission_atomic_rds(reference, file.path(output, "reference.rds"))
candidate <- run_one_replication_v2(".", config, 1L, base_seed = 202610272L,
  final = TRUE, allow_missing_benchmarks = FALSE, run_id = "submission_sqr_hybrid_check")
submission_atomic_rds(candidate, file.path(output, "candidate.rds"))
validate_paired_failed_coordinate_accounting_v2(reference$metrics, reference$coordinates)
validate_paired_failed_coordinate_accounting_v2(candidate$metrics, candidate$coordinates)
if (!all(reference$metrics$converged) || !all(candidate$metrics$converged))
  stop("Both complete SQR procedures must converge for numerical comparison.")
keys <- c("coordinate", "beta_tilde", "estimated_se", "ci_lower", "ci_upper")
first <- reference$coordinates[, keys]
second <- candidate$coordinates[, keys]
second <- second[match(first$coordinate, second$coordinate), ]
if (anyNA(second)) stop("The SQR coordinate comparison is incomplete.")
difference <- max(abs(as.matrix(first[, -1L]) - as.matrix(second[, -1L])))
first_tuning <- reference$method_diagnostics[["SQR-DEBIASED-IID"]]$tuning
second_tuning <- candidate$method_diagnostics[["SQR-DEBIASED-IID"]]$tuning
if (difference > 2e-5) stop("The complete SQR estimators differ beyond numerical accuracy.")
selection <- c("selected_multiplier", "selected_gamma_multiplier", "selected_gamma")
if (!all(selection %in% names(first_tuning)) || !all(selection %in% names(second_tuning)))
  stop("The complete SQR tuning records are unavailable.")
tuning_difference <- max(abs(unlist(first_tuning[selection]) - unlist(second_tuning[selection])))
if (tuning_difference > 1e-12) stop("The complete SQR tuning decisions differ.")
record <- data.frame(maximum_difference = difference, tuning_difference = tuning_difference,
  reference_seconds = reference$metrics$runtime_sec, candidate_seconds = candidate$metrics$runtime_sec,
  reference_solver = "full_HiGHS_complete_CV", candidate_solver = "hybrid_certified_exact_CV_pruning")
write.csv(record, file.path(output, "checks.csv"), row.names = FALSE)
print(record)

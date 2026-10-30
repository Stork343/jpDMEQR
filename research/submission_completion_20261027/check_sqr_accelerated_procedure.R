source("research/submission_completion_20261027/load_comparators.R")
source("research/submission_completion_20261027/load_compiled_clime.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Provide a recorded setting, replication and validation label.")
id <- args[1L]
replicate <- as.integer(args[2L])
label <- args[3L]
path <- file.path("results/raw/simulation/submission_20261028_main_sqr/replications",
  id, sprintf("replicate_%05d/result.rds", replicate))
reference <- readRDS(path)
output <- file.path("results/scratch/sqr_acceleration_20261030", paste0("full_", id, "_", replicate, "_", label))
if (dir.exists(output)) stop("The full procedure validation already exists.")
dir.create(output, recursive = TRUE)
globalCallingHandlers(error = function(condition) saveRDS(condition, file.path(output, "failure.rds")))
runner_environment <- new.env(parent = .GlobalEnv)
runner_environment$fit_benchmark_frozen_tuning_v2 <- function(method_id, train, tau, target_coords,
  tuning, seed, context = list(), control = list()) {
  if (method_id != "SQR-DEBIASED-IID") stop("This validation requires the registered SQR procedure.")
  control$solver_opts <- list(algorithm = "compiled_certified", native_guard = TRUE, column_loss_bounds = TRUE)
  control$exact_cv_pruning <- TRUE
  submission_comparator_reference(method_id, train, tau, target_coords, tuning, seed, context, control)
}
runner <- run_one_replication_v2
environment(runner) <- runner_environment
candidate <- runner(".", reference$config, replicate, base_seed = 20261027L,
  allow_missing_benchmarks = FALSE, final = TRUE,
  run_id = paste0("sqr_acceleration_", label))
saveRDS(candidate, file.path(output, "result.rds"), compress = "xz")
if (candidate$seed != reference$seed || !all(candidate$metrics$converged) ||
    !all(reference$metrics$converged)) stop("The matched procedure did not produce both completed fits.")
coordinates <- c("coordinate", "beta_tilde", "estimated_se", "ci_lower", "ci_upper")
first <- reference$coordinates[, coordinates]
second <- candidate$coordinates[, coordinates]
second <- second[match(first$coordinate, second$coordinate), ]
if (anyNA(second)) stop("The requested SQR intervals are incomplete.")
interval_difference <- max(abs(as.matrix(first[, -1L]) - as.matrix(second[, -1L])))
first_diagnostic <- reference$method_diagnostics[["SQR-DEBIASED-IID"]]
second_diagnostic <- candidate$method_diagnostics[["SQR-DEBIASED-IID"]]
selection <- c("selected_multiplier", "selected_gamma_multiplier", "selected_gamma")
tuning_difference <- max(abs(unlist(first_diagnostic$tuning[selection]) - unlist(second_diagnostic$tuning[selection])))
coefficient_difference <- max(abs(first_diagnostic$beta_hat - second_diagnostic$beta_hat))
if (interval_difference > 2e-5 || tuning_difference > 1e-12 || coefficient_difference > 2e-5)
  stop("The complete estimator, intervals or selected tuning changed.")
record <- data.frame(experiment_id = id, replicate = replicate,
  reference_seconds = reference$metrics$runtime_sec,
  candidate_seconds = candidate$metrics$runtime_sec,
  speedup = reference$metrics$runtime_sec / candidate$metrics$runtime_sec,
  interval_difference = interval_difference, tuning_difference = tuning_difference,
  coefficient_difference = coefficient_difference)
write.csv(record, file.path(output, "checks.csv"), row.names = FALSE)
print(record)

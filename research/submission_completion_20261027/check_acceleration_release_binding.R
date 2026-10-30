source("research/submission_completion_20261027/load_comparators.R")
reference <- readRDS("results/raw/simulation/submission_20261028_main_sqr/replications/A01/replicate_00001/result.rds")
candidate <- run_one_replication_v2(".", reference$config, 1L, base_seed = 20261027L,
  allow_missing_benchmarks = FALSE, final = TRUE, run_id = "sqr_release_binding_20261030")
if (!all(candidate$metrics$converged)) stop("The released SQR execution did not converge.")
fields <- c("beta_tilde", "estimated_se", "ci_lower", "ci_upper")
if (!identical(candidate$coordinates$coordinate, reference$coordinates$coordinate))
  stop("The released target coordinates differ.")
difference <- max(abs(as.matrix(candidate$coordinates[, fields]) - as.matrix(reference$coordinates[, fields])))
first <- reference$method_diagnostics[[1L]]$tuning
second <- candidate$method_diagnostics[[1L]]$tuning
selection <- c("selected_multiplier", "selected_gamma_multiplier", "selected_gamma")
tuning_difference <- max(abs(unlist(first[selection]) - unlist(second[selection])))
if (difference > 2e-5 || tuning_difference > 1e-12 ||
    second$gamma_cv$execution != "certified_column_and_fold_loss_bounds")
  stop("The released SQR calculation differs from its accepted definition.")
jsonlite::write_json(list(passed = TRUE,
  implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  interval_difference = difference, tuning_difference = tuning_difference,
  runtime_seconds = candidate$metrics$runtime_sec,
  numerical_specification = submission_numerical_specification$version),
  "results/preflight/sqr_acceleration_20261030/release_binding.json", auto_unbox = TRUE, pretty = TRUE)
print(c(interval_difference = difference, tuning_difference = tuning_difference,
  runtime_seconds = candidate$metrics$runtime_sec))

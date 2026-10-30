source("research/submission_completion_20261027/load_comparators.R")
source("research/submission_completion_20261027/load_compiled_clime.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Provide a registered SQR setting and replication.")
id <- args[1L]
replicate <- as.integer(args[2L])
registry <- submission_registry("main")
config <- config_row_to_list_v2(registry[registry$experiment_id == id, ])
if (!"SQR-DEBIASED-IID" %in% config$methods) stop("The setting does not register SQR.")
config$methods <- "SQR-DEBIASED-IID"
output <- file.path("results/scratch/sqr_acceleration_20261030", paste0("registered_", id, "_", replicate))
if (dir.exists(output)) stop("The registered-case validation already exists.")
dir.create(output, recursive = TRUE)
globalCallingHandlers(error = function(condition) saveRDS(condition, file.path(output, "failure.rds")))
scope <- new.env(parent = .GlobalEnv)
scope$fit_benchmark_frozen_tuning_v2 <- function(method_id, train, tau, target_coords, tuning,
  seed, context = list(), control = list()) {
  control$solver_opts <- list(algorithm = "compiled_certified", native_guard = TRUE, column_loss_bounds = TRUE)
  control$exact_cv_pruning <- TRUE
  submission_comparator_reference(method_id, train, tau, target_coords, tuning, seed, context, control)
}
runner <- run_one_replication_v2
environment(runner) <- scope
result <- runner(".", config, replicate, base_seed = 20261027L,
  allow_missing_benchmarks = FALSE, final = TRUE, run_id = "sqr_registered_validation_20261030")
saveRDS(result, file.path(output, "result.rds"), compress = "xz")
validate_paired_failed_coordinate_accounting_v2(result$metrics, result$coordinates)
if (!all(result$metrics$converged)) stop("The registered SQR fit did not converge.")
diagnostic <- result$method_diagnostics[[1L]]
data <- prepare_analysis_data_v2(simulate_from_config_v2(config, result$seed), config, result$seed)$data
data$X <- cbind(`(Intercept)` = 1, data$X)
data$p <- ncol(data$X)
cv <- diagnostic$tuning$gamma_cv
checks <- list()
for (fold in seq_len(5L)) {
  training <- subset_clusters_v2(data, cv$roles[[fold]]$training_subjects)
  validation <- subset_clusters_v2(data, cv$roles[[fold]]$validation_subjects)
  fit <- cv$fold_initial_estimators[[fold]]
  cell <- cv$cells[cv$cells$fold == fold & cv$cells$candidate == cv$selected_candidate, ]
  if (nrow(cell) != 1L || cell$status != "ok") stop("A selected validation fold is incomplete.")
  H <- sqr_iid_hessian_fn_v2(training$y, training$X, fit$beta, config$tau, fit$h)
  V <- sqr_iid_hessian_fn_v2(validation$y, validation$X, fit$beta, config$tau, fit$h)
  reference <- solve_clime_matrix_iid_v2(H, cell$gamma, solver_opts = list(algorithm = "hybrid_certified"))
  if (!reference$feasible) stop("The independent selected-fold precision is unavailable.")
  loss <- sum((profile_dense_crossprod_v2(V, reference$omega_symmetric) - diag(nrow(H)))^2) / nrow(H)
  difference <- abs(loss - cell$loss)
  if (difference > 2e-6) stop("The selected full-matrix validation loss differs.")
  checks[[fold]] <- data.frame(fold = fold, gamma = cell$gamma, loss_difference = difference)
}
fit <- diagnostic$tuning$full_training$fit
H <- sqr_iid_hessian_fn_v2(data$y, data$X, fit$beta, config$tau, fit$h)
reference <- solve_clime_matrix_iid_v2(H, diagnostic$tuning$selected_gamma,
  solver_opts = list(algorithm = "hybrid_certified"))
if (!reference$feasible) stop("The independent full-training precision is unavailable.")
coordinates <- match(result$coordinates$coordinate, colnames(data$X))
omega <- reference$omega_symmetric[coordinates, , drop = FALSE]
score <- sqr_iid_score_fn_v2(data$y, data$X, fit$beta, config$tau, fit$h)
estimates <- fit$beta[coordinates] - drop(omega %*% score)
standard_errors <- sqrt(config$tau * (1 - config$tau) * colSums((data$X %*% t(omega))^2)) / nrow(data$X)
difference <- max(abs(c(estimates - result$coordinates$beta_tilde,
  standard_errors - result$coordinates$estimated_se)))
if (difference > 2e-5) stop("The independent coefficient correction or variance differs.")
write.csv(do.call(rbind, checks), file.path(output, "fold_checks.csv"), row.names = FALSE)
record <- data.frame(experiment_id = id, replicate = replicate, p = config$p,
  n_clusters = config$n_clusters, candidate_seconds = result$metrics$runtime_sec,
  interval_difference = difference, maximum_validation_loss_difference = max(vapply(checks, function(x) x$loss_difference, numeric(1))))
write.csv(record, file.path(output, "checks.csv"), row.names = FALSE)
print(record)

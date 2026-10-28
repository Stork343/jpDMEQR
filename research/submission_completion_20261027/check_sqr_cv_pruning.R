source("research/submission_completion_20261027/load.R")
source("research/submission_completion_20261027/clime_certified.R")
source("research/submission_completion_20261027/clime_constraint_generation.R")
source("research/submission_completion_20261027/sqr_cv_pruning.R")
registry <- submission_registry("main")
records <- list()
for (tau in c(.25, .5, .75)) {
  config <- config_row_to_list_v2(registry[registry$experiment_id == "A01", ])
  config$n_clusters <- 40L
  config$p <- 12L
  config$s <- 3L
  config$tau <- tau
  seed <- 202610280L + as.integer(tau * 100)
  data <- simulate_from_config_v2(config, seed)
  coordinates <- resolve_target_coordinates_v2(config, data)
  tuning <- build_tuning_for_analysis_v2(config, data)
  control <- list(solver_opts = list(algorithm = "constraint_generation", subproblem = "clarabel"),
    fit_control = list(optimizer = "matrix_free_fista", max_iter = 10000L, max_backtrack = 50L, beta_tol = 1e-7))
  reference <- fit_benchmark_frozen_tuning_v2("SQR-DEBIASED-IID", data, tau, coordinates, tuning, seed,
    context = list(config = config), control = control)
  control$exact_cv_pruning <- TRUE
  candidate <- fit_benchmark_frozen_tuning_v2("SQR-DEBIASED-IID", data, tau, coordinates, tuning, seed,
    context = list(config = config), control = control)
  difference <- max(abs(c(reference$beta_tilde, reference$se) - c(candidate$beta_tilde, candidate$se)))
  if (reference$status != candidate$status || difference > 2e-5 ||
      reference$tuning_diagnostics$selected_gamma != candidate$tuning_diagnostics$selected_gamma ||
      reference$tuning_diagnostics$selected_multiplier != candidate$tuning_diagnostics$selected_multiplier)
    stop("Certified SQR tuning pruning changed the selected estimator.")
  old <- reference$tuning_diagnostics$gamma_cv$candidates
  new <- candidate$tuning_diagnostics$gamma_cv$candidates
  complete <- which(new$eligible)
  if (any(abs(old$loss[complete] - new$loss[complete]) > 1e-8)) stop("An evaluated SQR fold loss changed.")
  record <- data.frame(tau = tau, maximum_difference = difference,
    reference_seconds = reference$runtime_sec, candidate_seconds = candidate$runtime_sec,
    pruned_candidates = sum(new$certified_pruning),
    selected_gamma = candidate$tuning_diagnostics$selected_gamma)
  records[[length(records) + 1L]] <- record
  write.csv(do.call(rbind, records), "results/scratch/submission_completion_20261027/sqr_pruning_checks.csv", row.names = FALSE)
  print(record)
}

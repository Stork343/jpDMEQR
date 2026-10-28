source("research/submission_completion_20261027/load.R")
source("research/submission_completion_20261027/clime_primal.R")
source("research/submission_completion_20261027/clime_parametric.R")
registry <- read_simulation_registry_v2("config/simulation_main.csv")
config <- config_row_to_list_v2(registry[registry$experiment_id == "A01", ])
config$method_version <- "profile_pair_20261022"
output <- "results/scratch/submission_comparator_probes_20261027/SQR_A01_parametric_warm"
if (dir.exists(output)) stop("The complete SQR numerical comparison already exists.")
dir.create(output, recursive = TRUE)
seed <- seed_from_id_v2(20261027L, "A01", 1L, "numerical_probe")
dat <- prepare_analysis_data_v2(simulate_from_config_v2(config, seed), config, seed)$data
coordinates <- resolve_target_coordinates_v2(config, dat)
started <- proc.time()[[3L]]
globalCallingHandlers(error = function(condition) saveRDS(condition, file.path(output, "failure.rds")))
fit <- fit_benchmark_frozen_tuning_v2("SQR-DEBIASED-IID", dat, config$tau, coordinates,
  build_tuning_for_analysis_v2(config, dat), seed, context = list(config = config),
  control = list(solver_opts = list(algorithm = "parametric_warm"),
    fit_control = list(optimizer = "matrix_free_fista", max_iter = 10000L, max_backtrack = 50L, beta_tol = 1e-7)))
elapsed <- proc.time()[[3L]] - started
saveRDS(list(fit = fit, elapsed = elapsed), file.path(output, "result.rds"))
reference <- readRDS("results/scratch/submission_comparator_probes_20261027/SQR-DEBIASED-IID_A01_matrix_free_fista/result.rds")
error <- max(abs(c(fit$beta_tilde, fit$se) - c(reference$fit$beta_tilde, reference$fit$se)))
if (fit$status != reference$fit$status || error > 2e-5 ||
    fit$tuning_diagnostics$selected_multiplier != reference$fit$tuning_diagnostics$selected_multiplier ||
    abs(fit$inference_object$precision$gamma - reference$fit$inference_object$precision$gamma) > 1e-12)
  stop("The complete SQR procedure changed across LP formulations.")
record <- data.frame(maximum_difference = error, reference_seconds = reference$elapsed,
  candidate_seconds = elapsed, gamma = fit$inference_object$precision$gamma,
  multiplier = fit$tuning_diagnostics$selected_multiplier)
write.csv(record, file.path(output, "comparison.csv"), row.names = FALSE)
print(record)

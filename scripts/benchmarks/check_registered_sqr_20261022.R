args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
id <- if (length(args) >= 2L) args[2] else "P01"
mode <- if (length(args) >= 3L) args[3] else "fixed"
if (!mode %in% c("fixed", "cv")) stop("Choose fixed numerical tuning or the complete source CV procedure.")
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, environment())
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
registry <- read_simulation_registry_v2(file.path(root, "config", "simulation_pair_pilot_20261022.yml"))
row <- registry[registry$experiment_id == id, , drop = FALSE]
if (nrow(row) != 1L) stop("Supply an existing frozen pilot experiment.")
config <- config_row_to_list_v2(row)
if (!"SQR-DEBIASED-IID" %in% config$methods) stop("The frozen experiment does not include SQR.")
seed <- seed_from_id_v2(20260817L, id, 1L, "training")
dat <- prepare_analysis_data_v2(simulate_from_config_v2(config, seed), config, seed)$data
tuning <- build_tuning_for_analysis_v2(config, dat)
coordinates <- resolve_target_coordinates_v2(config, dat)
output <- file.path(root, "results", "scratch", "profile_pair_benchmark_20261022",
                     paste0("registered_sqr_", mode, "_", id, "_", format(Sys.time(), "%Y%m%dT%H%M%OS6")))
if (dir.exists(output)) stop("The registered numerical check directory already exists.")
dir.create(output, recursive = TRUE)
request <- list(config = config, seed = seed, replicate = 1L, tuning = tuning, tuning_mode = mode,
                 coordinates = coordinates, n_rows = nrow(dat$X), p = ncol(dat$X),
                 implementation_commit = current_commit_v2(root), purpose = "numerical_acceptance_not_pilot_results")
saveRDS(request, file.path(output, "request.rds"))
started <- proc.time()[[3L]]
withCallingHandlers({
  result <- if (mode == "cv") {
    fit_benchmark_frozen_tuning_v2(
      "SQR-DEBIASED-IID", dat, config$tau, coordinates, tuning,
      seed_from_id_v2(20260817L, id, 1L, "method-SQR-DEBIASED-IID"), context = list(config = config),
      control = list(fit_control = list(max_iter = 2000L, max_backtrack = 50L, beta_tol = 1e-7))
    )
  } else fit_benchmark_sqr_debiased_iid_v2(dat, config$tau, coordinates, tuning,
                                          seed_from_id_v2(20260817L, id, 1L, "SQR-DEBIASED-IID"))
  saveRDS(result, file.path(output, "result.rds"))
  precision <- result$inference_object$precision
  jsonlite::write_json(list(experiment_id = id, status = result$status,
                            runtime_sec = proc.time()[[3L]] - started,
                            n_clusters = dat$n_clusters, p = ncol(dat$X), n_rows = nrow(dat$X),
                            gamma = precision$gamma, solver = precision$backend,
                            actual_LP_calculations = nrow(precision$matrix_grid_record$attempts),
                            model_constructions = precision$matrix_grid_record$model_constructions,
                            requested_multiplier_grid = config$dantzig_multipliers,
                            tuning_mode = mode, method_version = result$method_version,
                            selected_lambda_multiplier = result$tuning_diagnostics$selected_multiplier,
                            gamma_candidates = result$tuning_diagnostics$gamma_cv$candidates),
                       file.path(output, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
}, error = function(condition) {
  saveRDS(condition, file.path(output, "failure.rds"))
  jsonlite::write_json(list(experiment_id = id, status = "failed", message = conditionMessage(condition),
                            failure_stage = condition$failure_stage %||% "SQR_numerical_acceptance",
                            runtime_sec = proc.time()[[3L]] - started),
                       file.path(output, "failed_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
})

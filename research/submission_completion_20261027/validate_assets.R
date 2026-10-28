source("research/submission_completion_20261027/load.R")
registry <- submission_registry("main")
records <- list()
output <- "results/preflight/submission_population_20261027"
if (dir.exists(output)) stop("The population validation already exists.")
dir.create(output, recursive = TRUE)
for (i in which(registry$target_mode == "profile_mc")) {
  config <- config_row_to_list_v2(registry[i, ])
  path <- submission_target_path(config)
  asset <- readRDS(path)
  if (config$experiment_id %in% c("B16", "B19")) {
    if (!identical(asset$design_identity$design, profile_pair_design_identity_v2(config)) ||
        asset$development || asset$n_population < 100000L || asset$repeats < 4L || !asset$all_converged ||
        anyDuplicated(asset$seeds) || length(asset$repeat_assets) != asset$repeats ||
        !isTRUE(all.equal(colMeans(asset$estimates_by_repeat), asset$beta_star_mc, tolerance = 1e-12))) stop("Conditional-scale population identity failed.")
    for (record in asset$repeat_assets) for (role in c("target", "refinement", "analysis")) {
      fit <- record[[role]]
      if (record$n_population != asset$n_population || !fit$converged || !fit$nuisance_converged ||
          fit$kkt_residual > asset$acceptance$score_tol ||
          fit$max_nuisance_gradient > asset$acceptance$nuisance_tol ||
          fit$profile_identity_error > asset$acceptance$identity_tol) stop("Conditional-scale population numerical verification failed.")
    }
  } else validate_profile_pair_target_asset_v2(asset, config, final = TRUE)
  records[[length(records) + 1L]] <- data.frame(experiment_id = config$experiment_id, type = "target",
    path = path, clusters = asset$n_population, repetitions = asset$repeats, passed = TRUE)
  write.csv(do.call(rbind, records), file.path(output, "checks.csv"), row.names = FALSE)
  cat("Validated target", config$experiment_id, "\n")
}
for (i in which(grepl("JOINT-PROFILE-POP-H", registry$methods, fixed = TRUE))) {
  config <- config_row_to_list_v2(registry[i, ])
  asset <- submission_direction(config)
  validate_profile_pair_direction_asset_v2(asset, config, final = TRUE)
  records[[length(records) + 1L]] <- data.frame(experiment_id = config$experiment_id, type = "direction",
    path = profile_pair_direction_asset_path_v2(".", config$experiment_id),
    clusters = asset$n_population, repetitions = asset$repeats, passed = TRUE)
  write.csv(do.call(rbind, records), file.path(output, "checks.csv"), row.names = FALSE)
  cat("Validated direction", config$experiment_id, "\n")
}
jsonlite::write_json(list(passed = TRUE, checked_assets = length(records),
  minimum_clusters = min(vapply(records, `[[`, numeric(1), "clusters")),
  minimum_repetitions = min(vapply(records, `[[`, numeric(1), "repetitions")),
  completed = as.character(Sys.time())), file.path(output, "completed.json"), auto_unbox = TRUE, pretty = TRUE)

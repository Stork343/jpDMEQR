#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/04_build_profile_targets.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
config_path <- cli$config %||% file.path(root, "config", "simulation_main.csv")
n_population <- as.integer(cli$`n-population` %||% 10000L)
repeats <- as.integer(cli$repeats %||% 2L)
base_seed <- as.integer(cli$seed %||% 87001L)
only_id <- cli$`experiment-id` %||% ""
if (n_population < 1000L) stop("n-population must be at least 1000 clusters.")
if (repeats < 2L) stop("At least two independent population runs are required.")

cfg_df <- read_simulation_registry_v2(config_path)
cfg_df <- cfg_df[cfg_df$target_mode == "profile_mc", , drop = FALSE]
if (nzchar(only_id)) cfg_df <- cfg_df[cfg_df$experiment_id == only_id, , drop = FALSE]
if (!nrow(cfg_df)) stop("No profile_mc configurations selected.")

out_dir <- file.path(root, "results", "intermediate", "targets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
all_diagnostics <- list(); diag_index <- 0L
for (ii in seq_len(nrow(cfg_df))) {
  cfg <- config_row_to_list_v2(cfg_df[ii, , drop = FALSE])
  target_columns <- sort(unique(c(seq_len(cfg$s),
                                  seq.int(cfg$s + 1L, min(cfg$p, cfg$s + 3L)))))
  estimates <- matrix(NA_real_, nrow = repeats, ncol = cfg$p,
                      dimnames = list(NULL, sprintf("x%05d", seq_len(cfg$p))))
  fits <- vector("list", repeats)
  for (rr in seq_len(repeats)) {
    seed <- seed_from_id_v2(base_seed, paste0(cfg$experiment_id, "_target"), rr)
    message(sprintf("Approximating %s: repeat %d/%d with %d population clusters",
                    cfg$experiment_id, rr, repeats, n_population))
    fit_obj <- approximate_profile_target_v2(
      config = cfg,
      target_columns = target_columns,
      n_population = n_population,
      seed = seed,
      fit_control = list(
        max_iter = 350L, beta_tol = 2e-7, kkt_tol = 2e-6,
        nuisance_control = list(reltol = 1e-10, maxit = 300L, grad_tol = 1e-7)
      )
    )
    full <- setNames(rep(0, cfg$p), sprintf("x%05d", seq_len(cfg$p)))
    full[names(fit_obj$beta_star_mc)] <- fit_obj$beta_star_mc
    estimates[rr, ] <- full
    fits[[rr]] <- fit_obj
    diag_index <- diag_index + 1L
    all_diagnostics[[diag_index]] <- data.frame(
      experiment_id = cfg$experiment_id, repeat = rr, seed = seed,
      n_population = n_population, target_dimension = length(target_columns),
      converged = fit_obj$converged, kkt_residual = fit_obj$kkt_residual,
      profile_identity_error = fit_obj$fit$components$profile_identity_error,
      max_nuisance_gradient = fit_obj$fit$components$max_nuisance_gradient,
      stringsAsFactors = FALSE
    )
  }
  beta_star_mc <- colMeans(estimates)
  se_by_coordinate <- apply(estimates, 2, stats::sd) / sqrt(repeats)
  target_mc_se <- max(se_by_coordinate[target_columns], na.rm = TRUE)
  between_run_max_difference <- max(apply(estimates[, target_columns, drop = FALSE], 2,
                                          function(x) max(x) - min(x)), na.rm = TRUE)
  obj <- list(
    experiment_id = cfg$experiment_id,
    beta_star_mc = beta_star_mc,
    target_mc_se = target_mc_se,
    target_mc_se_by_coordinate = se_by_coordinate,
    estimates_by_repeat = estimates,
    target_columns = target_columns,
    n_population = n_population,
    repeats = repeats,
    between_run_max_difference = between_run_max_difference,
    config = cfg,
    fits = fits,
    implementation_commit = current_commit_v2(root),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  saveRDS(obj, file.path(out_dir, paste0(cfg$experiment_id, ".rds")), compress = "xz")
  write.csv(data.frame(
    coordinate = names(beta_star_mc), beta_star_mc = beta_star_mc,
    target_mc_se = se_by_coordinate, is_optimised_coordinate = seq_len(cfg$p) %in% target_columns,
    stringsAsFactors = FALSE
  ), file.path(out_dir, paste0(cfg$experiment_id, "_coordinates.csv")), row.names = FALSE)
}
write.csv(do.call(rbind, all_diagnostics),
          file.path(out_dir, "profile_target_diagnostics.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Profile target approximations written to: ", out_dir)
message("Final paper runs should use >=100000 population clusters and >=4 repeats after optimization.")

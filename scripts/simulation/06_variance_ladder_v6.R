#!/usr/bin/env Rscript
# Round-6 P06-v5 T=L+Q+R variance decomposition and four-level meat ladder
# (METHOD_SPECIFICATION_ROUND6_AMENDMENT.md sections 2-4; actions item 4-8).
# Usage:
#   Rscript scripts/simulation/06_variance_ladder_v6.R <run_id> <shard> <n_shards> <config>
# After all shards finish, the same script aggregates when shard=0 (or use the
# aggregate entry point with --aggregate).
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/06_variance_ladder_v6.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
setwd(root)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source_v2_module(root, "benchmark_adapters_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
run_id <- args[1] %||% "pilot_v6_ladder_P06"
shard <- as.integer(args[2] %||% "1")
n_shards <- as.integer(args[3] %||% "1")
config_path <- if (length(args) > 3 && nzchar(args[4])) {
  file.path(root, "config", args[4])
} else file.path(root, "config", "simulation_pilot_v5.csv")
B <- 50L
experiment <- "P06"
out_dir <- file.path(root, "results", "raw", "simulation", run_id)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cfg <- read.csv(config_path, stringsAsFactors = FALSE)
row <- cfg[cfg$experiment_id == experiment, ]
config <- config_row_to_list_v2(row)
tau <- config$tau; lambda_gamma <- config$lambda_gamma

# POP-H population direction asset (P06).
asset_path <- population_direction_asset_path_v2(root, experiment)
if (!file.exists(asset_path)) stop("POP-H asset missing: ", asset_path)
asset <- readRDS(asset_path)
n_per_shard <- ceiling(B / n_shards)
reps <- seq.int(from = shard, to = B, by = n_shards)
for (r in reps) {
  seed <- seed_from_id_v2(20260817L, experiment, r, "training")
  dat0 <- simulate_from_config_v2(config, seed)
  prepared <- prepare_analysis_data_v2(dat0, config, seed)
  dat <- prepared$data
  beta_target <- dat$beta0
  target_names <- intersect(resolve_target_coordinates_v2(config, dat), names(beta_target))
  tuning <- build_tuning_for_analysis_v2(config, dat)
  fit <- tryCatch(
    fit_benchmark_profile_dqr_v2(
      dat, tau, target_names, tuning, seed,
      control = list(
        ci_level = 0.95,
        fit_control = list(
          max_iter = 2000L, max_backtrack = 50L,
          beta_tol = 1e-7, kkt_normalized_tol = 1e-3,
          nuisance_control = list(reltol = 1e-11, maxit = 400L, grad_tol = 1e-8)
        )
      )
    ),
    error = function(e) e
  )
  if (inherits(fit, "error") || !identical(fit$status, "ok")) {
    row_out <- data.frame(
      replicate = r, seed = seed, coordinate = NA_character_,
      T = NA_real_, L = NA_real_, Q = NA_real_, R = NA_real_,
      V_fit_fit = NA_real_, V_fit_target = NA_real_,
      V_pop_fit = NA_real_, V_pop_target = NA_real_,
      status = if (inherits(fit, "error")) paste("error:", conditionMessage(fit))
      else fit$failure_stage %||% fit$status,
      stringsAsFactors = FALSE
    )
  } else {
    vals <- lapply(target_names, function(nm) {
      omega_pop <- setNames(rep(0, length(names(dat$beta0))), names(dat$beta0))
      pop <- asset$directions[[nm]]$omega
      common <- intersect(names(pop), names(omega_pop))
      omega_pop[common] <- pop[common]
      out <- tryCatch(
        variance_ladder_v6_one_v2(
          fit, dat, beta_target, omega_pop, nm, tuning$h_inf, lambda_gamma
        ),
        error = function(e) list(T = NA_real_, L = NA_real_, Q = NA_real_, R = NA_real_,
                                 V_fit_fit = NA_real_, V_fit_target = NA_real_,
                                 V_pop_fit = NA_real_, V_pop_target = NA_real_)
      )
      data.frame(
        replicate = r, seed = seed, coordinate = nm,
        T = out$T, L = out$L, Q = out$Q, R = out$R,
        V_fit_fit = out$V_fit_fit, V_fit_target = out$V_fit_target,
        V_pop_fit = out$V_pop_fit, V_pop_target = out$V_pop_target,
        status = "ok", stringsAsFactors = FALSE
      )
    })
    row_out <- do.call(rbind, vals)
  }
  shard_file <- file.path(out_dir, sprintf("ladder_%s_s%d.csv", run_id, shard))
  write.table(row_out, shard_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(shard_file), append = file.exists(shard_file))
  cat(sprintf("shard %d: done rep %d\n", shard, r))
}
cat(sprintf("SHARD DONE %d/%d reps=%s\n", shard, n_shards, paste(reps, collapse = ",")))
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
source_v2_module(root, "metrics_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
config_path <- cli$config %||% file.path(root, "config", "simulation_main.csv")
development <- as_bool_cli_v2(cli$development, FALSE)
n_population <- as.integer(cli$`n-population` %||% if (development) 10000L else 100000L)
repeats <- as.integer(cli$repeats %||% if (development) 2L else 4L)
direction_repeats <- as.integer(cli$`direction-repeats` %||% if (development) 2L else 4L)
base_seed <- as.integer(cli$seed %||% 87001L)
only_id <- cli$`experiment-id` %||% ""
build_pop_h <- as_bool_cli_v2(cli$`build-pop-h`, TRUE)
audit_structural <- as_bool_cli_v2(cli$`audit-structural`, TRUE)

if (!development) {
  if (n_population < 100000L) stop("Final target construction requires at least 100000 clusters.")
  if (repeats < 4L) stop("Final target construction requires at least four repeats.")
  if (direction_repeats < 4L) stop("Final population directions require at least four repeats.")
} else {
  if (n_population < 1000L || repeats < 2L) {
    stop("Development target construction requires >=1000 clusters and >=2 repeats.")
  }
  message("DEVELOPMENT TARGET MODE: assets cannot satisfy the final freeze gate.")
}

cfg_df <- read_simulation_registry_v2(config_path)
if (nzchar(only_id)) {
  ids <- trimws(strsplit(only_id, ",", fixed = TRUE)[[1]])
  cfg_df <- cfg_df[cfg_df$experiment_id %in% ids, , drop = FALSE]
}
if (!nrow(cfg_df)) stop("No registry configurations selected.")
commit <- current_commit_v2(root)
config_sha <- file_checksum_v2(config_path)

target_dir <- file.path(root, "results", "intermediate", "targets")
direction_dir <- file.path(root, "results", "intermediate", "population_directions")
audit_dir <- file.path(root, "results", "intermediate", "target_audits")
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(direction_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
all_diagnostics <- list(); dd <- 0L

build_target_for_config <- function(cfg, repeats_local, label = "profile_target") {
  target_columns <- population_target_columns_v2(cfg)
  run_one <- function(rr) {
    seed <- seed_from_id_v2(base_seed, paste0(cfg$experiment_id, "_", label), rr, "population-target")
    message(sprintf(
      "%s %s repeat %d/%d with %d population clusters",
      label, cfg$experiment_id, rr, repeats_local, n_population
    ))
    approximate_profile_target_v2(
      config = cfg,
      target_columns = target_columns,
      n_population = n_population,
      seed = seed,
      fit_control = list(
        max_iter = 350L,
        beta_tol = 2e-7,
        kkt_tol = 2e-6,
        nuisance_control = list(reltol = 1e-10, maxit = 300L, grad_tol = 1e-7)
      ),
      retain_fit = FALSE
    )
  }
  # Parallel repeats: each repeat is an independent deterministic run with
  # its own seed. mclapply is used on non-Windows; Windows falls back to a
  # sequential loop (correctness identical, just slower).
  repeat_cores <- getOption("jpDMEQR.repeat_cores", 1L)
  fits <- if (repeat_cores > 1L && .Platform$OS.type != "windows" &&
              repeats_local > 1L) {
    parallel::mclapply(seq_len(repeats_local), run_one,
                       mc.cores = min(repeat_cores, repeats_local),
                       mc.preschedule = TRUE)
  } else {
    lapply(seq_len(repeats_local), run_one)
  }
  obj <- summarise_profile_target_repeats_v2(
    fits, cfg, target_columns, implementation_commit = commit
  )
  obj$config_sha256 <- config_sha
  obj$development <- development
  obj$label <- label
  obj
}

# Profile-MC targets.
profile_rows <- cfg_df[cfg_df$target_mode == "profile_mc", , drop = FALSE]
for (ii in seq_len(nrow(profile_rows))) {
  cfg <- config_row_to_list_v2(profile_rows[ii, , drop = FALSE])
  obj <- build_target_for_config(cfg, repeats, "profile_target")
  saveRDS(obj, target_asset_path_v2(root, cfg$experiment_id), compress = "xz")
  coord <- data.frame(
    coordinate = names(obj$beta_star_mc),
    beta_star_mc = as.numeric(obj$beta_star_mc),
    target_mc_sd = as.numeric(obj$target_mc_sd_by_coordinate),
    target_mc_se = as.numeric(obj$target_mc_se_by_coordinate),
    between_repeat_range = as.numeric(obj$between_run_range_by_coordinate),
    stringsAsFactors = FALSE
  )
  write_atomic_csv_v2(coord, file.path(target_dir, paste0(cfg$experiment_id, "_coordinates.csv")))
  for (rr in seq_along(obj$fits)) {
    dd <- dd + 1L
    fit <- obj$fits[[rr]]
    all_diagnostics[[dd]] <- data.frame(
      experiment_id = cfg$experiment_id,
      asset_type = "profile_target",
      repeat_id = rr,
      seed = fit$seed,
      n_population = fit$n_population,
      converged = fit$converged,
      kkt_residual = fit$kkt_residual,
      profile_identity_error = fit$profile_identity_error,
      max_nuisance_gradient = fit$max_nuisance_gradient,
      nuisance_converged = fit$nuisance_converged,
      stringsAsFactors = FALSE
    )
  }
}

# Prospective structural-target audit for the nuisance-penalty sensitivity rows.
if (audit_structural) {
  audit_rows <- cfg_df[
    cfg_df$module == "E_tuning" & cfg_df$target_mode == "structural" &
      cfg_df$lambda_gamma %in% c(0.25, 0.5, 2, 4), , drop = FALSE
  ]
  for (ii in seq_len(nrow(audit_rows))) {
    cfg <- config_row_to_list_v2(audit_rows[ii, , drop = FALSE])
    obj <- build_target_for_config(cfg, repeats, "structural_target_audit")
    beta0 <- make_beta_v2(max(obj$target_columns), cfg$s, cfg$signal)
    beta0 <- beta0[names(obj$beta_star_mc)]
    obj$beta0 <- beta0
    obj$target_displacement <- obj$beta_star_mc - beta0
    obj$max_absolute_displacement <- max(abs(obj$target_displacement))
    saveRDS(obj, file.path(audit_dir, paste0(cfg$experiment_id, ".rds")), compress = "xz")
    write_atomic_csv_v2(data.frame(
      coordinate = names(obj$beta_star_mc),
      beta0 = as.numeric(beta0),
      beta_star_mc = as.numeric(obj$beta_star_mc),
      displacement = as.numeric(obj$target_displacement),
      target_mc_se = as.numeric(obj$target_mc_se_by_coordinate),
      stringsAsFactors = FALSE
    ), file.path(audit_dir, paste0(cfg$experiment_id, "_coordinates.csv")))
  }
}

# Population effective-Hessian/direction assets.
if (build_pop_h) {
  pop_rows <- cfg_df[vapply(cfg_df$methods, function(x) {
    "PROFILE-DQR-POP-H" %in% parse_pipe_character_v2(x)
  }, logical(1)), , drop = FALSE]
  for (ii in seq_len(nrow(pop_rows))) {
    cfg <- config_row_to_list_v2(pop_rows[ii, , drop = FALSE])
    target_columns <- population_target_columns_v2(cfg)
    target_names <- sprintf("x%05d", target_columns)
    if (cfg$target_mode == "profile_mc") {
      target_obj <- readRDS(target_asset_path_v2(root, cfg$experiment_id))
      beta_target <- target_obj$beta_star_mc[target_names]
    } else {
      beta_target <- make_beta_v2(max(target_columns), cfg$s, cfg$signal)[target_names]
    }
    requested <- unique(c(
      head(target_names[target_columns <= cfg$s], ceiling(cfg$target_coordinate_count / 2)),
      head(target_names[target_columns > cfg$s], floor(cfg$target_coordinate_count / 2))
    ))
    direction_fits <- vector("list", direction_repeats)
    for (rr in seq_len(direction_repeats)) {
      seed <- seed_from_id_v2(base_seed + 100000L, paste0(cfg$experiment_id, "_popH"),
                              rr, "population-direction")
      message(sprintf(
        "population direction %s repeat %d/%d with %d clusters",
        cfg$experiment_id, rr, direction_repeats, n_population
      ))
      direction_fits[[rr]] <- approximate_population_direction_v2(
        config = cfg,
        beta_target = beta_target,
        target_columns = target_columns,
        requested_coordinates = requested,
        n_population = n_population,
        seed = seed
      )
    }
    H_mean <- Reduce(`+`, lapply(direction_fits, `[[`, "H_population")) /
      length(direction_fits)
    inv <- stable_inverse_symmetric_v2(H_mean)
    Omega <- inv$inverse
    colnames(Omega) <- rownames(Omega) <- target_names
    directions <- lapply(requested, function(nm) {
      k <- match(nm, target_names)
      omega <- Omega[, k]
      e <- numeric(length(target_names)); e[k] <- 1
      list(
        coordinate = nm,
        omega = setNames(as.numeric(omega), target_names),
        residual = max(abs(as.numeric(H_mean %*% omega - e))),
        l1_norm = sum(abs(omega)),
        l2_norm = sqrt(sum(omega^2))
      )
    })
    names(directions) <- requested
    asset <- list(
      experiment_id = cfg$experiment_id,
      H_population = H_mean,
      Omega_population = Omega,
      directions = directions,
      repeat_assets = direction_fits,
      target_columns = target_columns,
      target_names = target_names,
      beta_target = beta_target,
      n_population = n_population,
      repeats = direction_repeats,
      condition_number = inv$condition_number,
      raw_eigenvalues = inv$eigenvalues,
      eigen_floor = inv$eigen_floor,
      nuisance_converged = all(vapply(direction_fits, `[[`, logical(1), "nuisance_converged")),
      profile_identity_error = max(vapply(direction_fits, `[[`, numeric(1), "profile_identity_error")),
      max_nuisance_gradient = max(vapply(direction_fits, `[[`, numeric(1), "max_nuisance_gradient")),
      implementation_commit = commit,
      config_sha256 = config_sha,
      dependency_hash = pop_h_dependency_hash_v2(
        cfg, n_population, length(direction_fits),
        direction_fits[[1]]$n_analysis, direction_fits[[1]]$h_analysis
      ),
      config = cfg,
      development = development,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    saveRDS(asset, population_direction_asset_path_v2(root, cfg$experiment_id), compress = "xz")
    write_atomic_csv_v2(data.frame(
      coordinate = names(directions),
      residual = vapply(directions, `[[`, numeric(1), "residual"),
      omega_l1 = vapply(directions, `[[`, numeric(1), "l1_norm"),
      omega_l2 = vapply(directions, `[[`, numeric(1), "l2_norm"),
      stringsAsFactors = FALSE
    ), file.path(direction_dir, paste0(cfg$experiment_id, "_directions.csv")))
  }
}

if (length(all_diagnostics)) {
  write_atomic_csv_v2(
    do.call(rbind, all_diagnostics),
    file.path(target_dir, "profile_target_diagnostics.csv")
  )
}
write_json_v2(list(
  development = development,
  commit_sha = commit,
  config_sha256 = config_sha,
  n_population = n_population,
  repeats = repeats,
  direction_repeats = direction_repeats,
  profile_target_ids = profile_rows$experiment_id,
  population_direction_ids = if (exists("pop_rows")) pop_rows$experiment_id else character(),
  structural_audit_ids = if (exists("audit_rows")) audit_rows$experiment_id else character(),
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(target_dir, "build_manifest.json"))
capture.output(sessionInfo(), file = file.path(target_dir, "sessionInfo.txt"))
message("Target and population-direction assets written under results/intermediate/.")

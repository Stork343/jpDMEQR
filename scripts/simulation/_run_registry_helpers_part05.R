profile_pair_gate_configs_v2 <- function(cfg_df) {
  if (!nrow(cfg_df) || !"method_version" %in% names(cfg_df) ||
      anyNA(cfg_df$method_version) ||
      any(cfg_df$method_version != "profile_pair_20261022") ||
      anyDuplicated(cfg_df$experiment_id)) {
    stop("Paired target gates require a nonempty paired registry with unique experiment identities.")
  }
  configs <- lapply(seq_len(nrow(cfg_df)), function(i) config_row_to_list_v2(cfg_df[i, , drop = FALSE]))
  names(configs) <- cfg_df$experiment_id
  configs
}

profile_pair_structural_audit_path_v2 <- function(root, experiment_id) {
  file.path(root, "results", "intermediate", "profile_pair_20261022", "construction_v3",
            "structural_audits", paste0(experiment_id, ".rds"))
}

profile_pair_gate_audit_designs_v2 <- function(configs) {
  selected <- configs[vapply(configs, function(config) {
    identical(config$module, "E_tuning") && identical(config$target_mode, "structural") &&
      config$lambda_gamma %in% c(0.25, 0.5, 2, 4)
  }, logical(1))]
  if (length(selected)) {
    baseline <- if ("E03" %in% names(selected)) selected[["E03"]] else selected[[1]]
    baseline$experiment_id <- "LAMBDA1_BASELINE"
    baseline$lambda_gamma <- 1
    selected[[baseline$experiment_id]] <- baseline
  }
  selected
}

profile_pair_gate_requirements_v2 <- function(root, cfg_df) {
  configs <- profile_pair_gate_configs_v2(cfg_df)
  rows <- list()
  add <- function(id, kind, path) {
    rows[[length(rows) + 1L]] <<- data.frame(experiment_id = id, asset = kind,
                                           path = path, stringsAsFactors = FALSE)
  }
  for (config in configs) {
    if (identical(config$target_mode, "profile_mc")) {
      add(config$experiment_id, "profile_target", profile_pair_target_asset_path_v2(root, config$experiment_id))
    }
  }
  for (config in configs) {
    if ("PROFILE-DQR-POP-H" %in% config$methods) {
      add(config$experiment_id, "population_direction", profile_pair_direction_asset_path_v2(root, config$experiment_id))
    }
  }
  for (id in names(profile_pair_gate_audit_designs_v2(configs))) {
    add(id, "structural_target_audit", profile_pair_structural_audit_path_v2(root, id))
  }
  if (!length(rows)) return(data.frame(experiment_id = character(), asset = character(), path = character()))
  do.call(rbind, rows)
}

validate_profile_pair_structural_audit_gate_v2 <- function(object, registered_design, root,
                                                           final = TRUE, expected_commit = NULL) {
  previous_root <- getOption("jpDMEQR.source_root")
  on.exit(options(jpDMEQR.source_root = previous_root), add = TRUE)
  options(jpDMEQR.source_root = normalizePath(root, mustWork = TRUE))
  version <- profile_pair_asset_version_v2()
  if (!identical(object$asset_kind, "paired_structural_target_audit") ||
      !identical(object$audit_construction_version, "paired-structural-audit-20261022-v1") ||
      any(!vapply(names(version), function(field) identical(object[[field]], version[[field]]), logical(1)))) {
    stop("The structural audit has an incompatible paired definition or construction identity.")
  }
  columns <- population_target_columns_v2(registered_design)
  expected <- registered_design
  expected$p <- max(columns)
  if (!identical(object$experiment_id, registered_design$experiment_id) ||
      !identical(object$audit_scope, "all_active_coefficients_and_three_null_controls") ||
      !identical(object$audit_columns, columns) ||
      object$registered_p != registered_design$p || object$audit_p != expected$p ||
      !identical(profile_pair_design_identity_v2(object$audit_design),
                 profile_pair_design_identity_v2(registered_design)) ||
      !identical(profile_pair_design_identity_v2(object$config), profile_pair_design_identity_v2(expected))) {
    stop("The paired structural audit does not match the registered marginal design.")
  }
  source <- object$registered_config
  if (identical(object$experiment_id, "LAMBDA1_BASELINE")) {
    if (!identical(object$registered_experiment_id, "E03") ||
        !identical(source$experiment_id, "E03") ||
        !identical(object$audit_design_origin, "frozen_E_design_with_requested_Lambda1_reference")) {
      stop("The Lambda1 audit does not identify its frozen E-design source.")
    }
    source$experiment_id <- "LAMBDA1_BASELINE"
    source$lambda_gamma <- 1
  } else if (!identical(object$audit_design_origin, "registered_E_sensitivity") ||
             !identical(object$registered_experiment_id, registered_design$experiment_id)) {
    stop("The paired structural audit has an inconsistent registered source.")
  }
  if (!identical(profile_pair_design_identity_v2(source), profile_pair_design_identity_v2(registered_design))) {
    stop("The paired structural audit source differs from the required nuisance-penalty design.")
  }
  target <- object$paired_target_asset
  validate_profile_pair_target_asset_v2(target, expected, final = final, expected_commit = expected_commit)
  beta0 <- make_beta_v2(expected$p, expected$s, expected$signal)
  same <- identical(object$beta0, beta0) &&
    identical(object$beta_star_mc, target$beta_star_mc) &&
    identical(object$target_displacement, target$beta_star_mc - beta0) &&
    identical(object$target_mc_se_by_coordinate, target$target_mc_se_by_coordinate) &&
    identical(object$estimates_by_repeat, target$estimates_by_repeat) &&
    identical(object$target_displacement_by_repeat, sweep(target$estimates_by_repeat, 2, beta0)) &&
    identical(object$implementation_commit, target$implementation_commit) &&
    identical(object$n_population, target$n_population) && identical(object$repeats, target$repeats) &&
    identical(object$all_converged, target$all_converged) &&
    identical(object$target_mc_se, target$target_mc_se) &&
    identical(object$h_target, target$h_target) && identical(object$h_analysis, target$h_analysis) &&
    identical(object$max_kkt_residual, target$max_kkt_residual) &&
    identical(object$max_nuisance_gradient, target$max_nuisance_gradient) &&
    identical(object$profile_identity_error, target$profile_identity_error)
  if (!same) stop("The structural audit wrapper differs from its independently fitted paired roots.")
  validation <- validate_structural_target_audit_v2(object, final = final)
  if (!validation$valid) stop(paste(validation$problems, collapse = "; "))
  invisible(TRUE)
}

validate_profile_pair_target_gate_v2 <- function(root, cfg_df, final = TRUE) {
  previous_root <- getOption("jpDMEQR.source_root")
  on.exit(options(jpDMEQR.source_root = previous_root), add = TRUE)
  options(jpDMEQR.source_root = normalizePath(root, mustWork = TRUE))
  configs <- profile_pair_gate_configs_v2(cfg_df)
  audit_designs <- profile_pair_gate_audit_designs_v2(configs)
  required <- profile_pair_gate_requirements_v2(root, cfg_df)
  expected_commit <- if (isTRUE(final)) current_commit_v2(root) else NULL
  targets <- list()
  rows <- vector("list", nrow(required))
  for (i in seq_len(nrow(required))) {
    record <- required[i, , drop = FALSE]
    id <- record$experiment_id
    pass <- file.exists(record$path)
    reason <- if (pass) "accepted" else "missing"
    if (pass) {
      object <- readRDS(record$path)
      if (record$asset == "profile_target") {
        validate_profile_pair_target_asset_v2(object, configs[[id]], final, expected_commit)
        targets[[id]] <- list(beta = object$beta_star_mc, seeds = object$seeds,
                              commit = object$implementation_commit,
                              construction_version = object$construction_version,
                              h_target = object$h_target)
      } else if (record$asset == "population_direction") {
        validate_profile_pair_direction_asset_v2(object, configs[[id]], final, expected_commit)
        if (identical(configs[[id]]$target_mode, "profile_mc")) {
          target <- targets[[id]]
          if (is.null(target)) {
            pass <- FALSE
            reason <- "matching paired target is missing"
          } else if (!identical(object$beta_target, target$beta) ||
                     !identical(object$target_asset_identity$beta_star_mc, target$beta) ||
                     !identical(object$target_asset_identity$seeds, target$seeds) ||
                     !identical(object$target_asset_identity$implementation_commit, target$commit) ||
                     !identical(object$target_asset_identity$construction_version, target$construction_version) ||
                     !identical(object$target_asset_identity$h_target, target$h_target) ||
                     !identical(object$h_target, target$h_target)) {
            stop("Population direction and paired target identities differ: ", id)
          }
        }
      } else {
        validate_profile_pair_structural_audit_gate_v2(object, audit_designs[[id]], root,
                                                        final, expected_commit)
      }
      rm(object)
      invisible(gc(FALSE))
    }
    rows[[i]] <- data.frame(experiment_id = id, asset = record$asset,
                            pass = pass, reason = reason, stringsAsFactors = FALSE)
  }
  table <- if (length(rows)) do.call(rbind, rows) else data.frame(
    experiment_id = character(), asset = character(), pass = logical(), reason = character())
  list(pass = all(table$pass), table = table)
}

profile_pair_freeze_asset_paths_v2 <- function(root, cfg_df) {
  required <- profile_pair_gate_requirements_v2(root, cfg_df)
  prefix <- c(profile_target = "profile_target:", population_direction = "population_direction:",
              structural_target_audit = "structural_audit:")
  data.frame(label = paste0(unname(prefix[required$asset]), required$experiment_id),
              path = required$path, stringsAsFactors = FALSE)
}

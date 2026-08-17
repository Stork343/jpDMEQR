# One-replication execution and replication-level result construction.

build_tuning_for_analysis_v2 <- function(config, dat) {
  values <- reference_tuning_values_v2(config, p = ncol(dat$X), n = dat$n_clusters)
  lambda_multiplier <- config$lambda_beta_multipliers[
    which.min(abs(config$lambda_beta_multipliers - 1))
  ]
  list(
    h = config$h_multiplier * dat$n_clusters^(-1 / 3),
    lambda_beta = lambda_multiplier * sqrt(log(max(ncol(dat$X), 2)) / dat$n_clusters),
    lambda_gamma = config$lambda_gamma,
    mu_grid = config$dantzig_multipliers *
      (sqrt(log(max(ncol(dat$X), 2)) /
              (dat$n_clusters * (config$h_multiplier * dat$n_clusters^(-1 / 3)))) +
         (config$h_multiplier * dat$n_clusters^(-1 / 3))^2),
    full_reference = values
  )
}

method_requires_low_dimensional_subset_v2 <- function(method_id) {
  method_id %in% c("LQMM", "ORACLE-LQMM", "BIAS-ADJ-LQMM")
}

extract_selected_global_v2 <- function(fit, method_train, analysis_dat, screening) {
  selected_raw <- fit$selected %||% integer()
  if (!length(selected_raw)) return(integer())
  if (is.character(selected_raw)) {
    selected_names <- selected_raw
  } else if (!is.null(fit$beta_hat) && length(fit$beta_hat) >= max(selected_raw, 0)) {
    selected_names <- names(fit$beta_hat)[selected_raw]
  } else {
    selected_names <- colnames(method_train$X)[selected_raw]
  }
  selected_names <- selected_names[!is.na(selected_names)]
  match(selected_names, names(analysis_dat$beta0), nomatch = 0L)
}

inference_direction_table_v2 <- function(fit) {
  inf <- fit$inference_object %||% NULL
  if (is.null(inf) || is.null(inf$table)) return(data.frame())
  inf$table
}

method_theory_diagnostics_v2 <- function(config, replicate, seed, method_id,
                                         fit, target_obj, population_asset,
                                         analysis_dat) {
  fit_obj <- fit$fit_object %||% NULL
  H <- if (!is.null(fit_obj$components$hessian)) fit_obj$components$hessian else NULL
  eig <- if (!is.null(H) && all(is.finite(H))) {
    eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  } else numeric()
  hessian_max_error <- NA_real_
  precision_row_error <- NA_real_
  if (!is.null(H) && !is.null(population_asset$H_population)) {
    pop_names <- colnames(population_asset$H_population)
    fit_names <- names(fit_obj$beta)
    common <- intersect(pop_names, fit_names)
    if (length(common)) {
      hessian_max_error <- max(abs(
        H[match(common, fit_names), match(common, fit_names), drop = FALSE] -
          population_asset$H_population[common, common, drop = FALSE]
      ))
    }
  }
  inf <- inference_direction_table_v2(fit)
  if (nrow(inf) && !is.null(population_asset$directions)) {
    errs <- numeric()
    for (nm in intersect(inf$coordinate, names(population_asset$directions))) {
      direction <- fit$inference_object$directions[[which(inf$coordinate == nm)[1]]]
      if (!is.null(direction$omega)) {
        pop <- population_asset$directions[[nm]]$omega
        fit_names <- names(fit_obj$beta)
        full <- setNames(rep(0, length(fit_names)), fit_names)
        full[intersect(names(pop), fit_names)] <- pop[intersect(names(pop), fit_names)]
        errs <- c(errs, sqrt(sum((direction$omega - full)^2)))
      }
    }
    if (length(errs)) precision_row_error <- max(errs)
  }
  correction <- if (nrow(inf)) max(abs(inf$beta_tilde - inf$beta_hat), na.rm = TRUE) else NA_real_
  if (!is.finite(correction)) correction <- NA_real_
  data.frame(
    experiment_id = config$experiment_id,
    replicate = replicate,
    seed = seed,
    method_id = method_id,
    target_type = target_obj$target_type,
    target_mc_se = target_obj$target_mc_se,
    hessian_max_error = hessian_max_error,
    hessian_operator_error = NA_real_,
    hessian_min_eigenvalue = if (length(eig)) min(eig) else NA_real_,
    hessian_condition_number = if (length(eig) && min(abs(eig)) > 0) {
      max(abs(eig)) / min(abs(eig))
    } else NA_real_,
    profile_identity_error = fit$profile_identity_error %||%
      fit_obj$components$profile_identity_error %||% NA_real_,
    max_nuisance_gradient = fit$max_nuisance_gradient %||%
      fit_obj$components$max_nuisance_gradient %||% NA_real_,
    kkt_residual = fit$kkt_residual %||% NA_real_,
    max_dantzig_residual = if (nrow(inf) && "dantzig_residual" %in% names(inf) &&
                                any(is.finite(inf$dantzig_residual))) {
      max(inf$dantzig_residual[is.finite(inf$dantzig_residual)])
    } else NA_real_,
    max_omega_l1 = if (nrow(inf) && "omega_l1" %in% names(inf) &&
                       any(is.finite(inf$omega_l1))) {
      max(inf$omega_l1[is.finite(inf$omega_l1)])
    } else NA_real_,
    precision_row_error = precision_row_error,
    one_step_correction_max = correction,
    scaled_bahadur_remainder = NA_real_,
    score_skewness = NA_real_,
    score_excess_kurtosis = NA_real_,
    n_clusters = analysis_dat$n_clusters,
    p = ncol(analysis_dat$X),
    stringsAsFactors = FALSE
  )
}

run_one_replication_v2 <- function(root, config, replicate,
                                   base_seed = 20260817L,
                                   allow_missing_benchmarks = TRUE,
                                   final = FALSE) {
  seed <- seed_from_id_v2(base_seed, config$experiment_id, replicate, "training")
  commit <- current_commit_v2(root)
  dat0 <- simulate_from_config_v2(config, seed)
  prepared <- prepare_analysis_data_v2(dat0, config, seed)
  dat <- prepared$data
  target_obj <- load_profile_target_v2(root, config, dat, final = final)
  beta_target <- target_obj$beta
  target_names <- intersect(resolve_target_coordinates_v2(config, dat), names(beta_target))
  if (!length(target_names)) {
    target_names <- head(names(beta_target), min(2L, length(beta_target)))
  }
  tuning <- build_tuning_for_analysis_v2(config, dat)
  need_pop_h <- "PROFILE-DQR-POP-H" %in% config$methods
  population_asset <- load_population_direction_asset_v2(
    root, config, required = need_pop_h, final = final
  )

  test_seed <- seed_from_id_v2(base_seed, config$experiment_id, replicate, "test")
  test_dat0 <- simulate_from_config_v2(config, test_seed)
  test_dat <- apply_feature_map_to_test_v2(
    test_dat0, prepared$screening, config$fit_random_effects
  )
  if (!identical(names(test_dat$beta0), names(dat$beta0))) {
    stop("Training and test feature maps are inconsistent.")
  }

  method_rows <- list()
  coordinate_rows <- list()
  theory_rows <- list()
  mcount <- 0L
  ccount <- 0L
  tcount <- 0L

  for (method_id in config$methods) {
    method_train <- dat
    method_targets <- target_names
    context <- list(
      active = dat$active,
      screening = prepared$screening,
      population_direction_asset = population_asset,
      config = config,
      beta_target = beta_target
    )
    if (method_requires_low_dimensional_subset_v2(method_id)) {
      low <- common_lqmm_subset_v2(
        dat, target_names, config$tau,
        max_features = 20L,
        active = if (method_id == "ORACLE-LQMM") dat$active else NULL
      )
      method_train <- low$data
      method_targets <- intersect(target_names, colnames(method_train$X))
      context$common_screen_keep <- low$keep
    }

    start <- proc.time()[[3]]
    fit <- tryCatch(
      fit_benchmark_by_id_v2(
        method_id = method_id,
        train = method_train,
        tau = config$tau,
        target_coords = method_targets,
        tuning = tuning,
        seed = seed_from_id_v2(base_seed, config$experiment_id, replicate,
                               paste0("method-", method_id)),
        context = context,
        control = list(
          ci_level = 0.95,
          fit_control = list(
            max_iter = 180L,
            beta_tol = 1e-6,
            kkt_tol = 1e-5,
            nuisance_control = list(
              reltol = 1e-11,
              maxit = 400L,
              grad_tol = 1e-8
            )
          ),
          max_features = 20L,
          screen_fraction = config$screen_fraction,
          screen_dim = config$screen_dim,
          # Method-specific pilot settings: keep the expensive tuning
          # procedures small in pilot/development runs while final runs use
          # the full frozen grids. These are debugging thresholds, not
          # calibration. Final execution must pass lambda = NULL so that
          # QGEE-SCAD uses its full HBIC grid and DOUBLE-PEN-QLMM the full
          # frozen grid; the development runner below sets pilot defaults.
          lambda = if (isTRUE(final)) NULL else 0.5,
          lambda_beta_grid = c(0.5, 1, 2),
          lambda_alpha_grid = c(0.5, 1),
          B = 100L
        )
      ),
      error = function(e) benchmark_failure_v2(
        method_id, "unknown", conditionMessage(e),
        proc.time()[[3]] - start
      )
    )

    if (identical(fit$status, "not_implemented") && !allow_missing_benchmarks) {
      stop("Required benchmark is not implemented: ", method_id,
           " -- ", fit$failure_message)
    }

    beta_full <- if (!is.null(fit$beta_hat)) {
      map_beta_to_full_v2(fit$beta_hat, names(dat$beta0))
    } else {
      setNames(rep(NA_real_, length(dat$beta0)), names(dat$beta0))
    }
    est <- if (all(is.finite(beta_full))) {
      estimation_metrics_v2(beta_full, beta_target, dat$active)
    } else {
      data.frame(
        l1_error = NA_real_, l2_error = NA_real_, max_error = NA_real_,
        max_active_error = NA_real_, active_bias_mean = NA_real_,
        null_bias_mean = NA_real_
      )
    }
    selected <- extract_selected_global_v2(
      fit, method_train, dat, prepared$screening
    )
    sel <- selection_metrics_v2(selected, dat$active, length(dat$beta0))
    pred <- if (all(is.finite(beta_full))) {
      as.numeric(test_dat$X %*% beta_full)
    } else rep(NA_real_, length(test_dat$y))
    loss <- if (all(is.finite(pred))) {
      pinball_loss_v2(test_dat$y, pred, config$tau)
    } else NA_real_

    mcount <- mcount + 1L
    method_rows[[mcount]] <- data.frame(
      experiment_id = config$experiment_id,
      module = config$module,
      replicate = replicate,
      seed = seed,
      method_id = method_id,
      status = fit$status,
      failure_stage = fit$failure_stage %||% "",
      failure_message = fit$failure_message %||% "",
      warning_messages = paste(fit$warning_messages %||% character(), collapse = " | "),
      reference_identifier = fit$reference_identifier %||% NA_character_,
      adapter_fidelity_status = fit$adapter_fidelity_status %||% NA_character_,
      n_clusters = dat$n_clusters,
      total_rows = dat$total_rows,
      p = ncol(dat$X),
      p_original = config$p,
      s = config$s,
      tau = config$tau,
      q = ncol(dat$Z),
      m_rule = config$m_rule,
      m_min = dat$m_summary$m_min,
      m_mean = dat$m_summary$m_mean,
      m_median = dat$m_summary$m_median,
      m_max = dat$m_summary$m_max,
      error_dist = config$error_dist,
      error_dependence = config$error_dependence,
      design_type = config$design_type,
      random_effect_dist = config$random_effect_dist,
      x_b_corr = config$x_b_corr,
      informative_size = config$informative_size,
      nonlinear_re_strength = config$nonlinear_re_strength,
      screen_method = prepared$screening$method,
      screen_dimension = prepared$screening$screen_dimension_realised,
      screen_sure_inclusion = prepared$screening$sure_inclusion,
      forced_target_count = length(prepared$screening$forced_target_names),
      h = tuning$h,
      lambda_beta = tuning$lambda_beta,
      lambda_gamma = tuning$lambda_gamma,
      selected_size = sel$selected_size,
      tpr = sel$tpr,
      fdp = sel$fdp,
      fpr = sel$fpr,
      exact_support = sel$exact_support,
      l1_error = est$l1_error,
      l2_error = est$l2_error,
      max_active_error = est$max_active_error,
      pinball_loss = loss,
      quantile_calibration = if (all(is.finite(pred))) mean(test_dat$y <= pred) else NA_real_,
      runtime_sec = fit$runtime_sec %||% (proc.time()[[3]] - start),
      converged = fit$converged %||% FALSE,
      kkt_residual = fit$kkt_residual %||% NA_real_,
      profile_identity_error = fit$profile_identity_error %||% NA_real_,
      max_nuisance_gradient = fit$max_nuisance_gradient %||% NA_real_,
      target_type = target_obj$target_type,
      target_mc_se = target_obj$target_mc_se,
      implementation_commit = commit,
      requested_workers = config$workers,
      stringsAsFactors = FALSE
    )

    inf_tab <- inference_direction_table_v2(fit)
    if (nrow(inf_tab)) {
      coord_names <- intersect(inf_tab$coordinate, names(beta_target))
      for (nm in coord_names) {
        rr <- which(inf_tab$coordinate == nm)[1]
        se <- as.numeric(inf_tab$se[rr])
        bt <- as.numeric(inf_tab$beta_tilde[rr])
        lower <- as.numeric(inf_tab$ci_lower[rr])
        upper <- as.numeric(inf_tab$ci_upper[rr])
        target_se <- if (!is.null(target_obj$target_mc_se_by_coordinate) &&
                         nm %in% names(target_obj$target_mc_se_by_coordinate)) {
          target_obj$target_mc_se_by_coordinate[[nm]]
        } else target_obj$target_mc_se
        ccount <- ccount + 1L
        coordinate_rows[[ccount]] <- data.frame(
          experiment_id = config$experiment_id,
          module = config$module,
          replicate = replicate,
          seed = seed,
          method_id = method_id,
          coordinate = nm,
          coordinate_type = if (match(nm, names(dat$beta0)) %in% dat$active) "active" else "null",
          target_type = target_obj$target_type,
          target_value = as.numeric(beta_target[nm]),
          target_mc_se = as.numeric(target_se),
          beta_hat = as.numeric(beta_full[nm]),
          beta_tilde = bt,
          estimated_se = se,
          ci_level = 0.95,
          ci_lower = lower,
          ci_upper = upper,
          covered = is.finite(lower) && is.finite(upper) &&
            lower <= beta_target[nm] && upper >= beta_target[nm],
          covered_target_minus_1.96mcse = is.finite(lower) && is.finite(upper) &&
            lower <= beta_target[nm] - 1.96 * target_se &&
            upper >= beta_target[nm] - 1.96 * target_se,
          covered_target_plus_1.96mcse = is.finite(lower) && is.finite(upper) &&
            lower <= beta_target[nm] + 1.96 * target_se &&
            upper >= beta_target[nm] + 1.96 * target_se,
          wald_z = if (is.finite(se) && se > 0) (bt - beta_target[nm]) / se else NA_real_,
          interval_length = upper - lower,
          feasible = inf_tab$feasible[rr] %||% is.finite(se),
          dantzig_mu = inf_tab$mu[rr] %||% NA_real_,
          dantzig_residual = inf_tab$dantzig_residual[rr] %||% NA_real_,
          omega_l1 = inf_tab$omega_l1[rr] %||% NA_real_,
          omega_l2 = inf_tab$omega_l2[rr] %||% NA_real_,
          status = fit$status,
          implementation_commit = commit,
          stringsAsFactors = FALSE
        )
      }
    }

    tcount <- tcount + 1L
    theory_rows[[tcount]] <- method_theory_diagnostics_v2(
      config, replicate, seed, method_id, fit, target_obj,
      population_asset, dat
    )
  }

  list(
    metrics = do.call(rbind, method_rows),
    coordinates = if (length(coordinate_rows)) do.call(rbind, coordinate_rows) else data.frame(),
    theory = if (length(theory_rows)) do.call(rbind, theory_rows) else data.frame(),
    screening = prepared$screening,
    config = config,
    seed = seed
  )
}

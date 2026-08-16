run_one_replication_v2 <- function(root, config, replicate, base_seed = 20260817L,
                                   allow_missing_benchmarks = TRUE) {
  seed <- seed_from_id_v2(base_seed, config$experiment_id, replicate)
  commit <- current_commit_v2(root)
  dat0 <- simulate_from_config_v2(config, seed)
  prepared <- prepare_analysis_data_v2(dat0, config, seed)
  dat <- prepared$data
  target_obj <- load_profile_target_v2(root, config, dat)
  beta_target <- target_obj$beta

  target_global <- unique(c(seq_len(min(2L, config$s)),
                            seq.int(config$s + 1L, min(config$p, config$s + 2L))))
  target_names <- intersect(names(dat$beta0), sprintf("x%05d", target_global))
  if (!length(target_names)) target_names <- head(names(dat$beta0), min(2L, length(dat$beta0)))

  tuning_values <- reference_tuning_values_v2(config)
  lambda_multiplier <- config$lambda_beta_multipliers[
    which.min(abs(config$lambda_beta_multipliers - 1))
  ]
  lambda_beta <- lambda_multiplier * sqrt(log(max(ncol(dat$X), 2)) / dat$n_clusters)
  tuning <- list(
    h = config$h_multiplier * dat$n_clusters^(-1 / 3),
    lambda_beta = lambda_beta,
    lambda_gamma = config$lambda_gamma,
    mu_grid = config$dantzig_multipliers *
      (sqrt(log(max(ncol(dat$X), 2)) /
              (dat$n_clusters * (config$h_multiplier * dat$n_clusters^(-1 / 3)))) +
         (config$h_multiplier * dat$n_clusters^(-1 / 3))^2)
  )

  test_dat <- simulate_from_config_v2(config, seed + 1000003L)
  if (ncol(test_dat$Z) > ncol(dat$Z)) test_dat$Z <- test_dat$Z[, seq_len(ncol(dat$Z)), drop = FALSE]
  if (!identical(prepared$screening$method, "none")) {
    test_dat$X <- test_dat$X[, prepared$screening$selected_global, drop = FALSE]
    test_dat$beta0 <- test_dat$beta0[prepared$screening$selected_global]
  }

  method_rows <- list()
  coordinate_rows <- list()
  mcount <- 0L
  ccount <- 0L

  for (method_id in config$methods) {
    method_train <- dat
    method_targets <- target_names
    context <- list(active = dat$active)
    if (method_id %in% c("LQMM", "ORACLE-LQMM")) {
      low <- common_lqmm_subset_v2(dat, target_names, config$tau,
                                   max_features = 20L)
      method_train <- low$data
      method_targets <- intersect(target_names, colnames(method_train$X))
    }

    start <- proc.time()[[3]]
    fit <- tryCatch(
      fit_benchmark_by_id_v2(
        method_id = method_id,
        train = method_train,
        tau = config$tau,
        target_coords = method_targets,
        tuning = tuning,
        seed = seed + match(method_id, config$methods) * 7919L,
        context = context,
        control = list(
          ci_level = 0.95,
          fit_control = list(max_iter = 180L, beta_tol = 1e-6,
                             kkt_tol = 1e-5,
                             nuisance_control = list(reltol = 1e-11,
                                                     maxit = 400L,
                                                     grad_tol = 1e-8)),
          max_features = 20L
        )
      ),
      error = function(e) benchmark_failure_v2(method_id, "unknown", conditionMessage(e),
                                                proc.time()[[3]] - start)
    )

    if (identical(fit$status, "not_implemented") && !allow_missing_benchmarks) {
      stop("Required benchmark is not implemented: ", method_id, " -- ", fit$failure_message)
    }

    beta_full <- if (!is.null(fit$beta_hat)) map_beta_to_full_v2(fit$beta_hat, names(dat$beta0)) else
      setNames(rep(NA_real_, length(dat$beta0)), names(dat$beta0))
    est <- if (all(is.finite(beta_full))) estimation_metrics_v2(beta_full, beta_target, dat$active) else
      data.frame(l1_error = NA_real_, l2_error = NA_real_, max_error = NA_real_,
                 max_active_error = NA_real_, active_bias_mean = NA_real_, null_bias_mean = NA_real_)
    selected_raw <- fit$selected %||% integer()
    if (!is.null(fit$beta_hat) && length(selected_raw) &&
        length(fit$beta_hat) != length(dat$beta0)) {
      selected_names <- names(fit$beta_hat)[selected_raw]
      selected_raw <- match(selected_names, names(dat$beta0))
      selected_raw <- selected_raw[!is.na(selected_raw)]
    }
    sel <- selection_metrics_v2(selected_raw, dat$active, length(dat$beta0))
    pred <- if (all(is.finite(beta_full))) as.numeric(test_dat$X %*% beta_full) else rep(NA_real_, length(test_dat$y))
    loss <- if (all(is.finite(pred))) pinball_loss_v2(test_dat$y, pred, config$tau) else NA_real_

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
      n_clusters = dat$n_clusters,
      total_rows = dat$total_rows,
      p = ncol(dat$X),
      s = config$s,
      tau = config$tau,
      q = ncol(dat$Z),
      error_dist = config$error_dist,
      design_type = config$design_type,
      random_effect_dist = config$random_effect_dist,
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
      warning_count = length(fit$warning_messages %||% character()),
      target_type = target_obj$target_type,
      implementation_commit = commit,
      stringsAsFactors = FALSE
    )

    if (!is.null(fit$beta_tilde) && !is.null(fit$se)) {
      coord_names <- intersect(names(fit$beta_tilde), names(beta_target))
      for (nm in coord_names) {
        se <- as.numeric(fit$se[nm])
        bt <- as.numeric(fit$beta_tilde[nm])
        lower <- as.numeric(fit$ci_lower[nm] %||% (bt - stats::qnorm(.975) * se))
        upper <- as.numeric(fit$ci_upper[nm] %||% (bt + stats::qnorm(.975) * se))
        ccount <- ccount + 1L
        coordinate_rows[[ccount]] <- data.frame(
          experiment_id = config$experiment_id,
          replicate = replicate,
          seed = seed,
          method_id = method_id,
          coordinate = nm,
          coordinate_type = if (match(nm, names(dat$beta0)) %in% dat$active) "active" else "null",
          target_type = target_obj$target_type,
          target_value = as.numeric(beta_target[nm]),
          target_mc_se = target_obj$target_mc_se,
          beta_hat = as.numeric(beta_full[nm]),
          beta_tilde = bt,
          estimated_se = se,
          ci_level = 0.95,
          ci_lower = lower,
          ci_upper = upper,
          covered = lower <= beta_target[nm] && upper >= beta_target[nm],
          wald_z = (bt - beta_target[nm]) / se,
          interval_length = upper - lower,
          status = fit$status,
          implementation_commit = commit,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    metrics = do.call(rbind, method_rows),
    coordinates = if (length(coordinate_rows)) do.call(rbind, coordinate_rows) else data.frame(),
    screening = prepared$screening,
    config = config,
    seed = seed
  )
}

# One-replication execution and replication-level result construction.

build_tuning_for_analysis_v2 <- function(config, dat) {
  values <- reference_tuning_values_v2(config, p = ncol(dat$X), n = dat$n_clusters)
  # Round-3 first-stage lambda: cluster self-normalised profile-score penalty
  # (METHOD_SPECIFICATION_ROUND3_AMENDMENT.md section 1). lambda_0,n is the
  # calibrated base; the coordinate penalty is lambda_0,n * ell_j_final
  # computed inside fit_profile_lasso_calibrated_v2. lambda_beta (the old
  # sqrt(log p / n) anchor) is retained for benchmark comparators whose
  # first stage is unchanged (e.g. SQR-DEBIASED-IID).
  lambda_meta <- lambda_0_n_v2(dat$n_clusters, ncol(dat$X))
  # Primary inferential bandwidth h = c_h n^{-3/10} (theory decision:
  # PILOT_GATE_THEORY_DECISIONS.md). Gives sqrt(n) h^2 -> 0 and n h^3 -> Inf.
  h <- config$h_multiplier * dat$n_clusters^(-3 / 10)
  list(
    h = h,
    lambda_beta = 1 * sqrt(log(max(ncol(dat$X), 2)) / dat$n_clusters),
    lambda_0_n = lambda_meta$lambda_0,
    lambda_alpha = lambda_meta$alpha,
    lambda_normal_quantile = lambda_meta$normal_quantile,
    lambda_safety_constant = lambda_meta$safety_constant,
    lambda_gamma = config$lambda_gamma,
    mu_grid = config$dantzig_multipliers *
      (sqrt(log(max(ncol(dat$X), 2)) / (dat$n_clusters * h)) + h^2),
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

# POP-H row-accuracy and inverse-defect diagnostics (theory decision section
# 3.5): E2_k, E1_k, cosine similarity, D_k, and the scaled Bahadur remainder.
# Diagnostics/gates only; never used to select mu.
pop_direction_diagnostics_v2 <- function(fit, population_asset, nm,
                                         delta_l1, n_clusters, tau,
                                         beta_star = NA_real_) {
  out <- list(E1 = NA_real_, E2 = NA_real_, cosine_pop = NA_real_,
              D_k = NA_real_, A_k = NA_real_, bahadur_scaled = NA_real_)
  if (is.null(population_asset) || is.null(population_asset$directions) ||
      !nm %in% names(population_asset$directions)) return(out)
  inf <- fit$inference_object$table
  rr <- which(inf$coordinate == nm)[1]
  if (!length(rr) || !isTRUE(inf$feasible[rr])) return(out)
  fit_names <- names(fit$beta_hat)
  pop_full <- setNames(rep(0, length(fit_names)), fit_names)
  pop <- population_asset$directions[[nm]]$omega
  common <- intersect(names(pop), fit_names)
  pop_full[common] <- pop[common]
  omega_hat <- as.numeric(fit$inference_object$directions[[rr]]$omega %||%
                            rep(NA_real_, length(fit_names)))
  if (length(omega_hat) != length(fit_names) || any(!is.finite(omega_hat))) return(out)
  diffv <- omega_hat - pop_full
  n_pop <- sqrt(sum(pop_full^2))
  out$E2 <- sqrt(sum(diffv^2)) / max(n_pop, 1e-12)
  out$E1 <- sum(abs(diffv)) / (1 + sum(abs(pop_full)))
  denom <- sqrt(sum(omega_hat^2)) * n_pop
  out$cosine_pop <- if (denom > 0) drop(crossprod(omega_hat, pop_full)) / denom else NA_real_
  sigma0 <- population_asset$sigma0_pop[[nm]] %||% NA_real_
  delta_k <- as.numeric(inf$dantzig_residual[rr])
  if (is.finite(sigma0) && sigma0 > 0 && is.finite(delta_k) && is.finite(delta_l1)) {
    out$D_k <- sqrt(n_clusters) * delta_k * delta_l1 / sigma0
  }
  # A_k: the ACTUAL normalized inverse-defect inner product (round-2, B2):
  #   A_k = sqrt(n) |(e_k - H omega_hat)' Delta| / sigma0_pop,
  # Delta = beta_hat - beta_star (full vector set by caller as fit$delta_full).
  Hf <- fit$fit_object$components$hessian
  delta_vec <- fit$delta_full
  if (is.finite(sigma0) && sigma0 > 0 && !is.null(Hf) && length(delta_vec) == length(fit_names) &&
      all(is.finite(delta_vec))) {
    e <- numeric(length(fit_names)); e[match(nm, fit_names)] <- 1
    rv <- as.numeric(e - Hf %*% omega_hat)
    out$A_k <- sqrt(n_clusters) * abs(sum(rv * delta_vec)) / sigma0
  }
  # Scaled Bahadur remainder using the oracle (population) direction:
  #   R_B,k = sqrt(n)(beta_tilde_k - beta*_k) - n^{-1/2} sum_i xi_{ik}^{oracle},
  #   xi_{ik}^{oracle} = -omega_k^{pop}' g_i^{(0)} (unsmoothed sample scores),
  # scaled by the population unsmoothed meat SD sigma0.
  bt <- as.numeric(inf$beta_tilde[rr])
  comp <- fit$fit_object$components
  G0 <- comp$unsmoothed_cluster_scores
  if (is.finite(bt) && is.finite(beta_star) && !is.null(G0) &&
      ncol(G0) == length(fit_names) && is.finite(sigma0) && sigma0 > 0) {
    proj <- as.numeric(G0 %*% pop_full)
    xi_sum <- sum(proj)  # sum_i omega_pop' g_i^(0); xi = -omega'g so -sum(xi) = +xi_sum
    rB <- sqrt(n_clusters) * (bt - beta_star) + (1 / sqrt(n_clusters)) * xi_sum
    out$bahadur_scaled <- rB / sigma0
  }
  out
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
  # POP-H row-accuracy / inverse-defect aggregates (theory decision 3.5).
  Dvals <- numeric(); E1vals <- numeric(); E2vals <- numeric()
  Bvals <- numeric(); Avals <- numeric()
  if (nrow(inf) && !is.null(population_asset$directions)) {
    delta_l1 <- if (!is.null(fit_obj$beta) && !is.null(target_obj$beta)) {
      fitn <- names(fit_obj$beta)
      tgt <- target_obj$beta
      common <- intersect(fitn, names(tgt))
      if (length(common)) {
        fit$delta_full <- setNames(rep(0, length(fitn)), fitn)
        fit$delta_full[common] <- fit_obj$beta[common] - tgt[common]
        sum(abs(fit$delta_full))
      } else NA_real_
    } else NA_real_
    for (nm in intersect(inf$coordinate, names(population_asset$directions))) {
      pd <- pop_direction_diagnostics_v2(
        fit, population_asset, nm, delta_l1, analysis_dat$n_clusters,
        config$tau, beta_star = as.numeric(target_obj$beta[nm] %||% NA_real_)
      )
      Dvals <- c(Dvals, pd$D_k); E1vals <- c(E1vals, pd$E1)
      E2vals <- c(E2vals, pd$E2); Bvals <- c(Bvals, pd$bahadur_scaled)
      Avals <- c(Avals, pd$A_k)
    }
  }
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
    max_E1_pop = if (length(E1vals)) max(E1vals[is.finite(E1vals)], na.rm = TRUE) else NA_real_,
    max_E2_pop = if (length(E2vals)) max(E2vals[is.finite(E2vals)], na.rm = TRUE) else NA_real_,
    max_D_k = if (length(Dvals)) max(Dvals[is.finite(Dvals)], na.rm = TRUE) else NA_real_,
    median_D_k = if (length(Dvals)) stats::median(Dvals[is.finite(Dvals)], na.rm = TRUE) else NA_real_,
    q90_D_k = if (length(Dvals)) stats::quantile(Dvals[is.finite(Dvals)], 0.9, na.rm = TRUE, names = FALSE) else NA_real_,
    max_A_k = if (length(Avals)) max(Avals[is.finite(Avals)], na.rm = TRUE) else NA_real_,
    median_A_k = if (length(Avals)) stats::median(Avals[is.finite(Avals)], na.rm = TRUE) else NA_real_,
    max_bahadur_scaled = if (length(Bvals)) max(abs(Bvals[is.finite(Bvals)]), na.rm = TRUE) else NA_real_,
    scaled_bahadur_remainder = if (length(Bvals)) max(Bvals[is.finite(Bvals)], na.rm = TRUE) else NA_real_,
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
          # CLARABEL remains the production Dantzig solver: on the dense
          # profile-Hessian LPs it solves a full mu-grid row in ~26s while
          # HIGHS grinds >30min on the same problem (infeasibility proofs).
          # HIGHS/ECOS/SCS stay available as audit solvers (round-2 parity
          # suite: HIGHS results match CLARABEL where it terminates).
          solver_preference = c("CLARABEL", "ECOS", "SCS"),
          fit_control = list(
            # Round-3 reference solver controls
            # (METHOD_SPECIFICATION_ROUND3_AMENDMENT.md section 3):
            # max_iter >= 2000, max_backtrack >= 50, beta_tol = 1e-7.
            max_iter = 2000L,
            max_backtrack = 50L,
            beta_tol = 1e-7,
            kkt_normalized_tol = 1e-3,
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
          B = if (isTRUE(final)) 100L else 25L,
          nK = if (isTRUE(final)) 15L else 7L
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
      lambda_0_n = tuning$lambda_0_n,
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
      lambda_rule = fit$lambda_rule %||% NA_character_,
      lambda_alpha = fit$lambda_alpha %||% NA_real_,
      lambda_normal_quantile = fit$lambda_normal_quantile %||% NA_real_,
      lambda_safety_constant = fit$lambda_safety_constant %||% NA_real_,
      lambda_base = fit$lambda_base %||% NA_real_,
      lambda_loading_pass0_min = fit$lambda_loading_pass0_min %||% NA_real_,
      lambda_loading_pass0_median = fit$lambda_loading_pass0_median %||% NA_real_,
      lambda_loading_pass0_max = fit$lambda_loading_pass0_max %||% NA_real_,
      lambda_loading_pass1_min = fit$lambda_loading_pass1_min %||% NA_real_,
      lambda_loading_pass1_median = fit$lambda_loading_pass1_median %||% NA_real_,
      lambda_loading_pass1_max = fit$lambda_loading_pass1_max %||% NA_real_,
      lambda_coordinate_min = fit$lambda_coordinate_min %||% NA_real_,
      lambda_coordinate_median = fit$lambda_coordinate_median %||% NA_real_,
      lambda_coordinate_max = fit$lambda_coordinate_max %||% NA_real_,
      zero_profile_score_max = fit$zero_profile_score_max %||% NA_real_,
      zero_kkt_ratio_max = fit$zero_kkt_ratio_max %||% NA_real_,
      preliminary_kkt_normalized = fit$preliminary_kkt_normalized %||% NA_real_,
      final_kkt_normalized = fit$final_kkt_normalized %||% NA_real_,
      final_kkt_absolute = fit$final_kkt_absolute %||% NA_real_,
      first_stage_iterations = fit$first_stage_iterations %||% NA_integer_,
      first_stage_beta_change = fit$first_stage_beta_change %||% NA_real_,
      first_stage_nonzero_count = fit$first_stage_nonzero_count %||% NA_integer_,
      target_type = target_obj$target_type,
      target_mc_se = target_obj$target_mc_se,
      implementation_commit = commit,
      requested_workers = config$workers,
      stringsAsFactors = FALSE
    )

    inf_tab <- inference_direction_table_v2(fit)
    if (nrow(inf_tab)) {
      if (all(is.finite(beta_full)) && all(is.finite(beta_target))) {
        fit$delta_full <- beta_full - beta_target
        delta_l1 <- sum(abs(fit$delta_full))
      } else {
        fit$delta_full <- NULL
        delta_l1 <- NA_real_
      }
      coord_names <- intersect(inf_tab$coordinate, names(beta_target))
      # Oracle variance ladder (round-2 A1): target unsmoothed scores at beta*
      # computed once per replication and shared across coordinates. Only for
      # reference-profile methods where the POP-H asset and exact rows exist.
      ladder_Gs <- NULL
      # Ladder is computed on the FULL analysis design (dat), so only for the
      # full-design reference methods; TRUE-SUPPORT operates on the oracle
      # sub-design and its target-scores need a separate sub-design pass.
      ladder_wanted <- method_id %in% c("PROFILE-DQR", "PROFILE-DQR-POP-H") &&
        !is.null(population_asset) && !is.null(population_asset$directions) &&
        length(coord_names) > 0 && exists("target_residual_scores_v2", mode = "function")
      if (ladder_wanted) {
        ladder_Gs <- tryCatch(
          target_residual_scores_v2(
            dat, beta_target, config$tau, tuning$lambda_gamma,
            h = tuning$h,
            nuisance_control = list(reltol = 1e-10, maxit = 300L, grad_tol = 1e-7)
          ),
          error = function(e) NULL
        )
      }
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
        popd <- pop_direction_diagnostics_v2(
          fit, population_asset, nm, delta_l1, dat$n_clusters, config$tau,
          beta_star = as.numeric(beta_target[nm])
        )
        ladd <- if (!is.null(ladder_Gs) && exists("variance_ladder_one_v2", mode = "function")) {
          tryCatch(
            variance_ladder_one_v2(
              dat, fit, beta_target, population_asset, nm, config$tau,
              tuning$lambda_gamma, Gs = ladder_Gs, h = tuning$h
            ),
            error = function(e) NULL
          )
        } else NULL
        ncl <- dat$n_clusters
        lse1 <- if (!is.null(ladd)) sqrt(max(ladd$meat1, 0) / ncl) else NA_real_
        lse2 <- if (!is.null(ladd)) sqrt(max(ladd$meat2, 0) / ncl) else NA_real_
        lse3 <- if (!is.null(ladd)) sqrt(max(ladd$meat3, 0) / ncl) else NA_real_
        lse4 <- if (!is.null(ladd)) sqrt(max(ladd$meat4, 0) / ncl) else NA_real_
        Lk <- if (!is.null(ladd)) ladd$L_k else NA_real_
        # Finite-sample sandwich diagnostics (round-2 A1): CR0/CR3/KC are
        # cheap per-coordinate; only meaningful where the exact row is close
        # (TRUE-SUPPORT / low-dim), but computed for all reference-profile rows.
        fsd <- if (exists("finite_sample_sandwich_diagnostics_v2", mode = "function") &&
                   !is.null(fit$fit_object$components$unsmoothed_cluster_scores) &&
                   isTRUE(inf_tab$feasible[rr])) {
          tryCatch(finite_sample_sandwich_diagnostics_v2(
            fit$fit_object$components$unsmoothed_cluster_scores,
            as.numeric(fit$inference_object$directions[[rr]]$omega %||% rep(NA_real_, ncol(dat$X)))
          ), error = function(e) NULL)
        } else NULL
        fs_cr0 <- if (!is.null(fsd)) fsd$se_cr0 else NA_real_
        fs_cr3 <- if (!is.null(fsd)) fsd$se_cr3 else NA_real_
        fs_kc <- if (!is.null(fsd)) fsd$se_kc else NA_real_
        fs_lev <- if (!is.null(fsd)) fsd$leverage_max else NA_real_
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
          estimated_se_smoothed = as.numeric(inf_tab$se_smoothed[rr] %||% NA_real_),
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
          adjacent_stability = inf_tab$adjacent_stability[rr] %||% NA_real_,
          mu_cv_min_defect = inf_tab$mu_cv_min_defect[rr] %||% NA_real_,
          mu_cv_defect = inf_tab$mu_cv_defect[rr] %||% NA_real_,
          mu_cv_quad = inf_tab$mu_cv_quad[rr] %||% NA_real_,
          E1_pop = popd$E1,
          E2_pop = popd$E2,
          cosine_pop = popd$cosine_pop,
          D_k = popd$D_k,
          A_k = popd$A_k,
          bahadur_scaled = popd$bahadur_scaled,
          ladder_L_k = Lk,
          ladder_se1 = lse1, ladder_se2 = lse2, ladder_se3 = lse3, ladder_se4 = lse4,
          fs_se_cr0 = fs_cr0, fs_se_cr3 = fs_cr3, fs_se_kc = fs_kc, fs_leverage_max = fs_lev,
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

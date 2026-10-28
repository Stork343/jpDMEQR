# One-replication execution and replication-level result construction.

build_tuning_for_analysis_v2 <- function(config, dat) {
  values <- reference_tuning_values_v2(config, p = ncol(dat$X), n = dat$n_clusters)
  # Round-4 dual bandwidths (METHOD_SPECIFICATION_ROUND4_AMENDMENT.md section 1):
  #   h_est = (log(p_P)/n)^{1/4}, c_est = 1 -- first-stage estimation only;
  #   h_inf = c_h n^{-3/10}                  -- inference only (unchanged).
  # p_P = number of penalised coordinates in the fitted design (all-X default;
  # the calibrated wrapper recomputes h_est from its own fit design).
  p_P <- ncol(dat$X)
  h_est <- (log(max(p_P, 2)) / dat$n_clusters)^(1 / 4)
  # Round-3 first-stage lambda: cluster self-normalised profile-score penalty
  # (METHOD_SPECIFICATION_ROUND3_AMENDMENT.md section 1). lambda_0,n is the
  # calibrated base; the coordinate penalty is lambda_0,n * ell_j_final
  # computed inside fit_profile_lasso_calibrated_v2. lambda_beta (the old
  # sqrt(log p / n) anchor) is retained for benchmark comparators whose
  # first stage is unchanged (e.g. SQR-DEBIASED-IID).
  lambda_meta <- lambda_0_n_v2(dat$n_clusters, ncol(dat$X))
  # Primary inferential bandwidth h_inf = c_h n^{-3/10} (theory decision:
  # PILOT_GATE_THEORY_DECISIONS.md). Gives sqrt(n) h^2 -> 0 and n h^3 -> Inf.
  h_inf <- config$h_multiplier * dat$n_clusters^(-3 / 10)
  list(
    method_version = config$method_version %||% "registered_profile_v6",
    mean_score_multiplier = as.numeric(config$mean_score_multiplier %||% 1),
    lambda_multiplier = if (length(config$lambda_beta_multipliers) == 1L)
      config$lambda_beta_multipliers else 1,
    h = h_inf,
    h_est = h_est,
    h_inf = h_inf,
    lambda_beta = 1 * sqrt(log(max(ncol(dat$X), 2)) / dat$n_clusters),
    lambda_0_n = lambda_meta$lambda_0,
    lambda_alpha = lambda_meta$alpha,
    lambda_normal_quantile = lambda_meta$normal_quantile,
    lambda_safety_constant = lambda_meta$safety_constant,
    lambda_gamma = config$lambda_gamma,
    mu_grid = config$dantzig_multipliers *
      (sqrt(log(max(ncol(dat$X), 2)) / (dat$n_clusters * h_inf)) + h_inf^2),
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

failed_paired_coordinate_rows_v2 <- function(run_id, config, replicate, seed, method_id,
                                             target_names, target_obj, full_names,
                                             active, implementation_commit) {
  if (!length(target_names) || !all(target_names %in% full_names) ||
      !all(target_names %in% names(target_obj$beta))) {
    stop("Failed paired coordinates must retain every requested target name.")
  }
  rows <- lapply(target_names, function(nm) {
    target_se <- if (nm %in% names(target_obj$target_mc_se_by_coordinate)) {
      target_obj$target_mc_se_by_coordinate[[nm]]
    } else target_obj$target_mc_se
    data.frame(
      run_id = run_id, experiment_id = config$experiment_id, module = config$module,
      replicate = replicate, seed = seed, method_id = method_id,
      method_version = "profile_pair_20261022",
      score_definition = "profile_bandwidth_pair_2_half_minus_full",
      coordinate = nm,
      coordinate_type = if (match(nm, full_names) %in% active) "active" else "null",
      target_type = target_obj$target_type, target_value = as.numeric(target_obj$beta[nm]),
      target_mc_se = as.numeric(target_se), beta_hat = NA_real_, beta_tilde = NA_real_,
      estimated_se = NA_real_, estimated_se_smoothed = NA_real_, ci_level = 0.95,
      ci_lower = NA_real_, ci_upper = NA_real_, covered = NA,
      ci_level_secondary = 0.90, ci_lower_secondary = NA_real_, ci_upper_secondary = NA_real_,
      covered_secondary = NA, interval_length_secondary = NA_real_,
      covered_target_minus_1.96mcse = NA, covered_target_plus_1.96mcse = NA,
      wald_z = NA_real_, p_value = NA_real_, interval_length = NA_real_, feasible = FALSE,
      dantzig_mu = NA_real_, mu_source = NA_character_, dantzig_residual = NA_real_,
      omega_l1 = NA_real_, omega_l2 = NA_real_, adjacent_stability = NA_real_,
      mu_cv_min_defect = NA_real_, mu_cv_defect = NA_real_, mu_cv_quad = NA_real_,
      mu_cv_folds_used = NA_integer_, E1_pop = NA_real_, E2_pop = NA_real_, cosine_pop = NA_real_,
      D_k = NA_real_, A_k = NA_real_, A_k_lasso_start = NA_real_,
      bahadur_scaled = NA_real_, bahadur_refit_scaled = NA_real_, ladder_L_k = NA_real_,
      ladder_se1 = NA_real_, ladder_se2 = NA_real_, ladder_se3 = NA_real_, ladder_se4 = NA_real_,
      fs_se_cr0 = NA_real_, fs_se_cr3 = NA_real_, fs_se_kc = NA_real_, fs_leverage_max = NA_real_,
      status = "failed", implementation_commit = implementation_commit,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

coordinate_wald_statistics_v2 <- function(estimate, se, target) {
  if (!is.finite(estimate) || !is.finite(se) || se <= 0 || !is.finite(target)) {
    return(list(wald_z = NA_real_, p_value = NA_real_))
  }
  list(wald_z = (estimate - target) / se,
       p_value = 2 * stats::pnorm(-abs(estimate / se)))
}

paired_population_score_centering_v2 <- function(asset, beta_target, h, config) {
  result <- list(matched = FALSE, reason = "population_asset_missing", mean = NULL)
  if (is.null(asset)) return(result)
  if (!identical(asset$score_definition, "profile_bandwidth_pair_2_half_minus_full") ||
      !identical(asset$design_identity, profile_pair_design_identity_v2(config))) {
    result$reason <- "population_score_design_mismatch"
    return(result)
  }
  nm <- names(beta_target)
  if (is.null(nm) || !setequal(nm, names(asset$beta_target)) ||
      !setequal(nm, names(asset$score_mean_population)) ||
      any(!is.finite(asset$score_mean_population))) {
    result$reason <- "population_score_mean_or_coordinates_missing"
    return(result)
  }
  if (!isTRUE(all.equal(as.numeric(asset$beta_target[nm]), as.numeric(beta_target),
                       tolerance = 0))) {
    result$reason <- "population_score_beta_mismatch"
    return(result)
  }
  if (!identical(as.numeric(asset$h_analysis), as.numeric(h))) {
    result$reason <- "population_score_bandwidth_mismatch"
    return(result)
  }
  list(matched = TRUE, reason = "matched_population_beta_bandwidth_design",
       mean = asset$score_mean_population[nm], beta_target = asset$beta_target[nm],
       h_analysis = asset$h_analysis, experiment_id = asset$experiment_id,
       construction_version = asset$construction_version)
}

coordinate_coverage_summary_v2 <- function(data, column = "covered") {
  covered <- data[[column]]
  failed <- data$status %in% c("failed", "not_implemented")
  primary <- covered
  primary[failed & is.na(primary)] <- FALSE
  conditional <- summarise_coverage_v2(covered, scheduled = nrow(data))
  output <- summarise_coverage_v2(primary, scheduled = nrow(data))
  output$evaluable <- sum(!is.na(covered))
  output$evaluable_fraction <- output$evaluable / nrow(data)
  output$failed_inference_coordinates <- sum(failed)
  output$coverage_denominator <- sum(!is.na(primary))
  output$coverage_successful <- conditional$coverage
  output$mcse_successful <- conditional$mcse
  output$binom_lower_successful <- conditional$binom_lower
  output$binom_upper_successful <- conditional$binom_upper
  output
}

coordinate_replicate_coverage_v2 <- function(data, column = "covered") {
  covered <- data[[column]]
  failed <- data$status %in% c("failed", "not_implemented")
  covered[failed & is.na(covered)] <- FALSE
  replicate_coverage_summary_v2(covered, data$replicate)
}

validate_paired_failed_coordinate_accounting_v2 <- function(metrics, coordinates) {
  failed <- metrics$status %in% c("failed", "not_implemented") &
    (metrics$score_definition %in% "profile_bandwidth_pair_2_half_minus_full" |
       startsWith(metrics$method_id, "JOINT-PROFILE") |
       metrics$method_id %in% c("SQR-DEBIASED-IID", "BIAS-ADJ-LQMM"))
  if (!any(failed)) return(invisible(TRUE))
  keys <- c("experiment_id", "replicate", "method_id")
  for (i in which(failed)) {
    keep <- if (nrow(coordinates)) {
      coordinates$experiment_id == metrics$experiment_id[i] &
        coordinates$replicate == metrics$replicate[i] &
        coordinates$method_id == metrics$method_id[i]
    } else logical()
    d <- coordinates[keep, , drop = FALSE]
    if (!nrow(d) || anyDuplicated(d$coordinate) || any(d$status != "failed") ||
        any(d$feasible) || any(!is.na(d$beta_tilde)) || any(!is.na(d$estimated_se)) ||
        any(!is.na(d$ci_lower)) || any(!is.na(d$ci_upper))) {
      stop("Failed paired inference requires explicit unavailable rows for every requested coordinate: ",
           paste(unlist(metrics[i, keys]), collapse = "/"))
    }
  }
  invisible(TRUE)
}

# POP-H row-accuracy and inverse-defect diagnostics (theory decision section
# 3.5): E2_k, E1_k, cosine similarity, D_k, and the scaled Bahadur remainder.
# Diagnostics/gates only; never used to select mu.
pop_direction_diagnostics_v2 <- function(fit, population_asset, nm,
                                         delta_l1, n_clusters, tau,
                                         beta_star = NA_real_,
                                         delta_override = NULL) {
  out <- list(E1 = NA_real_, E2 = NA_real_, cosine_pop = NA_real_,
              D_k = NA_real_, A_k = NA_real_, bahadur_scaled = NA_real_,
              score_skewness = NA_real_, score_excess_kurtosis = NA_real_)
  paired <- identical(fit$score_definition, "profile_bandwidth_pair_2_half_minus_full")
  if (identical(population_asset$score_definition,
                 "profile_bandwidth_pair_2_half_minus_full") &&
      (!paired || !fit$method_id %in% c("PROFILE-DQR", "PROFILE-DQR-POP-H",
                                       "POSTREFIT-EXACT-H"))) return(out)
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
  # Delta = beta_hat - beta_star (full vector set by caller as fit$delta_full;
  # round-5: delta_override lets the caller plug the L1-start delta for the
  # A_k_lasso_start diagnostic).
  Hf <- fit$fit_object$components$hessian
  delta_vec <- delta_override %||% fit$delta_full
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
  G0 <- if (paired) fit$target_cluster_scores else comp$unsmoothed_cluster_scores
  if (paired && is.null(G0)) stop("Paired Bahadur diagnostics require paired target scores.")
  if (is.finite(bt) && is.finite(beta_star) && !is.null(G0) &&
      ncol(G0) == length(fit_names) && is.finite(sigma0) && sigma0 > 0) {
    proj <- as.numeric(G0 %*% pop_full)
    centered <- proj - mean(proj)
    sd_proj <- stats::sd(proj)
    if (is.finite(sd_proj) && sd_proj > 0) {
      out$score_skewness <- mean(centered^3) / sd_proj^3
      out$score_excess_kurtosis <- mean(centered^4) / sd_proj^4 - 3
    }
    centering <- fit$population_score_centering
    if (!paired || isTRUE(centering$matched)) {
      population_mean <- if (paired) sum(centering$mean[fit_names] * pop_full) else 0
      xi_sum <- sum(proj - population_mean)
      rB <- sqrt(n_clusters) * (bt - beta_star) + (1 / sqrt(n_clusters)) * xi_sum
      out$bahadur_scaled <- rB / sigma0
    }
  }
  out
}

method_theory_diagnostics_v2 <- function(config, replicate, seed, method_id,
                                         fit, target_obj, population_asset,
                                         analysis_dat, run_id = NA_character_) {
  if (identical(population_asset$score_definition,
                 "profile_bandwidth_pair_2_half_minus_full") &&
      (!identical(fit$score_definition, population_asset$score_definition) ||
       !method_id %in% c("PROFILE-DQR", "PROFILE-DQR-POP-H", "POSTREFIT-EXACT-H"))) {
    population_asset <- NULL
  }
  fit_obj <- fit$fit_object %||% NULL
  H <- if (!is.null(fit_obj$components$hessian)) fit_obj$components$hessian else NULL
  eig <- if (!is.null(H) && all(is.finite(H))) {
    eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  } else numeric()
  hessian_max_error <- NA_real_
  hessian_operator_error <- NA_real_
  precision_row_error <- NA_real_
  if (!is.null(H) && !is.null(population_asset$H_population)) {
    pop_names <- colnames(population_asset$H_population)
    fit_names <- names(fit_obj$beta)
    common <- intersect(pop_names, fit_names)
    if (length(common)) {
      difference <- H[match(common, fit_names), match(common, fit_names), drop = FALSE] -
        population_asset$H_population[common, common, drop = FALSE]
      hessian_max_error <- max(abs(difference))
      hessian_operator_error <- max(abs(
        eigen((difference + t(difference)) / 2, symmetric = TRUE,
              only.values = TRUE)$values
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
  correction <- if (nrow(inf)) max_or_na_v2(abs(inf$beta_tilde - inf$beta_hat)) else NA_real_
  # POP-H row-accuracy / inverse-defect aggregates (theory decision 3.5).
  Dvals <- numeric(); E1vals <- numeric(); E2vals <- numeric()
  Bvals <- numeric(); Avals <- numeric(); Svals <- numeric(); Kvals <- numeric()
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
      Svals <- c(Svals, pd$score_skewness)
      Kvals <- c(Kvals, pd$score_excess_kurtosis)
    }
  }
  data.frame(
    run_id = run_id,
    experiment_id = config$experiment_id,
    replicate = replicate,
    seed = seed,
    method_id = method_id,
    method_version = fit$method_version %||% fit$implementation_version %||% NA_character_,
    score_definition = fit$score_definition %||% NA_character_,
    target_type = target_obj$target_type,
    target_mc_se = target_obj$target_mc_se,
    hessian_max_error = hessian_max_error,
    hessian_operator_error = hessian_operator_error,
    hessian_min_eigenvalue = if (length(eig)) min(eig) else NA_real_,
    hessian_condition_number = if (length(eig) && min(abs(eig)) > 0) {
      max(abs(eig)) / min(abs(eig))
    } else NA_real_,
    beta_h_target_distance = target_obj$smoothing_target_distance %||% NA_real_,
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
    max_E1_pop = max_or_na_v2(E1vals),
    max_E2_pop = max_or_na_v2(E2vals),
    max_D_k = max_or_na_v2(Dvals),
    median_D_k = if (length(Dvals)) stats::median(Dvals[is.finite(Dvals)], na.rm = TRUE) else NA_real_,
    q90_D_k = if (length(Dvals)) stats::quantile(Dvals[is.finite(Dvals)], 0.9, na.rm = TRUE, names = FALSE) else NA_real_,
    max_A_k = max_or_na_v2(Avals),
    median_A_k = if (length(Avals)) stats::median(Avals[is.finite(Avals)], na.rm = TRUE) else NA_real_,
    max_bahadur_scaled = max_or_na_v2(abs(Bvals)),
    scaled_bahadur_remainder = max_or_na_v2(Bvals),
    score_skewness = mean_or_na_v2(Svals),
    score_excess_kurtosis = mean_or_na_v2(Kvals),
    n_clusters = analysis_dat$n_clusters,
    p = ncol(analysis_dat$X),
    stringsAsFactors = FALSE
  )
}

run_one_replication_v2 <- function(root, config, replicate,
                                   base_seed = 20260817L,
                                   allow_missing_benchmarks = TRUE,
                                   final = FALSE,
                                   seed_experiment = NULL,
                                   run_id = NA_character_, implementation_commit = NULL) {
  # seed_experiment overrides the experiment_id used in seed derivation. Used
  # by the round-4 paired probe to nest P06 clusters inside P05 (common random
  # numbers, unchanged marginal DGPs).
  seed_exp <- seed_experiment %||% config$seed_family %||% config$experiment_id
  seed <- seed_from_id_v2(base_seed, seed_exp, replicate, "training")
  commit <- implementation_commit %||% current_commit_v2(root)
  dat0 <- simulate_from_config_v2(config, seed)
  prepared <- prepare_analysis_data_v2(dat0, config, seed)
  dat <- prepared$data
  target_obj <- load_profile_target_v2(root, config, dat, final = final)
  beta_target <- target_obj$beta
  joint_execution <- identical(config$method_version, "joint_profile_design_20261026")
  evaluation_data <- if (joint_execution) dat0 else dat
  target_names <- intersect(resolve_target_coordinates_v2(config, dat), names(beta_target))
  if (!length(target_names)) {
    target_names <- head(names(beta_target), min(2L, length(beta_target)))
  }
  tuning <- build_tuning_for_analysis_v2(config, dat)
  # Round-4 first-stage RSC diagnostics context (shared across methods):
  # target components at h_est computed once per replication; per-cell
  # population restricted curvature from a diagnostic asset (truth-permitted,
  # never used in tuning).
  rsc_ctx <- NULL
  if (any(c("PROFILE-DQR", "PROFILE-DQR-POP-H") %in% config$methods) &&
      exists("first_stage_rsc_diagnostics_v2", mode = "function")) {
    rsc_ctx <- list(
      target_components = tryCatch(
        profile_components_v2(
          dat$y, dat$X, dat$Z, dat$cluster_id, beta_target, config$tau,
          tuning$h_est, tuning$lambda_gamma, need_hessian = TRUE,
          nuisance_control = list(reltol = 1e-11, maxit = 400L, grad_tol = 1e-8)
        ),
        error = function(e) NULL
      ),
      pop_restricted_ev = {
        pop_path <- file.path(root, "results", "preflight",
                              paste0("rsc_pop_curvature_", config$experiment_id, ".rds"))
        if (file.exists(pop_path)) {
          tryCatch(readRDS(pop_path)$lambda_min, error = function(e) NA_real_)
        } else NA_real_
      }
    )
  }
  need_pop_h <- "PROFILE-DQR-POP-H" %in% config$methods
  population_asset <- load_population_direction_asset_v2(
    root, config, required = need_pop_h, final = final
  )

  test_seed <- seed_from_id_v2(base_seed, seed_exp, replicate, "test")
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
  method_diagnostics <- list()
  mcount <- 0L
  ccount <- 0L
  tcount <- 0L

  # Lever 2 缓存：重复级环境，校准第一阶段请求相同的方法共享同一份拟合。
  # 调度器只在 JPDMEQR_SHARED_FIRST_STAGE 打开时把它注入 control。
  shared_first_stage <- new.env(parent = emptyenv())

  # 预测基线：训练样本 tau 分位数作为截距型分位数预测，与合并分位数回归基线
  # 一起构成相对技巧的配对基线。
  intercept_quantile <- as.numeric(stats::quantile(dat$y, probs = config$tau, names = FALSE))
  intercept_loss <- pinball_loss_v2(
    test_dat$y, rep(intercept_quantile, length(test_dat$y)), config$tau
  )
  prediction_data <- if (joint_execution) test_dat0 else test_dat
  if (joint_execution) intercept_loss <- submission_pinball(prediction_data$y,
    rep(intercept_quantile, length(prediction_data$y)), prediction_data$cluster_id, config$tau)

  for (method_id in config$methods) {
    inference_target_obj <- if (joint_execution) submission_interval_target(method_id, config, target_obj) else target_obj
    method_train <- dat
    method_targets <- target_names
    context <- list(
      active = dat$active,
      raw_training = dat0,
      screening = prepared$screening,
      population_direction_asset = population_asset,
      config = config,
      beta_target = beta_target,
      rsc_ctx = rsc_ctx,
      shared_first_stage = shared_first_stage
    )
    if (method_requires_low_dimensional_subset_v2(method_id) ||
        (joint_execution && method_id == "DOUBLE-PEN-QLMM")) {
      feature_cap <- 20L
      low <- common_lqmm_subset_v2(
        dat, target_names, config$tau,
        max_features = feature_cap,
        active = if (method_id == "ORACLE-LQMM") dat$active else NULL,
        seed = seed_from_id_v2(base_seed, config$experiment_id, replicate, "common-lowdim-roles"),
        screen_fraction = config$screen_fraction
      )
      method_train <- low$data
      method_targets <- intersect(target_names, colnames(method_train$X))
      context$common_screen_keep <- low$keep
      context$common_screen_roles <- low$roles
    }

    mem_before <- memory_snapshot_v2()
    start <- proc.time()[[3]]
    native_warning_records <- list()
    fit <- tryCatch(
      withCallingHandlers(fit_benchmark_frozen_tuning_v2(
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
          final_assets = final,
          implementation_commit = commit,
          solver_preference = c("CLARABEL", "ECOS", "SCS"),
          fit_control = list(
            # Round-3 reference solver controls
            # (METHOD_SPECIFICATION_ROUND3_AMENDMENT.md section 3):
            # max_iter >= 2000, max_backtrack >= 50, beta_tol = 1e-7.
            optimizer = if (identical(config$method_version, "joint_profile_design_20261026")) "matrix_free_fista" else "reference",
            max_iter = if (identical(config$method_version, "joint_profile_design_20261026")) 10000L else 2000L,
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
          maxit = if (identical(config$method_version, "joint_profile_design_20261026")) 1000L else 100L,
          screen_fraction = config$screen_fraction,
          screen_dim = config$screen_dim,
          # 倍率配置使用 registry 中的完整数值；其余适配器控制按运行类型设定。
          lambda = if (isTRUE(final)) NULL else 0.5,
          lambda_beta_grid = config$lambda_beta_multipliers,
          lambda_alpha_grid = c(0.5, 1),
          B = if (isTRUE(final)) 100L else 25L,
          nK = if (isTRUE(final)) 15L else 7L
        )
      ), warning = function(condition) {
        native_warning_records[[length(native_warning_records) + 1L]] <<- list(
          sequence = length(native_warning_records) + 1L,
          message = conditionMessage(condition),
          call = paste(deparse(conditionCall(condition)), collapse = " "),
          condition_class = class(condition),
          stage = "benchmark_fitting_and_tuning"
        )
      }),
      error = function(e) {
        paired_method <- identical(tuning$method_version, "profile_pair_20261022") &&
          method_id %in% c("PROFILE-DQR", "PROFILE-DQR-TRUE-SUPPORT",
                            "PROFILE-DQR-TRUE-NUISANCE", "PROFILE-DQR-POP-H",
                            "POSTREFIT-EXACT-H", "PROFILE-DQR-SPLIT")
        failed <- benchmark_failure_v2(
          method_id, e$failure_stage %||% "unknown", conditionMessage(e),
          proc.time()[[3]] - start
        )
        failed$method_version <- e$method_version %||% if (
          identical(tuning$method_version, "profile_pair_20261022") && !paired_method
        ) NA_character_ else tuning$method_version
        failed$score_definition <- e$score_definition %||% if (paired_method) {
          "profile_bandwidth_pair_2_half_minus_full"
        } else NA_character_
        failed$fit_object <- e$fit
        failed$failure_diagnostics <- e$diagnostics
        failed
      }
    )
    fit$native_warning_records <- native_warning_records
    native_messages <- vapply(native_warning_records, function(record) record$message, character(1))
    fit$warning_messages <- c(fit$warning_messages %||% character(), native_messages)
    if (length(native_warning_records) && identical(fit$status, "ok")) fit$status <- "warning"
    mem_after_fit <- memory_snapshot_v2()

    if (identical(fit$status, "not_implemented") && !allow_missing_benchmarks) {
      stop("Required benchmark is not implemented: ", method_id,
           " -- ", fit$failure_message)
    }

    beta_full <- if (!is.null(fit$beta_hat)) {
      map_beta_to_full_v2(fit$beta_hat, names(evaluation_data$beta0))
    } else {
      setNames(rep(NA_real_, length(evaluation_data$beta0)), names(evaluation_data$beta0))
    }
    est <- if (all(is.finite(beta_full))) {
      estimation_metrics_v2(beta_full, beta_target, evaluation_data$active)
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
    if (joint_execution) selected <- if (method_id == "JOINT-PROFILE-TRUE-SUPPORT") evaluation_data$active else
      match(colnames(dat$X)[selected], colnames(evaluation_data$X))
    sel <- selection_metrics_v2(selected, evaluation_data$active, length(evaluation_data$beta0))
    pred <- if (all(is.finite(beta_full))) {
      as.numeric(prediction_data$X %*% beta_full) + (fit$prediction_intercept %||% 0)
    } else rep(NA_real_, length(prediction_data$y))
    loss <- if (all(is.finite(pred))) {
      if (joint_execution) submission_pinball(prediction_data$y, pred, prediction_data$cluster_id, config$tau) else
        pinball_loss_v2(test_dat$y, pred, config$tau)
    } else NA_real_
    cluster_loss <- if (all(is.finite(pred))) {
      cluster_loss_quantiles_v2(prediction_data$y, pred, prediction_data$cluster_id, config$tau)
    } else c(q10 = NA_real_, q50 = NA_real_, q90 = NA_real_)

    mcount <- mcount + 1L
    method_rows[[mcount]] <- data.frame(
      run_id = run_id,
      experiment_id = config$experiment_id,
      module = config$module,
      replicate = replicate,
      seed = seed,
      method_id = method_id,
      method_version = fit$method_version %||% fit$implementation_version %||% NA_character_,
      score_definition = fit$score_definition %||% NA_character_,
      status = fit$status,
      failure_stage = fit$failure_stage %||% "",
      failure_message = fit$failure_message %||% "",
      warning_messages = paste(fit$warning_messages %||% character(), collapse = " | "),
      reference_identifier = fit$reference_identifier %||% NA_character_,
      adapter_fidelity_status = fit$adapter_fidelity_status %||% NA_character_,
      n_clusters = method_train$n_clusters,
      n_clusters_original = config$n_clusters,
      total_rows = method_train$total_rows,
      p = fit$working_dimension %||% ncol(method_train$X),
      p_original = config$p,
      s = config$s,
      tau = config$tau,
      q = ncol(dat$Z),
      m_rule = config$m_rule,
      m_min = min(table(method_train$cluster_id)),
      m_mean = mean(table(method_train$cluster_id)),
      m_median = stats::median(as.numeric(table(method_train$cluster_id))),
      m_max = max(table(method_train$cluster_id)),
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
      h = fit$tuning_values$h %||% tuning$h,
      lambda_beta = fit$tuning_values$lambda_beta %||% tuning$lambda_beta,
      lambda_0_n = fit$lambda_base %||% tuning$lambda_0_n,
      lambda_stochastic = fit$lambda_stochastic %||% NA_real_,
      lambda_mean = fit$lambda_mean %||% NA_real_,
      lambda_mean_multiplier = fit$lambda_mean_multiplier %||% NA_real_,
      lambda_gamma = tuning$lambda_gamma,
      selected_size = sel$selected_size,
      selected_names = {
        keep_selected <- selected[selected > 0L]
        if (length(keep_selected)) {
          paste(names(evaluation_data$beta0)[keep_selected], collapse = "|")
        } else ""
      },
      tpr = sel$tpr,
      fdp = sel$fdp,
      fpr = sel$fpr,
      exact_support = sel$exact_support,
      l1_error = est$l1_error,
      l2_error = est$l2_error,
      max_active_error = est$max_active_error,
      pinball_loss = loss,
      pinball_loss_intercept = intercept_loss,
      pinball_skill_intercept = if (is.finite(loss) && is.finite(intercept_loss) &&
                                    intercept_loss > 0) {
        1 - loss / intercept_loss
      } else NA_real_,
      pinball_skill_pooled = NA_real_,
      cluster_loss_q10 = unname(cluster_loss[["q10"]]),
      cluster_loss_q50 = unname(cluster_loss[["q50"]]),
      cluster_loss_q90 = unname(cluster_loss[["q90"]]),
      quantile_calibration = if (all(is.finite(pred))) {
        if (joint_execution) mean(tapply(prediction_data$y <= pred, prediction_data$cluster_id, mean)) else mean(test_dat$y <= pred)
      } else NA_real_,
      objective = fit$profile_objective %||% fit$objective %||% NA_real_,
      runtime_sec = fit$runtime_sec %||% (proc.time()[[3]] - start),
      first_stage_sec = fit$first_stage_sec %||% NA_real_,
      refit_sec = fit$refit_sec %||% NA_real_,
      inference_sec = fit$inference_sec %||% NA_real_,
      component_timing_available = all(is.finite(c(
        fit$first_stage_sec %||% NA_real_,
        fit$refit_sec %||% NA_real_,
        fit$inference_sec %||% NA_real_
      ))),
      peak_memory_mb = max(mem_before[["rss_mb"]], mem_after_fit[["rss_mb"]],
                           na.rm = TRUE),
      process_peak_wset_mb = mem_after_fit[["peak_wset_mb"]],
      process_peak_rss_mb = mem_after_fit[["peak_rss_mb"]],
      outer_iterations = {
        iteration_counts <- c(fit$first_stage_iterations %||% NA_integer_,
                              fit$post_refit_iterations %||% NA_integer_)
        if (all(is.na(iteration_counts))) {
          NA_integer_
        } else {
          sum(iteration_counts, na.rm = TRUE)
        }
      },
      nuisance_iterations = fit$nuisance_iterations %||% NA_integer_,
      dantzig_feasible_fraction = NA_real_,
      warning_count = length(fit$warning_messages %||% character()),
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
      # Round-5 post-refit audit fields (METHOD_SPECIFICATION_ROUND5_AMENDMENT.md
      # section 10); l1/l2 errors are simulation-only.
      selected_support_size = fit$selected_support_size %||% NA_integer_,
      refit_set_size = fit$refit_set_size %||% NA_integer_,
      refit_contains_targets = fit$refit_contains_targets %||% NA,
      post_refit_status = fit$post_refit_status %||% NA_character_,
      post_refit_iterations = fit$post_refit_iterations %||% NA_integer_,
      post_refit_gradient_max = fit$post_refit_gradient_max %||% NA_real_,
      post_refit_nuisance_gradient_max = fit$post_refit_nuisance_gradient_max %||% NA_real_,
      post_refit_beta_change = fit$post_refit_beta_change %||% NA_real_,
      post_refit_hessian_min_eigenvalue = fit$post_refit_hessian_min_eigenvalue %||% NA_real_,
      post_refit_hessian_condition_number = fit$post_refit_hessian_condition_number %||% NA_real_,
      l1_error_lasso = if (!is.null(fit$fit_object$beta_l1) && all(is.finite(beta_full))) {
        l1v <- fit$fit_object$beta_l1
        common <- intersect(names(l1v), names(beta_target))
        if (length(common)) sum(abs(l1v[common] - beta_target[common])) else NA_real_
      } else NA_real_,
      l1_error_refit = if (!is.null(fit$fit_object$beta_refit) && all(is.finite(beta_full))) {
        l1v <- fit$fit_object$beta_refit
        common <- intersect(names(l1v), names(beta_target))
        if (length(common)) sum(abs(l1v[common] - beta_target[common])) else NA_real_
      } else NA_real_,
      l2_error_lasso = if (!is.null(fit$fit_object$beta_l1) && all(is.finite(beta_full))) {
        l1v <- fit$fit_object$beta_l1
        common <- intersect(names(l1v), names(beta_target))
        if (length(common)) sqrt(sum((l1v[common] - beta_target[common])^2)) else NA_real_
      } else NA_real_,
      l2_error_refit = if (!is.null(fit$fit_object$beta_refit) && all(is.finite(beta_full))) {
        l1v <- fit$fit_object$beta_refit
        common <- intersect(names(l1v), names(beta_target))
        if (length(common)) sqrt(sum((l1v[common] - beta_target[common])^2)) else NA_real_
      } else NA_real_,
      omitted_active_count = if (!is.null(fit$fit_object$refit_set)) {
        sum(!(dat$active %in% fit$fit_object$refit_set))
      } else NA_integer_,
      target_type = target_obj$target_type,
      target_mc_se = target_obj$target_mc_se,
      implementation_commit = commit,
      requested_workers = config$workers,
      stringsAsFactors = FALSE
    )

    # Round-4 first-stage RSC diagnostics (simulation-only; never tune).
    rsc_cols <- if (method_id %in% c("PROFILE-DQR", "PROFILE-DQR-POP-H") &&
                    !is.null(context$rsc_ctx) &&
                    isTRUE(fit$fit_object$first_stage_calibrated %||% FALSE) &&
                    exists("first_stage_rsc_diagnostics_v2", mode = "function")) {
      tryCatch(
        first_stage_rsc_diagnostics_v2(
          fit$fit_object, dat, beta_target, tuning$h_est, tuning$h_inf,
          rsc_ss = dat$active,
          tgt_components = context$rsc_ctx$target_components,
          pop_restricted_ev = context$rsc_ctx$pop_restricted_ev,
          cone_seed = seed
        ),
        error = function(e) NULL
      )
    } else NULL
    if (is.null(rsc_cols)) {
      rsc_cols <- data.frame(
        h_est = NA_real_, h_inf = NA_real_, rsc_order_lower = NA_real_,
        h_est_over_rsc_order = NA_real_, lambda_min_Hest_SS_population = NA_real_,
        lambda_min_Hest_SS_target_sample = NA_real_,
        lambda_min_Hest_SS_fitted_sample = NA_real_,
        condition_Hest_SS_population = NA_real_,
        condition_Hest_SS_target_sample = NA_real_,
        condition_Hest_SS_fitted_sample = NA_real_,
        cone_curvature_proxy_min = NA_real_
      )
    }
    method_rows[[mcount]] <- cbind(method_rows[[mcount]], rsc_cols)

    inf_tab <- inference_direction_table_v2(fit)
    failed_paired <- identical(fit$score_definition, "profile_bandwidth_pair_2_half_minus_full") &&
      !fit$status %in% c("ok", "warning")
    failed_source_interval <- (startsWith(method_id, "JOINT-PROFILE") || method_id %in% c(
      "SQR-DEBIASED-IID", "BIAS-ADJ-LQMM"
    )) && !fit$status %in% c("ok", "warning")
    if (failed_paired || failed_source_interval) {
      missing <- failed_paired_coordinate_rows_v2(
        run_id, config, replicate, seed, method_id, method_targets, inference_target_obj,
        names(evaluation_data$beta0), evaluation_data$active, commit
      )
      if (failed_source_interval) {
        missing$method_version <- fit$method_version %||% fit$implementation_version %||% NA_character_
        missing$score_definition <- fit$score_definition %||% NA_character_
      }
      for (j in seq_len(nrow(missing))) {
        ccount <- ccount + 1L
        coordinate_rows[[ccount]] <- missing[j, , drop = FALSE]
      }
    } else if (nrow(inf_tab)) {
      if (all(is.finite(beta_full)) && all(is.finite(beta_target))) {
        fit$delta_full <- beta_full - beta_target
        delta_l1 <- sum(abs(fit$delta_full))
      } else {
        fit$delta_full <- NULL
        delta_l1 <- NA_real_
      }
      coord_names <- intersect(inf_tab$coordinate, names(beta_target))
      if (identical(fit$score_definition, "profile_bandwidth_pair_2_half_minus_full")) {
        fit$population_score_centering <- paired_population_score_centering_v2(
          population_asset, beta_target, tuning$h, config
        )
      }
      # Oracle variance ladder (round-2 A1): target unsmoothed scores at beta*
      # computed once per replication and shared across coordinates. Only for
      # reference-profile methods where the POP-H asset and exact rows exist.
      ladder_Gs <- NULL
      # Ladder is computed on the FULL analysis design (dat), so only for the
      # full-design reference methods; TRUE-SUPPORT operates on the oracle
      # sub-design and its target-scores need a separate sub-design pass.
      ladder_wanted <- method_id %in% c("PROFILE-DQR", "PROFILE-DQR-POP-H",
                                         "POSTREFIT-EXACT-H") &&
        !is.null(population_asset) && !is.null(population_asset$directions) &&
        length(coord_names) > 0 && exists("target_residual_scores_v2", mode = "function")
      if (ladder_wanted && identical(fit$score_definition,
                                      "profile_bandwidth_pair_2_half_minus_full")) {
        ladder_Gs <- target_residual_scores_v2(
          dat, beta_target, config$tau, tuning$lambda_gamma, h = tuning$h,
          nuisance_control = list(reltol = 1e-10, maxit = 300L, grad_tol = 1e-7),
          score_definition = fit$score_definition
        )
        fit$target_cluster_scores <- ladder_Gs
      } else if (ladder_wanted) {
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
        coordinate_target <- as.numeric(inference_target_obj$beta[nm])
        target_se <- if (!is.null(inference_target_obj$target_mc_se_by_coordinate) &&
                         nm %in% names(inference_target_obj$target_mc_se_by_coordinate)) {
          inference_target_obj$target_mc_se_by_coordinate[[nm]]
        } else inference_target_obj$target_mc_se
        popd <- pop_direction_diagnostics_v2(
          fit, population_asset, nm, delta_l1, dat$n_clusters, config$tau,
          beta_star = as.numeric(beta_target[nm])
        )
        # Round-5 A_k from the L1 starting estimator (selection-start delta);
        # the primary A_k above uses the post-refit start.
        popd_l1 <- NULL
        if (!is.null(fit$fit_object$beta_l1)) {
          l1v <- fit$fit_object$beta_l1
          full_l1 <- setNames(rep(0, length(beta_full)), names(beta_full))
          common <- intersect(names(l1v), names(beta_full))
          full_l1[common] <- l1v[common]
          dl1 <- full_l1 - beta_target
          popd_l1 <- tryCatch(pop_direction_diagnostics_v2(
            fit, population_asset, nm, sum(abs(dl1)), dat$n_clusters, config$tau,
            beta_star = as.numeric(beta_target[nm]), delta_override = dl1
          ), error = function(e) NULL)
        }
        A_k_lasso <- if (!is.null(popd_l1)) popd_l1$A_k else NA_real_
        ladd <- if (!is.null(ladder_Gs) && exists("variance_ladder_one_v2", mode = "function")) {
          tryCatch(
            variance_ladder_one_v2(
              dat, fit, beta_target, population_asset, nm, config$tau,
              tuning$lambda_gamma, Gs = if (
                identical(fit$score_definition, "profile_bandwidth_pair_2_half_minus_full") &&
                isTRUE(fit$population_score_centering$matched)
              ) sweep(ladder_Gs, 2, fit$population_score_centering$mean[colnames(ladder_Gs)]) else ladder_Gs,
              h = tuning$h
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
        if (identical(fit$score_definition, "profile_bandwidth_pair_2_half_minus_full") &&
            !isTRUE(fit$population_score_centering$matched)) Lk <- NA_real_
        # Finite-sample sandwich diagnostics (round-2 A1): CR0/CR3/KC are
        # cheap per-coordinate; only meaningful where the exact row is close
        # (TRUE-SUPPORT / low-dim), but computed for all reference-profile rows.
        fsd <- if (exists("finite_sample_sandwich_diagnostics_v2", mode = "function") &&
                   !is.null(inferential_cluster_scores_v2(fit$fit_object$components)) &&
                   isTRUE(inf_tab$feasible[rr])) {
          tryCatch(finite_sample_sandwich_diagnostics_v2(
            inferential_cluster_scores_v2(fit$fit_object$components),
            as.numeric(fit$inference_object$directions[[rr]]$omega %||% rep(NA_real_, ncol(dat$X)))
          ), error = function(e) NULL)
        } else NULL
        fs_cr0 <- if (!is.null(fsd)) fsd$se_cr0 else NA_real_
        fs_cr3 <- if (!is.null(fsd)) fsd$se_cr3 else NA_real_
        fs_kc <- if (!is.null(fsd)) fsd$se_kc else NA_real_
        fs_lev <- if (!is.null(fsd)) fsd$leverage_max else NA_real_
        bounds_secondary <- interval_bounds_v2(bt, se, 0.90)
        wald <- coordinate_wald_statistics_v2(bt, se, coordinate_target)
        ccount <- ccount + 1L
        coordinate_rows[[ccount]] <- data.frame(
          run_id = run_id,
          experiment_id = config$experiment_id,
          module = config$module,
          replicate = replicate,
          seed = seed,
          method_id = method_id,
          method_version = fit$method_version %||% fit$implementation_version %||% NA_character_,
          score_definition = fit$score_definition %||% NA_character_,
          coordinate = nm,
          coordinate_type = if (match(nm, names(evaluation_data$beta0)) %in% evaluation_data$active) "active" else "null",
          target_type = inference_target_obj$target_type,
          target_value = coordinate_target,
          target_mc_se = as.numeric(target_se),
          beta_hat = as.numeric(beta_full[nm]),
          beta_tilde = bt,
          estimated_se = se,
          estimated_se_smoothed = as.numeric(inf_tab$se_smoothed[rr] %||% NA_real_),
          ci_level = 0.95,
          ci_lower = lower,
          ci_upper = upper,
          covered = coverage_indicator_v2(lower, upper, coordinate_target),
          ci_level_secondary = 0.90,
          ci_lower_secondary = unname(bounds_secondary[["lower"]]),
          ci_upper_secondary = unname(bounds_secondary[["upper"]]),
          covered_secondary = coverage_indicator_v2(
            unname(bounds_secondary[["lower"]]), unname(bounds_secondary[["upper"]]),
            coordinate_target
          ),
          interval_length_secondary = unname(
            bounds_secondary[["upper"]] - bounds_secondary[["lower"]]
          ),
          covered_target_minus_1.96mcse = coverage_indicator_v2(
            lower, upper, coordinate_target - 1.96 * target_se
          ),
          covered_target_plus_1.96mcse = coverage_indicator_v2(
            lower, upper, coordinate_target + 1.96 * target_se
          ),
          wald_z = wald$wald_z,
          p_value = wald$p_value,
          interval_length = upper - lower,
          feasible = inf_tab$feasible[rr] %||% is.finite(se),
          dantzig_mu = inf_tab$mu[rr] %||% NA_real_,
          mu_source = inf_tab$mu_source[rr] %||% NA_character_,
          dantzig_residual = inf_tab$dantzig_residual[rr] %||% NA_real_,
          omega_l1 = inf_tab$omega_l1[rr] %||% NA_real_,
          omega_l2 = inf_tab$omega_l2[rr] %||% NA_real_,
          adjacent_stability = inf_tab$adjacent_stability[rr] %||% NA_real_,
          mu_cv_min_defect = inf_tab$mu_cv_min_defect[rr] %||% NA_real_,
          mu_cv_defect = inf_tab$mu_cv_defect[rr] %||% NA_real_,
          mu_cv_quad = inf_tab$mu_cv_quad[rr] %||% NA_real_,
          mu_cv_folds_used = inf_tab$mu_cv_folds_used[rr] %||% NA_integer_,
          E1_pop = popd$E1,
          E2_pop = popd$E2,
          cosine_pop = popd$cosine_pop,
          D_k = popd$D_k,
          A_k = popd$A_k,
          A_k_lasso_start = A_k_lasso,
          bahadur_scaled = popd$bahadur_scaled,
          bahadur_refit_scaled = popd$bahadur_scaled,
          ladder_L_k = Lk,
          ladder_se1 = lse1, ladder_se2 = lse2, ladder_se3 = lse3, ladder_se4 = lse4,
          fs_se_cr0 = fs_cr0, fs_se_cr3 = fs_cr3, fs_se_kc = fs_kc, fs_leverage_max = fs_lev,
          status = fit$status,
          implementation_commit = commit,
          stringsAsFactors = FALSE
        )
      }
    }

    mem_end <- memory_snapshot_v2()
    memory_samples <- c(mem_before[["rss_mb"]], mem_after_fit[["rss_mb"]],
                        mem_end[["rss_mb"]])
    memory_samples <- memory_samples[is.finite(memory_samples)]
    method_rows[[mcount]]$peak_memory_mb <- if (length(memory_samples)) {
      max(memory_samples)
    } else NA_real_
    method_rows[[mcount]]$process_peak_wset_mb <- mem_end[["peak_wset_mb"]]
    method_rows[[mcount]]$process_peak_rss_mb <- mem_end[["peak_rss_mb"]]
    # 逐坐标精度程序只属于解同一行式 Dantzig 问题的路径；对照方法使用各自的
    # 区间公式，其 feasible 字段没有可行性的含义，因此按缺失记录。
    method_rows[[mcount]]$dantzig_feasible_fraction <- if (nrow(inf_tab) &&
                                                           any(is.finite(inf_tab$mu))) {
      mean(as.logical(inf_tab$feasible))
    } else NA_real_

    tcount <- tcount + 1L
    theory_rows[[tcount]] <- method_theory_diagnostics_v2(
      config, replicate, seed, method_id, fit, target_obj,
      population_asset, dat, run_id = run_id
    )
    method_diagnostics[[method_id]] <- list(
      status = fit$status,
      failure_stage = fit$failure_stage,
      failure_message = fit$failure_message,
      method_version = fit$method_version,
      score_definition = fit$score_definition,
      native_warning_records = fit$native_warning_records,
      beta_hat = fit$beta_hat %||% fit$fit_object$beta,
      beta_l1 = fit$fit_object$beta_l1,
      beta_refit = fit$fit_object$beta_refit,
      score = fit$fit_object$score %||% fit$fit_object$components$score,
      score_at_h = fit$fit_object$components$components_h$score,
      score_at_half_h = fit$fit_object$components$components_half_h$score,
      lambda_components = fit$fit_object$lambda_components,
      selected_support = fit$fit_object$selected_support,
      refit_set = fit$fit_object$refit_set,
      diagnostics = fit$fit_object$diagnostics,
      tuning = fit$tuning_diagnostics %||% fit$failure_diagnostics$tuning,
      cv = fit$inference_object$cv,
      directions = fit$inference_object$directions,
      source_precision = if (!is.null(fit$inference_object$precision)) list(
        definition = fit$inference_object$precision$precision_definition,
        common_gamma = fit$inference_object$precision$gamma,
        gamma_grid_requested = fit$inference_object$precision$gamma_grid_requested,
        excluded_zero_direction = fit$inference_object$precision$excluded_zero_direction,
        selected_matrix_attempts = fit$inference_object$precision$attempts,
        full_grid_attempts = fit$inference_object$precision$matrix_grid_record$attempts,
        full_matrix_feasibility_certificate =
          fit$inference_object$precision$matrix_grid_record$feasibility_certificate
      ) else NULL,
      requested_target_coordinates = method_targets,
      submission = fit$submission_diagnostics,
      prediction_intercept = fit$prediction_intercept,
      population_score_centering = fit$population_score_centering,
      failure_diagnostics = fit$failure_diagnostics,
      partial_fit = if (!fit$status %in% c("ok", "warning")) fit$fit_object else NULL,
      training_design = list(n_clusters = method_train$n_clusters,
                                total_rows = method_train$total_rows,
                                p = ncol(method_train$X), q = ncol(method_train$Z),
                                common_screen_roles = context$common_screen_roles)
    )
  }

  metrics <- do.call(rbind, method_rows)
  pooled_rows <- which(metrics$method_id == "POOLED-QR-LASSO")
  pooled_baseline_loss <- if (length(pooled_rows) == 1L) {
    metrics$pinball_loss[pooled_rows]
  } else NA_real_
  if (is.finite(pooled_baseline_loss) && pooled_baseline_loss > 0) {
    usable <- is.finite(metrics$pinball_loss)
    metrics$pinball_skill_pooled[usable] <-
      1 - metrics$pinball_loss[usable] / pooled_baseline_loss
  }
  list(
    metrics = metrics,
    coordinates = if (length(coordinate_rows)) do.call(rbind, coordinate_rows) else data.frame(),
    theory = if (length(theory_rows)) do.call(rbind, theory_rows) else data.frame(),
    method_diagnostics = method_diagnostics,
    screening = prepared$screening,
    config = config,
    seed = seed,
    run_id = run_id
  )
}

# Oracle and split-screening benchmark adapters plus the authoritative dispatcher.

benchmark_add_metadata_v2 <- function(ans,
                                      reference_identifier,
                                      fidelity_status = "pending") {
  ans$reference_identifier <- reference_identifier
  ans$adapter_fidelity_status <- fidelity_status
  ans
}

fit_benchmark_profile_true_nuisance_v2 <- function(train,
                                                   tau,
                                                   target_coords,
                                                   tuning,
                                                   seed,
                                                   control = list()) {
  if (is.null(train$random_effects) || is.null(rownames(train$random_effects))) {
    return(benchmark_failure_v2(
      "PROFILE-DQR-TRUE-NUISANCE", "oracle_context",
      "Generated random effects with cluster row names are required."
    ))
  }
  cluster_match <- match(as.character(train$cluster_id), rownames(train$random_effects))
  if (anyNA(cluster_match)) {
    return(benchmark_failure_v2(
      "PROFILE-DQR-TRUE-NUISANCE", "oracle_context",
      "Could not align generated random effects with cluster identifiers."
    ))
  }
  q <- ncol(train$Z)
  if (ncol(train$random_effects) < q) {
    return(benchmark_failure_v2(
      "PROFILE-DQR-TRUE-NUISANCE", "oracle_context",
      "Generated nuisance dimension is smaller than the fitted nuisance design."
    ))
  }
  b_rows <- train$random_effects[cluster_match, seq_len(q), drop = FALSE]
  nuisance_truth <- rowSums(train$Z * b_rows)
  oracle <- train
  oracle$y <- as.numeric(train$y - nuisance_truth)
  oracle$Z <- matrix(0, nrow = length(oracle$y), ncol = 1L,
                     dimnames = list(NULL, "known_nuisance_removed"))

  ans <- fit_benchmark_profile_dqr_v2(
    oracle, tau, target_coords, tuning, seed,
    control = modifyList(control, list(
      penalty_factor = control$penalty_factor %||% rep(1, ncol(oracle$X))
    ))
  )
  ans$method_id <- "PROFILE-DQR-TRUE-NUISANCE"
  ans$implementation_version <- "profile-v2-known-generated-nuisance"
  ans$oracle_nuisance_removed <- nuisance_truth
  ans$target_scope <- "structural_known_nuisance"
  benchmark_add_metadata_v2(
    ans,
    reference_identifier = "docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md#profile-dqr-true-nuisance",
    fidelity_status = "implementation_present_acceptance_pending"
  )
}

population_direction_for_fit_v2 <- function(asset, fit_names, coordinate) {
  if (is.null(asset$directions[[coordinate]])) {
    stop("Population-direction asset does not contain coordinate: ", coordinate)
  }
  source_omega <- asset$directions[[coordinate]]$omega
  omega <- setNames(rep(0, length(fit_names)), fit_names)
  common <- intersect(names(source_omega), fit_names)
  omega[common] <- source_omega[common]
  if (!coordinate %in% names(omega)) stop("Coordinate is not in fitted design: ", coordinate)
  list(
    omega = omega,
    source_residual = asset$directions[[coordinate]]$residual,
    source_l1 = asset$directions[[coordinate]]$l1_norm,
    source_l2 = asset$directions[[coordinate]]$l2_norm
  )
}

fit_benchmark_profile_pop_h_v2 <- function(train,
                                           tau,
                                           target_coords,
                                           tuning,
                                           seed,
                                           population_direction_asset,
                                           control = list()) {
  if (is.null(population_direction_asset)) {
    return(benchmark_failure_v2(
      "PROFILE-DQR-POP-H", "oracle_asset",
      "A scenario-specific population-direction asset is required."
    ))
  }
  set.seed(seed)
  start <- proc.time()[[3]]
  fit <- tryCatch(
    fit_profile_lasso_v2(
      y = train$y,
      X = train$X,
      Z = train$Z,
      cluster_id = train$cluster_id,
      tau = tau,
      h = tuning$h,
      lambda_beta = tuning$lambda_beta,
      lambda_gamma = tuning$lambda_gamma,
      penalty_factor = control$penalty_factor %||% rep(1, ncol(train$X)),
      control = control$fit_control %||% list()
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(benchmark_failure_v2(
      "PROFILE-DQR-POP-H", "penalised_fit", conditionMessage(fit),
      proc.time()[[3]] - start
    ))
  }

  coords <- if (is.character(target_coords)) target_coords else names(fit$beta)[as.integer(target_coords)]
  alpha <- 1 - (control$ci_level %||% 0.95)
  zcrit <- stats::qnorm(1 - alpha / 2)
  rows <- vector("list", length(coords))
  directions <- vector("list", length(coords))
  names(directions) <- coords

  for (ii in seq_along(coords)) {
    nm <- coords[ii]
    direction <- tryCatch(
      population_direction_for_fit_v2(population_direction_asset, names(fit$beta), nm),
      error = function(e) e
    )
    if (inherits(direction, "error")) {
      return(benchmark_failure_v2(
        "PROFILE-DQR-POP-H", "oracle_asset", conditionMessage(direction),
        proc.time()[[3]] - start
      ))
    }
    omega <- direction$omega
    k <- match(nm, names(fit$beta))
    beta_tilde <- as.numeric(fit$beta[k] - sum(omega * fit$components$score))
    sand <- cluster_sandwich_coordinate_v2(
      fit$components$cluster_scores, omega,
      df_correction = isTRUE(control$df_correction)
    )
    empirical_residual <- max(abs(as.numeric(fit$components$hessian %*% omega -
                                                replace(numeric(length(omega)), k, 1))))
    rows[[ii]] <- data.frame(
      coordinate = nm,
      index = k,
      feasible = TRUE,
      beta_hat = as.numeric(fit$beta[k]),
      beta_tilde = beta_tilde,
      se = sand$se,
      ci_lower = beta_tilde - zcrit * sand$se,
      ci_upper = beta_tilde + zcrit * sand$se,
      mu = NA_real_,
      dantzig_residual = empirical_residual,
      population_residual = direction$source_residual,
      omega_l1 = sum(abs(omega)),
      omega_l2 = sqrt(sum(omega^2)),
      stringsAsFactors = FALSE
    )
    directions[[ii]] <- list(
      omega = omega,
      residual = empirical_residual,
      population_residual = direction$source_residual,
      feasible = TRUE
    )
  }

  tab <- do.call(rbind, rows)
  elapsed <- proc.time()[[3]] - start
  ans <- list(
    method_id = "PROFILE-DQR-POP-H",
    status = if (fit$converged) "ok" else "warning",
    failure_stage = "",
    failure_message = "",
    beta_hat = fit$beta,
    beta_tilde = setNames(tab$beta_tilde, tab$coordinate),
    se = setNames(tab$se, tab$coordinate),
    ci_lower = setNames(tab$ci_lower, tab$coordinate),
    ci_upper = setNames(tab$ci_upper, tab$coordinate),
    selected = which(abs(fit$beta) > (control$selection_tol %||% 1e-8)),
    fit_object = fit,
    inference_object = list(table = tab, directions = directions),
    population_direction_asset = population_direction_asset,
    runtime_sec = elapsed,
    converged = fit$converged,
    kkt_residual = fit$kkt_residual,
    profile_identity_error = fit$components$profile_identity_error,
    max_nuisance_gradient = fit$components$max_nuisance_gradient,
    warning_messages = if (fit$converged) character() else
      "Penalised profile fit did not meet all stopping rules.",
    implementation_version = "profile-v2-population-direction",
    target_scope = "regularised_profile_population_direction"
  )
  benchmark_add_metadata_v2(
    ans,
    reference_identifier = "docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md#profile-dqr-pop-h",
    fidelity_status = "implementation_present_acceptance_pending"
  )
}

fit_benchmark_profile_split_v2 <- function(train,
                                           tau,
                                           target_coords,
                                           tuning,
                                           seed,
                                           context = list(),
                                           control = list()) {
  screening <- context$screening %||% NULL
  if (!is.null(screening) && identical(screening$method, "split_quantile_score")) {
    if (length(intersect(screening$screen_clusters, screening$fit_clusters))) {
      return(benchmark_failure_v2(
        "PROFILE-DQR-SPLIT", "screening",
        "Screening and inference cluster sets overlap."
      ))
    }
    ans <- fit_benchmark_profile_dqr_v2(train, tau, target_coords, tuning, seed, control)
    ans$method_id <- "PROFILE-DQR-SPLIT"
    ans$implementation_version <- "profile-v2-independent-cluster-screening"
    ans$screening <- screening
    ans$target_scope <- "regularised_profile_independent_split"
    return(benchmark_add_metadata_v2(
      ans,
      reference_identifier = "docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md#profile-dqr-split",
      fidelity_status = "implementation_present_acceptance_pending"
    ))
  }

  fraction <- control$screen_fraction %||% 0.5
  d <- control$screen_dim %||% min(300L, ncol(train$X))
  split <- split_clusters_v2(train$cluster_id, fraction = fraction, seed = seed + 17L)
  screen_dat <- subset_clusters_v2(train, split$a)
  fit_dat <- subset_clusters_v2(train, split$b)
  scr <- quantile_score_screen_v2(
    screen_dat$y, screen_dat$X, screen_dat$cluster_id, tau, d
  )
  target_idx <- match(target_coords, colnames(train$X))
  target_idx <- target_idx[!is.na(target_idx)]
  selected <- sort(unique(c(scr$selected, target_idx)))
  fit_dat <- subset_features_v2(fit_dat, selected)
  fit_targets <- intersect(as.character(target_coords), colnames(fit_dat$X))
  ans <- fit_benchmark_profile_dqr_v2(fit_dat, tau, fit_targets, tuning, seed, control)
  full_beta <- setNames(rep(0, ncol(train$X)), colnames(train$X))
  full_beta[names(ans$beta_hat)] <- ans$beta_hat
  ans$beta_hat <- full_beta
  ans$selected <- selected
  ans$method_id <- "PROFILE-DQR-SPLIT"
  ans$implementation_version <- "profile-v2-independent-cluster-screening"
  ans$screening <- list(
    method = "split_quantile_score",
    screen_clusters = split$a,
    fit_clusters = split$b,
    selected_global = selected,
    scores = scr$scores
  )
  ans$target_scope <- "regularised_profile_independent_split"
  benchmark_add_metadata_v2(
    ans,
    reference_identifier = "docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md#profile-dqr-split",
    fidelity_status = "implementation_present_acceptance_pending"
  )
}

# Authoritative dispatcher. Later files may replace individual pending adapters,
# but method IDs and target restrictions must remain governed by the acceptance
# contract.
fit_benchmark_by_id_v2 <- function(method_id, train, tau, target_coords,
                                   tuning, seed, context = list(), control = list()) {
  ans <- switch(
    method_id,
    `PROFILE-DQR` = fit_benchmark_profile_dqr_v2(
      train, tau, target_coords, tuning, seed, control
    ),
    `POOLED-QR-LASSO` = fit_benchmark_pooled_qr_lasso_v2(
      train, tau, target_coords, tuning, seed, control
    ),
    `SQR-DEBIASED-IID` = fit_benchmark_sqr_debiased_iid_v2(
      train, tau, target_coords, tuning, seed, control
    ),
    `PROFILE-DQR-TRUE-SUPPORT` = fit_benchmark_profile_true_support_v2(
      train, tau, target_coords, tuning, seed,
      active = context$active %||% stop("context$active is required"),
      control = control
    ),
    `PROFILE-DQR-TRUE-NUISANCE` = fit_benchmark_profile_true_nuisance_v2(
      train, tau, target_coords, tuning, seed, control
    ),
    `PROFILE-DQR-POP-H` = fit_benchmark_profile_pop_h_v2(
      train, tau, target_coords, tuning, seed,
      population_direction_asset = context$population_direction_asset,
      control = control
    ),
    `PROFILE-DQR-SPLIT` = fit_benchmark_profile_split_v2(
      train, tau, target_coords, tuning, seed,
      context = context, control = control
    ),
    `LQMM` = fit_benchmark_lqmm_v2(
      train, tau, target_coords, tuning, seed,
      random_slope = ncol(train$Z) >= 2L, control = control
    ),
    `ORACLE-LQMM` = fit_benchmark_lqmm_v2(
      train, tau, target_coords, tuning, seed,
      random_slope = ncol(train$Z) >= 2L, control = control
    ),
    `BIAS-ADJ-LQMM` = benchmark_not_implemented_v2(
      method_id,
      "Faithful Battagliola et al. bias-adjusted implementation and fidelity evidence are required."
    ),
    `DOUBLE-PEN-QLMM` = benchmark_not_implemented_v2(
      method_id,
      "Faithful Li, Liu and Luo double-penalised QLMM implementation and fidelity evidence are required."
    ),
    `QGEE-SCAD` = benchmark_not_implemented_v2(
      method_id,
      "Faithful Zu et al. quantile GEE/SCAD implementation and fidelity evidence are required."
    ),
    `QIF-SEE` = benchmark_not_implemented_v2(
      method_id,
      "Faithful Bhattacharya et al. QIF smooth-threshold implementation and fidelity evidence are required."
    ),
    `BAYES-MIXED-LASSO` = benchmark_not_implemented_v2(
      method_id,
      "Excluded from final quantile benchmarks; mean-prediction supplement only."
    ),
    stop("Unknown benchmark method_id: ", method_id)
  )

  if (is.null(ans$reference_identifier)) {
    ans$reference_identifier <- switch(
      method_id,
      `PROFILE-DQR` = "docs/METHOD_SPECIFICATION.md",
      `POOLED-QR-LASSO` = "quantreg::rq.fit.lasso",
      `SQR-DEBIASED-IID` = "Yan-Wang-Zhang-2023-JMLR-22-1217",
      `LQMM` = "Geraci-Bottai-2014",
      `ORACLE-LQMM` = "Geraci-Bottai-2014",
      "pending"
    )
  }
  if (is.null(ans$adapter_fidelity_status)) {
    ans$adapter_fidelity_status <- if (identical(ans$status, "not_implemented")) {
      "not_implemented"
    } else if (method_id == "SQR-DEBIASED-IID") {
      "formula_fidelity_pending"
    } else {
      "acceptance_pending"
    }
  }
  ans
}

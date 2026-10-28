profile_pair_asset_version_v2 <- function() {
  list(method_version = "profile_pair_20261022",
       score_definition = "profile_bandwidth_pair_2_half_minus_full",
       construction_version = "profile-pair-population-assets-20261022-v3")
}

profile_pair_asset_commit_v2 <- function(implementation_commit, development) {
  if (!is.character(implementation_commit) || length(implementation_commit) != 1L ||
      is.na(implementation_commit) || !nzchar(trimws(implementation_commit))) {
    stop("A nonempty implementation commit identity is required.")
  }
  commit <- trimws(implementation_commit)
  is_commit <- grepl("^([0-9a-fA-F]{40}|[0-9a-fA-F]{64})$", commit)
  if (!isTRUE(development) && !is_commit) {
    stop("Final paired assets require a real implementation commit identity.")
  }
  if (isTRUE(development) && !is_commit && !identical(commit, "working_tree")) {
    stop("Development commit identity must be working_tree or a real implementation commit.")
  }
  if (is_commit) tolower(commit) else commit
}

profile_pair_asset_control_v2 <- function(control = list()) {
  value <- utils::modifyList(list(
    max_iter = 1000L, reltol = 1e-12, max_root_iter = 100L,
    max_backtrack = 40L, score_tol = 1e-7, nuisance_tol = 1e-7,
    identity_tol = 1e-8, eigen_min = 1e-8, condition_max = 1e10,
    population_engine = "reference", batch_clusters = 1000L
  ), control)
  positive <- c("max_iter", "reltol", "max_root_iter", "max_backtrack",
                "score_tol", "nuisance_tol", "identity_tol", "eigen_min", "condition_max")
  if (any(!vapply(value[positive], function(x) {
    length(x) == 1L && is.finite(x) && x > 0
  }, logical(1)))) stop("Population control values must be positive finite scalars.")
  if (value$score_tol > 1e-7 || value$nuisance_tol > 1e-7 ||
      value$identity_tol > 1e-8 || value$eigen_min < 1e-8 ||
      value$condition_max > 1e10) stop("Population control exceeds the declared acceptance limits.")
  if (!value$population_engine %in% c("reference", "streamed", "conditional_gaussian") ||
      length(value$batch_clusters) != 1L || !is.finite(value$batch_clusters) ||
      value$batch_clusters < 1L || value$batch_clusters != as.integer(value$batch_clusters)) {
    stop("Population engine or batch size is invalid.")
  }
  value
}

profile_pair_design_identity_v2 <- function(config) {
  fields <- c("p", "s", "tau", "q", "m_values", "m_rule", "m_multiset",
              "design_type", "rho_x", "signal", "error_dist", "heteroskedastic",
              "error_dependence", "copula_rho", "random_effect_dist", "sigma_b0",
              "sigma_b1", "rho_b", "x_b_corr", "informative_size",
              "nonlinear_re_strength", "fit_random_effects", "lambda_gamma", "target_mode")
  list(experiment_id = as.character(config$experiment_id),
       design = stats::setNames(lapply(fields, function(x) config[[x]]), fields),
       target_columns = seq_len(as.integer(config$p)),
       design_scope = "full_registered_fixed_effect_design")
}

profile_pair_asset_scale_v2 <- function(n_population, repeats, development) {
  if (length(n_population) != 1L || !is.finite(n_population) ||
      n_population != as.integer(n_population) || n_population < 2L ||
      length(repeats) != 1L || !is.finite(repeats) ||
      repeats != as.integer(repeats) || repeats < 2L) {
    stop("Population construction needs integer n_population >= 2 and repeats >= 2.")
  }
  if (!is.logical(development) || length(development) != 1L || is.na(development)) {
    stop("development must be TRUE or FALSE.")
  }
  if (!development && (n_population < 100000L || repeats < 4L)) {
    stop("Final paired assets require n_population >= 100000 and repeats >= 4.")
  }
  invisible(TRUE)
}

profile_pair_curvature_v2 <- function(H, control = list()) {
  ctl <- profile_pair_asset_control_v2(control)
  H <- assert_numeric_matrix_v2(H, "paired_population_hessian")
  if (nrow(H) != ncol(H) || nrow(H) < 1L) stop("Paired Hessian must be nonempty and square.")
  if (max(abs(H - t(H))) > 1e-8) stop("Paired Hessian is not symmetric.")
  values <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  condition <- max(values) / min(values)
  if (min(values) <= ctl$eigen_min || !is.finite(condition) ||
      condition >= ctl$condition_max) {
    stop("Paired population curvature failed: min eigenvalue = ", signif(min(values), 6),
         ", condition = ", signif(condition, 6), ".")
  }
  list(min_eigenvalue = min(values), max_eigenvalue = max(values),
       condition_number = condition, eigenvalues = values)
}

profile_pair_population_data_v2 <- function(config, n_population, seed) {
  if (length(config$p) != 1L || !is.finite(config$p) || config$p < 1L) {
    stop("config$p must specify the complete positive fixed-effect dimension.")
  }
  if (!is.null(config$m_multiset) && length(config$m_multiset) != n_population) {
    stop("The frozen visit multiset must match the population size; no resampling rule was declared.")
  }
  cfg <- config
  cfg$n_clusters <- as.integer(n_population)
  dat <- simulate_from_config_v2(cfg, seed = seed)
  fit_working_nuisance_design_v2(dat, config$fit_random_effects)
}

profile_pair_moments_v2 <- function(components) {
  expected <- profile_pair_asset_version_v2()$score_definition
  if (!identical(components$score_definition, expected)) stop("Paired score components are required.")
  if (!is.null(components$population_moments)) return(components$population_moments)
  half <- components$components_half_h$cluster_scores
  full <- components$components_h$cluster_scores
  G <- components$cluster_scores
  means <- list(score = colMeans(G), half = colMeans(half), full = colMeans(full))
  centered <- sweep(G, 2, means$score)
  centered_half <- sweep(half, 2, means$half)
  centered_full <- sweep(full, 2, means$full)
  n <- nrow(G)
  covariance <- crossprod(centered) / n
  cov_half <- crossprod(centered_half) / n
  cov_full <- crossprod(centered_full) / n
  cross <- crossprod(centered_half, centered_full) / n
  identity <- max(abs(covariance - (4 * cov_half + cov_full - 2 * cross - 2 * t(cross))))
  if (!is.finite(identity) || identity > 1e-8) stop("Cross-bandwidth covariance identity failed.")
  list(n_clusters = n, score_mean = means$score,
       score_mean_half = means$half, score_mean_full = means$full,
       Sigma_population = covariance, Sigma_half = cov_half, Sigma_full = cov_full,
       cross_covariance_half_full = cross, covariance_identity_error = identity,
       H_population = components$hessian,
       H_half = components$components_half_h$hessian,
       H_full = components$components_h$hessian)
}

solve_profile_pair_target_v2 <- function(dat, tau, h, lambda_gamma,
                                         beta_start = NULL, control = list()) {
  ctl <- profile_pair_asset_control_v2(control)
  validate_tau_h_v2(tau, h)
  X <- assert_numeric_matrix_v2(dat$X, "X")
  if (is.null(colnames(X)) || anyDuplicated(colnames(X))) stop("X must have unique coordinate names.")
  initialisation_source <- if (is.null(beta_start)) {
    "generated_parameter_initialisation_only"
  } else "supplied_initialisation"
  if (is.null(beta_start)) beta_start <- dat$beta0
  if (length(beta_start) != ncol(X) || any(!is.finite(beta_start))) {
    stop("A finite full-design starting vector is required.")
  }
  cache <- new.env(parent = emptyenv())
  cache$beta <- NULL
  cache$components <- NULL
  cache$level <- -1L
  evaluations <- 0L
  evaluate <- function(beta, level = 0L) {
    if (!is.null(cache$beta) && identical(as.numeric(beta), cache$beta) && cache$level >= level) {
      return(cache$components)
    }
    components <- profile_pair_population_evaluate_v2(dat, beta, tau, h, lambda_gamma, ctl, level)
    cache$beta <- as.numeric(beta)
    cache$components <- components
    cache$level <- if (identical(ctl$population_engine, "reference")) 2L else level
    evaluations <<- evaluations + 1L
    components
  }
  started <- proc.time()
  optimum <- stats::optim(
    as.numeric(beta_start), function(beta) evaluate(beta)$objective,
    function(beta) as.numeric(evaluate(beta)$score), method = "BFGS",
    control = list(maxit = ctl$max_iter, reltol = ctl$reltol)
  )
  if (optimum$convergence != 0L) {
    stop("Population BFGS did not converge, status ", optimum$convergence, ".")
  }
  beta <- as.numeric(optimum$par)
  root_iterations <- 0L
  last_step <- 0
  for (iteration in seq_len(ctl$max_root_iter)) {
    components <- evaluate(beta, 1L)
    score <- as.numeric(components$score)
    if (max(abs(score)) <= ctl$score_tol) break
    curvature <- profile_pair_curvature_v2(components$hessian, ctl)
    direction <- -as.numeric(solve(components$hessian, score))
    merit <- sum(score^2)
    accepted <- FALSE
    for (backtrack in 0:ctl$max_backtrack) {
      step <- 2^(-backtrack)
      candidate <- beta + step * direction
      trial <- evaluate(candidate)
      if (sum(trial$score^2) <= (1 - 1e-4 * step) * merit) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) stop("Population score-root line search failed.")
    last_step <- max(abs(candidate - beta))
    beta <- candidate
    root_iterations <- iteration
  }
  components <- evaluate(beta, 2L)
  curvature <- profile_pair_curvature_v2(components$hessian, ctl)
  kkt <- max(abs(components$score))
  if (!is.finite(kkt) || kkt > ctl$score_tol ||
      components$max_nuisance_gradient > ctl$nuisance_tol ||
      components$profile_identity_error > ctl$identity_tol) {
    stop("Population root failed score, nuisance, or profile-identity acceptance.")
  }
  names(beta) <- colnames(X)
  elapsed <- proc.time() - started
  list(beta = beta, components = components, moments = profile_pair_moments_v2(components),
       converged = TRUE, kkt_residual = kkt,
       max_nuisance_gradient = components$max_nuisance_gradient,
       nuisance_converged = components$nuisance_converged,
       profile_identity_error = components$profile_identity_error,
       curvature = curvature, bandwidths = c(h / 2, h),
       optimizer = "BFGS_then_Newton_score_root", optimizer_status = optimum$convergence,
       initial_beta = stats::setNames(as.numeric(beta_start), colnames(X)),
       initialisation_source = initialisation_source,
       optimizer_counts = optimum$counts, root_iterations = root_iterations,
       population_engine = ctl$population_engine, batch_clusters = ctl$batch_clusters,
       population_matrix_backend = getOption("jpDMEQR.population_matrix_backend", "R_BLAS"),
       population_nuisance_backend = getOption("jpDMEQR.population_nuisance_backend", "vectorised_R"),
       last_root_step = last_step, component_evaluations = evaluations,
       runtime_sec = unname(elapsed[["elapsed"]]),
       cpu_sec = unname(elapsed[["user.self"]] + elapsed[["sys.self"]]))
}

profile_pair_checkpoint_v2 <- function(checkpoint, label, value) {
  if (!is.null(checkpoint)) {
    if (!is.function(checkpoint)) stop("checkpoint must be a function or NULL.")
    checkpoint(label, value)
  }
  invisible(NULL)
}

profile_pair_fit_record_v2 <- function(fit) {
  fit$components <- NULL
  fit
}

build_profile_pair_target_asset_v2 <- function(config, n_population = 100000L,
                                               repeats = 4L, seed = 87001L,
                                               h_target = NULL, development = FALSE,
                                               implementation_commit = "working_tree",
                                               control = list(), checkpoint = NULL) {
  implementation_commit <- profile_pair_asset_commit_v2(implementation_commit, development)
  profile_pair_asset_scale_v2(n_population, repeats, development)
  ctl <- profile_pair_asset_control_v2(control)
  if (identical(ctl$population_engine, "conditional_gaussian")) {
    stop("Population target roots require the full reference or streamed empirical score.")
  }
  h_analysis <- config$h_multiplier * config$n_clusters^(-3 / 10)
  if (is.null(h_target)) h_target <- min(h_analysis / 4, n_population^(-1 / 3))
  if (length(h_target) != 1L || !is.finite(h_target) || h_target <= 0 || h_target >= h_analysis) {
    stop("h_target must be positive and smaller than the analysis bandwidth.")
  }
  fits <- vector("list", repeats)
  seeds <- vapply(seq_len(repeats), function(rr) {
    seed_from_id_v2(seed, config$experiment_id, rr, "paired_population_target")
  }, integer(1))
  for (rr in seq_len(repeats)) {
    dat <- profile_pair_population_data_v2(config, n_population, seeds[rr])
    if (identical(ctl$population_engine, "streamed")) dat$population_layout <- profile_pair_population_layout_v2(dat)
    base <- solve_profile_pair_target_v2(dat, config$tau, h_target,
                                        config$lambda_gamma, control = ctl)
    profile_pair_checkpoint_v2(checkpoint, sprintf("target_%02d_base", rr), base)
    refinement <- solve_profile_pair_target_v2(dat, config$tau, h_target / 2,
                                              config$lambda_gamma, base$beta, ctl)
    profile_pair_checkpoint_v2(checkpoint, sprintf("target_%02d_refinement", rr), refinement)
    analysis <- solve_profile_pair_target_v2(dat, config$tau, h_analysis,
                                            config$lambda_gamma, base$beta, ctl)
    profile_pair_checkpoint_v2(checkpoint, sprintf("target_%02d_analysis", rr), analysis)
    at_analysis <- profile_pair_population_evaluate_v2(
      dat, base$beta, config$tau, h_analysis, config$lambda_gamma, ctl, 2L)
    fits[[rr]] <- list(seed = seeds[rr], n_population = as.integer(n_population),
                       target = profile_pair_fit_record_v2(base),
                       refinement = profile_pair_fit_record_v2(refinement),
                       analysis = profile_pair_fit_record_v2(analysis),
                       analysis_score_at_target = c(profile_pair_moments_v2(at_analysis), list(
                         nuisance_converged = at_analysis$nuisance_converged,
                         max_nuisance_gradient = at_analysis$max_nuisance_gradient,
                         profile_identity_error = at_analysis$profile_identity_error)))
    rm(dat, base, refinement, analysis, at_analysis)
    invisible(gc(FALSE))
  }
  estimates <- do.call(rbind, lapply(fits, function(x) x$target$beta))
  refined <- do.call(rbind, lapply(fits, function(x) x$refinement$beta))
  analysis_estimates <- do.call(rbind, lapply(fits, function(x) x$analysis$beta))
  version <- profile_pair_asset_version_v2()
  asset <- c(version, list(
    asset_kind = "paired_profile_target", development = development,
    target_type = "numerical_approximation_unsmoothed_regularised_profile_target",
    target_scope = "full_registered_fixed_effect_design",
    experiment_id = config$experiment_id, design_identity = profile_pair_design_identity_v2(config),
    target_columns = seq_len(config$p), target_names = colnames(estimates),
    beta_star_mc = colMeans(estimates), estimates_by_repeat = estimates,
    target_source = "paired_population_target_mc",
    evaluation_point = colMeans(estimates),
    target_mc_sd_by_coordinate = apply(estimates, 2, stats::sd),
    target_mc_se_by_coordinate = apply(estimates, 2, stats::sd) / sqrt(repeats),
    target_mc_se = max(apply(estimates, 2, stats::sd) / sqrt(repeats)),
    between_run_range_by_coordinate = apply(estimates, 2, function(x) diff(range(x))),
    between_run_max_difference = max(apply(estimates, 2, function(x) diff(range(x)))),
    beta_refinement_mc = colMeans(refined), refinement_estimates_by_repeat = refined,
    refinement_mc_se_by_coordinate = apply(refined, 2, stats::sd) / sqrt(repeats),
    beta_analysis_bandwidth_mc = colMeans(analysis_estimates),
    analysis_estimates_by_repeat = analysis_estimates,
    analysis_bandwidth_mc_se_by_coordinate = apply(analysis_estimates, 2, stats::sd) / sqrt(repeats),
    target_bandwidth_displacement_by_repeat = refined - estimates,
    smoothing_target_displacement_by_repeat = analysis_estimates - estimates,
    smoothing_target_distance = mean(sqrt(rowSums((analysis_estimates - estimates)^2))),
    smoothing_target_distance_l1 = mean(rowSums(abs(analysis_estimates - estimates))),
    smoothing_target_distance_by_repeat = sqrt(rowSums((analysis_estimates - estimates)^2)),
    target_refinement_distance_by_repeat = sqrt(rowSums((refined - estimates)^2)),
    h_target = h_target, bandwidths_target = c(h_target / 2, h_target),
    bandwidths_refinement = c(h_target / 4, h_target / 2),
    h_analysis = h_analysis, bandwidths_analysis = c(h_analysis / 2, h_analysis),
    h_analysis_rule = h_analysis, h_population_rule = h_target,
    n_analysis = config$n_clusters, n_population = as.integer(n_population),
    repeats = as.integer(repeats), seeds = seeds, repeat_assets = fits,
    all_converged = TRUE, acceptance = ctl, solver_configuration = solver_configuration_v2(),
    population_matrix_backend = getOption("jpDMEQR.population_matrix_backend", "R_BLAS"),
    max_kkt_residual = max(vapply(fits, function(x) x$target$kkt_residual, numeric(1))),
    max_nuisance_gradient = max(vapply(fits, function(x) {
      max(x$target$max_nuisance_gradient, x$refinement$max_nuisance_gradient,
          x$analysis$max_nuisance_gradient, x$analysis_score_at_target$max_nuisance_gradient)
    }, numeric(1))),
    profile_identity_error = max(vapply(fits, function(x) {
      max(x$target$profile_identity_error, x$refinement$profile_identity_error,
          x$analysis$profile_identity_error, x$analysis_score_at_target$profile_identity_error)
    }, numeric(1))), nuisance_converged = TRUE,
    implementation_commit = implementation_commit, config = config,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ))
  validate_profile_pair_asset_v2(asset, "paired_profile_target", config, final = !development)
  asset
}

profile_pair_pool_covariance_v2 <- function(fits, mean_name, covariance_name) {
  means <- lapply(fits, `[[`, mean_name)
  pooled_mean <- Reduce(`+`, means) / length(means)
  second <- Map(function(x, mu) x[[covariance_name]] + tcrossprod(mu), fits, means)
  Reduce(`+`, second) / length(second) - tcrossprod(pooled_mean)
}

build_profile_pair_direction_asset_v2 <- function(config, target_asset = NULL,
                                                  n_population = 100000L,
                                                  repeats = 4L, seed = 97001L,
                                                  requested_coordinates = NULL,
                                                  development = FALSE,
                                                  implementation_commit = "working_tree",
                                                  control = list(), checkpoint = NULL) {
  implementation_commit <- profile_pair_asset_commit_v2(implementation_commit, development)
  profile_pair_asset_scale_v2(n_population, repeats, development)
  if (!is.null(target_asset)) {
    validate_profile_pair_asset_v2(target_asset, "paired_profile_target", config,
                                  final = !development)
  } else if (!identical(config$target_mode, "structural")) {
    stop("A profile_mc population direction requires its independent target asset.")
  }
  ctl <- profile_pair_asset_control_v2(control)
  target_mode <- config$target_mode
  if (!target_mode %in% c("structural", "profile_mc")) {
    stop("Paired population directions require structural or profile_mc target_mode.")
  }
  if (identical(ctl$population_engine, "conditional_gaussian")) profile_pair_gaussian_row_conditions_v2(config)
  target_source <- if (identical(target_mode, "structural")) {
    "generated_structural_coefficient"
  } else "paired_population_target_mc"
  beta_target <- if (identical(target_mode, "profile_mc")) target_asset$beta_star_mc else NULL
  h_analysis <- config$h_multiplier * config$n_clusters^(-3 / 10)
  h_target <- target_asset$h_target %||% min(h_analysis / 4, n_population^(-1 / 3))
  seeds <- vapply(seq_len(repeats), function(rr) {
    seed_from_id_v2(seed, config$experiment_id, rr, "paired_population_direction")
  }, integer(1))
  if (any(seeds %in% target_asset$seeds)) stop("Target and direction population seeds must be disjoint.")
  fits <- vector("list", repeats)
  for (rr in seq_len(repeats)) {
    dat <- if (identical(ctl$population_engine, "conditional_gaussian")) {
      profile_pair_conditional_gaussian_data_v2(config, n_population, seeds[rr])
    } else profile_pair_population_data_v2(config, n_population, seeds[rr])
    if (identical(ctl$population_engine, "streamed")) dat$population_layout <- profile_pair_population_layout_v2(dat)
    generated_beta0 <- stats::setNames(as.numeric(dat$beta0), names(dat$beta0))
    if (identical(target_mode, "structural")) {
      if (rr == 1L) beta_target <- generated_beta0
      if (!identical(generated_beta0, beta_target)) {
        stop("Structural coefficients differ across population direction repeats.")
      }
    }
    comp <- profile_pair_population_evaluate_v2(
      dat, beta_target, config$tau, h_analysis, config$lambda_gamma, ctl, 2L)
    curvature <- profile_pair_curvature_v2(comp$hessian, ctl)
    fits[[rr]] <- c(profile_pair_moments_v2(comp), list(
      seed = seeds[rr], curvature = curvature, nuisance_converged = comp$nuisance_converged,
      target_source = target_source, evaluation_point = beta_target,
      generated_beta0 = generated_beta0,
      population_engine = comp$population_engine %||% "reference",
      population_matrix_backend = comp$population_matrix_backend %||%
        getOption("jpDMEQR.population_matrix_backend", "R_BLAS"),
      population_nuisance_backend = comp$population_nuisance_backend %||%
        getOption("jpDMEQR.population_nuisance_backend", "vectorised_R"),
      population_data_definition = dat$population_data_definition %||% "full_reference_DGP_observations",
      max_nuisance_gradient = comp$max_nuisance_gradient,
      profile_identity_error = comp$profile_identity_error,
      Omega_population = solve(comp$hessian)
    ))
    profile_pair_checkpoint_v2(checkpoint, sprintf("direction_%02d", rr), fits[[rr]])
    rm(dat, comp)
    invisible(gc(FALSE))
  }
  H <- Reduce(`+`, lapply(fits, `[[`, "H_population")) / repeats
  curvature <- profile_pair_curvature_v2(H, ctl)
  Omega <- solve(H)
  Sigma <- profile_pair_pool_covariance_v2(fits, "score_mean", "Sigma_population")
  Sigma_half <- profile_pair_pool_covariance_v2(fits, "score_mean_half", "Sigma_half")
  Sigma_full <- profile_pair_pool_covariance_v2(fits, "score_mean_full", "Sigma_full")
  mean_half <- Reduce(`+`, lapply(fits, `[[`, "score_mean_half")) / repeats
  mean_full <- Reduce(`+`, lapply(fits, `[[`, "score_mean_full")) / repeats
  cross <- Reduce(`+`, lapply(fits, function(x) {
    x$cross_covariance_half_full + tcrossprod(x$score_mean_half, x$score_mean_full)
  })) / repeats - tcrossprod(mean_half, mean_full)
  covariance_identity_error <- max(abs(Sigma - (4 * Sigma_half + Sigma_full - 2 * cross - 2 * t(cross))))
  if (covariance_identity_error > ctl$identity_tol) stop("Pooled covariance identity failed.")
  names_target <- names(beta_target)
  dimnames(Omega) <- dimnames(H) <- dimnames(Sigma) <- list(names_target, names_target)
  if (is.null(requested_coordinates)) requested_coordinates <- names_target
  if (anyDuplicated(requested_coordinates) || !length(requested_coordinates) ||
      !all(requested_coordinates %in% names_target)) stop("Invalid requested population coordinates.")
  directions <- lapply(requested_coordinates, function(nm) {
    k <- match(nm, names_target)
    omega <- Omega[, k]
    e <- replace(numeric(length(omega)), k, 1)
    list(coordinate = nm, omega = stats::setNames(as.numeric(omega), names_target),
         residual = max(abs(H %*% omega - e)),
         l1_norm = sum(abs(omega)), l2_norm = sqrt(sum(omega^2)))
  })
  names(directions) <- requested_coordinates
  repeat_matrices <- function(field) {
    array(unlist(lapply(fits, `[[`, field), use.names = FALSE),
          dim = c(length(names_target), length(names_target), repeats),
          dimnames = list(names_target, names_target, NULL))
  }
  omega_repeats <- repeat_matrices("Omega_population")
  hessian_repeats <- repeat_matrices("H_population")
  covariance_repeats <- repeat_matrices("Sigma_population")
  sigma <- sqrt(colSums(Omega * profile_pair_population_crossprod_v2(Sigma, Omega)))
  score_means_by_repeat <- do.call(rbind, lapply(fits, `[[`, "score_mean"))
  sigma_by_repeat <- do.call(rbind, lapply(fits, function(fit) {
    if (!is.null(fit$conditional_integration)) {
      return(sqrt(diag(fit$Omega_population) *
        fit$conditional_integration$covariance_coefficient /
        fit$conditional_integration$hessian_coefficient))
    }
    sqrt(colSums(fit$Omega_population * profile_pair_population_crossprod_v2(
      fit$Sigma_population, fit$Omega_population)))
  }))
  colnames(sigma_by_repeat) <- names_target
  if (any(!is.finite(sigma)) || any(sigma <= 0)) stop("Population projected variances must be positive.")
  version <- profile_pair_asset_version_v2()
  asset <- c(version, list(
    asset_kind = "paired_population_direction", development = development,
    experiment_id = config$experiment_id, design_identity = profile_pair_design_identity_v2(config),
    target_columns = seq_len(config$p), target_names = names_target,
    beta_target = beta_target, evaluation_point = beta_target,
    target_source = target_source, target_mode = target_mode,
    target_asset_identity = list(
      reference_source = if (is.null(target_asset)) "structural_definition" else "independent_population_target",
      construction_version = target_asset$construction_version %||% version$construction_version,
      implementation_commit = target_asset$implementation_commit %||% implementation_commit,
      seeds = target_asset$seeds %||% integer(), h_target = h_target,
      beta_star_mc = target_asset$beta_star_mc),
    H_population = H, Omega_population = Omega, Sigma_population = Sigma,
    Sigma_half = Sigma_half, Sigma_full = Sigma_full,
    cross_covariance_half_full = cross, covariance_identity_error = covariance_identity_error,
    score_mean_population = Reduce(`+`, lapply(fits, `[[`, "score_mean")) / repeats,
    score_mean_by_repeat = score_means_by_repeat,
    score_mean_mc_se = apply(score_means_by_repeat, 2, stats::sd) / sqrt(repeats),
    sigma_population = stats::setNames(sigma, names_target), directions = directions,
    sigma_by_repeat = sigma_by_repeat,
    sigma_mc_se = apply(sigma_by_repeat, 2, stats::sd) / sqrt(repeats),
    sigma0_pop = stats::setNames(sigma, names_target),
    H_mc_se = apply(hessian_repeats, c(1, 2), stats::sd) / sqrt(repeats),
    Omega_mc_se = apply(omega_repeats, c(1, 2), stats::sd) / sqrt(repeats),
    Sigma_mc_se = apply(covariance_repeats, c(1, 2), stats::sd) / sqrt(repeats),
    repeat_assets = fits, n_population = as.integer(n_population), repeats = as.integer(repeats),
    n_analysis = config$n_clusters, h_analysis = h_analysis,
    bandwidths_analysis = c(h_analysis / 2, h_analysis),
    h_target = h_target, bandwidths_target = c(h_target / 2, h_target),
    seeds = seeds, all_converged = TRUE, curvature = curvature,
    condition_number = curvature$condition_number, raw_eigenvalues = curvature$eigenvalues,
    nuisance_converged = TRUE,
    max_nuisance_gradient = max(vapply(fits, `[[`, numeric(1), "max_nuisance_gradient")),
    profile_identity_error = max(vapply(fits, `[[`, numeric(1), "profile_identity_error")),
    acceptance = ctl, solver_configuration = solver_configuration_v2(),
    population_matrix_backend = getOption("jpDMEQR.population_matrix_backend", "R_BLAS"),
    implementation_commit = implementation_commit, config = config,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ))
  validate_profile_pair_asset_v2(asset, "paired_population_direction", config, final = !development)
  asset
}

validate_profile_pair_asset_v2 <- function(asset, kind, config = NULL,
                                           final = TRUE, expected_commit = NULL) {
  if (!is.logical(final) || length(final) != 1L || is.na(final)) {
    stop("final must be TRUE or FALSE.")
  }
  profile_pair_asset_commit_v2(asset$implementation_commit, development = !final)
  version <- profile_pair_asset_version_v2()
  for (nm in names(version)) {
    if (!identical(asset[[nm]], version[[nm]])) {
      stop("Paired asset rejected: incompatible ", nm, ".")
    }
  }
  if (!kind %in% c("paired_profile_target", "paired_population_direction") ||
      !identical(asset$asset_kind, kind)) stop("Paired asset rejected: asset_kind mismatch.")
  profile_pair_asset_scale_v2(asset$n_population, asset$repeats, asset$development)
  if (final && (isTRUE(asset$development) || asset$n_population < 100000L || asset$repeats < 4L)) {
    stop("Development paired assets cannot pass the final population gate.")
  }
  if (!is.null(config) && !identical(asset$design_identity, profile_pair_design_identity_v2(config))) {
    stop("Paired asset rejected: scenario or design identity mismatch.")
  }
  if (!identical(asset$design_identity, profile_pair_design_identity_v2(asset$config)) ||
      !identical(as.character(asset$experiment_id), as.character(asset$config$experiment_id))) {
    stop("Paired asset rejected: stored scenario or working-design identity is inconsistent.")
  }
  if (!is.null(config) && (!identical(as.integer(asset$n_analysis), as.integer(config$n_clusters)) ||
      !isTRUE(all.equal(asset$h_analysis, config$h_multiplier * config$n_clusters^(-3 / 10),
                        tolerance = 1e-14)))) {
    stop("Paired asset rejected: analysis bandwidth identity mismatch.")
  }
  if (!is.null(expected_commit) &&
      !profile_pair_population_commit_matches_v2(
        asset$implementation_commit, expected_commit,
        engine = asset$acceptance$population_engine %||% "reference")) {
    stop("Paired asset rejected: implementation commit mismatch.")
  }
  if (!isTRUE(asset$all_converged) || length(asset$repeat_assets) != asset$repeats ||
      length(asset$seeds) != asset$repeats || anyDuplicated(asset$seeds)) {
    stop("Paired asset rejected: incomplete or non-independent repeat records.")
  }
  expected_names <- sprintf("x%05d", asset$target_columns)
  if (!identical(asset$target_names, expected_names) ||
      !identical(asset$target_columns, seq_len(asset$config$p))) {
    stop("Paired asset rejected: full-design coordinate coverage is required.")
  }
  if (!is.finite(asset$h_target) || asset$h_target <= 0 || asset$h_target >= asset$h_analysis ||
      !identical(asset$bandwidths_target, c(asset$h_target / 2, asset$h_target)) ||
      !identical(asset$bandwidths_analysis, c(asset$h_analysis / 2, asset$h_analysis))) {
    stop("Paired asset rejected: target or analysis bandwidth records are inconsistent.")
  }
  ctl <- profile_pair_asset_control_v2(asset$acceptance)
  if (identical(ctl$population_engine, "conditional_gaussian")) {
    if (kind != "paired_population_direction") stop("Conditional Gaussian assets must be population directions.")
    profile_pair_gaussian_row_conditions_v2(asset$config)
  }
  if (kind == "paired_profile_target") {
    if (!identical(asset$target_source, "paired_population_target_mc") ||
        !identical(asset$evaluation_point, asset$beta_star_mc)) {
      stop("Paired target source or evaluation point is inconsistent.")
    }
    if (!identical(names(asset$beta_star_mc), expected_names) ||
        any(!is.finite(asset$beta_star_mc)) ||
        any(!is.finite(asset$target_mc_se_by_coordinate)) ||
        any(!is.finite(asset$target_bandwidth_displacement_by_repeat))) {
      stop("Paired target estimate or uncertainty is missing or non-finite.")
    }
    if (!is.matrix(asset$estimates_by_repeat) ||
        !identical(dim(asset$estimates_by_repeat), c(asset$repeats, length(expected_names))) ||
        !identical(colnames(asset$estimates_by_repeat), expected_names) ||
        any(!is.finite(asset$estimates_by_repeat)) ||
        !isTRUE(all.equal(colMeans(asset$estimates_by_repeat), asset$beta_star_mc, tolerance = 1e-12))) {
      stop("Paired target repeat estimates are incomplete or inconsistent.")
    }
    for (rr in seq_along(asset$repeat_assets)) {
      record <- asset$repeat_assets[[rr]]
      if (!identical(as.integer(record$seed), as.integer(asset$seeds[rr])) ||
          record$n_population != asset$n_population) {
        stop("Paired target repeat identity or population size is inconsistent.")
      }
      for (role in c("target", "refinement", "analysis")) {
        fit <- record[[role]]
        if (!isTRUE(fit$converged) || !isTRUE(fit$nuisance_converged) ||
            fit$moments$n_clusters != asset$n_population ||
            fit$kkt_residual > ctl$score_tol || fit$max_nuisance_gradient > ctl$nuisance_tol ||
            fit$profile_identity_error > ctl$identity_tol) stop("Paired target diagnostics failed.")
        profile_pair_curvature_v2(fit$moments$H_population, ctl)
      }
      diagnostic <- record$analysis_score_at_target
      if (!isTRUE(diagnostic$nuisance_converged) ||
          diagnostic$max_nuisance_gradient > ctl$nuisance_tol ||
          diagnostic$profile_identity_error > ctl$identity_tol) {
        stop("Analysis-bandwidth target-score diagnostics failed.")
      }
    }
  } else {
    target_mode <- asset$config$target_mode
    if (!target_mode %in% c("structural", "profile_mc") ||
        !identical(asset$target_mode, target_mode) ||
        !identical(names(asset$beta_target), expected_names) ||
        any(!is.finite(asset$beta_target)) ||
        !identical(asset$evaluation_point, asset$beta_target)) {
      stop("Paired direction target mode or evaluation point is inconsistent.")
    }
    expected_source <- if (identical(target_mode, "structural")) {
      "generated_structural_coefficient"
    } else "paired_population_target_mc"
    expected_point <- if (identical(target_mode, "structural")) {
      make_beta_v2(asset$config$p, asset$config$s, asset$config$signal)
    } else asset$target_asset_identity$beta_star_mc
    if (!identical(asset$target_source, expected_source) ||
        !identical(asset$beta_target, expected_point) ||
        !identical(asset$target_asset_identity$construction_version, version$construction_version)) {
      stop("Paired direction target source or evaluation point does not match its definition.")
    }
    profile_pair_curvature_v2(asset$H_population, ctl)
    if (any(!is.finite(asset$Sigma_population)) ||
        max(abs(asset$H_population %*% asset$Omega_population - diag(length(expected_names)))) > 1e-7 ||
        asset$covariance_identity_error > ctl$identity_tol ||
        any(!is.finite(asset$sigma_population)) || any(asset$sigma_population <= 0)) {
      stop("Paired population direction or covariance failed acceptance.")
    }
    covariance_error <- max(abs(asset$Sigma_population -
      (4 * asset$Sigma_half + asset$Sigma_full -
         2 * asset$cross_covariance_half_full - 2 * t(asset$cross_covariance_half_full))))
    if (!is.finite(covariance_error) || covariance_error > ctl$identity_tol ||
        !identical(asset$sigma0_pop, asset$sigma_population) ||
        any(asset$seeds %in% asset$target_asset_identity$seeds)) {
      stop("Paired direction covariance or target-sample independence failed.")
    }
    if (max(abs(asset$Sigma_population - t(asset$Sigma_population))) > ctl$identity_tol ||
        min(eigen(asset$Sigma_population, symmetric = TRUE, only.values = TRUE)$values) < -ctl$identity_tol) {
      stop("Paired population covariance must be symmetric positive semidefinite.")
    }
    for (nm in names(asset$directions)) {
      direction <- asset$directions[[nm]]
      k <- match(nm, expected_names)
      if (is.na(k) || !identical(names(direction$omega), expected_names) ||
          any(!is.finite(direction$omega)) ||
          max(abs(asset$H_population %*% direction$omega -
                    replace(numeric(length(expected_names)), k, 1))) > 1e-7) {
        stop("A paired population direction has inconsistent coordinates or inverse defect.")
      }
    }
    for (rr in seq_along(asset$repeat_assets)) {
      record <- asset$repeat_assets[[rr]]
      if (!identical(as.integer(record$seed), as.integer(asset$seeds[rr])) ||
          record$n_clusters != asset$n_population ||
          !identical(record$target_source, asset$target_source) ||
          !identical(record$evaluation_point, asset$beta_target) ||
          (identical(target_mode, "structural") &&
             !identical(record$generated_beta0, asset$beta_target))) {
        stop("Direction repeat identity or population size is inconsistent.")
      }
      profile_pair_curvature_v2(record$H_population, ctl)
      if (!isTRUE(record$nuisance_converged) ||
          record$max_nuisance_gradient > ctl$nuisance_tol ||
          record$profile_identity_error > ctl$identity_tol) stop("Direction repeat diagnostics failed.")
      if (identical(ctl$population_engine, "conditional_gaussian") &&
          (!identical(record$population_data_definition, "conditional_gaussian_X_integrated") ||
           !identical(record$conditional_integration$score_mean_definition, "conditional_expectation_exactly_zero") ||
           any(record$score_mean != 0))) {
        stop("Conditional Gaussian population direction identity is incomplete.")
      }
    }
  }
  invisible(list(valid = TRUE, problems = character()))
}

profile_pair_target_asset_path_v2 <- function(root, experiment_id) {
  file.path(root, "results", "intermediate", "profile_pair_20261022", "construction_v3", "targets",
            paste0(experiment_id, ".rds"))
}

profile_pair_direction_asset_path_v2 <- function(root, experiment_id) {
  file.path(root, "results", "intermediate", "profile_pair_20261022", "construction_v3", "directions",
            paste0(experiment_id, ".rds"))
}

validate_profile_pair_target_asset_v2 <- function(asset, config = NULL, final = TRUE,
                                                  expected_commit = NULL) {
  validate_profile_pair_asset_v2(asset, "paired_profile_target", config, final, expected_commit)
}

validate_profile_pair_direction_asset_v2 <- function(asset, config = NULL, final = TRUE,
                                                     expected_commit = NULL) {
  validate_profile_pair_asset_v2(asset, "paired_population_direction", config, final, expected_commit)
}

profile_pair_target_scores_v2 <- function(dat, target_asset, h_analysis = target_asset$h_analysis,
                                          nuisance_control = list(grad_tol = 1e-7)) {
  validate_profile_pair_target_asset_v2(target_asset, final = FALSE)
  dat <- fit_working_nuisance_design_v2(dat, target_asset$config$fit_random_effects)
  if (!identical(colnames(dat$X), target_asset$target_names)) {
    stop("Target-score diagnostics require the same full fixed-effect design.")
  }
  comp <- profile_bandwidth_pair_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, target_asset$beta_star_mc,
    target_asset$config$tau, h_analysis, target_asset$config$lambda_gamma,
    nuisance_control = nuisance_control, identity_check = TRUE, need_fold_sums = FALSE
  )
  list(cluster_scores = comp$cluster_scores, score = comp$score, hessian = comp$hessian,
       score_definition = comp$score_definition, bandwidths = comp$bandwidths,
       max_nuisance_gradient = comp$max_nuisance_gradient,
       profile_identity_error = comp$profile_identity_error)
}

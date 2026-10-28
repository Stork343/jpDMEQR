profile_pair_fit_from_data_v2 <- function(train, tau, target_coords, tuning, control) {
  cache <- control$shared_first_stage
  request <- list(
    y = train$y, X = train$X, Z = train$Z, cluster_id = train$cluster_id,
    tau = tau, target_coords = target_coords, h_inf = tuning$h_inf %||% tuning$h,
    h_est = tuning$h_est, lambda_gamma = tuning$lambda_gamma,
    mean_score_multiplier = tuning$mean_score_multiplier %||% 1,
    lambda_multiplier = tuning$lambda_multiplier %||% 1,
    penalty_factor = control$penalty_factor,
    fit_control = control$fit_control, refit_control = control$refit_control,
    numerical_configuration = solver_configuration_v2()
  )
  if (!is.null(cache) && !is.environment(cache)) stop("Shared fit cache must be an environment.")
  if (!is.null(cache) && exists("paired_request", envir = cache, inherits = FALSE) &&
      identical(cache$paired_request, request)) {
    result <- cache$paired_fit
    result$diagnostics$cache <- list(reused = TRUE,
                                     original_first_stage_sec = result$first_stage_sec,
                                     original_refit_sec = result$refit_sec,
                                     original_reprofile_sec = result$paired_reprofile_sec)
    result$first_stage_sec <- result$refit_sec <- result$paired_reprofile_sec <- 0
    return(result)
  }
  result <- fit_profile_pair_v2(
    y = train$y, X = train$X, Z = train$Z,
    cluster_id = train$cluster_id, tau = tau,
    h_inf = tuning$h_inf %||% tuning$h, h_est = tuning$h_est,
    lambda_gamma = tuning$lambda_gamma,
    mean_score_multiplier = tuning$mean_score_multiplier %||% 1,
    lambda_multiplier = tuning$lambda_multiplier %||% 1,
    penalty_factor = control$penalty_factor,
    target_coordinates = target_coords,
    control = list(
      fit_control = control$fit_control %||% list(),
      refit_control = control$refit_control %||% list(),
      nuisance_control = control$fit_control$nuisance_control %||% list()
    )
  )
  result$diagnostics$cache <- list(reused = FALSE)
  if (!is.null(cache)) {
    cache$paired_request <- request
    cache$paired_fit <- result
  }
  result
}

profile_pair_prediction_v2 <- function(beta) {
  force(beta)
  function(newX) {
    newX <- as.matrix(newX)
    if (is.null(colnames(newX)) || !all(names(beta) %in% colnames(newX))) {
      stop("Prediction design does not contain every fitted coordinate.")
    }
    as.numeric(newX[, names(beta), drop = FALSE] %*% beta)
  }
}

profile_pair_benchmark_output_v2 <- function(fit, inf, method_id, elapsed) {
  tab <- inf$table
  if (!nrow(tab) || any(!tab$feasible) || any(!is.finite(tab$se))) {
    stop("Paired profile inference requires a valid interval for every target.")
  }
  output <- list(
    method_id = method_id, status = "ok", failure_stage = "", failure_message = "",
    beta_hat = fit$beta,
    beta_tilde = setNames(tab$beta_tilde, tab$coordinate),
    se = setNames(tab$se, tab$coordinate),
    ci_lower = setNames(tab$ci_lower, tab$coordinate),
    ci_upper = setNames(tab$ci_upper, tab$coordinate),
    selected = fit$selected_support,
    fit_object = fit, inference_object = inf,
    prediction_function = profile_pair_prediction_v2(fit$beta),
    runtime_sec = elapsed, converged = TRUE,
    kkt_residual = fit$kkt_residual,
    profile_identity_error = fit$components$profile_identity_error,
    max_nuisance_gradient = fit$components$max_nuisance_gradient,
    warning_messages = character(),
    implementation_version = "profile_pair_20261022",
    method_version = "profile_pair_20261022",
    score_definition = "profile_bandwidth_pair_2_half_minus_full",
    lambda_stochastic = fit$lambda_components$stochastic,
    lambda_mean = fit$lambda_components$mean,
    lambda_mean_multiplier = fit$lambda_components$mean_score_multiplier %||%
      fit$mean_score_multiplier,
    reference_identifier = "docs/METHOD_SPECIFICATION.md",
    adapter_fidelity_status = "paired_profile_checks_required"
  )
  fields <- c(
    "lambda_rule", "lambda_alpha", "lambda_normal_quantile", "lambda_safety_constant",
    "lambda_base", "lambda_loading_pass0_min", "lambda_loading_pass0_median",
    "lambda_loading_pass0_max", "lambda_loading_pass1_min", "lambda_loading_pass1_median",
    "lambda_loading_pass1_max", "lambda_coordinate_min", "lambda_coordinate_median",
    "lambda_coordinate_max", "zero_profile_score_max", "zero_kkt_ratio_max",
    "preliminary_kkt_normalized", "final_kkt_normalized", "final_kkt_absolute",
    "first_stage_iterations", "first_stage_beta_change", "first_stage_nonzero_count",
    "first_stage_sec", "refit_sec", "nuisance_iterations", "selected_support_size",
    "refit_set_size", "refit_contains_targets", "post_refit_status",
    "post_refit_iterations", "post_refit_gradient_max", "post_refit_nuisance_gradient_max",
    "post_refit_beta_change", "post_refit_hessian_min_eigenvalue",
    "post_refit_hessian_condition_number"
  )
  output[fields] <- lapply(fields, function(nm) fit[[nm]] %||% NA_real_)
  output
}

fit_benchmark_profile_pair_v2 <- function(train, tau, target_coords,
                                          tuning, seed, control = list()) {
  set.seed(seed)
  start <- proc.time()[[3]]
  fit <- profile_pair_fit_from_data_v2(train, tau, target_coords, tuning, control)
  inference_start <- proc.time()[[3]]
  inf <- infer_profile_pair_v2(
    fit, target_coords, tuning$mu_grid,
    ci_level = control$ci_level %||% 0.95,
    solver_preference = control$solver_preference %||% "CLARABEL"
  )
  output <- profile_pair_benchmark_output_v2(
    fit, inf, "PROFILE-DQR", proc.time()[[3]] - start
  )
  output$inference_sec <- proc.time()[[3]] - inference_start
  output
}

profile_pair_fixed_direction_inference_v2 <- function(fit, target_coords,
                                                       precision, source,
                                                       ci_level = 0.95) {
  coordinates <- if (is.character(target_coords)) {
    match(target_coords, names(fit$beta))
  } else as.integer(target_coords)
  if (anyNA(coordinates)) stop("Unknown paired-profile inference coordinate.")
  tab <- profile_bandwidth_pair_wald_v2(
    fit$components, fit$beta, precision, coordinates, alpha = 1 - ci_level
  )
  directions <- vector("list", length(coordinates))
  names(directions) <- tab$coordinate
  tab$beta_hat <- tab$beta_start
  tab$feasible <- TRUE
  tab$index <- coordinates
  tab$mu <- NA_real_
  tab$mu_source <- source
  tab$mu_cv_min_defect <- NA_real_
  tab$mu_cv_defect <- NA_real_
  tab$mu_cv_quad <- NA_real_
  tab$mu_cv_folds_used <- NA_integer_
  tab$dantzig_residual <- tab$inverse_defect
  tab$omega_l1 <- colSums(abs(precision))
  tab$omega_l2 <- sqrt(colSums(precision^2))
  tab$adjacent_stability <- NA_real_
  tab$solver <- source
  tab$method_version <- "profile_pair_20261022"
  tab$se_smoothed <- NA_real_
  for (i in seq_along(coordinates)) {
    omega <- precision[, i]
    tab$se_smoothed[i] <- cluster_sandwich_coordinate_v2(
      fit$components$components_h$cluster_scores, omega
    )$se
    directions[[i]] <- list(omega = omega, feasible = TRUE,
                            residual = tab$inverse_defect[i],
                            solver = source, mu = NA_real_)
  }
  list(table = tab, directions = directions, fit = fit,
       score_definition = fit$score_definition, method_version = fit$method_version)
}

fit_benchmark_profile_pair_pop_h_v2 <- function(train, tau, target_coords,
                                                tuning, seed, asset,
                                                control = list(), config = NULL) {
  if (is.null(asset)) stop("A paired population-direction asset is required.")
  if (is.null(config)) stop("A paired population direction requires context$config.")
  if (is.null(control$implementation_commit)) {
    stop("A paired population direction requires an implementation commit identity.")
  }
  validate_profile_pair_direction_asset_v2(
    asset, config = config, final = isTRUE(control$final_assets %||% TRUE),
    expected_commit = control$implementation_commit
  )
  set.seed(seed)
  start <- proc.time()[[3]]
  fit <- profile_pair_fit_from_data_v2(train, tau, target_coords, tuning, control)
  if (abs(asset$h_analysis - fit$h_inf) > 1e-12 * max(1, fit$h_inf)) {
    stop("Population direction and fitted paired bandwidth differ.")
  }
  coords <- if (is.character(target_coords)) target_coords else names(fit$beta)[target_coords]
  precision <- vapply(coords, function(nm) {
    omega <- asset$directions[[nm]]$omega
    if (is.null(omega) || !all(names(fit$beta) %in% names(omega))) {
      stop("Paired population direction does not cover the fitted design.")
    }
    as.numeric(omega[names(fit$beta)])
  }, numeric(length(fit$beta)))
  if (is.null(dim(precision))) precision <- matrix(precision, ncol = 1L)
  inf <- profile_pair_fixed_direction_inference_v2(
    fit, target_coords, precision, "population_paired_direction",
    ci_level = control$ci_level %||% 0.95
  )
  profile_pair_benchmark_output_v2(fit, inf, "PROFILE-DQR-POP-H",
                                   proc.time()[[3]] - start)
}

fit_benchmark_profile_pair_exact_h_v2 <- function(train, tau, target_coords,
                                                  tuning, seed, control = list()) {
  set.seed(seed)
  start <- proc.time()[[3]]
  fit <- profile_pair_fit_from_data_v2(train, tau, target_coords, tuning, control)
  coords <- if (is.character(target_coords)) match(target_coords, names(fit$beta)) else
    as.integer(target_coords)
  selected <- fit$refit_set
  if (is.character(selected)) selected <- match(selected, names(fit$beta))
  if (anyNA(coords) || anyNA(selected) || !all(coords %in% selected)) {
    stop("The paired refit design must contain every inference target.")
  }
  J <- fit$components$hessian[selected, selected, drop = FALSE]
  chol(J)
  canonical <- diag(length(selected))[, match(coords, selected), drop = FALSE]
  precision <- matrix(0, nrow = length(fit$beta), ncol = length(coords))
  precision[selected, ] <- solve(J, canonical)
  inf <- profile_pair_fixed_direction_inference_v2(
    fit, target_coords, precision, "selected_refit_paired_inverse",
    ci_level = control$ci_level %||% 0.95
  )
  profile_pair_benchmark_output_v2(fit, inf, "POSTREFIT-EXACT-H",
                                   proc.time()[[3]] - start)
}

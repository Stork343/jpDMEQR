conditional_scale_target_source_paths_v2 <- function(kind) {
  if (!kind %in% c("paired_conditional_scale_target", "sqr_marginal_qr_target")) stop("Unknown conditional target type.")
  unique(c(if (kind == "paired_conditional_scale_target") profile_pair_population_source_paths_v2("streamed") else
    c("R/v2_profile_v2_part01.R", paste0("R/v2_simulation_v2_part", sprintf("%02d", c(1, 2, 4, 5)), ".R")),
    "R/v2_simulation_v2_part07.R", "R/v2_simulation_v2_part08.R"))
}

conditional_scale_target_assert_build_source_v2 <- function(commit, kind,
    root = getOption("jpDMEQR.source_root", ".")) {
  paths <- conditional_scale_target_source_paths_v2(kind)
  tracked <- system2("git", c("-C", shQuote(root), "ls-files", "--error-unmatch", "--", shQuote(paths)), stdout = FALSE, stderr = FALSE)
  unchanged <- system2("git", c("-C", shQuote(root), "diff", "--quiet", commit, "--", shQuote(paths)), stdout = FALSE, stderr = FALSE)
  if (tracked != 0L || unchanged != 0L) stop("Final conditional target construction requires its committed source dependencies.")
  invisible(TRUE)
}

conditional_scale_target_commit_matches_v2 <- function(asset_commit, expected_commit, kind,
                                                       root = getOption("jpDMEQR.source_root", ".")) {
  if (identical(asset_commit, expected_commit)) return(TRUE)
  if (!all(grepl("^[0-9a-fA-F]{40}$", c(asset_commit, expected_commit)))) return(FALSE)
  ancestor <- system2("git", c("-C", shQuote(root), "merge-base", "--is-ancestor", asset_commit, expected_commit), stdout = FALSE, stderr = FALSE)
  if (ancestor != 0L) return(FALSE)
  unchanged <- system2("git", c("-C", shQuote(root), "diff", "--quiet", asset_commit, expected_commit,
    "--", shQuote(conditional_scale_target_source_paths_v2(kind))), stdout = FALSE, stderr = FALSE)
  identical(unchanged, 0L)
}

conditional_scale_target_conditions_v2 <- function(config, marginal = FALSE) {
  required <- c("experiment_id", "p", "s", "tau", "q", "n_clusters", "h_multiplier", "target_mode",
    "design_type", "rho_x", "heteroskedastic", "error_dependence", "x_b_corr", "informative_size",
    "nonlinear_re_strength", "fit_random_effects", "error_dist", "random_effect_dist", "lambda_gamma")
  if (!all(required %in% names(config))) stop("The complete conditional target configuration is required.")
  valid <- config$p >= config$s && config$s >= 1L && config$p == as.integer(config$p) &&
    config$s == as.integer(config$s) && config$tau > 0 && config$tau < 1 &&
    identical(config$target_mode, "profile_mc") && identical(config$design_type, "ar1") &&
    is.finite(config$rho_x) && abs(config$rho_x) < 1 && isTRUE(config$heteroskedastic) &&
    identical(config$error_dependence, "independent") && config$x_b_corr == 0 &&
    !isTRUE(config$informative_size) && config$nonlinear_re_strength == 0 &&
    config$error_dist %in% c("gaussian", "t3") &&
    ((config$q == 1L && config$fit_random_effects == "intercept") ||
      (config$q == 2L && config$fit_random_effects == "intercept+slope")) &&
    is.null(config$m_multiset)
  if (!isTRUE(valid)) stop("The conditional Gaussian scale target conditions are not satisfied.")
  if (marginal && (config$q != 1L || config$random_effect_dist != "normal")) {
    stop("Marginal scale quadrature requires the registered Normal random intercept.")
  }
  invisible(TRUE)
}

conditional_scale_target_identity_v2 <- function(config) {
  list(design = profile_pair_design_identity_v2(config),
    scale_definition = "exp(0.35*clip(X1,-3,3))", design_rows = "independent_Gaussian_AR1",
    integrated_coordinates = if (config$p > 1L) 2:config$p else integer(),
    signal_support = seq_len(config$s), finite_full_design_empirical_root = FALSE)
}

conditional_scale_target_native_v2 <- function() {
  old <- options(jpDMEQR.population_nuisance_backend = "compiled",
    jpDMEQR.population_matrix_backend = "Eigen_single_thread")
  profile_pair_population_eigen_init_v2(getOption("jpDMEQR.source_root", "."))
  old
}

conditional_scale_target_data_v2 <- function(config, n_population, seed) {
  conditional_scale_target_conditions_v2(config)
  cfg <- config
  cfg$n_clusters <- as.integer(n_population)
  cfg$p <- min(config$p, config$s + 3L)
  dat <- simulate_from_config_v2(cfg, seed)
  full_beta <- make_beta_v2(config$p, config$s, config$signal)
  if (!identical(dat$active, seq_len(config$s)) || any(full_beta[-seq_len(cfg$p)] != 0) ||
      !identical(as.numeric(dat$beta0), as.numeric(full_beta[seq_len(cfg$p)]))) stop("The complete generating signal was not retained.")
  scalar <- list(y = dat$y - drop(dat$X %*% dat$beta0), X = dat$X[, 1L, drop = FALSE],
    Z = dat$Z, cluster_id = dat$cluster_id, beta0 = setNames(0, "x00001"), tau = config$tau)
  scalar$population_layout <- profile_pair_population_layout_v2(scalar)
  list(scalar = scalar, full_beta0 = full_beta, generated_prefix_p = cfg$p,
    generated_support = dat$active, m_summary = dat$m_summary, total_rows = dat$total_rows)
}

conditional_scale_profile_geometry_v2 <- function(component, delta, scalar, config, h) {
  layout <- scalar$population_layout
  Lambda <- make_lambda_matrix_v2(config$lambda_gamma, ncol(scalar$Z))
  residual <- scalar$y - drop(scalar$X * delta)
  half <- profile_pair_population_nuisance_v2(residual, scalar$Z, layout, config$tau, h / 2, Lambda, 1e-9)
  full <- profile_pair_population_nuisance_v2(residual, scalar$Z, layout, config$tau, h, Lambda, 1e-9)
  ch <- profile_pair_conditional_hessian_coefficients_v2(half, scalar$Z, layout, Lambda)
  cf <- profile_pair_conditional_hessian_coefficients_v2(full, scalar$Z, layout, Lambda)
  coefficients <- c(alpha = component$hessian[1L, 1L], residual_trace = mean(2 * ch$coefficient - cf$coefficient))
  identity_coefficients <- c(alpha = component$hessian_identity[1L, 1L], residual_trace = mean(2 * ch$identity - cf$identity))
  a <- setNames(config$rho_x^(seq_len(config$p) - 1L), names(make_beta_v2(config$p, config$s, config$signal)))
  C <- ar1_covariance_v2(config$p, config$rho_x) - tcrossprod(a)
  difference <- coefficients - identity_coefficients
  identity_error <- max(abs(difference[1L] * tcrossprod(a) + difference[2L] * C))
  lower <- min(coefficients) * (1 - abs(config$rho_x)) / (1 + abs(config$rho_x))
  condition <- max(coefficients) / min(coefficients) * ((1 + abs(config$rho_x)) / (1 - abs(config$rho_x)))^2
  inv_m2 <- 1 / layout$sizes[layout$row_cluster]^2
  traces <- c(pair = sum((2 * half$psi - full$psi)^2 * inv_m2) / layout$n,
    half = sum(half$psi^2 * inv_m2) / layout$n, full = sum(full$psi^2 * inv_m2) / layout$n,
    cross = sum(half$psi * full$psi * inv_m2) / layout$n)
  list(delta = as.numeric(delta), n_clusters = layout$n, h = h,
    full_score = a * component$score[1L], scalar_score = as.numeric(component$score[1L]),
    hessian_coefficients = coefficients, hessian_identity_coefficients = identity_coefficients,
    profile_identity_error = identity_error, minimum_eigenvalue_lower_bound = lower,
    condition_number_upper_bound = condition,
    covariance_coefficients = c(var_g1 = component$population_moments$Sigma_population[1L, 1L], residual_trace = unname(traces[1L])),
    covariance_trace_components = traces,
    covariance_identity_error = abs(traces[1L] - (4 * traces[2L] + traces[3L] - 4 * traces[4L])),
    nuisance_converged = component$nuisance_converged,
    max_nuisance_gradient = max(component$max_nuisance_gradient, half$max_nuisance_gradient, full$max_nuisance_gradient))
}

conditional_scale_profile_fit_record_v2 <- function(fit, beta0, geometry) {
  delta <- as.numeric(fit$beta)
  value <- profile_pair_fit_record_v2(fit)
  value$scalar_beta <- delta
  value$scalar_moments <- value$moments
  value$moments <- NULL
  value$beta <- beta0
  value$beta[1L] <- value$beta[1L] + delta
  value$conditional_geometry <- geometry
  value$kkt_residual <- abs(geometry$scalar_score)
  value$profile_identity_error <- geometry$profile_identity_error
  value$max_nuisance_gradient <- geometry$max_nuisance_gradient
  value$population_engine <- "conditional_gaussian_scale_v1"
  value
}

build_profile_pair_conditional_scale_target_v2 <- function(config, n_population = 100000L,
    repeats = 4L, seed = 87001L, h_target = NULL, development = FALSE,
    implementation_commit = "working_tree", control = list(), checkpoint = NULL) {
  conditional_scale_target_conditions_v2(config)
  implementation_commit <- profile_pair_asset_commit_v2(implementation_commit, development)
  profile_pair_asset_scale_v2(n_population, repeats, development)
  if (!development) conditional_scale_target_assert_build_source_v2(implementation_commit, "paired_conditional_scale_target")
  ctl <- profile_pair_asset_control_v2(utils::modifyList(list(population_engine = "streamed", nuisance_tol = 1e-9), control))
  if (ctl$population_engine != "streamed") stop("Conditional target scalar roots use the streamed complete profile.")
  old <- conditional_scale_target_native_v2()
  on.exit(options(old), add = TRUE)
  h_analysis <- config$h_multiplier * config$n_clusters^(-3 / 10)
  if (is.null(h_target)) h_target <- min(h_analysis / 4, n_population^(-1 / 3))
  if (length(h_target) != 1L || !is.finite(h_target) || h_target <= 0 || h_target >= h_analysis) stop("Invalid conditional target bandwidth.")
  seeds <- vapply(seq_len(repeats), function(r) seed_from_id_v2(seed, config$experiment_id, r, "paired_population_target"), integer(1))
  records <- vector("list", repeats)
  for (r in seq_len(repeats)) {
    generated <- conditional_scale_target_data_v2(config, n_population, seeds[r])
    scalar <- generated$scalar
    base <- solve_profile_pair_target_v2(scalar, config$tau, h_target, config$lambda_gamma, beta_start = 0, control = ctl)
    profile_pair_checkpoint_v2(checkpoint, sprintf("conditional_scalar_%02d_target", r), profile_pair_fit_record_v2(base))
    refinement <- solve_profile_pair_target_v2(scalar, config$tau, h_target / 2, config$lambda_gamma, beta_start = base$beta, control = ctl)
    profile_pair_checkpoint_v2(checkpoint, sprintf("conditional_scalar_%02d_refinement", r), profile_pair_fit_record_v2(refinement))
    analysis <- solve_profile_pair_target_v2(scalar, config$tau, h_analysis, config$lambda_gamma, beta_start = base$beta, control = ctl)
    profile_pair_checkpoint_v2(checkpoint, sprintf("conditional_scalar_%02d_analysis", r), profile_pair_fit_record_v2(analysis))
    at_analysis <- profile_pair_population_evaluate_v2(scalar, base$beta, config$tau, h_analysis, config$lambda_gamma, ctl, 2L)
    make_record <- function(fit, h) conditional_scale_profile_fit_record_v2(fit, generated$full_beta0,
      conditional_scale_profile_geometry_v2(fit$components, fit$beta, scalar, config, h))
    records[[r]] <- list(seed = seeds[r], n_population = as.integer(n_population),
      target = make_record(base, h_target), refinement = make_record(refinement, h_target / 2),
      analysis = make_record(analysis, h_analysis),
      analysis_score_at_target = conditional_scale_profile_geometry_v2(at_analysis, base$beta, scalar, config, h_analysis),
      computational_prefix_p = generated$generated_prefix_p, generated_complete_support = generated$generated_support,
      m_summary = generated$m_summary, total_rows = generated$total_rows)
    profile_pair_checkpoint_v2(checkpoint, sprintf("conditional_target_repeat_%02d", r), records[[r]])
  }
  estimates <- do.call(rbind, lapply(records, function(x) x$target$beta))
  refined <- do.call(rbind, lapply(records, function(x) x$refinement$beta))
  analyses <- do.call(rbind, lapply(records, function(x) x$analysis$beta))
  beta0 <- make_beta_v2(config$p, config$s, config$signal)
  asset <- list(asset_kind = "paired_conditional_scale_target", population_engine = "conditional_gaussian_scale_v1",
    construction_version = "conditional-target-assets-20261023-v4", method_version = "profile_pair_20261022",
    score_definition = "profile_bandwidth_pair_2_half_minus_full", development = development,
    experiment_id = as.character(config$experiment_id), config = config,
    design_identity = conditional_scale_target_identity_v2(config), target_scope = "full_registered_fixed_effect_design",
    target_source = "conditional_expectation_integration_mc", finite_full_design_empirical_root = FALSE,
    target_columns = seq_len(config$p), target_names = names(beta0),
    beta_star_mc = colMeans(estimates), evaluation_point = colMeans(estimates), estimates_by_repeat = estimates,
    target_mc_sd_by_coordinate = apply(estimates, 2L, stats::sd),
    target_mc_se_by_coordinate = apply(estimates, 2L, stats::sd) / sqrt(repeats),
    target_mc_se = max(apply(estimates, 2L, stats::sd) / sqrt(repeats)),
    between_run_range_by_coordinate = apply(estimates, 2L, function(x) diff(range(x))),
    between_run_max_difference = max(apply(estimates, 2L, function(x) diff(range(x)))),
    beta_refinement_mc = colMeans(refined), refinement_estimates_by_repeat = refined,
    refinement_mc_se_by_coordinate = apply(refined, 2L, stats::sd) / sqrt(repeats),
    beta_analysis_bandwidth_mc = colMeans(analyses), analysis_estimates_by_repeat = analyses,
    target_bandwidth_displacement_by_repeat = refined - estimates,
    smoothing_target_displacement_by_repeat = analyses - estimates,
    smoothing_target_distance = mean(sqrt(rowSums((analyses - estimates)^2))),
    smoothing_target_distance_l1 = mean(rowSums(abs(analyses - estimates))),
    smoothing_target_distance_by_repeat = sqrt(rowSums((analyses - estimates)^2)),
    target_refinement_distance_by_repeat = sqrt(rowSums((refined - estimates)^2)),
    h_target = h_target, h_analysis = h_analysis, n_analysis = as.integer(config$n_clusters),
    bandwidths_target = c(h_target / 2, h_target), bandwidths_refinement = c(h_target / 4, h_target / 2),
    bandwidths_analysis = c(h_analysis / 2, h_analysis),
    n_population = as.integer(n_population), repeats = as.integer(repeats), seeds = seeds, repeat_assets = records,
    all_converged = TRUE, acceptance = ctl, implementation_commit = implementation_commit,
    max_kkt_residual = max(vapply(records, function(record) max(vapply(record[c("target", "refinement", "analysis")], function(fit) fit$kkt_residual, numeric(1))), numeric(1))),
    max_nuisance_gradient = max(vapply(records, function(record) max(c(vapply(record[c("target", "refinement", "analysis")], function(fit) fit$max_nuisance_gradient, numeric(1)), record$analysis_score_at_target$max_nuisance_gradient)), numeric(1))),
    profile_identity_error = max(vapply(records, function(record) max(c(vapply(record[c("target", "refinement", "analysis")], function(fit) fit$profile_identity_error, numeric(1)), record$analysis_score_at_target$profile_identity_error)), numeric(1))),
    source_paths = conditional_scale_target_source_paths_v2("paired_conditional_scale_target"),
    integration_definition = "X2 through Xp integrated conditionally on all cluster X1; beta_star=beta0+delta*e1",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  validate_profile_pair_conditional_scale_target_v2(asset, config, final = !development)
  asset
}

conditional_scale_target_validate_common_v2 <- function(asset, kind, config, final, expected_commit) {
  if (!is.logical(final) || length(final) != 1L || is.na(final)) stop("final must be one logical value.")
  if (!identical(asset$asset_kind, kind) || !identical(asset$construction_version, "conditional-target-assets-20261023-v4")) stop("Conditional target type or construction version mismatch.")
  profile_pair_asset_commit_v2(asset$implementation_commit, !final)
  profile_pair_asset_scale_v2(asset$n_population, asset$repeats, asset$development)
  if (final && (isTRUE(asset$development) || asset$n_population < 100000L || asset$repeats < 4L)) stop("Development conditional targets cannot pass the final gate.")
  conditional_scale_target_conditions_v2(asset$config, marginal = kind == "sqr_marginal_qr_target")
  if (!is.null(config) && !identical(asset$design_identity, conditional_scale_target_identity_v2(config))) stop("Conditional target scenario identity mismatch.")
  if (!identical(asset$design_identity, conditional_scale_target_identity_v2(asset$config)) ||
      !identical(asset$source_paths, conditional_scale_target_source_paths_v2(kind))) stop("Conditional target design or source identity mismatch.")
  if (!is.null(expected_commit) && !conditional_scale_target_commit_matches_v2(asset$implementation_commit, expected_commit, kind)) stop("Conditional target implementation commit mismatch.")
  expected <- sprintf("x%05d", seq_len(asset$config$p))
  if (!identical(asset$target_columns, seq_len(asset$config$p)) || !identical(asset$target_names, expected) ||
      !identical(names(asset$beta_star_mc), expected) || any(!is.finite(asset$beta_star_mc)) ||
      length(asset$seeds) != asset$repeats || anyDuplicated(asset$seeds) || length(asset$repeat_assets) != asset$repeats ||
      !identical(asset$finite_full_design_empirical_root, FALSE)) stop("Conditional target full-coordinate or independent-repeat records are incomplete.")
  if (length(asset$target_mc_se_by_coordinate) != length(expected) || any(!is.finite(asset$target_mc_se_by_coordinate)) ||
      any(asset$target_mc_se_by_coordinate < 0)) stop("Conditional target uncertainty is invalid.")
  invisible(TRUE)
}

validate_profile_pair_conditional_scale_target_v2 <- function(asset, config = NULL, final = TRUE, expected_commit = NULL) {
  conditional_scale_target_validate_common_v2(asset, "paired_conditional_scale_target", config, final, expected_commit)
  if (!identical(asset$population_engine, "conditional_gaussian_scale_v1") ||
      !identical(asset$method_version, "profile_pair_20261022") ||
      !identical(asset$score_definition, "profile_bandwidth_pair_2_half_minus_full") ||
      !identical(asset$target_source, "conditional_expectation_integration_mc") ||
      !identical(asset$evaluation_point, asset$beta_star_mc) || !isTRUE(asset$all_converged)) stop("Conditional profile target definition mismatch.")
  ctl <- profile_pair_asset_control_v2(asset$acceptance)
  cfg <- asset$config
  beta0 <- make_beta_v2(cfg$p, cfg$s, cfg$signal)
  a <- cfg$rho_x^(seq_len(cfg$p) - 1L)
  C <- ar1_covariance_v2(cfg$p, cfg$rho_x) - tcrossprod(a)
  expected_h <- cfg$h_multiplier * cfg$n_clusters^(-3 / 10)
  if (!isTRUE(all.equal(asset$h_analysis, expected_h, tolerance = 1e-14)) ||
      !identical(asset$bandwidths_target, c(asset$h_target / 2, asset$h_target)) ||
      !identical(asset$bandwidths_refinement, c(asset$h_target / 4, asset$h_target / 2)) ||
      !identical(asset$bandwidths_analysis, c(expected_h / 2, expected_h)) || asset$h_target <= 0 || asset$h_target >= expected_h) stop("Conditional target bandwidth identity mismatch.")
  matrices <- list(target = asset$estimates_by_repeat, refinement = asset$refinement_estimates_by_repeat, analysis = asset$analysis_estimates_by_repeat)
  for (role in names(matrices)) {
    M <- matrices[[role]]
    if (!is.matrix(M) || !identical(dim(M), c(asset$repeats, cfg$p)) ||
        !identical(colnames(M), names(beta0)) || any(!is.finite(M))) stop("Conditional target repeat matrix is invalid.")
    for (r in seq_len(asset$repeats)) {
      record <- asset$repeat_assets[[r]]
      fit <- record[[role]]
      g <- fit$conditional_geometry
      expected_beta <- beta0
      expected_beta[1L] <- expected_beta[1L] + fit$scalar_beta
      h <- switch(role, target = asset$h_target, refinement = asset$h_target / 2, analysis = asset$h_analysis)
      difference <- g$hessian_coefficients - g$hessian_identity_coefficients
      identity <- max(abs(difference[1L] * tcrossprod(a) + difference[2L] * C))
      lower <- min(g$hessian_coefficients) * (1 - abs(cfg$rho_x)) / (1 + abs(cfg$rho_x))
      condition <- max(g$hessian_coefficients) / min(g$hessian_coefficients) * ((1 + abs(cfg$rho_x)) / (1 - abs(cfg$rho_x)))^2
      traces <- g$covariance_trace_components
      covariance_identity <- abs(traces[1L] - (4 * traces[2L] + traces[3L] - 4 * traces[4L]))
      if (record$seed != asset$seeds[r] || record$n_population != asset$n_population ||
          !identical(record$generated_complete_support, seq_len(cfg$s)) ||
          record$computational_prefix_p < cfg$s || record$computational_prefix_p > cfg$p ||
          !identical(fit$beta, expected_beta) || !isTRUE(all.equal(as.numeric(M[r, ]), as.numeric(fit$beta), tolerance = 0)) ||
          !isTRUE(fit$converged) || !isTRUE(fit$nuisance_converged) || !isTRUE(g$nuisance_converged) ||
          g$n_clusters != asset$n_population || g$h != h ||
          !isTRUE(all.equal(as.numeric(g$full_score), a * g$scalar_score, tolerance = 1e-12)) ||
          abs(g$scalar_score) > ctl$score_tol || g$max_nuisance_gradient > ctl$nuisance_tol ||
          identity > ctl$identity_tol || g$profile_identity_error > ctl$identity_tol ||
          min(g$hessian_coefficients) <= 0 || lower <= ctl$eigen_min ||
          !isTRUE(all.equal(lower, g$minimum_eigenvalue_lower_bound, tolerance = 1e-12)) ||
          !is.finite(condition) || condition >= ctl$condition_max ||
          !isTRUE(all.equal(condition, g$condition_number_upper_bound, tolerance = 1e-12)) ||
          covariance_identity > ctl$identity_tol || g$covariance_identity_error > ctl$identity_tol) stop("Conditional target root, geometry, or nuisance diagnostics failed.")
    }
  }
  for (record in asset$repeat_assets) {
    g <- record$analysis_score_at_target
    if (!isTRUE(g$nuisance_converged) || g$n_clusters != asset$n_population || g$h != asset$h_analysis ||
        any(!is.finite(g$full_score)) || !isTRUE(all.equal(as.numeric(g$full_score), a * g$scalar_score, tolerance = 1e-12)) ||
        g$max_nuisance_gradient > ctl$nuisance_tol || g$profile_identity_error > ctl$identity_tol) stop("Conditional analysis-bandwidth score diagnostics failed.")
  }
  if (!isTRUE(all.equal(colMeans(asset$estimates_by_repeat), asset$beta_star_mc, tolerance = 1e-12)) ||
      !isTRUE(all.equal(apply(asset$estimates_by_repeat, 2L, stats::sd) / sqrt(asset$repeats), asset$target_mc_se_by_coordinate, tolerance = 1e-12)) ||
      !isTRUE(all.equal(asset$target_bandwidth_displacement_by_repeat,
        asset$refinement_estimates_by_repeat - asset$estimates_by_repeat, tolerance = 1e-12))) stop("Conditional target mean, MCSE, or bandwidth displacement is inconsistent.")
  invisible(list(valid = TRUE, problems = character()))
}

profile_pair_conditional_scale_target_asset_path_v2 <- function(root, experiment_id) {
  file.path(root, "results", "intermediate", "profile_pair_20261022", "construction_v4", "targets", paste0(experiment_id, ".rds"))
}

sqr_marginal_scale_target_asset_path_v2 <- function(root, experiment_id) {
  file.path(root, "results", "intermediate", "profile_pair_20261022", "construction_v4", "comparator_targets", "SQR-DEBIASED-IID", paste0(experiment_id, ".rds"))
}

sqr_marginal_scale_quadrature_v2 <- function(config, tolerance = 1e-10) {
  conditional_scale_target_conditions_v2(config, marginal = TRUE)
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 || tolerance > 1e-8) stop("The marginal quadrature tolerance is invalid.")
  quantile <- error_quantile_function_v2(config$tau, config$error_dist, config$tau)
  lower <- error_quantile_function_v2(1e-12, config$error_dist, config$tau)
  upper <- error_quantile_function_v2(1 - 1e-12, config$error_dist, config$tau)
  raw_cdf <- function(value) {
    result <- if (config$error_dist == "gaussian") stats::pnorm(value) else stats::pt(value * sqrt(3), df = 3)
    result[value < lower] <- 0
    result[value >= upper] <- 1
    result
  }
  max_inner_error <- 0
  equation <- function(delta) {
    conditional <- function(x) {
      scale <- exp(.35 * pmax(pmin(x, 3), -3))
      value <- stats::integrate(function(b) raw_cdf((delta * x - b) / scale + quantile) *
        stats::dnorm(b, sd = config$sigma_b0), -Inf, Inf, rel.tol = tolerance / 10,
        abs.tol = tolerance / 100, subdivisions = 1000L)
      max_inner_error <<- max(max_inner_error, value$abs.error)
      x * (value$value - config$tau) * stats::dnorm(x)
    }
    sum(vapply(list(c(-Inf, -3), c(-3, 3), c(3, Inf)), function(interval) stats::integrate(
      function(x) vapply(x, conditional, numeric(1)), interval[1L], interval[2L],
      rel.tol = tolerance, abs.tol = tolerance / 10, subdivisions = 1000L)$value, numeric(1)))
  }
  root <- stats::uniroot(equation, interval = c(-2, 2), tol = 1e-10)
  score <- equation(root$root)
  list(delta = root$root, scalar_score = score, score_residual = abs(score), iterations = root$iter,
    estimated_precision = root$estim.prec, tolerance = tolerance,
    maximum_inner_reported_error = max_inner_error,
    error_probability_clipping = c(1e-12, 1 - 1e-12))
}

build_sqr_marginal_scale_target_v2 <- function(config, n_population = 100000L, repeats = 4L,
    seed = 137001L, development = FALSE, implementation_commit = "working_tree",
    quadrature_tolerance = 1e-10, checkpoint = NULL) {
  conditional_scale_target_conditions_v2(config, marginal = TRUE)
  implementation_commit <- profile_pair_asset_commit_v2(implementation_commit, development)
  profile_pair_asset_scale_v2(n_population, repeats, development)
  if (!development) conditional_scale_target_assert_build_source_v2(implementation_commit, "sqr_marginal_qr_target")
  ordinary <- sqr_marginal_scale_quadrature_v2(config, min(1e-8, quadrature_tolerance * 100))
  strict <- sqr_marginal_scale_quadrature_v2(config, quadrature_tolerance)
  if (abs(ordinary$delta - strict$delta) > 1e-7 || strict$score_residual > 1e-7) stop("The original marginal QR equation failed quadrature agreement.")
  beta0 <- make_beta_v2(config$p, config$s, config$signal)
  beta <- beta0
  beta[1L] <- beta[1L] + strict$delta
  a <- config$rho_x^(seq_len(config$p) - 1L)
  seeds <- vapply(seq_len(repeats), function(r) seed_from_id_v2(seed, config$experiment_id, r, "sqr_marginal_population_score_audit"), integer(1))
  records <- vector("list", repeats)
  for (r in seq_len(repeats)) {
    generated <- conditional_scale_target_data_v2(config, n_population, seeds[r])
    scalar <- generated$scalar
    layout <- scalar$population_layout
    residual <- scalar$y - drop(scalar$X * strict$delta)
    contribution <- -drop(scalar$X) * (config$tau - as.numeric(residual < 0))
    totals <- as.numeric(rowsum(contribution, layout$row_cluster, reorder = FALSE))
    mean_score <- sum(totals) / sum(layout$sizes)
    influence <- (totals - mean_score * layout$sizes) / mean(layout$sizes)
    mcse <- sqrt(mean(influence^2) / layout$n)
    records[[r]] <- list(seed = seeds[r], n_population = as.integer(n_population),
      total_rows = generated$total_rows, generated_complete_support = generated$generated_support,
      computational_prefix_p = generated$generated_prefix_p,
      score_mean = setNames(a * mean_score, names(beta)), score_mcse = setNames(abs(a) * mcse, names(beta)),
      scalar_score_mean = mean_score, scalar_score_mcse = mcse,
      score_weighting = "pooled_visit_ratio_with_independent_cluster_MCSE",
      target = beta, m_summary = generated$m_summary)
    profile_pair_checkpoint_v2(checkpoint, sprintf("sqr_marginal_audit_repeat_%02d", r), records[[r]])
  }
  asset <- list(asset_kind = "sqr_marginal_qr_target", population_engine = "marginal_scale_quadrature_v1",
    construction_version = "conditional-target-assets-20261023-v4", method_version = "source_20261023_CV",
    score_definition = "iid_unsmoothed_quantile_score", development = development,
    experiment_id = as.character(config$experiment_id), config = config,
    design_identity = conditional_scale_target_identity_v2(config),
    target_scope = "own_pooled_visit_marginal_QR_population_projection",
    target_source = "deterministic_conditional_CDF_population_quadrature",
    target_columns = seq_len(config$p), target_names = names(beta), beta_star_mc = beta,
    beta_marginal = beta, beta_target = beta, scalar_delta = strict$delta,
    target_mc_se_by_coordinate = setNames(numeric(config$p), names(beta)), target_mc_se = 0,
    deterministic_target = TRUE, finite_full_design_empirical_root = FALSE,
    population_score_at_target = setNames(a * strict$scalar_score, names(beta)),
    max_kkt_residual = strict$score_residual,
    quadrature = strict, quadrature_comparison = ordinary,
    quadrature_root_difference = abs(strict$delta - ordinary$delta),
    n_population = as.integer(n_population), repeats = as.integer(repeats), seeds = seeds,
    repeat_assets = records, all_converged = TRUE, implementation_commit = implementation_commit,
    source_paths = conditional_scale_target_source_paths_v2("sqr_marginal_qr_target"),
    equation = "E[X1{F_{b+sigma(X1)(U-q_tau)|X1}(delta*X1)-tau}]=0",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  validate_sqr_marginal_scale_target_v2(asset, config, final = !development)
  asset
}

validate_sqr_marginal_scale_target_v2 <- function(asset, config = NULL, final = TRUE, expected_commit = NULL) {
  conditional_scale_target_validate_common_v2(asset, "sqr_marginal_qr_target", config, final, expected_commit)
  if (!identical(asset$population_engine, "marginal_scale_quadrature_v1") ||
      !identical(asset$method_version, "source_20261023_CV") ||
      !identical(asset$score_definition, "iid_unsmoothed_quantile_score") ||
      !identical(asset$target_source, "deterministic_conditional_CDF_population_quadrature") ||
      !identical(asset$deterministic_target, TRUE) || !isTRUE(asset$all_converged) ||
      !identical(asset$beta_marginal, asset$beta_star_mc) || !identical(asset$beta_target, asset$beta_star_mc) ||
      any(asset$target_mc_se_by_coordinate != 0)) stop("The own marginal QR target definition is inconsistent.")
  beta <- make_beta_v2(asset$config$p, asset$config$s, asset$config$signal)
  beta[1L] <- beta[1L] + asset$scalar_delta
  if (!identical(asset$beta_marginal, beta) || !is.finite(asset$quadrature$score_residual) ||
      asset$quadrature$score_residual > 1e-7 || asset$quadrature$tolerance > 1e-10 ||
      !identical(asset$quadrature$error_probability_clipping, c(1e-12, 1 - 1e-12)) ||
      abs(asset$quadrature$delta - asset$quadrature_comparison$delta) > 1e-7 ||
      abs(asset$scalar_delta - asset$quadrature$delta) > 1e-12) stop("The marginal QR population equation or precision failed.")
  a <- asset$config$rho_x^(seq_len(asset$config$p) - 1L)
  if (!identical(names(asset$population_score_at_target), names(beta)) ||
      !isTRUE(all.equal(as.numeric(asset$population_score_at_target), a * asset$quadrature$scalar_score, tolerance = 1e-12)) ||
      max(abs(asset$population_score_at_target)) > 1e-7) stop("The full marginal QR score identity failed.")
  for (r in seq_len(asset$repeats)) {
    record <- asset$repeat_assets[[r]]
    if (record$seed != asset$seeds[r] || record$n_population != asset$n_population ||
        !identical(record$generated_complete_support, seq_len(asset$config$s)) ||
        !identical(record$target, beta) || !is.finite(record$scalar_score_mean) ||
        !is.finite(record$scalar_score_mcse) || record$scalar_score_mcse <= 0 ||
        !isTRUE(all.equal(as.numeric(record$score_mean), a * record$scalar_score_mean, tolerance = 1e-12)) ||
        !isTRUE(all.equal(as.numeric(record$score_mcse), abs(a) * record$scalar_score_mcse, tolerance = 1e-12)) ||
        !identical(record$score_weighting, "pooled_visit_ratio_with_independent_cluster_MCSE")) stop("The independent marginal target audit is incomplete.")
  }
  invisible(list(valid = TRUE, problems = character()))
}

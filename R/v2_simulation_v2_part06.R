profile_pair_gaussian_row_conditions_v2 <- function(config) {
  valid <- identical(config$target_mode, "structural") &&
    identical(config$error_dist, "gaussian") && !isTRUE(config$heteroskedastic) &&
    identical(config$error_dependence, "independent") &&
    !isTRUE(config$informative_size) && identical(as.numeric(config$x_b_corr), 0) &&
    identical(as.numeric(config$nonlinear_re_strength), 0) &&
    config$design_type %in% c("ar1", "independent", "block", "sparse_precision", "factor", "dense_precision") &&
    ((config$q == 1L && identical(config$fit_random_effects, "intercept")) ||
       (config$q == 2L && identical(config$fit_random_effects, "intercept+slope")))
  if (!isTRUE(valid)) stop("Conditional Gaussian population integration conditions are not satisfied.")
  invisible(TRUE)
}

profile_pair_consume_design_rng_v2 <- function(n_rows, p, design_type) {
  consume <- function(count) {
    remaining <- as.double(count)
    while (remaining > 0) {
      take <- min(remaining, 2000000)
      invisible(stats::rnorm(take))
      remaining <- remaining - take
    }
  }
  if (design_type == "block") {
    blocks <- split(seq_len(p), ceiling(seq_len(p) / 25L))
    for (ids in blocks) {
      consume(n_rows)
      consume(as.double(n_rows) * length(ids))
    }
  } else if (design_type == "factor") {
    consume(as.double(n_rows) * 3L)
    consume(as.double(n_rows) * p)
  } else if (design_type == "dense_precision") {
    consume(n_rows)
    consume(as.double(n_rows) * p)
  } else consume(as.double(n_rows) * p)
  invisible(NULL)
}

profile_pair_conditional_gaussian_data_v2 <- function(config, n_population, seed) {
  profile_pair_gaussian_row_conditions_v2(config)
  set.seed(seed)
  p <- as.integer(config$p)
  q <- as.integer(config$q)
  sampler <- make_design_sampler_v2(p, config$design_type, config$rho_x)
  Psi <- make_random_effect_covariance_v2(q, config$sigma_b0, config$sigma_b1, config$rho_b)
  random_effects <- r_random_effects_v2(n_population, Psi, config$random_effect_dist)
  latent_size <- stats::rnorm(n_population)
  sizes <- draw_cluster_sizes_v2(
    n_population, config$m_values, config$m_rule,
    random_intercept = random_effects[, 1] / max(config$sigma_b0, 1e-8),
    latent_size = latent_size, m_multiset = config$m_multiset)
  n_rows <- sum(sizes)
  profile_pair_consume_design_rng_v2(n_rows, p, config$design_type)
  time <- unlist(lapply(sizes, function(m) sort(stats::runif(m, -1, 1))), use.names = FALSE)
  Z <- if (q == 1L) matrix(1, n_rows, 1L) else cbind(1, time)
  cluster_index <- rep(seq_len(n_population), sizes)
  epsilon <- numeric(n_rows)
  starts <- cumsum(c(1L, head(sizes, -1L)))
  for (i in seq_len(n_population)) {
    rows <- starts[i]:(starts[i] + sizes[i] - 1L)
    epsilon[rows] <- r_cluster_quantile_error_v2(
      sizes[i], config$tau, config$error_dist, x1 = rep(0, sizes[i]),
      random_intercept = rep(random_effects[i, 1], sizes[i]),
      heteroskedastic = FALSE, dependence = "independent", copula_rho = config$copula_rho)
  }
  residual_base <- rowSums(Z * random_effects[cluster_index, seq_len(q), drop = FALSE]) + epsilon
  covariance <- sampler$covariance %||%
    make_design_covariance_v2(p, config$design_type, config$rho_x)
  nm <- sprintf("x%05d", seq_len(p))
  dimnames(covariance) <- list(nm, nm)
  cluster_id <- rep(sprintf("c%05d", seq_len(n_population)), sizes)
  dat <- list(Z = Z, cluster_id = cluster_id, beta0 = make_beta_v2(p, config$s, config$signal),
              residual_base = residual_base, Sigma_x = covariance, epsilon = epsilon,
              time = time, random_effects = random_effects, m_i = sizes,
              n_clusters = n_population, total_rows = n_rows, p = p, q = q,
              population_data_definition = "conditional_gaussian_X_integrated",
              seed = seed)
  dat$population_layout <- profile_pair_population_layout_v2(dat)
  dat
}

profile_pair_conditional_hessian_coefficients_v2 <- function(profile, Z, layout, Lambda) {
  n <- layout$n
  q <- ncol(Z)
  coefficient <- identity <- numeric(n)
  for (group in layout$groups) {
    ids <- group$clusters
    b <- length(ids)
    m <- group$m
    rows <- as.vector(group$rows)
    w <- matrix(profile$weights[rows], b, m)
    z1 <- matrix(Z[rows, 1], b, m)
    Q11 <- rowSums(w * z1^2)
    V11 <- rowSums(w^2 * z1^2)
    b11 <- profile$B11[ids]
    if (q == 1L) {
      correction <- V11 / b11
      M11 <- V11 / b11^2
      identity[ids] <- rowSums(w) - 2 * correction + Q11 * M11 + Lambda[1, 1] * M11
    } else {
      z2 <- matrix(Z[rows, 2], b, m)
      Q12 <- rowSums(w * z1 * z2)
      Q22 <- rowSums(w * z2^2)
      V12 <- rowSums(w^2 * z1 * z2)
      V22 <- rowSums(w^2 * z2^2)
      b12 <- profile$B12[ids]
      b22 <- profile$B22[ids]
      determinant <- b11 * b22 - b12^2
      a <- b22 / determinant
      d <- b11 / determinant
      c <- -b12 / determinant
      correction <- a * V11 + 2 * c * V12 + d * V22
      M11 <- a^2 * V11 + 2 * a * c * V12 + c^2 * V22
      M12 <- a * c * V11 + (a * d + c^2) * V12 + c * d * V22
      M22 <- c^2 * V11 + 2 * c * d * V12 + d^2 * V22
      identity[ids] <- rowSums(w) - 2 * correction +
        Q11 * M11 + 2 * Q12 * M12 + Q22 * M22 +
        Lambda[1, 1] * M11 + 2 * Lambda[1, 2] * M12 + Lambda[2, 2] * M22
    }
    coefficient[ids] <- rowSums(w) - correction
  }
  list(coefficient = coefficient, identity = identity)
}

profile_pair_conditional_gaussian_components_v2 <- function(dat, beta, tau, h, lambda_gamma,
                                                            nuisance_tol = 1e-9) {
  if (!identical(as.numeric(beta), as.numeric(dat$beta0)) ||
      !identical(dat$population_data_definition, "conditional_gaussian_X_integrated")) {
    stop("Conditional Gaussian integration requires the generated structural coefficient.")
  }
  layout <- dat$population_layout
  Lambda <- make_lambda_matrix_v2(lambda_gamma, ncol(dat$Z))
  half <- profile_pair_population_nuisance_v2(dat$residual_base, dat$Z, layout, tau, h / 2,
                                              Lambda, nuisance_tol)
  full <- profile_pair_population_nuisance_v2(dat$residual_base, dat$Z, layout, tau, h,
                                              Lambda, nuisance_tol)
  ch <- profile_pair_conditional_hessian_coefficients_v2(half, dat$Z, layout, Lambda)
  cf <- profile_pair_conditional_hessian_coefficients_v2(full, dat$Z, layout, Lambda)
  inverse_m2 <- (1 / layout$sizes[layout$row_cluster])^2
  n <- layout$n
  Sigma <- dat$Sigma_x
  Hh <- mean(ch$coefficient) * Sigma
  Hf <- mean(cf$coefficient) * Sigma
  covariance_half <- sum(half$psi^2 * inverse_m2) / n * Sigma
  covariance_full <- sum(full$psi^2 * inverse_m2) / n * Sigma
  cross <- sum(half$psi * full$psi * inverse_m2) / n * Sigma
  covariance <- sum((2 * half$psi - full$psi)^2 * inverse_m2) / n * Sigma
  H <- 2 * Hh - Hf
  identity_H <- mean(2 * ch$identity - cf$identity) * Sigma
  covariance_identity <- max(abs(covariance -
    (4 * covariance_half + covariance_full - 2 * cross - 2 * t(cross))))
  zero <- stats::setNames(numeric(dat$p), names(dat$beta0))
  list(objective = 2 * half$objective - full$objective, score = zero,
       hessian = H, hessian_identity = identity_H,
       profile_identity_error = max(abs(H - identity_H)),
       max_nuisance_gradient = max(half$max_nuisance_gradient, full$max_nuisance_gradient),
       nuisance_converged = TRUE, n_clusters = n, h = h, bandwidths = c(h / 2, h),
       score_definition = "profile_bandwidth_pair_2_half_minus_full",
       population_engine = "conditional_gaussian_rows_rao_blackwell_v1",
       population_nuisance_backend = getOption("jpDMEQR.population_nuisance_backend", "vectorised_R"),
       population_moments = list(
         n_clusters = n, score_mean = zero, score_mean_half = zero, score_mean_full = zero,
         Sigma_population = covariance, Sigma_half = covariance_half, Sigma_full = covariance_full,
         cross_covariance_half_full = cross, covariance_identity_error = covariance_identity,
         H_population = H, H_half = Hh, H_full = Hf,
         conditional_integration = list(
           definition = "independent_gaussian_design_rows_integrated_given_nuisance_variables",
           hessian_coefficient = mean(2 * ch$coefficient - cf$coefficient),
           covariance_coefficient = sum((2 * half$psi - full$psi)^2 * inverse_m2) / n,
           score_mean_definition = "conditional_expectation_exactly_zero",
           rng_definition = "reference_design_normal_draws_consumed_in_original_order")))
}

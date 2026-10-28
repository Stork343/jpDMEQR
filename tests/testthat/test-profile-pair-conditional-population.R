conditional_root <- v2_require_root()
source(file.path(conditional_root, "scripts", "00_source_v2.R"), local = TRUE)
source_v2_module(conditional_root, "profile_v2", environment())
source_v2_module(conditional_root, "simulation_v2", environment())

testthat::test_that("conditional Gaussian population integration retains reference nuisance observations", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  for (design in c("ar1", "independent", "block", "sparse_precision", "factor", "dense_precision")) for (q in 1:2) {
    config <- list(
      experiment_id = paste0("CONDITIONAL-", design), n_clusters = 40L, p = 8L, s = 2L,
      tau = 0.25, q = q, m_values = 3:5, m_rule = "uniform", m_multiset = NULL,
      design_type = design, rho_x = 0.3, signal = 0.75, error_dist = "gaussian",
      heteroskedastic = FALSE, error_dependence = "independent", copula_rho = 0,
      random_effect_dist = "normal", sigma_b0 = 0.8, sigma_b1 = 0.5, rho_b = 0.4,
      x_b_corr = 0, informative_size = FALSE, nonlinear_re_strength = 0,
      fit_random_effects = if (q == 1L) "intercept" else "intercept+slope",
      lambda_gamma = if (q == 1L) matrix(0.8, 1L) else matrix(c(1, 0.25, 0.25, 0.8), 2L),
      h_multiplier = 1, target_mode = "structural")
    full <- profile_pair_population_data_v2(config, 40L, 20261203L)
    conditional <- profile_pair_conditional_gaussian_data_v2(config, 40L, 20261203L)
    testthat::expect_identical(conditional$epsilon, full$epsilon)
    testthat::expect_identical(conditional$time, full$time)
    testthat::expect_equal(unname(conditional$random_effects), unname(full$random_effects), tolerance = 0)
    testthat::expect_equal(unname(conditional$m_i), unname(full$m_i), tolerance = 0)
    testthat::expect_lt(max(abs(conditional$residual_base - (full$y - full$X %*% full$beta0))), 1e-12)
    reference <- profile_bandwidth_pair_v2(
      full$y, full$X, full$Z, full$cluster_id, full$beta0, config$tau, 0.2,
      config$lambda_gamma, nuisance_control = list(grad_tol = 1e-8), need_fold_sums = FALSE)
    observed <- profile_pair_conditional_gaussian_components_v2(
      conditional, conditional$beta0, config$tau, 0.2, config$lambda_gamma)
    ci <- cluster_index_v2(full$cluster_id)
    coefficient <- function(component) {
      mean(vapply(seq_len(ci$n), function(i) {
        Z <- full$Z[ci$rows[[i]], , drop = FALSE]
        w <- phi_smooth_v2(component$residuals[[i]], config$tau, component$h) / nrow(Z)
        W <- diag(w)
        L <- W - W %*% Z %*% solve(crossprod(Z, W %*% Z) + config$lambda_gamma) %*% t(Z) %*% W
        sum(diag(L))
      }, numeric(1)))
    }
    expected_h <- (2 * coefficient(reference$components_half_h) - coefficient(reference$components_h)) *
      conditional$Sigma_x
    score_coefficient <- mean(vapply(seq_len(ci$n), function(i) {
      pair <- 2 * psi_smooth_v2(reference$components_half_h$residuals[[i]], config$tau, 0.1) -
        psi_smooth_v2(reference$components_h$residuals[[i]], config$tau, 0.2)
      sum(pair^2) / length(pair)^2
    }, numeric(1)))
    testthat::expect_lt(max(abs(observed$hessian - expected_h)), 1e-6)
    testthat::expect_lt(max(abs(observed$population_moments$Sigma_population -
                                 score_coefficient * conditional$Sigma_x)), 1e-7)
    testthat::expect_lt(observed$profile_identity_error, 1e-10)
    testthat::expect_lt(observed$population_moments$covariance_identity_error, 1e-10)
    testthat::expect_true(all(observed$population_moments$score_mean == 0))
    incompatible <- config
    incompatible$x_b_corr <- 0.4
    testthat::expect_error(profile_pair_gaussian_row_conditions_v2(incompatible), "conditions")
  }
})

testthat::test_that("conditional structural directions retain full dimension and independent repeats", {
  config <- list(
    experiment_id = "CONDITIONAL-DIRECTION", n_clusters = 60L, p = 20L, s = 2L,
    tau = 0.5, q = 1L, m_values = 3:5, m_rule = "uniform", m_multiset = NULL,
    design_type = "factor", rho_x = 0.3, signal = 0.75, error_dist = "gaussian",
    heteroskedastic = FALSE, error_dependence = "independent", copula_rho = 0,
    random_effect_dist = "normal", sigma_b0 = 0.8, sigma_b1 = 0.5, rho_b = 0.4,
    x_b_corr = 0, informative_size = FALSE, nonlinear_re_strength = 0,
    fit_random_effects = "intercept", lambda_gamma = 1, h_multiplier = 1,
    target_mode = "structural")
  direction <- build_profile_pair_direction_asset_v2(
    config, n_population = 60L, repeats = 2L, development = TRUE,
    control = list(population_engine = "conditional_gaussian"))
  testthat::expect_true(validate_profile_pair_direction_asset_v2(direction, config, final = FALSE)$valid)
  testthat::expect_equal(dim(direction$H_population), c(20L, 20L))
  testthat::expect_equal(dim(direction$Sigma_population), c(20L, 20L))
  testthat::expect_equal(length(direction$directions), 20L)
  testthat::expect_true(all(direction$score_mean_population == 0))
  testthat::expect_identical(direction$repeat_assets[[1]]$population_data_definition,
                             "conditional_gaussian_X_integrated")
  testthat::expect_identical(direction$repeat_assets[[1]]$conditional_integration$score_mean_definition,
                             "conditional_expectation_exactly_zero")
  testthat::expect_true(all(is.finite(direction$H_mc_se)))
  testthat::expect_length(unique(direction$seeds), 2L)
})

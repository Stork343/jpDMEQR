streamed_root <- v2_require_root()
source(file.path(streamed_root, "scripts", "00_source_v2.R"), local = TRUE)
source_v2_module(streamed_root, "profile_v2", environment())
source_v2_module(streamed_root, "simulation_v2", environment())

testthat::test_that("streamed population components equal reference calculations on identical observations", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  for (q in 1:2) for (tau in c(0.1, 0.25, 0.5, 0.75)) {
    dat <- generate_profile_qr_data_v2(40L, 8L, 2L, tau, q = q,
                                       m_multiset = rep(c(2L, 3L, 5L, 4L), 10L),
                                       seed = 20261100L + q * 10L + round(tau * 100))
    beta <- dat$beta0 + seq(-0.03, 0.04, length.out = 8L)
    Lambda <- if (q == 1L) 0.8 else matrix(c(1, 0.25, 0.25, 0.8), 2L)
    reference <- profile_bandwidth_pair_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, 0.2, Lambda,
      nuisance_control = list(grad_tol = 1e-8), identity_check = TRUE,
      assembly = "sequential", need_fold_sums = FALSE
    )
    expected <- profile_pair_moments_v2(reference)
    for (batch in c(7L, 40L)) {
      observed <- profile_pair_population_streamed_v2(dat, beta, tau, 0.2, Lambda,
                                                       batch_clusters = batch)
      testthat::expect_equal(observed$objective, reference$objective, tolerance = 1e-9)
      testthat::expect_lt(max(abs(observed$score - reference$score)), 2e-8)
      testthat::expect_equal(observed$hessian, reference$hessian, tolerance = 2e-7)
      testthat::expect_lt(observed$profile_identity_error, 1e-10)
      testthat::expect_lt(observed$max_nuisance_gradient, 1e-8)
      for (field in c("Sigma_population", "Sigma_half", "Sigma_full", "cross_covariance_half_full")) {
        testthat::expect_equal(observed$population_moments[[field]], expected[[field]], tolerance = 2e-8)
      }
      score_only <- profile_pair_population_streamed_v2(
        dat, beta, tau, 0.2, Lambda, batch_clusters = batch,
        need_hessian = FALSE, need_moments = FALSE)
      testthat::expect_equal(score_only$score, observed$score, tolerance = 1e-12)
    }
  }
})

testthat::test_that("streamed roots and structural directions preserve the population definition", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  config <- list(
    experiment_id = "STREAMED-STRUCTURAL", n_clusters = 60L, p = 2L, s = 1L,
    tau = 0.5, q = 1L, m_values = 3:4, m_rule = "uniform", m_multiset = NULL,
    design_type = "ar1", rho_x = 0.3, signal = 0.75, error_dist = "gaussian",
    heteroskedastic = FALSE, error_dependence = "independent", copula_rho = 0,
    random_effect_dist = "normal", sigma_b0 = 0.8, sigma_b1 = 0.5, rho_b = 0.4,
    x_b_corr = 0, informative_size = FALSE, nonlinear_re_strength = 0,
    fit_random_effects = "intercept", lambda_gamma = 1, h_multiplier = 1,
    target_mode = "structural"
  )
  dat <- profile_pair_population_data_v2(config, 60L, 20261123L)
  reference <- solve_profile_pair_target_v2(dat, 0.5, 0.15, 1)
  streamed <- solve_profile_pair_target_v2(dat, 0.5, 0.15, 1,
                                           control = list(population_engine = "streamed", batch_clusters = 13L))
  testthat::expect_equal(streamed$beta, reference$beta, tolerance = 1e-6)
  testthat::expect_lt(streamed$kkt_residual, 1e-7)
  direction <- build_profile_pair_direction_asset_v2(
    config, n_population = 60L, repeats = 2L, development = TRUE,
    control = list(population_engine = "streamed", batch_clusters = 13L))
  testthat::expect_true(validate_profile_pair_direction_asset_v2(direction, config, final = FALSE)$valid)
  testthat::expect_identical(direction$target_asset_identity$reference_source, "structural_definition")
  testthat::expect_identical(direction$beta_target, dat$beta0)
  testthat::expect_identical(direction$acceptance$population_engine, "streamed")
  testthat::expect_lt(direction$covariance_identity_error, 1e-10)
})

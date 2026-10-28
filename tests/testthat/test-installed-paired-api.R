testthat::test_that("installed paired APIs retain the cluster score correction and sandwich", {
  withr::local_envvar(c(
    JPDMEQR_CLUSTER_SOLVER = "compiled", JPDMEQR_ASSEMBLY = "batch",
    JPDMEQR_MATRIX_BACKEND = "Eigen_single_thread", JPDMEQR_IDENTITY_DETAIL = "aggregate",
    JPDMEQR_UNPENALISED_OPTIMIZER = "newton_armijo", JPDMEQR_PROXIMAL_SCALING = "cluster_rms",
    JPDMEQR_DANTZIG_SOLVER = "dual_warm", JPDMEQR_HESSIAN_REFRESH = "20",
    JPDMEQR_SHARED_FIRST_STAGE = "TRUE", JPDMEQR_CV_CORES = "1"
  ))
  generate <- getFromNamespace("generate_profile_qr_data_v2", "jpDMEQR")
  dat <- generate(n_clusters = 40L, p = 6L, s = 2L, tau = 0.25, q = 2L,
                  m_values = 3:5, seed = 20261121L)
  Lambda <- matrix(c(1, 0.2, 0.2, 0.8), 2L)
  targets <- colnames(dat$X)[c(1L, 4L)]
  fit <- jpDMEQR::fit_profile_pair_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, 0.25,
    h_inf = 40^(-3 / 10), lambda_gamma = Lambda, target_coordinates = targets
  )
  inf <- jpDMEQR::infer_profile_pair_v2(fit, targets, c(0.1, 0.3, 0.6))
  direct <- jpDMEQR::profile_bandwidth_pair_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, fit$beta_refit, 0.25,
    h = fit$h_inf, lambda_gamma = Lambda
  )
  testthat::expect_true(fit$converged)
  testthat::expect_identical(fit$score_definition, "profile_bandwidth_pair_2_half_minus_full")
  testthat::expect_equal(fit$Lambda, Lambda)
  testthat::expect_true(all(inf$table$feasible))
  testthat::expect_true(all(is.finite(inf$table$se) & inf$table$se > 0))
  testthat::expect_equal(fit$components$score, direct$score, tolerance = 1e-10)
  testthat::expect_equal(fit$components$hessian, direct$hessian, tolerance = 1e-10)
  for (i in seq_along(targets)) {
    omega <- inf$directions[[i]]$omega
    projected <- as.numeric(direct$cluster_scores %*% omega)
    testthat::expect_equal(inf$table$beta_tilde[[i]],
                           unname(fit$beta_refit[targets[[i]]] - sum(omega * direct$score)),
                           tolerance = 1e-12)
    testthat::expect_equal(inf$table$se[[i]]^2 * dat$n_clusters,
                           mean((projected - mean(projected))^2), tolerance = 1e-12)
  }
})

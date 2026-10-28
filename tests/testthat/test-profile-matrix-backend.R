root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", .GlobalEnv)
source_v2_module(root, "simulation_v2", .GlobalEnv)

testthat::test_that("aggregate Eigen assembly preserves actual paired derivatives", {
  for (q in c(1L, 2L)) for (tau in c(0.25, 0.5, 0.75)) {
    dat <- generate_profile_qr_data_v2(
      n_clusters = 30L, p = 40L, s = 3L, q = q, tau = tau,
      seed = 20261300L + 10L * q + round(100 * tau)
    )
    beta <- dat$beta0 * 0.8
    withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "compiled",
                          JPDMEQR_ASSEMBLY = "batch",
                          JPDMEQR_MATRIX_BACKEND = "R_BLAS",
                          JPDMEQR_IDENTITY_DETAIL = "cluster"))
    reference <- profile_bandwidth_pair_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, 0.3, identity_check = TRUE
    )
    withr::local_envvar(c(JPDMEQR_MATRIX_BACKEND = "Eigen_single_thread",
                          JPDMEQR_IDENTITY_DETAIL = "aggregate"))
    optimized <- profile_bandwidth_pair_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, 0.3, identity_check = TRUE
    )
    testthat::expect_equal(optimized$score, reference$score, tolerance = 1e-12)
    testthat::expect_equal(optimized$hessian, reference$hessian, tolerance = 1e-11)
    testthat::expect_equal(optimized$cluster_scores, reference$cluster_scores, tolerance = 1e-12)
    testthat::expect_equal(optimized$hessian_fold_sums, reference$hessian_fold_sums, tolerance = 1e-10)
    testthat::expect_lte(optimized$profile_identity_error, 1e-8)
    testthat::expect_false(optimized$components_h$cluster_identity_checked)
    testthat::expect_true(all(is.na(optimized$components_h$cluster_identity_errors)))
  }
})

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", .GlobalEnv)
source_v2_module(root, "simulation_v2", .GlobalEnv)

testthat::test_that("compiled cluster profiles preserve the complete criterion", {
  testthat::expect_true(is.loaded("jpdmeqr_cluster_profiles", PACKAGE = "jpDMEQR"))
  for (q in c(1L, 2L)) for (tau in c(0.1, 0.25, 0.5, 0.75)) {
    dat <- generate_profile_qr_data_v2(
      n_clusters = 18L, p = 7L, s = 2L, tau = tau, q = q,
      m_multiset = rep(c(2L, 3L, 5L), 6L), seed = 20261100L + 10L * q + round(100 * tau)
    )
    Lambda <- if (q == 1L) matrix(0.8, 1L) else matrix(c(1, 0.2, 0.2, 0.7), 2L)
    beta <- dat$beta0 * 0.9
    withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
    reference <- profile_components_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, h = 0.3,
      lambda_gamma = Lambda, identity_check = TRUE, assembly = "batch"
    )
    withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "compiled"))
    compiled <- profile_components_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, h = 0.3,
      lambda_gamma = Lambda, identity_check = TRUE, assembly = "batch"
    )
    testthat::expect_equal(compiled$objective, reference$objective, tolerance = 1e-10)
    testthat::expect_equal(compiled$score, reference$score, tolerance = 1e-8)
    testthat::expect_equal(compiled$hessian, reference$hessian, tolerance = 1e-7)
    testthat::expect_equal(compiled$cluster_scores, reference$cluster_scores, tolerance = 1e-8)
    testthat::expect_equal(compiled$gamma, reference$gamma, tolerance = 1e-7)
    testthat::expect_true(compiled$nuisance_converged)
    testthat::expect_lte(compiled$max_nuisance_gradient, 1e-8)
    testthat::expect_lte(compiled$profile_identity_error, 1e-8)
    for (fold in seq_len(4L)) {
      testthat::expect_equal(compiled$hessian_fold_sums[[fold]],
                             reference$hessian_fold_sums[[fold]], tolerance = 1e-7)
    }
    warm <- profile_components_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, h = 0.3,
      lambda_gamma = Lambda, gamma_init = compiled$gamma, need_hessian = FALSE
    )
    testthat::expect_equal(warm$score, compiled$score, tolerance = 1e-12)
  }
})

testthat::test_that("compiled paired fits preserve inference on actual clustered data", {
  withr::local_envvar(c(JPDMEQR_UNPENALISED_OPTIMIZER = "proximal_gradient",
                        JPDMEQR_PROXIMAL_SCALING = "none",
                        JPDMEQR_HESSIAN_REFRESH = "1", JPDMEQR_IDENTITY_DETAIL = "cluster",
                        JPDMEQR_DANTZIG_SOLVER = "clarabel_native"))
  dat <- generate_profile_qr_data_v2(
    n_clusters = 40L, p = 10L, s = 2L, tau = 0.5, q = 2L, seed = 20261203L
  )
  targets <- colnames(dat$X)[c(1L, 4L)]
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton", JPDMEQR_ASSEMBLY = "batch"))
  reference <- fit_profile_pair_v2(dat$y, dat$X, dat$Z, dat$cluster_id,
                                  0.5, 40^(-3 / 10), target_coordinates = targets)
  reference_inference <- infer_profile_pair_v2(reference, targets, c(0.05, 0.15, 0.4))
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "compiled"))
  compiled <- fit_profile_pair_v2(dat$y, dat$X, dat$Z, dat$cluster_id,
                                 0.5, 40^(-3 / 10), target_coordinates = targets)
  compiled_inference <- infer_profile_pair_v2(compiled, targets, c(0.05, 0.15, 0.4))
  testthat::expect_equal(compiled$beta_l1, reference$beta_l1, tolerance = 1e-7)
  testthat::expect_equal(compiled$beta_refit, reference$beta_refit, tolerance = 1e-7)
  testthat::expect_equal(compiled$components$score, reference$components$score, tolerance = 1e-7)
  testthat::expect_equal(compiled$components$hessian, reference$components$hessian, tolerance = 1e-6)
  testthat::expect_equal(compiled_inference$table$mu, reference_inference$table$mu)
  testthat::expect_equal(compiled_inference$table$beta_tilde,
                         reference_inference$table$beta_tilde, tolerance = 1e-6)
  testthat::expect_equal(compiled_inference$table$se, reference_inference$table$se, tolerance = 1e-6)
})

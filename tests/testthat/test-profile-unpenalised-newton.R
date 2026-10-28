newton_root <- v2_require_root()
source(file.path(newton_root, "scripts", "00_source_v2.R"), local = TRUE)
source_v2_module(newton_root, "profile_v2", environment())
source_v2_module(newton_root, "simulation_v2", environment())

testthat::test_that("Newton 与 proximal 在完整凸 profile 上取得相同解", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "compiled"))
  testthat::expect_true(is.loaded("jpdmeqr_cluster_profiles", PACKAGE = "jpDMEQR"))
  for (q in 1:2) {
    dat <- generate_profile_qr_data_v2(
      n_clusters = 40L, p = 4L, s = 2L, tau = 0.25, q = q,
      m_values = 3:5, seed = 20261022L + q
    )
    Lambda <- if (q == 1L) 1 else matrix(c(1, .2, .2, .8), 2L)
    control <- list(max_iter = 5000L, beta_tol = 1e-9,
                    kkt_normalized_tol = 1e-8,
                    nuisance_control = list(grad_tol = 1e-9))
    for (tau in c(.25, .5, .75)) {
      reference <- fit_profile_lasso_v2(
        dat$y, dat$X, dat$Z, dat$cluster_id, tau, .4, 0,
        Lambda, rep(0, 4), control = modifyList(control, list(
          unpenalised_optimizer = "proximal_gradient"
        ))
      )
      result <- fit_profile_lasso_v2(
        dat$y, dat$X, dat$Z, dat$cluster_id, tau, .4, 0,
        Lambda, rep(0, 4), control = modifyList(control, list(
          unpenalised_optimizer = "newton_armijo"
        ))
      )
      testthat::expect_true(reference$converged)
      testthat::expect_true(result$converged)
      testthat::expect_lte(result$kkt_normalized, 1e-8)
      testthat::expect_lte(result$components$max_nuisance_gradient, 1e-9)
      testthat::expect_equal(result$beta, reference$beta, tolerance = 2e-7)
      testthat::expect_equal(result$profile_objective, reference$profile_objective,
                             tolerance = 1e-12)
      testthat::expect_equal(result$components$score, reference$components$score,
                             tolerance = 2e-8)
      testthat::expect_true(all(diff(result$history$profile_objective) <= 1e-12))
      testthat::expect_true(isTRUE(result$components$identity_checked))
      testthat::expect_lt(result$components$profile_identity_error, 1e-10)
    }
  }
})

testthat::test_that("显式 Newton 配置保持 L1 路径和无效输入检验", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "compiled"))
  dat <- generate_profile_qr_data_v2(
    n_clusters = 30L, p = 4L, s = 2L, tau = .5, q = 1L,
    m_values = 3:5, seed = 20261029L
  )
  a <- fit_profile_lasso_v2(dat$y, dat$X, dat$Z, dat$cluster_id, .5, .4, .1,
                            control = list(unpenalised_optimizer = "proximal_gradient"))
  b <- fit_profile_lasso_v2(dat$y, dat$X, dat$Z, dat$cluster_id, .5, .4, .1,
                            control = list(unpenalised_optimizer = "newton_armijo"))
  testthat::expect_identical(a$beta, b$beta)
  testthat::expect_identical(a$history, b$history)
  testthat::expect_error(fit_profile_unpenalised_newton_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, .5, .4,
    lambda_beta = .1, penalty_factor = rep(1, 4)), "unpenalised")
  testthat::expect_error(fit_profile_lasso_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, .5, .4, 0,
    control = list(unpenalised_optimizer = "unknown")), "Unknown")
})

scaling_root <- v2_require_root()
source(file.path(scaling_root, "scripts", "00_source_v2.R"), local = TRUE)
source_v2_module(scaling_root, "profile_v2", environment())
source_v2_module(scaling_root, "simulation_v2", environment())

testthat::test_that("cluster RMS proximal 更新保持同一 weighted L1 criterion", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton",
                        JPDMEQR_UNPENALISED_OPTIMIZER = "proximal_gradient"))
  for (q in 1:2) {
    dat <- generate_profile_qr_data_v2(
      n_clusters = 40L, p = 4L, s = 2L, tau = .5, q = q,
      m_values = 3:5, seed = 20261102L + q
    )
    X <- sweep(dat$X, 2, c(1, 2, .5, 1.5), `*`)
    pf <- c(0, 1, 1.5, .5)
    Lambda <- if (q == 1L) 1 else matrix(c(1, .2, .2, .8), 2L)
    control <- list(max_iter = 5000L, beta_tol = 1e-9,
                    kkt_normalized_tol = 1e-7,
                    nuisance_control = list(grad_tol = 1e-9))
    for (tau in c(.25, .5, .75)) {
      reference <- fit_profile_lasso_v2(
        dat$y, X, dat$Z, dat$cluster_id, tau, .4, .1, Lambda, pf,
        control = modifyList(control, list(proximal_scaling = "none"))
      )
      result <- fit_profile_lasso_v2(
        dat$y, X, dat$Z, dat$cluster_id, tau, .4, .1, Lambda, pf,
        control = modifyList(control, list(proximal_scaling = "cluster_rms"))
      )
      testthat::expect_true(reference$converged)
      testthat::expect_true(result$converged)
      testthat::expect_equal(result$beta, reference$beta, tolerance = 2e-6)
      testthat::expect_equal(result$composite_objective, reference$composite_objective,
                             tolerance = 1e-11)
      testthat::expect_equal(result$composite_objective,
                             result$components$objective + .1 * sum(pf * abs(result$beta)),
                             tolerance = 1e-12)
      testthat::expect_lte(result$kkt_normalized, 1e-7)
      testthat::expect_equal(result$kkt_normalized,
                             kkt_residual_normalized_v2(result$beta, result$components$score, .1 * pf))
      cluster <- cluster_index_v2(dat$cluster_id)
      rms <- sqrt(Reduce(`+`, lapply(cluster$rows, function(index) {
        colMeans(X[index, , drop = FALSE]^2)
      })) / cluster$n)
      testthat::expect_equal(unname(result$coordinate_scale), unname(rms), tolerance = 1e-14)
      testthat::expect_true(all(diff(result$history$composite_objective) <= 1e-12))
    }
  }
})

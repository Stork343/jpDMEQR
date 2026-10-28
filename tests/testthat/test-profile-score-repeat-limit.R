root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)

testthat::test_that("repeated observations satisfy the exact profile geometry", {
  for (tau in c(0.25, 0.5, 0.75)) {
    for (lambda in c(0.5, 1, 2)) {
      for (h in c(0.2, 0.05)) {
        for (t in c(-0.7, -0.2, 0.4, 0.8)) {
          s <- tau - 1 + kernel_epan_cdf_v2(t)
          w <- s / lambda + h * t
          comp <- profile_components_v2(
            y = rep(w, 2), X = matrix(1, 2, 1), Z = matrix(1, 2, 1),
            cluster_id = c(1, 1), beta = 0, tau = tau, h = h,
            lambda_gamma = lambda, identity_check = TRUE,
            nuisance_control = list(grad_tol = 1e-11, reltol = 1e-13)
          )
          density <- kernel_epan_pdf_v2(t) / h
          H <- lambda * density / (lambda + density)
          testthat::expect_lt(abs(comp$score + s), 2e-8)
          testthat::expect_lt(abs(comp$hessian[1, 1] - H), 2e-7)
          testthat::expect_lt(comp$profile_identity_error, 1e-12)
          testthat::expect_lt(abs(as.numeric(comp$gamma) - s / lambda), 2e-8)
        }
      }
    }
  }
})

testthat::test_that("the profile limit retains the nuisance subgradient", {
  tau <- 0.25
  lambda <- 1
  w <- -0.2
  h <- 0.001
  comp <- profile_components_v2(
    y = rep(w, 2), X = matrix(1, 2, 1), Z = matrix(1, 2, 1),
    cluster_id = c(1, 1), beta = 0, tau = tau, h = h,
    lambda_gamma = lambda,
    nuisance_control = list(grad_tol = 1e-10, reltol = 1e-13)
  )
  profile_limit <- -lambda * w
  ordinary_limit <- -tau
  testthat::expect_equal(comp$score, profile_limit, tolerance = 0.002)
  testthat::expect_equal(as.numeric(comp$unsmoothed_cluster_scores[1, 1]), ordinary_limit,
                        tolerance = 1e-12)
  testthat::expect_gt(abs(comp$score - comp$unsmoothed_cluster_scores[1, 1]), 0.44)
})

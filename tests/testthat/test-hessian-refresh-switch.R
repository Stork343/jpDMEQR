root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)

refresh_fixture <- function(seed = 20260930L) {
  generate_profile_qr_data_v2(24L, 8L, 3L, 0.5, q = 2L, m_values = 3:5,
                              seed = seed)
}

refresh_fit <- function(d, control = list()) {
  fit_profile_lasso_v2(d$y, d$X, d$Z, d$cluster_id, 0.5, 0.4,
                       lambda_beta = 0.05, control = control)
}

no_refresh_env <- function(code) {
  withr::with_envvar(c(JPDMEQR_HESSIAN_REFRESH = NA), code)
}

testthat::test_that("the frozen default is five and the environment selects another value", {
  testthat::expect_identical(no_refresh_env(hessian_refresh_v2()), 5L)
  testthat::expect_identical(
    withr::with_envvar(c(JPDMEQR_HESSIAN_REFRESH = "20"), hessian_refresh_v2()), 20L
  )
  testthat::expect_identical(
    no_refresh_env(withr::with_options(list(jpDMEQR.hessian_refresh = 9L),
                                       hessian_refresh_v2())),
    9L
  )
  testthat::expect_identical(
    no_refresh_env(withr::with_options(
      list(jpDMEQR.hessian_refresh = 9L),
      withr::with_envvar(c(JPDMEQR_HESSIAN_REFRESH = "7"), hessian_refresh_v2())
    )),
    7L
  )
})

testthat::test_that("the plain fit reads the control field and keeps the frozen default", {
  d <- refresh_fixture()
  base <- no_refresh_env(refresh_fit(d))
  same <- no_refresh_env(refresh_fit(d, control = list(hessian_refresh = 5L)))
  testthat::expect_identical(base$beta, same$beta)

  every <- refresh_fit(d, control = list(hessian_refresh = 1L))
  testthat::expect_true(every$converged)
  testthat::expect_lte(every$kkt_normalized, every$kkt_normalized_tol)

  wide <- refresh_fit(d, control = list(hessian_refresh = 40L))
  testthat::expect_true(wide$converged)
  testthat::expect_lte(wide$kkt_normalized, wide$kkt_normalized_tol)
  # 两次拟合停在同一个接受水平内的不同点，系数差别以该接受水平为界。
  testthat::expect_lte(max(abs(wide$beta - base$beta)),
                       base$kkt_normalized_tol * max(1, max(abs(base$beta))))
  testthat::expect_gt(wide$iterations, base$iterations)
})

testthat::test_that("the plain fit ignores the environment", {
  d <- refresh_fixture()
  without <- no_refresh_env(refresh_fit(d))
  with_env <- withr::with_envvar(c(JPDMEQR_HESSIAN_REFRESH = "40"), refresh_fit(d))
  testthat::expect_identical(with_env$beta, without$beta)
  testthat::expect_identical(with_env$iterations, without$iterations)
})

testthat::test_that("the calibrated first stage takes the environment value", {
  set.seed(20261001L)
  n_cluster <- 12L
  cluster_id <- rep(seq_len(n_cluster), each = 3L)
  N <- length(cluster_id)
  X <- matrix(stats::rnorm(N * 6L), N, 6L, dimnames = list(NULL, paste0("x", 1:6)))
  Z <- matrix(1, N, 1L, dimnames = list(NULL, "intercept"))
  beta0 <- c(0.6, -0.6, rep(0, 4))
  b0 <- stats::rnorm(n_cluster)
  y <- as.numeric(X %*% beta0 + b0[match(cluster_id, unique(cluster_id))] +
                    stats::rnorm(N))
  lam0 <- lambda_0_n_v2(n_cluster, ncol(X))$lambda_0
  calibrate <- function(control) {
    fit_profile_lasso_calibrated_v2(
      y, X, Z, cluster_id, tau = 0.5, h_est = 0.45, h_inf = 0.35,
      lambda_gamma = 1, lambda_0_n = lam0, control = control
    )
  }

  from_env <- withr::with_envvar(c(JPDMEQR_HESSIAN_REFRESH = "3"),
                                 calibrate(list()))
  from_control <- no_refresh_env(
    calibrate(list(fit_control = list(hessian_refresh = 3L)))
  )
  testthat::expect_identical(from_env$iterations, from_control$iterations)
  testthat::expect_identical(from_env$beta, from_control$beta)
})

testthat::test_that("a refresh interval that is not a positive integer stops the fit", {
  d <- refresh_fixture()
  testthat::expect_error(refresh_fit(d, control = list(hessian_refresh = 0L)),
                         "hessian_refresh")
  testthat::expect_error(refresh_fit(d, control = list(hessian_refresh = 2.5)),
                         "hessian_refresh")
  testthat::expect_error(
    withr::with_envvar(c(JPDMEQR_HESSIAN_REFRESH = "0"), hessian_refresh_v2()),
    "JPDMEQR_HESSIAN_REFRESH"
  )
  testthat::expect_error(
    withr::with_envvar(c(JPDMEQR_HESSIAN_REFRESH = "many"), hessian_refresh_v2()),
    "JPDMEQR_HESSIAN_REFRESH"
  )
})
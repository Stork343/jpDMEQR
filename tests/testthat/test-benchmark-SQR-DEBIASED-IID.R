# Deterministic acceptance tests for SQR-DEBIASED-IID (JMLR 22-1217 fidelity).

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
source_v2_module(root, "benchmark_adapters_v2", envir = .GlobalEnv)

make_iid_train <- function(n = 90, p = 25, s = 3, tau = 0.5, seed = 11L) {
  dat <- generate_profile_qr_data_v2(
    n_clusters = n, p = p, s = s, tau = tau, q = 1L,
    m_values = 1L, m_rule = "uniform",
    error_dist = "gaussian", seed = seed
  )
  list(
    y = dat$y,
    X = dat$X,
    Z = matrix(0, nrow = length(dat$y), ncol = 1L),
    cluster_id = dat$cluster_id,
    beta0 = dat$beta0,
    active = dat$active
  )
}

testthat::test_that("SQR-DEBIASED-IID adapter is deterministic and schema-conformant", {
  train <- make_iid_train(seed = 202)
  tuning <- list(h = 0.5, lambda_beta = 0.15, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2), lambda_beta_multiplier = 1)
  coords <- colnames(train$X)[1:5]
  a1 <- fit_benchmark_sqr_debiased_iid_v2(train, 0.5, coords, tuning, seed = 77L)
  a2 <- fit_benchmark_sqr_debiased_iid_v2(train, 0.5, coords, tuning, seed = 77L)

  testthat::expect_equal(a1$status, "ok")
  testthat::expect_identical(a1$beta_hat, a2$beta_hat)
  testthat::expect_identical(a1$beta_tilde, a2$beta_tilde)
  testthat::expect_true(all(c("method_id", "status", "beta_hat", "beta_tilde",
                              "se", "ci_lower", "ci_upper", "selected",
                              "runtime_sec", "converged", "kkt_residual",
                              "warning_messages", "implementation_version",
                              "reference_identifier", "adapter_fidelity_status",
                              "inference_object") %in% names(a1)))
  testthat::expect_equal(a1$method_id, "SQR-DEBIASED-IID")
  testthat::expect_equal(a1$adapter_fidelity_status, "formula_level_implemented")
  testthat::expect_true(all(names(a1$beta_hat) == colnames(train$X)))
  inf <- a1$inference_object
  testthat::expect_true(isTRUE(inf$iid))
  testthat::expect_true(all(c("coordinate", "beta_tilde", "se", "ci_lower",
                              "ci_upper", "feasible", "mu", "dantzig_residual",
                              "omega_l1", "omega_l2") %in% names(inf$table)))
  testthat::expect_true(all(inf$table$feasible))
  testthat::expect_true(all(is.finite(inf$table$se)))
  testthat::expect_true(all(is.finite(a1$beta_tilde)))
})

testthat::test_that("SQR-DEBIASED-IID uses observation-level sample size and iid variance", {
  train <- make_iid_train(seed = 303)
  tuning <- list(h = 0.5, lambda_beta = 0.15, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2), lambda_beta_multiplier = 1)
  coords <- colnames(train$X)[1:2]
  ans <- fit_benchmark_sqr_debiased_iid_v2(train, 0.5, coords, tuning, seed = 5L)
  inf <- ans$inference_object
  # variance = tau(1-tau) * omega' Sigma_hat omega / N
  N <- length(train$y)
  Sigma_hat <- crossprod(train$X) / N
  for (ii in seq_len(nrow(inf$table))) {
    nm <- inf$table$coordinate[ii]
    k <- match(nm, colnames(train$X))
    omega <- inf$directions[[ii]]$omega
    expected_var <- 0.5 * 0.5 * drop(crossprod(omega, Sigma_hat %*% omega)) / N
    testthat::expect_equal(inf$table$se[ii]^2, expected_var, tolerance = 1e-12)
  }
  testthat::expect_equal(inf$bandwidth,
                         max(sqrt(0.25) * sqrt(log(ncol(train$X))) * N^(-3/10), 0.05))
})

testthat::test_that("SQR-DEBIASED-IID limiting case: no penalty approaches smoothed QR", {
  train <- make_iid_train(seed = 404)
  tuning <- list(h = 0.5, lambda_beta = 1e-10, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2), lambda_beta_multiplier = 1e-10)
  coords <- colnames(train$X)[1:3]
  ans <- fit_benchmark_sqr_debiased_iid_v2(train, 0.5, coords, tuning, seed = 9L)
  # With lambda ~ 0, the one-step correction vanishes (KKT score ~ 0 at the
  # unpenalised optimum), so beta_tilde should be close to beta_hat.
  diff <- ans$beta_tilde - ans$beta_hat[names(ans$beta_tilde)]
  testthat::expect_true(all(is.finite(diff)))
  testthat::expect_true(max(abs(diff)) < 0.05)
})

testthat::test_that("SQR-DEBIASED-IID smooth loss matches analytic derivatives", {
  u <- seq(-2, 2, length.out = 11)
  h <- 0.4; tau <- 0.3
  eps <- 1e-6
  for (uu in u) {
    fd <- (sqr_iid_smooth_loss_v2(uu + eps, tau, h) -
             sqr_iid_smooth_loss_v2(uu - eps, tau, h)) / (2 * eps)
    analytic <- tau - stats::pnorm(-uu / h)
    testthat::expect_equal(fd, analytic, tolerance = 1e-4)
  }
  # Gaussian-kernel identity: l(u) = u*(tau - Phi(-u/h)) + h*phi(u/h)
  testthat::expect_equal(
    sqr_iid_smooth_loss_v2(1.3, tau, h),
    1.3 * (tau - stats::pnorm(-1.3 / h)) + h * stats::dnorm(1.3 / h),
    tolerance = 1e-12
  )
})

testthat::test_that("SQR-DEBIASED-IID score and Hessian match finite differences", {
  train <- make_iid_train(n = 60, p = 8, s = 2, seed = 505)
  beta <- rnorm(8) * 0.2
  h <- 0.6; tau <- 0.5
  score <- sqr_iid_score_fn_v2(train$y, train$X, beta, tau, h)
  obj <- function(b) mean(sqr_iid_smooth_loss_v2(train$y - as.numeric(train$X %*% b), tau, h))
  fd_grad <- finite_difference_gradient_v2(obj, beta, eps = 1e-5)
  testthat::expect_equal(as.numeric(score), as.numeric(fd_grad), tolerance = 2e-4)
  H <- sqr_iid_hessian_fn_v2(train$y, train$X, beta, tau, h)
  score_fn <- function(b) sqr_iid_score_fn_v2(train$y, train$X, b, tau, h)
  fd_hess <- finite_difference_jacobian_v2(score_fn, beta, eps = 1e-5)
  fd_hess <- (fd_hess + t(fd_hess)) / 2
  testthat::expect_equal(as.numeric(H), as.numeric(fd_hess), tolerance = 3e-3)
})

testthat::test_that("SQR-DEBIASED-IID ignores clustering and is labelled naive_iid", {
  # Two clusters with strongly different means: the iid adapter must not use
  # cluster structure (fixed-effect design only), matching a plain iid fit.
  dat <- generate_profile_qr_data_v2(
    n_clusters = 40, p = 15, s = 2, tau = 0.5, q = 1L,
    m_values = 2:4, m_rule = "uniform", sigma_b0 = 3,
    error_dist = "gaussian", seed = 606
  )
  train <- list(y = dat$y, X = dat$X,
                Z = matrix(0, nrow = length(dat$y), ncol = 1L),
                cluster_id = dat$cluster_id, beta0 = dat$beta0, active = dat$active)
  tuning <- list(h = 0.5, lambda_beta = 0.2, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2), lambda_beta_multiplier = 1)
  ans <- fit_benchmark_sqr_debiased_iid_v2(train, 0.5, colnames(train$X)[1:2],
                                           tuning, seed = 1L)
  testthat::expect_equal(ans$variance_scope, "naive_iid")
  testthat::expect_equal(ans$target_scope, "iid_smoothed_quantile")
})

testthat::test_that("SQR-DEBIASED-IID complete CLIME matrix and symmetrisation match independent solvers", {
  train <- make_iid_train(n = 80L, p = 6L, s = 2L, seed = 722L)
  beta <- fit_sqr_iid_lasso_v2(train$y, train$X, .35, .45, .05)
  grid <- c(.05, .15, .4, 1.1)
  dual <- solve_clime_matrix_iid_v2(beta$hessian, grid, solver_opts = list(backend = "dual_warm"))
  reference <- solve_clime_matrix_iid_v2(beta$hessian, grid, solver_opts = list(backend = "cvxr"))
  testthat::expect_true(dual$feasible && reference$feasible)
  testthat::expect_equal(dual$gamma, reference$gamma, tolerance = 1e-12)
  testthat::expect_equal(dual$omega_raw, reference$omega_raw, tolerance = 1e-5)
  testthat::expect_equal(dual$omega_symmetric, reference$omega_symmetric, tolerance = 1e-5)
  testthat::expect_equal(dual$omega_symmetric, t(dual$omega_symmetric), tolerance = 0)
  testthat::expect_true(all(abs(dual$omega_symmetric) <= abs(dual$omega_raw) + 1e-14))
  testthat::expect_lte(max(abs(beta$hessian %*% t(dual$omega_raw) - diag(6L))), dual$gamma + 1e-6)
  testthat::expect_identical(dual$excluded_zero_direction$status, "excluded_zero_direction")
  testthat::expect_length(dual$raw_grid_records, 6L)
  testthat::expect_true(any(dual$omega_raw == 0))
  testthat::expect_equal(dual$matrix_grid_record$model_constructions, 1L)
  certificate <- dual$matrix_grid_record$feasibility_certificate
  testthat::expect_true(certificate$verified)
  testthat::expect_lte(max(abs(beta$hessian %*% certificate$direction_matrix - diag(6L))), min(grid))
  testthat::expect_equal(nrow(dual$matrix_grid_record$attempts), 6L)
  testthat::expect_true(all(dual$raw_grid_records[[1L]]$status[-1L] ==
                            "feasible_by_same_full_matrix_certificate"))
})

testthat::test_that("SQR-DEBIASED-IID common grid preserves genuine rank deficient matrix failures", {
  train <- make_iid_train(n = 8L, p = 18L, s = 2L, seed = 723L)
  H <- sqr_iid_hessian_fn_v2(train$y, train$X, rep(0, 18L), .5, .4)
  grid <- c(.05, .2, .6, .9)
  dual <- solve_clime_matrix_iid_v2(H, grid, solver_opts = list(backend = "dual_warm"))
  reference <- solve_clime_matrix_iid_v2(H, grid, solver_opts = list(backend = "cvxr"))
  testthat::expect_identical(dual$feasible, reference$feasible)
  if (dual$feasible) testthat::expect_equal(dual$gamma, reference$gamma, tolerance = 1e-12)
  testthat::expect_equal(dual$matrix_grid_record$model_constructions, 1L)
  testthat::expect_true(any(!dual$matrix_grid_record$attempts$feasible))
  cutoff <- dual$matrix_grid_record$global_infeasible_at
  testthat::expect_true(is.finite(cutoff))
  testthat::expect_true(all(dual$matrix_grid_record$attempts$gamma >= cutoff))
})

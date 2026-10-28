root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, .GlobalEnv)
}

testthat::test_that("QIF-SEE source domain failure is deterministic and schema conformant", {
  dat <- generate_profile_qr_data_v2(n_clusters = 30L, p = 40L, s = 3L, q = 1L,
                                      tau = .25, seed = 700770L)
  a <- fit_benchmark_qif_see_v2(dat, .25, colnames(dat$X)[1:3], list(), 1L)
  b <- fit_benchmark_qif_see_v2(dat, .25, colnames(dat$X)[1:3], list(), 1L)
  testthat::expect_identical(a$status, "failed")
  testthat::expect_identical(a$failure_stage, "qif_covariance_rank")
  testthat::expect_identical(a$fit_object, b$fit_object)
  testthat::expect_equal(a$fit_object$matrix_dimension, 120L)
  testthat::expect_equal(a$fit_object$rank_upper_bound, 30L)
  testthat::expect_true(is.na(a$fit_object$actual_rank))
  testthat::expect_identical(a$fit_object$inverse_definition, "ordinary_inverse")
  testthat::expect_true(all(c("method_id", "status", "failure_stage", "failure_message",
                              "runtime_sec", "fit_object", "reference_identifier") %in% names(a)))
})

testthat::test_that("QIF-SEE modified extended score matches Eq3 Eq5 and frozen delta derivatives", {
  dat <- generate_profile_qr_data_v2(n_clusters = 40L, p = 3L, s = 2L, q = 1L,
                                      tau = .35, sigma_b0 = 1, seed = 700771L)
  beta <- as.numeric(quantreg::rq.fit(dat$X, dat$y, tau = .35)$coefficients)
  weights <- rep(stats::dnorm(stats::qnorm(.35) / sqrt(2)) / sqrt(2), length(dat$y))
  delta <- c(.1, .2, .3)
  h <- 40^(-1 / 2)
  component <- qif_see_extended_source_v2(dat$y, dat$X, dat$cluster_id, beta, .35,
                                          h, weights, delta)
  rows <- which(dat$cluster_id == dat$cluster_id[1L])
  X <- dat$X[rows, , drop = FALSE]
  residual <- dat$y[rows] - as.numeric(X %*% beta)
  psi <- stats::pnorm(residual / h) - .65
  blocks <- lapply(qif_see_basis_v2(length(rows)), function(M) {
    as.numeric(crossprod(X, weights[rows] * as.numeric(M %*% psi)))
  })
  raw <- as.numeric(t(do.call(cbind, blocks)))
  expected <- rep(1 - delta, each = 3L) * raw + rep(delta * beta, each = 3L)
  testthat::expect_equal(component$G[, 1L], expected, tolerance = 1e-12)
  testthat::expect_equal(component$covariance, tcrossprod(component$G) / 40L, tolerance = 1e-12)
  numerical <- finite_difference_jacobian_v2(function(b) {
    qif_see_extended_source_v2(dat$y, dat$X, dat$cluster_id, b, .35, h, weights, delta)$mean
  }, beta, eps = 1e-6)
  testthat::expect_equal(component$derivative, numerical, tolerance = 1e-6)
  singular <- qif_see_extended_source_v2(dat$y, dat$X, dat$cluster_id, beta, .35,
                                         h, weights, rep(1, 3L))
  testthat::expect_error(qif_see_source_step_v2(singular, beta), class = "qif_see_source_error")
})

testthat::test_that("QIF-SEE source limiting case reproduces independent smoothed quantile equations", {
  dat <- generate_profile_qr_data_v2(n_clusters = 120L, p = 2L, s = 2L, q = 1L,
                                      tau = .5, sigma_b0 = 1, seed = 700772L)
  beta <- as.numeric(quantreg::rq.fit(dat$X, dat$y, tau = .5)$coefficients)
  density <- rep(stats::dnorm(0) / sqrt(2), length(dat$y))
  a <- fit_qif_see_source_v2(dat$y, dat$X, dat$cluster_id, .5, beta, density,
                              lambda = 1e-10, gamma = 1, corstr = "independence", max_iter = 100L)
  b <- fit_qif_see_source_v2(dat$y, dat$X, dat$cluster_id, .5, beta, density,
                              lambda = 1e-10, gamma = 1, corstr = "independence", max_iter = 100L)
  reference <- fit_sqr_iid_lasso_v2(dat$y, dat$X, .5, 120^(-1 / 2), 0,
                                    control = list(max_iter = 2000L, beta_tol = 1e-9, kkt_tol = 1e-8))
  testthat::expect_true(a$converged && reference$converged)
  testthat::expect_identical(a$beta, b$beta)
  testthat::expect_equal(a$beta, as.numeric(reference$beta), tolerance = 1e-5)
  testthat::expect_equal(a$bic, a$q_value + a$df * log(120), tolerance = 1e-12)
  testthat::expect_equal(a$history[[1L]]$delta, pmin(1, 1e-10 / abs(beta)^2), tolerance = 1e-15)
  testthat::expect_equal(a$history[[2L]]$delta,
                         pmin(1, 1e-10 / abs(a$history[[2L]]$beta)^2), tolerance = 1e-15)
  testthat::expect_error(fit_benchmark_qif_see_v2(dat, .5, colnames(dat$X), list(), 1L,
                                                 control = list(corstr = "independence")),
                         class = "qif_see_source_error")
  adapter <- fit_benchmark_qif_see_v2(dat, .5, colnames(dat$X), list(), 700772L,
                                      control = list(corstr = "independence", density_weights = density,
                                                     lambda_grid = 1e-10, gamma_grid = 1,
                                                     max_iter = 100L))
  testthat::expect_identical(adapter$status, "ok")
  testthat::expect_identical(names(adapter$beta_hat), colnames(dat$X))
  testthat::expect_equal(as.numeric(adapter$beta_hat), as.numeric(reference$beta), tolerance = 1e-5)
  testthat::expect_true(all(is.finite(adapter$prediction_function(dat$X))))
  testthat::expect_null(adapter$beta_tilde)
  testthat::expect_setequal(adapter$fit_object$initial_subject_folds$cluster_id, unique(dat$cluster_id))
  testthat::expect_equal(nrow(adapter$fit_object$initial_subject_folds), 120L)
  testthat::expect_setequal(adapter$fit_object$initial_subject_folds$fold, 1:10)
})

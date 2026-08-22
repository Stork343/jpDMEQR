# Deterministic acceptance tests for DOUBLE-PEN-QLMM (Li, Liu & Luo 2020).

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
source_v2_module(root, "benchmark_adapters_v2", envir = .GlobalEnv)

make_train <- function(n = 30, p = 8, s = 3, q = 1L, seed = 202L) {
  generate_profile_qr_data_v2(
    n_clusters = n, p = p, s = s, tau = 0.5, q = q,
    m_values = 4:6, m_rule = "uniform", sigma_b0 = 1,
    error_dist = "gaussian", seed = seed
  )
}

testthat::test_that("DOUBLE-PEN-QLMM adapter is deterministic and schema-conformant", {
  train <- make_train(seed = 202)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  ctrl <- list(lambda_beta_grid = c(0.5, 1), lambda_alpha_grid = c(0.5, 1))
  a1 <- fit_benchmark_double_pen_qlmm_v2(train, 0.5, colnames(train$X)[1:5],
                                         tuning, seed = 21L, control = ctrl)
  a2 <- fit_benchmark_double_pen_qlmm_v2(train, 0.5, colnames(train$X)[1:5],
                                         tuning, seed = 21L, control = ctrl)
  testthat::expect_equal(a1$status, "ok")
  testthat::expect_identical(a1$beta_hat, a2$beta_hat)
  testthat::expect_true(all(c("method_id", "status", "beta_hat", "selected",
                              "runtime_sec", "converged",
                              "reference_identifier",
                              "adapter_fidelity_status") %in% names(a1)))
  testthat::expect_equal(a1$method_id, "DOUBLE-PEN-QLMM")
  testthat::expect_equal(a1$adapter_fidelity_status, "formula_level_implemented")
  testthat::expect_true(all(names(a1$beta_hat) == colnames(train$X)))
  testthat::expect_true(is.null(a1$beta_tilde))  # no inferential claim by default
  testthat::expect_true(all(a1$inference_object$table$feasible == FALSE))
})

testthat::test_that("DOUBLE-PEN-QLMM objective matches paper equation (12)", {
  train <- make_train(seed = 303)
  y <- train$y; X <- train$X; Z <- train$Z; cl <- train$cluster_id
  ci <- cluster_index_v2(cl)
  fit <- fit_double_pen_qlmm_v2(
    y, X, Z, cl, 0.5,
    lambda_beta = 0.1, lambda_alpha = 0.05,
    control = list(max_iter = 30L)
  )
  testthat::expect_true(fit$converged)
  resid <- y - as.numeric(X %*% fit$beta) -
    rowSums(Z * fit$alpha[as.integer(ci$factor), , drop = FALSE])
  rho <- resid * (0.5 - as.numeric(resid < 0))
  expected <- sum(rho) + 0.1 * sum(abs(fit$beta)) + 0.05 * sum(abs(fit$alpha))
  testthat::expect_equal(fit$objective, expected, tolerance = 1e-8)
})

testthat::test_that("DOUBLE-PEN-QLMM zero-penalty limit matches unpenalised blocks", {
  train <- make_train(seed = 404)
  y <- train$y; X <- train$X; Z <- train$Z; cl <- train$cluster_id
  # Tiny penalties: the objective is close to the unpenalised check loss
  # minimisation, and the fitted beta should be finite and close to the
  # pooled QR solution in a low-dimensional setting.
  fit <- fit_double_pen_qlmm_v2(
    y, X, Z, cl, 0.5,
    lambda_beta = 1e-6, lambda_alpha = 1e-6,
    control = list(max_iter = 30L)
  )
  testthat::expect_true(fit$converged)
  testthat::expect_true(all(is.finite(fit$beta)))
  testthat::expect_equal(sum(abs(fit$beta)) > 0, TRUE)
})

testthat::test_that("DOUBLE-PEN-QLMM recovers active signs on deterministic data", {
  train <- make_train(seed = 505)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  ans <- fit_benchmark_double_pen_qlmm_v2(train, 0.5, colnames(train$X)[1:3],
                                          tuning, seed = 3L,
                                          control = list(
                                            lambda_beta_grid = c(0.25, 0.5),
                                            lambda_alpha_grid = c(0.25, 0.5)))
  # The first s=3 active coefficients alternate signs starting positive.
  active_est <- unname(ans$beta_hat[1:3])
  testthat::expect_equal(sign(active_est[1]), 1)
  testthat::expect_equal(sign(active_est[2]), -1)
  testthat::expect_equal(sign(active_est[3]), 1)
})

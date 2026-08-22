# Acceptance tests for the patch-provided adapters:
# PROFILE-DQR, POOLED-QR-LASSO, LQMM, PROFILE-DQR-TRUE-SUPPORT.

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
source_v2_module(root, "benchmark_adapters_v2", envir = .GlobalEnv)

make_train <- function(n = 40, p = 12, s = 3, q = 1L, sigma_b0 = 1,
                       seed = 202L) {
  dat <- generate_profile_qr_data_v2(
    n_clusters = n, p = p, s = s, tau = 0.5, q = q,
    m_values = 3:5, m_rule = "uniform", sigma_b0 = sigma_b0,
    error_dist = "gaussian", seed = seed
  )
  dat
}

tuning <- list(h = 0.5, lambda_beta = 0.15, lambda_gamma = 1,
               mu_grid = c(0.3, 0.6, 1.2))

# ---------------------------------------------------------------------------
# PROFILE-DQR (proposed practical method)
# ---------------------------------------------------------------------------

testthat::test_that("PROFILE-DQR adapter is deterministic and schema-conformant", {
  train <- make_train(seed = 202)
  coords <- colnames(train$X)[1:4]
  a1 <- fit_benchmark_profile_dqr_v2(train, 0.5, coords, tuning, seed = 7L)
  a2 <- fit_benchmark_profile_dqr_v2(train, 0.5, coords, tuning, seed = 7L)
  testthat::expect_equal(a1$status, "ok")
  testthat::expect_identical(a1$beta_hat, a2$beta_hat)
  testthat::expect_identical(a1$beta_tilde, a2$beta_tilde)
  testthat::expect_true(all(c("method_id", "status", "beta_hat", "beta_tilde",
                              "se", "ci_lower", "ci_upper", "selected",
                              "inference_object", "runtime_sec",
                              "converged", "kkt_residual",
                              "profile_identity_error") %in% names(a1)))
  testthat::expect_equal(a1$method_id, "PROFILE-DQR")
  inf <- a1$inference_object$table
  testthat::expect_true(all(inf$feasible))
  testthat::expect_true(all(is.finite(inf$se)))
  # One-step correction is a small finite adjustment to beta_hat.
  testthat::expect_true(max(abs(inf$beta_tilde - inf$beta_hat)) < 2)
  # Cluster sandwich variance is positive.
  testthat::expect_true(all(inf$se > 0))
})

testthat::test_that("PROFILE-DQR matches its reference fit on a small design", {
  train <- make_train(seed = 303)
  coords <- colnames(train$X)[1:3]
  ans <- fit_benchmark_profile_dqr_v2(train, 0.5, coords, tuning, seed = 3L)
  # Round-3 reference: cluster self-normalised calibrated first stage with the
  # same fit-design lambda_0,n (METHOD_SPECIFICATION_ROUND3_AMENDMENT.md 1-3).
  lam0 <- lambda_0_n_v2(length(unique(train$cluster_id)), ncol(train$X))$lambda_0
  ref <- fit_profile_lasso_calibrated_v2(
    y = train$y, X = train$X, Z = train$Z,
    cluster_id = train$cluster_id, tau = 0.5,
    h = tuning$h, lambda_0_n = lam0,
    lambda_gamma = tuning$lambda_gamma,
    base_penalty_factor = rep(1, ncol(train$X)),
    control = list(fit_control = list(max_iter = 2000L, max_backtrack = 50L,
                                      beta_tol = 1e-7, kkt_normalized_tol = 1e-3))
  )
  testthat::expect_equal(as.numeric(ans$beta_hat),
                         as.numeric(ref$beta), tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# POOLED-QR-LASSO (classical negative control)
# ---------------------------------------------------------------------------

testthat::test_that("POOLED-QR-LASSO is deterministic and pooled", {
  train <- make_train(seed = 404)
  coords <- colnames(train$X)[1:4]
  a1 <- fit_benchmark_pooled_qr_lasso_v2(train, 0.5, coords, tuning, seed = 5L)
  a2 <- fit_benchmark_pooled_qr_lasso_v2(train, 0.5, coords, tuning, seed = 5L)
  testthat::expect_equal(a1$status, "ok")
  testthat::expect_identical(a1$beta_hat, a2$beta_hat)
  testthat::expect_true(all(c("method_id", "status", "beta_hat", "selected",
                              "prediction_function", "runtime_sec") %in% names(a1)))
  testthat::expect_equal(a1$method_id, "POOLED-QR-LASSO")
  # Cluster ids are ignored: permutation of rows within clusters must not
  # change the estimate.
  ord <- sample(seq_along(train$y))
  train2 <- train
  train2$y <- train$y[ord]
  train2$X <- train$X[ord, , drop = FALSE]
  train2$cluster_id <- train$cluster_id[ord]
  a3 <- fit_benchmark_pooled_qr_lasso_v2(train2, 0.5, coords, tuning, seed = 5L)
  testthat::expect_equal(unname(a1$beta_hat), unname(a3$beta_hat),
                         tolerance = 1e-8)
  # Prediction function works on new rows.
  pred <- a1$prediction_function(train$X[1:5, , drop = FALSE])
  testthat::expect_length(pred, 5)
  testthat::expect_true(all(is.finite(pred)))
})

testthat::test_that("POOLED-QR-LASSO has no naive interval by default", {
  train <- make_train(seed = 505)
  ans <- fit_benchmark_pooled_qr_lasso_v2(train, 0.5, colnames(train$X)[1:3],
                                          tuning, seed = 5L)
  testthat::expect_null(ans$beta_tilde)
  testthat::expect_null(ans$se)
})

# ---------------------------------------------------------------------------
# LQMM (classical direct comparator)
# ---------------------------------------------------------------------------

testthat::test_that("LQMM adapter is deterministic on a low-dimensional set", {
  train <- make_train(seed = 606)
  # LQMM is a low-dimensional method: restrict the design.
  keep <- c(1:4, 7)  # 5 features including active and null
  sub <- train
  sub$X <- train$X[, keep, drop = FALSE]
  sub$beta0 <- train$beta0[keep]
  coords <- colnames(sub$X)[1:3]
  a1 <- fit_benchmark_lqmm_v2(sub, 0.5, coords, tuning, seed = 9L,
                              random_slope = FALSE)
  a2 <- fit_benchmark_lqmm_v2(sub, 0.5, coords, tuning, seed = 9L,
                              random_slope = FALSE)
  testthat::expect_equal(a1$status, "ok")
  testthat::expect_identical(a1$beta_hat, a2$beta_hat)
  testthat::expect_true(all(c("method_id", "status", "beta_hat", "selected",
                              "runtime_sec", "converged") %in% names(a1)))
  testthat::expect_equal(a1$method_id, "LQMM")
  testthat::expect_true(all(is.finite(a1$beta_hat)))
})

testthat::test_that("LQMM limiting case: negligible random-effect variance approaches pooled QR", {
  train <- make_train(seed = 707, sigma_b0 = 0.02)
  keep <- 1:5
  sub <- train
  sub$X <- train$X[, keep, drop = FALSE]
  sub$beta0 <- train$beta0[keep]
  coords <- colnames(sub$X)[1:3]
  ans <- fit_benchmark_lqmm_v2(sub, 0.5, coords, tuning, seed = 9L,
                               random_slope = FALSE)
  testthat::expect_equal(ans$status, "ok")
  # With negligible cluster variance the LQMM fixed effects stay finite and
  # reasonably close to a pooled quantile fit.
  pooled <- quantreg::rq.fit(x = sub$X, y = sub$y, tau = 0.5)
  testthat::expect_true(max(abs(unname(ans$beta_hat[coords]) -
                                  as.numeric(pooled$coefficients[coords]))) < 1.5)
})

# ---------------------------------------------------------------------------
# PROFILE-DQR-TRUE-SUPPORT (oracle selection)
# ---------------------------------------------------------------------------

testthat::test_that("TRUE-SUPPORT uses the oracle support and maps back to the full design", {
  train <- make_train(seed = 808)
  coords <- colnames(train$X)[1:4]
  ans <- fit_benchmark_profile_true_support_v2(
    train, 0.5, coords, tuning, seed = 11L,
    active = train$active
  )
  testthat::expect_equal(ans$method_id, "PROFILE-DQR-TRUE-SUPPORT")
  testthat::expect_equal(ans$status, "ok")
  # Full-design vector returned with zeros outside the oracle support.
  testthat::expect_true(all(names(ans$beta_hat) == colnames(train$X)))
  testthat::expect_true(all(is.finite(ans$beta_tilde)))
  # All active coordinates are in the support.
  testthat::expect_true(all(train$active %in% ans$selected))
  # Oracle label is mandatory in metadata.
  testthat::expect_true("PROFILE-DQR-TRUE-SUPPORT" %in% ans$method_id)
})

testthat::test_that("TRUE-SUPPORT requires the active set", {
  train <- make_train(seed = 909)
  testthat::expect_error(
    fit_benchmark_profile_true_support_v2(train, 0.5,
                                          colnames(train$X)[1:3],
                                          tuning, seed = 11L),
    "active"
  )
})

# Deterministic acceptance tests for QIF-SEE (Bhattacharya et al. 2026).

root <- if (file.exists(file.path("scripts", "00_source_v2.R"))) "." else
  normalizePath(file.path(testthat::test_path(), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
source_v2_module(root, "benchmark_adapters_v2", envir = .GlobalEnv)

make_train <- function(n = 40, p = 10, s = 3, seed = 202L) {
  generate_profile_qr_data_v2(
    n_clusters = n, p = p, s = s, tau = 0.5, q = 1L,
    m_values = 4:6, m_rule = "uniform", sigma_b0 = 1,
    error_dist = "gaussian", seed = seed
  )
}

testthat::test_that("QIF-SEE adapter is deterministic and schema-conformant", {
  train <- make_train(seed = 202)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  ctrl <- list(lambda_grid = c(0.1, 0.3, 0.5), gamma = 1)
  a1 <- fit_benchmark_qif_see_v2(train, 0.5, colnames(train$X)[1:5],
                                 tuning, seed = 21L, control = ctrl)
  a2 <- fit_benchmark_qif_see_v2(train, 0.5, colnames(train$X)[1:5],
                                 tuning, seed = 21L, control = ctrl)
  testthat::expect_equal(a1$status, "ok")
  testthat::expect_identical(a1$beta_hat, a2$beta_hat)
  testthat::expect_identical(a1$selected, a2$selected)
  testthat::expect_true(all(c("method_id", "status", "beta_hat", "selected",
                              "runtime_sec", "converged",
                              "reference_identifier",
                              "adapter_fidelity_status") %in% names(a1)))
  testthat::expect_equal(a1$method_id, "QIF-SEE")
  testthat::expect_equal(a1$adapter_fidelity_status,
                         "formula_level_implemented")
  testthat::expect_true(all(names(a1$beta_hat) == colnames(train$X)))
  testthat::expect_true(is.null(a1$beta_tilde))  # no coverage claim by default
  testthat::expect_true(all(a1$inference_object$table$feasible == FALSE))
})

testthat::test_that("QIF-SEE selects active variables on deterministic data", {
  train <- make_train(seed = 303)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  ans <- fit_benchmark_qif_see_v2(train, 0.5, colnames(train$X)[1:3],
                                  tuning, seed = 3L,
                                  control = list(lambda_grid = c(0.1, 0.3),
                                                 gamma = 1))
  # The SEE variable-selection mechanism should retain the true active set
  # on a clean deterministic DGP (possibly with extra noise coordinates).
  testthat::expect_true(all(1:3 %in% ans$selected))
  active_est <- unname(ans$beta_hat[1:3])
  for (kk in which(active_est != 0)) {
    expected_sign <- if (kk %% 2L == 1L) 1 else -1
    testthat::expect_equal(sign(active_est[kk]), expected_sign)
  }
})

testthat::test_that("QIF-SEE is invariant to cluster permutation", {
  train <- make_train(seed = 404)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  ctrl <- list(lambda_grid = c(0.1, 0.3), gamma = 1)
  a1 <- fit_benchmark_qif_see_v2(train, 0.5, colnames(train$X)[1:3],
                                 tuning, seed = 9L, control = ctrl)
  # Permute the order of clusters (rows and identifiers together).
  ord <- order(runif(length(train$cluster_id)))
  train2 <- train
  train2$y <- train$y[ord]
  train2$X <- train$X[ord, , drop = FALSE]
  train2$Z <- train$Z[ord, , drop = FALSE]
  train2$cluster_id <- train$cluster_id[ord]
  a2 <- fit_benchmark_qif_see_v2(train2, 0.5, colnames(train2$X)[1:3],
                                 tuning, seed = 9L, control = ctrl)
  # The QIF extended-score construction is permutation invariant, but the
  # rqPen initial estimator draws internal CV folds from the row order, so
  # the initial beta_tilde can differ slightly under permutation. Require
  # the same selected set and matching signs/direction of the retained
  # active coefficients rather than exact numerical equality.
  testthat::expect_identical(a1$selected, a2$selected)
  testthat::expect_equal(sign(a1$beta_hat[a1$selected]),
                         sign(a2$beta_hat[a2$selected]))
})

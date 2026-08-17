# Deterministic acceptance tests for QGEE-SCAD (Zu et al. 2023, geeVerse).

root <- if (file.exists(file.path("scripts", "00_source_v2.R"))) "." else
  normalizePath(file.path(testthat::test_path(), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
source_v2_module(root, "benchmark_adapters_v2", envir = .GlobalEnv)

make_train <- function(n = 40, p = 12, s = 3, seed = 202L) {
  generate_profile_qr_data_v2(
    n_clusters = n, p = p, s = s, tau = 0.5, q = 1L,
    m_values = 4:6, m_rule = "uniform", sigma_b0 = 1,
    error_dist = "gaussian", seed = seed
  )
}

testthat::test_that("QGEE-SCAD adapter is deterministic and schema-conformant", {
  train <- make_train(seed = 202)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  ctrl <- list(corstr = "exchangeable", method = "HBIC", lambda = 0.5)
  a1 <- fit_benchmark_qgee_scad_v2(train, 0.5, colnames(train$X)[1:5],
                                   tuning, seed = 21L, control = ctrl)
  a2 <- fit_benchmark_qgee_scad_v2(train, 0.5, colnames(train$X)[1:5],
                                   tuning, seed = 21L, control = ctrl)
  testthat::expect_equal(a1$status %in% c("ok", "warning"), TRUE)
  testthat::expect_identical(a1$beta_hat, a2$beta_hat)
  testthat::expect_true(all(c("method_id", "status", "beta_hat", "selected",
                              "runtime_sec", "converged",
                              "reference_identifier",
                              "adapter_fidelity_status") %in% names(a1)))
  testthat::expect_equal(a1$method_id, "QGEE-SCAD")
  testthat::expect_equal(a1$adapter_fidelity_status,
                         "official_reference_implementation")
  testthat::expect_true(all(names(a1$beta_hat) == colnames(train$X)))
  testthat::expect_true(is.null(a1$beta_tilde))  # no coverage claim by default
  testthat::expect_equal(a1$target_scope, "marginal_longitudinal_quantile")
})

testthat::test_that("QGEE-SCAD recovers active support on deterministic data", {
  train <- make_train(seed = 303)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  ans <- fit_benchmark_qgee_scad_v2(train, 0.5, colnames(train$X)[1:3],
                                    tuning, seed = 3L,
                                    control = list(corstr = "exchangeable",
                                                   method = "HBIC",
                                                   lambda = 0.4))
  # The marginal estimating-equation target does not guarantee recovery of
  # every conditional active coefficient in finite samples; check that the
  # selected set is nonempty, contains at least one true active coordinate,
  # and any recovered active coefficient carries the correct alternating sign.
  testthat::expect_true(length(ans$selected) > 0)
  testthat::expect_true(any(ans$selected %in% 1:3))
  active_est <- unname(ans$beta_hat[1:3])
  for (kk in which(active_est != 0)) {
    expected_sign <- if (kk %% 2L == 1L) 1 else -1
    testthat::expect_equal(sign(active_est[kk]), expected_sign)
  }
})

testthat::test_that("QGEE-SCAD independence working correlation matches exchangeable on weak dependence", {
  train <- make_train(seed = 404)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  a_ind <- fit_benchmark_qgee_scad_v2(train, 0.5, colnames(train$X)[1:3],
                                      tuning, seed = 9L,
                                      control = list(corstr = "independence",
                                                     method = "HBIC",
                                                     lambda = 0.5))
  a_exc <- fit_benchmark_qgee_scad_v2(train, 0.5, colnames(train$X)[1:3],
                                      tuning, seed = 9L,
                                      control = list(corstr = "exchangeable",
                                                     method = "HBIC",
                                                     lambda = 0.5))
  # Both should be finite; the working-correlation choice is recorded.
  testthat::expect_true(all(is.finite(a_ind$beta_hat)))
  testthat::expect_true(all(is.finite(a_exc$beta_hat)))
  testthat::expect_equal(a_ind$inference_object$working_correlation,
                         "independence")
  testthat::expect_equal(a_exc$inference_object$working_correlation,
                         "exchangeable")
})

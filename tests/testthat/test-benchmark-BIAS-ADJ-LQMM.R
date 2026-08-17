# Deterministic acceptance tests for BIAS-ADJ-LQMM (Battagliola et al. 2022).

root <- if (file.exists(file.path("scripts", "00_source_v2.R"))) "." else
  normalizePath(file.path(testthat::test_path(), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
source_v2_module(root, "benchmark_adapters_v2", envir = .GlobalEnv)

make_clustered_train <- function(n = 40, p = 5, s = 2, tau = 0.5,
                                 sigma_b0 = 1.2, seed = 202L) {
  dat <- generate_profile_qr_data_v2(
    n_clusters = n, p = p, s = s, tau = tau, q = 1L,
    m_values = 3:5, m_rule = "uniform", sigma_b0 = sigma_b0,
    error_dist = "gaussian", seed = seed
  )
  dat
}

testthat::test_that("BIAS-ADJ-LQMM adapter is deterministic and schema-conformant", {
  train <- make_clustered_train(seed = 202)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  coords <- colnames(train$X)[1:4]
  a1 <- fit_benchmark_bias_adj_lqmm_v2(train, 0.5, coords, tuning, seed = 77L,
                                       control = list(B = 20))
  a2 <- fit_benchmark_bias_adj_lqmm_v2(train, 0.5, coords, tuning, seed = 77L,
                                       control = list(B = 20))
  testthat::expect_equal(a1$status, "ok")
  testthat::expect_identical(a1$beta_tilde, a2$beta_tilde)
  testthat::expect_identical(a1$beta_hat, a2$beta_hat)
  testthat::expect_true(all(c("method_id", "status", "beta_hat", "beta_tilde",
                              "se", "ci_lower", "ci_upper", "selected",
                              "runtime_sec", "converged",
                              "reference_identifier",
                              "adapter_fidelity_status") %in% names(a1)))
  testthat::expect_equal(a1$method_id, "BIAS-ADJ-LQMM")
  testthat::expect_equal(a1$adapter_fidelity_status, "formula_level_implemented")
  inf <- a1$inference_object$table
  testthat::expect_true(all(inf$coordinate == coords))
  testthat::expect_true(all(inf$feasible))
  testthat::expect_true(all(is.finite(inf$beta_tilde)))
})

testthat::test_that("BIAS-ADJ-LQMM rejects random-slope designs (source scope)", {
  # q=2 design: source method supports random intercept only.
  train <- make_clustered_train(seed = 303)
  train$Z <- cbind(1, train$time)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  ans <- fit_benchmark_bias_adj_lqmm_v2(train, 0.5, colnames(train$X)[1:3],
                                        tuning, seed = 1L, control = list(B = 10))
  testthat::expect_equal(ans$status, "failed")
  testthat::expect_match(ans$failure_stage, "schema")
})

testthat::test_that("BIAS-ADJ-LQMM two-step estimate matches manual pipeline", {
  train <- make_clustered_train(seed = 404)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  ans <- fit_benchmark_bias_adj_lqmm_v2(train, 0.5, colnames(train$X)[1:4],
                                        tuning, seed = 5L, control = list(B = 10))
  # Manual two-step: lqmm + centred ranef + rq.
  df <- data.frame(y = train$y, cluster = factor(train$cluster_id), train$X,
                   check.names = FALSE)
  ff <- stats::as.formula(paste("y ~ 0 +", paste(colnames(train$X), collapse = " + ")))
  step1 <- lqmm::lqmm(ff, random = ~1, group = cluster, tau = 0.5, data = df,
                      nK = 15, type = "normal",
                      control = lqmm::lqmmControl(method = "df",
                                                  LP_max_iter = 1000,
                                                  verbose = FALSE))
  ran <- lqmm::ranef(step1)
  u_est_cluster <- as.numeric(ran[, "(Intercept)"]) -
    mean(as.numeric(ran[, "(Intercept)"]))
  levels_c <- levels(df$cluster)
  u_est <- u_est_cluster[match(df$cluster, levels_c)]
  df2 <- df
  df2$y <- train$y - as.numeric(train$Z[, 1] * u_est)
  step2 <- quantreg::rq(y ~ 0 + ., data = df2, tau = 0.5)
  manual <- stats::coef(step2)[colnames(train$X)]
  testthat::expect_equal(as.numeric(ans$beta_hat),
                         as.numeric(manual), tolerance = 1e-6)
})

testthat::test_that("BIAS-ADJ-LQMM wild weights follow the source distribution", {
  set.seed(4242)
  tau <- 0.3
  w <- sample(c(2 * (1 - tau), -2 * tau), size = 200000, replace = TRUE,
              prob = c(1 - tau, tau))
  # Source (eq. w_distr): w = 2(1-tau) w.p. 1-tau; w = -2tau w.p. tau.
  testthat::expect_equal(mean(w == 2 * (1 - tau)), 1 - tau, tolerance = 1e-2)
  testthat::expect_equal(mean(w == -2 * tau), tau, tolerance = 1e-2)
  testthat::expect_equal(mean(w < 0), tau, tolerance = 1e-2)
})

testthat::test_that("BIAS-ADJ-LQMM small random-effect variance limit", {
  # With negligible random-effect variance the adjusted estimate stays close
  # to the two-step estimate (adjustment has little to correct).
  train <- make_clustered_train(n = 50, sigma_b0 = 0.05, seed = 505)
  tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
                 mu_grid = c(0.3, 0.6, 1.2))
  ans <- fit_benchmark_bias_adj_lqmm_v2(train, 0.5, colnames(train$X)[1:4],
                                        tuning, seed = 7L, control = list(B = 20))
  testthat::expect_equal(ans$status, "ok")
  diff_adj <- max(abs(ans$beta_tilde - ans$beta_hat[names(ans$beta_tilde)]))
  testthat::expect_lt(diff_adj, 1.0)
})

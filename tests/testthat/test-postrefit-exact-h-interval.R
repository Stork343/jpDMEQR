# Regression test: POSTREFIT-EXACT-H interval width must be ~3.92*se.

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
source_v2_module(root, "benchmark_adapters_v2", envir = .GlobalEnv)

make_tiny_data <- function() {
  set.seed(20260825)
  n_cluster <- 20L
  m <- rep(4L, n_cluster)
  cluster_id <- rep(sprintf("c%02d", seq_len(n_cluster)), m)
  N <- length(cluster_id)
  X <- matrix(stats::rnorm(N * 10L), N, 10L, dimnames = list(NULL, paste0("x", 1:10)))
  Z <- matrix(1, N, 1L, dimnames = list(NULL, "intercept"))
  beta0 <- c(0.75, -0.75, 0.5, rep(0, 7))
  b0 <- stats::rnorm(n_cluster, sd = 1)
  y <- as.numeric(X %*% beta0 + b0[match(cluster_id, unique(cluster_id))] +
                    stats::rnorm(N))
  list(y = y, X = X, Z = Z, cluster_id = cluster_id, beta0 = beta0)
}

testthat::test_that("POSTREFIT-EXACT-H intervals satisfy width = 3.92*se", {
  train <- make_tiny_data()
  n <- length(unique(train$cluster_id))
  lam0 <- lambda_0_n_v2(n, ncol(train$X))$lambda_0
  tuning <- list(h = 0.5 / n^(0) * n^(-1/3), h_est = 0.4, lambda_0_n = lam0,
                 lambda_gamma = 1, mu_grid = c(0.2, 0.5, 1))
  tuning$h <- 0.35
  ans <- fit_benchmark_postrefit_exact_h_v2(
    train, 0.5, colnames(train$X)[1:3], tuning, seed = 7L
  )
  testthat::expect_equal(ans$status, "ok")
  tab <- ans$inference_object$table
  testthat::expect_true(all(tab$feasible))
  testthat::expect_equal(tab$ci_upper - tab$ci_lower,
                         2 * stats::qnorm(0.975) * tab$se, tolerance = 1e-12)
})
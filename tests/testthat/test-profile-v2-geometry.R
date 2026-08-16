root <- if (file.exists(file.path("scripts", "00_source_v2.R"))) "." else normalizePath(file.path(testthat::test_path(), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)

testthat::test_that("v2 profile score and Hessian match finite differences", {
  set.seed(20260817)
  n_cluster <- 12L
  m <- rep(4L, n_cluster)
  cluster_id <- rep(sprintf("c%02d", seq_len(n_cluster)), m)
  N <- length(cluster_id)
  X <- matrix(stats::rnorm(N * 5L), N, 5L)
  Z <- cbind(1, rep(seq(-1, 1, length.out = 4L), n_cluster))
  beta <- c(0.25, -0.15, 0.10, 0, 0)
  b <- matrix(stats::rnorm(n_cluster * 2L, sd = 0.4), n_cluster, 2L)
  y <- as.numeric(X %*% beta + rowSums(Z * b[match(cluster_id, unique(cluster_id)), ]) +
                    stats::rnorm(N))
  out <- validate_profile_geometry_v2(
    y, X, Z, cluster_id, beta = beta + 0.02,
    tau = 0.5, h = 0.45, lambda_gamma = diag(c(0.8, 1.2)), eps = 2e-5
  )
  testthat::expect_true(out$nuisance_converged)
  testthat::expect_lt(out$gradient_max_error, 2e-4)
  testthat::expect_lt(out$hessian_max_error, 3e-3)
  testthat::expect_lt(out$profile_identity_error, 1e-9)
  testthat::expect_lt(out$max_nuisance_gradient, 1e-7)
})

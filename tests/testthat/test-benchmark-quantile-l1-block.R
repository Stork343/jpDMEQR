l1_block_root <- v2_require_root()
source(file.path(l1_block_root, "scripts", "00_source_v2.R"), local = TRUE)
for (module in c("profile_v2", "simulation_v2", "benchmark_adapters_v2")) {
  source_v2_module(l1_block_root, module, environment())
}

testthat::test_that("quantile L1 blocks match an independent explicit convex programme", {
  dat <- generate_profile_qr_data_v2(n_clusters = 25L, p = 3L, s = 2L,
                                     tau = 0.3, q = 1L, m_values = 3:4,
                                     seed = 20261115L)
  penalty <- c(0.7, 1.1, 0.4)
  fit <- fit_quantile_l1_block_v2(dat$X, dat$y, 0.3, penalty)
  beta <- CVXR::Variable(ncol(dat$X))
  residual <- dat$y - dat$X %*% beta
  objective <- CVXR::sum_entries(0.3 * CVXR::pos(residual) +
                                  0.7 * CVXR::pos(-residual)) +
    CVXR::p_norm(diag(penalty) %*% beta, 1)
  problem <- CVXR::Problem(CVXR::Minimize(objective))
  reference <- CVXR::psolve(problem, solver = "CLARABEL", abstol = 1e-9,
                            reltol = 1e-9, feastol = 1e-9)
  reference_beta <- as.numeric(reference$getValue(beta))
  actual_residual <- dat$y - as.numeric(dat$X %*% fit$coefficients)
  actual_objective <- sum(actual_residual * (0.3 - (actual_residual < 0))) +
    sum(penalty * abs(fit$coefficients))
  testthat::expect_equal(actual_objective, reference$value, tolerance = 1e-6)
  testthat::expect_equal(as.numeric(fit$coefficients), reference_beta, tolerance = 1e-5)
})

testthat::test_that("scalar random intercept blocks retain their declared L1 penalty", {
  set.seed(20261116L)
  response <- 2 + stats::rnorm(41L)
  design <- matrix(1, 41L, 1L, dimnames = list(NULL, "intercept"))
  penalty <- 1.37
  fit <- fit_quantile_l1_block_v2(design, response, 0.4, penalty)
  expected <- sort(response)[ceiling(length(response) * 0.4 - penalty)]
  testthat::expect_gt(expected, 0)
  testthat::expect_equal(as.numeric(fit$coefficients), expected, tolerance = 1e-5)
  plain <- fit_quantile_l1_block_v2(design, response, 0.4, 0)
  testthat::expect_equal(as.numeric(plain$coefficients),
                         sort(response)[ceiling(length(response) * 0.4)], tolerance = 1e-5)
  testthat::expect_lt(fit$coefficients, plain$coefficients)
})

testthat::test_that("quantreg vector loading convention matches the declared objective", {
  dat <- generate_profile_qr_data_v2(n_clusters = 25L, p = 2L, s = 1L,
                                     tau = 0.5, q = 1L, m_values = 3:4,
                                     seed = 20261117L)
  penalty <- c(0.7, 1.2)
  fit <- fit_quantile_l1_block_v2(dat$X, dat$y, 0.5, penalty)
  reference <- quantreg::rq.fit.lasso(dat$X, dat$y, tau = 0.5, lambda = 2 * penalty)
  testthat::expect_equal(as.numeric(fit$coefficients),
                         as.numeric(reference$coefficients), tolerance = 1e-5)
})

paired_profile_root <- v2_require_root()
source(file.path(paired_profile_root, "scripts", "00_source_v2.R"), local = TRUE)
source_v2_module(paired_profile_root, "profile_v2", environment())

paired_profile_problem <- function() {
  values <- c(-1.37, -0.91, -0.52, -0.29, -0.08, 0.11,
              0.26, 0.43, 0.67, 0.94, 1.31)
  list(
    y = rep(values, each = 2L),
    X = matrix(1, nrow = 2L * length(values), ncol = 1L,
               dimnames = list(NULL, "intercept")),
    Z = matrix(1, nrow = 2L * length(values), ncol = 1L),
    cluster_id = rep(seq_along(values), each = 2L),
    beta = 0.03,
    tau = 0.5,
    h = 0.24
  )
}

paired_profile_evaluate <- function(problem, beta = problem$beta) {
  profile_bandwidth_pair_v2(
    problem$y, problem$X, problem$Z, problem$cluster_id,
    beta, problem$tau, problem$h,
    identity_check = TRUE
  )
}

testthat::test_that("paired score and Hessian differentiate the paired objective", {
  problem <- paired_profile_problem()
  result <- paired_profile_evaluate(problem)
  step <- 1e-4
  plus <- paired_profile_evaluate(problem, problem$beta + step)
  minus <- paired_profile_evaluate(problem, problem$beta - step)
  objective_derivative <- (plus$objective - minus$objective) / (2 * step)
  score_derivative <- (plus$score - minus$score) / (2 * step)

  # 使用既有几何检查的数值容差。
  testthat::expect_lt(abs(objective_derivative - result$score), 2e-4)
  testthat::expect_lt(abs(score_derivative - result$hessian[1, 1]), 3e-3)
  testthat::expect_lt(result$profile_identity_error, 1e-10)
  total <- Reduce(`+`, result$hessian_fold_sums)
  testthat::expect_equal(total / result$n_clusters, result$hessian,
                         tolerance = 1e-12)
})

testthat::test_that("Wald variance includes covariance between bandwidth scores", {
  result <- paired_profile_evaluate(paired_profile_problem())
  precision <- solve(result$hessian)
  interval <- profile_bandwidth_pair_wald_v2(result, 0.03, precision, 1L)
  small <- as.numeric(result$components_half_h$cluster_scores %*% precision)
  large <- as.numeric(result$components_h$cluster_scores %*% precision)
  small <- small - mean(small)
  large <- large - mean(large)
  exact_variance <- mean((2 * small - large)^2)
  separated_variance <- 4 * mean(small^2) + mean(large^2)

  testthat::expect_equal(interval$se^2 * result$n_clusters, exact_variance,
                         tolerance = 1e-12)
  testthat::expect_gt(abs(exact_variance - separated_variance), 0.01)
  testthat::expect_lt(interval$inverse_defect, 1e-12)
  testthat::expect_equal(interval$beta_tilde,
                         0.03 - as.numeric(precision * result$score))
})

testthat::test_that("paired inference rejects incomplete coordinates", {
  result <- paired_profile_evaluate(paired_profile_problem())
  testthat::expect_error(
    profile_bandwidth_pair_wald_v2(result, 0.03, matrix(1, 1, 1), 1.5),
    "coordinates"
  )
  testthat::expect_error(
    profile_bandwidth_pair_wald_v2(result, 0.03, matrix(1, 1, 1), 1L, 1),
    "alpha"
  )
  testthat::expect_error(
    profile_bandwidth_pair_wald_v2(result, 0.03, matrix(1, 1, 1), Inf),
    "coordinates"
  )
})

testthat::test_that("paired calculus retains the random slope profile", {
  source_v2_module(paired_profile_root, "simulation_v2", environment())
  previous <- Sys.getenv("JPDMEQR_CLUSTER_SOLVER", unset = NA_character_)
  on.exit({
    if (is.na(previous)) Sys.unsetenv("JPDMEQR_CLUSTER_SOLVER") else
      Sys.setenv(JPDMEQR_CLUSTER_SOLVER = previous)
  }, add = TRUE)
  Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "newton")
  dat <- generate_profile_qr_data_v2(
    n_clusters = 20L, p = 4L, s = 2L, tau = 0.25,
    q = 2L, m_values = 3:5, seed = 20261022L
  )
  beta <- c(0.15, -0.12, 0.03, 0.07)
  evaluate <- function(b) profile_bandwidth_pair_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, b, 0.25, 0.3,
    lambda_gamma = matrix(c(1, 0.2, 0.2, 0.8), 2, 2),
    identity_check = TRUE
  )
  result <- evaluate(beta)
  step <- 1e-4
  plus <- beta
  minus <- beta
  plus[2] <- plus[2] + step
  minus[2] <- minus[2] - step
  forward <- evaluate(plus)
  backward <- evaluate(minus)

  testthat::expect_lt(
    abs((forward$objective - backward$objective) / (2 * step) -
          result$score[2]), 2e-4
  )
  testthat::expect_lt(
    max(abs((forward$score - backward$score) / (2 * step) -
              result$hessian[, 2])), 3e-3
  )
  testthat::expect_lt(result$profile_identity_error, 1e-8)
  testthat::expect_equal(result$hessian, t(result$hessian))
  testthat::expect_identical(nrow(result$components_h$gamma), 20L)
  testthat::expect_identical(ncol(result$components_half_h$gamma), 2L)
})

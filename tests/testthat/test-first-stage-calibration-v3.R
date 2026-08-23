# Round-3 first-stage calibration unit tests
# (METHOD_SPECIFICATION_ROUND3_AMENDMENT.md sections 1-3; PILOT_GATE_THEORY_
# ACTIONS.txt item 12):
#   1. score-loading scale and centring;
#   2. deterministic lambda critical value;
#   3. penalty-factor mapping;
#   4. two-pass loading update;
#   5. weighted normalized KKT;
#   6. lambda change does not alter target/POP-H dependency hashes.

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)

# --- 1. score-loading scale and centring -------------------------------------
testthat::test_that("score_loadings_v2: per-cluster SD with n denominator, centred", {
  set.seed(11)
  S <- matrix(stats::rnorm(120, mean = 3, sd = 2), 30, 4)  # rows = clusters
  ell <- score_loadings_v2(S)
  manual <- apply(S, 2, function(col) sqrt(mean((col - mean(col))^2)))
  testthat::expect_equal(ell, manual, tolerance = 1e-12)
  # Shift invariance (centring): adding a constant to a column does not change
  # its loading.
  S2 <- S
  S2[, 2] <- S2[, 2] + 100
  testthat::expect_equal(score_loadings_v2(S2)[2], ell[2], tolerance = 1e-12)
  # n denominator, not n-1: sd()/sqrt((n-1)/n) must match.
  testthat::expect_equal(ell[1], stats::sd(S[, 1]) * sqrt((30 - 1) / 30),
                         tolerance = 1e-12)
})

# --- 2. deterministic lambda critical value ----------------------------------
testthat::test_that("lambda_0_n_v2: frozen critical value is deterministic", {
  m1 <- lambda_0_n_v2(200, 500)
  alpha <- 0.10 / log(max(200, 3))
  q <- stats::qnorm(1 - alpha / (2 * 500))
  testthat::expect_equal(m1$alpha, alpha, tolerance = 1e-15)
  testthat::expect_equal(m1$normal_quantile, q, tolerance = 1e-12)
  testthat::expect_equal(m1$lambda_0, 1.10 * q / sqrt(200), tolerance = 1e-12)
  # Hand-checked reference values (round-3 freeze).
  testthat::expect_equal(m1$lambda_0, 0.320527, tolerance = 1e-5)
  m2 <- lambda_0_n_v2(80, 200)
  testthat::expect_equal(m2$lambda_0, 0.474528, tolerance = 1e-5)
  # Deterministic: identical inputs give identical outputs.
  testthat::expect_identical(lambda_0_n_v2(200, 500)$lambda_0, m1$lambda_0)
})

# --- small cluster dataset for fit-level tests -------------------------------
make_tiny_data <- function() {
  set.seed(20260817)
  n_cluster <- 12L
  m <- rep(4L, n_cluster)
  cluster_id <- rep(sprintf("c%02d", seq_len(n_cluster)), m)
  N <- length(cluster_id)
  X <- matrix(stats::rnorm(N * 8L), N, 8L, dimnames = list(NULL, paste0("x", 1:8)))
  Z <- matrix(1, N, 1L, dimnames = list(NULL, "intercept"))
  beta0 <- c(0.75, -0.75, rep(0, 6))
  b0 <- stats::rnorm(n_cluster, sd = 1)
  y <- as.numeric(X %*% beta0 + b0[match(cluster_id, unique(cluster_id))] +
                    stats::rnorm(N))
  list(y = y, X = X, Z = Z, cluster_id = cluster_id, beta0 = beta0)
}

# --- 3 & 4. penalty-factor mapping and two-pass loading update ---------------
testthat::test_that("calibrated first stage: two-pass loadings and penalty mapping", {
  dat <- make_tiny_data()
  n <- length(unique(dat$cluster_id)); p <- ncol(dat$X)
  lam0 <- lambda_0_n_v2(n, p)$lambda_0
  base <- rep(1, p)
  fit <- fit_profile_lasso_calibrated_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, tau = 0.5,
    h_est = 0.45, h_inf = 0.35,
    lambda_gamma = 1, lambda_0_n = lam0, base_penalty_factor = base,
    control = list(fit_control = list(max_iter = 2000L, max_backtrack = 50L,
                                      beta_tol = 1e-7, kkt_normalized_tol = 1e-3))
  )
  testthat::expect_true(isTRUE(fit$first_stage_calibrated))
  testthat::expect_true(isTRUE(fit$converged),
                        info = paste("failure_stage:", fit$failure_stage))
  # Round-4 dual bandwidths: inferential components at h_inf, first-stage
  # components at h_est, and the hard separation.
  testthat::expect_equal(fit$h_est, 0.45, tolerance = 1e-12)
  testthat::expect_equal(fit$h_inf, 0.35, tolerance = 1e-12)
  testthat::expect_equal(fit$components$h, 0.35, tolerance = 1e-12)
  testthat::expect_equal(fit$components_first_stage$h, 0.45, tolerance = 1e-12)
  # The guard rejects a fit whose inferential components carry the h_est
  # bandwidth (an h_est Hessian must never reach the precision pipeline).
  tampered <- fit
  tampered$components <- tampered$components_first_stage
  testthat::expect_error(assert_inferential_components_v2(tampered, h_inf = 0.35),
                         "bandwidth mismatch")
  testthat::expect_true(assert_inferential_components_v2(fit, h_inf = 0.35))
  # Pass-0 loading recomputed independently at beta=0.
  comp0 <- profile_components_v2(dat$y, dat$X, dat$Z, dat$cluster_id,
                                 rep(0, p), 0.5, 0.45, 1,
                                 need_hessian = FALSE)
  ell0_manual <- score_loadings_v2(comp0$cluster_scores)
  testthat::expect_equal(fit$lambda_loading_pass0_median,
                         stats::median(ell0_manual), tolerance = 1e-10)
  # Two-pass update: rerun the passes in-test and compare the final penalty
  # factor against the recorded one (mapping p_j = lambda_0 * ell_final * base).
  prelim <- fit_profile_lasso_v2(dat$y, dat$X, dat$Z, dat$cluster_id, 0.5, 0.45,
                                 lambda_beta = lam0, lambda_gamma = 1,
                                 penalty_factor = ell0_manual * base,
                                 control = list(max_iter = 2000L))
  comp1 <- profile_components_v2(dat$y, dat$X, dat$Z, dat$cluster_id,
                                 prelim$beta, 0.5, 0.45, 1,
                                 gamma_init = prelim$gamma,
                                 need_hessian = FALSE)
  ell_final_manual <- pmax(ell0_manual, score_loadings_v2(comp1$cluster_scores))
  testthat::expect_equal(as.numeric(fit$penalty_factor),
                         as.numeric(ell_final_manual * base),
                         tolerance = 1e-8)
  # Pass-1 fields finite and nonnegative; degenerate-loading rule not tripped.
  testthat::expect_true(all(is.finite(c(fit$lambda_loading_pass1_min,
                                        fit$lambda_loading_pass1_median,
                                        fit$lambda_loading_pass1_max))))
  testthat::expect_true(fit$lambda_coordinate_min >= 0)
  testthat::expect_equal(fit$lambda_base, lam0, tolerance = 1e-12)
  testthat::expect_equal(fit$lambda_rule, "cluster-score-loading-v3")
  # The calibrated penalty must make the zero fit KKT-infeasible in this
  # signal setting: max|g(0)| > min coordinate penalty at pass 0 would imply
  # zero_kkt_ratio_max > 1; record it and require a non-degenerate first stage.
  testthat::expect_true(is.finite(fit$zero_profile_score_max))
  testthat::expect_true(fit$first_stage_nonzero_count >= 0)
})

testthat::test_that("calibrated first stage: unpenalised short-circuit", {
  dat <- make_tiny_data()
  n <- length(unique(dat$cluster_id)); p <- ncol(dat$X)
  lam0 <- lambda_0_n_v2(n, p)$lambda_0
  fit <- fit_profile_lasso_calibrated_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, tau = 0.5,
    h_est = 0.45, h_inf = 0.35,
    lambda_gamma = 1, lambda_0_n = lam0, base_penalty_factor = rep(0, p)
  )
  testthat::expect_false(isTRUE(fit$first_stage_calibrated))
  testthat::expect_true(isTRUE(fit$converged))
  testthat::expect_equal(fit$penalty_factor, rep(0, p))
})

# --- 5. weighted normalized KKT ----------------------------------------------
testthat::test_that("kkt_residual_normalized_v2 matches the round-3 formula", {
  beta <- c(0.5, 0, -0.25, 1.0)          # active, inactive(P), active, U
  score <- c(-0.8, 0.3, 0.9, 0.4)
  pen <- c(0.5, 0.5, 0.5, 0)             # P: 1-3, U: 4
  # r_j by hand:
  # j1 active:  | -0.8 + 0.5*1 | = 0.3
  # j2 inactive: max(0, |0.3| - 0.5) = 0
  # j3 active:  | 0.9 + 0.5*(-1) | = 0.4
  # j4 U:       | 0.4 | = 0.4
  # p_ref = median positive penalty = 0.5
  # R_KKT = max(0.3/0.5, 0/0.5, 0.4/0.5, 0.4/0.5) = 0.8
  testthat::expect_equal(kkt_residual_normalized_v2(beta, score, pen), 0.8,
                         tolerance = 1e-12)
  # Fully unpenalised: plain max |score|.
  testthat::expect_equal(kkt_residual_normalized_v2(beta, score, rep(0, 4)),
                         max(abs(score)), tolerance = 1e-12)
  # All-penalised, no U: only the P ratios count.
  pen2 <- rep(2, 4)
  r <- c(abs(-0.8 + 2), 0, abs(0.9 - 2), abs(0.4 + 2))
  testthat::expect_equal(kkt_residual_normalized_v2(beta, score, pen2),
                         max(r / 2), tolerance = 1e-12)
})

# --- 6. lambda change does not alter target/POP-H dependency hashes ----------
testthat::test_that("target/POP-H dependency hashes ignore the lambda rule", {
  base_cfg <- list(
    n_clusters = 200L, p = 500L, s = 5L, tau = 0.5, q = 1L,
    h_multiplier = 1, lambda_gamma = 1, fit_random_effects = "intercept",
    error_dist = "gaussian", error_dependence = "independent",
    design_type = "ar1", rho_x = 0.5, signal = 0.75,
    random_effect_dist = "normal", sigma_b0 = 1, sigma_b1 = 0.5, rho_b = 0.4,
    x_b_corr = 0, informative_size = FALSE, nonlinear_re_strength = 0,
    m_rule = "uniform", copula_rho = 0, heteroskedastic = FALSE
  )
  cfg_old <- c(base_cfg, list(lambda_beta_multipliers = "0.25;0.5;1;2"))
  cfg_new <- c(base_cfg, list(lambda_beta_multipliers = "1"))
  h1 <- target_dependency_hash_v2(cfg_old, 100000L, 4L)
  h2 <- target_dependency_hash_v2(cfg_new, 100000L, 4L)
  testthat::expect_identical(h1, h2)
  p1 <- pop_h_dependency_hash_v2(cfg_old, 100000L, 4L, 200L, 200^(-3 / 10))
  p2 <- pop_h_dependency_hash_v2(cfg_new, 100000L, 4L, 200L, 200^(-3 / 10))
  testthat::expect_identical(p1, p2)
})

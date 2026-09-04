# Round-5 post-refit unit tests (METHOD_SPECIFICATION_ROUND5_AMENDMENT.md;
# PILOT_GATE_THEORY_ACTIONS.txt item 13):
#   S_R construction and forced target inclusion;
#   refit dimension gate;
#   zero-L1 reduced profile refit at h_inf;
#   no h_est/reduced Hessian leakage into primary full-p inference;
#   post-refit numerical acceptance and explicit failure paths;
#   target/POP-H dependency hashes unchanged by the practical post-refit rule.

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)

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

# --- 1 & 2. S_R construction, forced target inclusion, dimension gate --------
testthat::test_that("post-refit: S_R contains targets, gate is explicit", {
  dat <- make_tiny_data()
  n <- length(unique(dat$cluster_id)); p <- ncol(dat$X)
  lam0 <- lambda_0_n_v2(n, p)$lambda_0
  fit <- fit_profile_lasso_calibrated_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, tau = 0.5,
    h_est = 0.45, h_inf = 0.35, lambda_gamma = 1, lambda_0_n = lam0,
    target_coordinates = c("x1", "x6")
  )
  testthat::expect_true(isTRUE(fit$converged),
                        info = paste("failure_stage:", fit$failure_stage))
  testthat::expect_true(all(c("x1", "x6") %in% fit$refit_set_names))
  testthat::expect_true(isTRUE(fit$refit_contains_targets))
  testthat::expect_true(fit$refit_set_size >= 2L)
  testthat::expect_equal(fit$post_refit_status, "ok")
  # beta_refit is embedded in the full p-vector; beta_l1 retained separately.
  testthat::expect_length(fit$beta_refit, p)
  testthat::expect_length(fit$beta_l1, p)
  testthat::expect_equal(fit$beta, fit$beta_refit)
  # Dimension gate: an absurd target set that pushes d_R over floor(n/log n).
  big_targets <- colnames(dat$X)   # all 8 coordinates forced in
  fit2 <- fit_profile_lasso_calibrated_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, tau = 0.5,
    h_est = 0.45, h_inf = 0.35, lambda_gamma = 1, lambda_0_n = lam0,
    target_coordinates = big_targets
  )
  gate <- floor(n / log(max(n, 3)))
  if (8L > gate) {  # n=12 -> gate = floor(12/log(12)) = floor(4.83) = 4 < 8
    testthat::expect_false(isTRUE(fit2$converged))
    testthat::expect_equal(fit2$failure_stage, "post_refit_dimension")
    testthat::expect_equal(fit2$post_refit_status, "dimension_gate_failed")
  }
})

# --- 3. zero-L1 reduced refit at h_inf: unpenalised optimum ------------------
testthat::test_that("post-refit: reduced model is the unpenalised h_inf optimum", {
  dat <- make_tiny_data()
  n <- length(unique(dat$cluster_id)); p <- ncol(dat$X)
  lam0 <- lambda_0_n_v2(n, p)$lambda_0
  fit <- fit_profile_lasso_calibrated_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, tau = 0.5,
    h_est = 0.45, h_inf = 0.35, lambda_gamma = 1, lambda_0_n = lam0,
    target_coordinates = c("x1", "x6")
  )
  testthat::expect_true(isTRUE(fit$converged))
  S_R <- fit$refit_set
  # Independent zero-L1 refit on the same set reproduces beta_refit[S_R].
  ref <- fit_profile_lasso_v2(
    dat$y, dat$X[, S_R, drop = FALSE], dat$Z, dat$cluster_id, 0.5, 0.35,
    lambda_beta = 0, lambda_gamma = 1, penalty_factor = rep(0, length(S_R)),
    control = list(max_iter = 1000L, max_backtrack = 50L, beta_tol = 1e-8,
                   kkt_normalized_tol = 1e-7)
  )
  testthat::expect_equal(as.numeric(fit$beta[S_R]), as.numeric(ref$beta),
                         tolerance = 1e-6)
  testthat::expect_true(all(fit$beta[-S_R] == 0))
})

# --- 4. no h_est / reduced Hessian leakage into primary full-p inference -----
testthat::test_that("post-refit: full-p h_inf components are authoritative", {
  dat <- make_tiny_data()
  n <- length(unique(dat$cluster_id)); p <- ncol(dat$X)
  lam0 <- lambda_0_n_v2(n, p)$lambda_0
  fit <- fit_profile_lasso_calibrated_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, tau = 0.5,
    h_est = 0.45, h_inf = 0.35, lambda_gamma = 1, lambda_0_n = lam0,
    target_coordinates = c("x1", "x6")
  )
  testthat::expect_true(isTRUE(fit$converged))
  # Primary components: full p x p Hessian at h_inf at beta_refit.
  testthat::expect_equal(fit$components$h, 0.35, tolerance = 1e-12)
  testthat::expect_equal(dim(fit$components$hessian), c(p, p))
  # First-stage diagnostic components: h_est.
  testthat::expect_equal(fit$components_first_stage$h, 0.45, tolerance = 1e-12)
  # Guard rejects h_est components in the inferential pipeline.
  tampered <- fit
  tampered$components <- tampered$components_first_stage
  testthat::expect_error(assert_inferential_components_v2(tampered, 0.35),
                         "bandwidth mismatch")
  testthat::expect_true(assert_inferential_components_v2(fit, 0.35))
})

# --- 5. post-refit acceptance failure path -----------------------------------
testthat::test_that("post-refit: acceptance failure is explicit, no fallback", {
  dat <- make_tiny_data()
  n <- length(unique(dat$cluster_id)); p <- ncol(dat$X)
  lam0 <- lambda_0_n_v2(n, p)$lambda_0
  # Force a post-refit acceptance failure by requiring an unattainable
  # gradient tolerance (the reduced profile gradient cannot reach 1e-12 at
  # this n in the available budget).
  fit <- fit_profile_lasso_calibrated_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, tau = 0.5,
    h_est = 0.45, h_inf = 0.35, lambda_gamma = 1, lambda_0_n = lam0,
    target_coordinates = c("x1", "x6"),
    control = list(refit_control = list(max_iter = 2L, kkt_normalized_tol = 1e-12))
  )
  testthat::expect_false(isTRUE(fit$converged))
  testthat::expect_equal(fit$failure_stage, "post_refit")
  testthat::expect_true(fit$post_refit_status %in%
                          c("acceptance_failed", "error"))
})

# --- 6. target/POP-H dependency hashes unchanged by the post-refit rule ------
testthat::test_that("dependency hashes ignore the practical starting-estimator rule", {
  base_cfg <- list(
    n_clusters = 200L, p = 500L, s = 5L, tau = 0.5, q = 1L,
    h_multiplier = 1, lambda_gamma = 1, fit_random_effects = "intercept",
    error_dist = "gaussian", error_dependence = "independent",
    design_type = "ar1", rho_x = 0.5, signal = 0.75,
    random_effect_dist = "normal", sigma_b0 = 1, sigma_b1 = 0.5, rho_b = 0.4,
    x_b_corr = 0, informative_size = FALSE, nonlinear_re_strength = 0,
    m_rule = "uniform", copula_rho = 0, heteroskedastic = FALSE,
    methods = "PROFILE-DQR|POSTREFIT-EXACT-H"
  )
  cfg_alt <- base_cfg
  cfg_alt$methods <- "PROFILE-DQR"   # post-refit rule is not a hash input
  testthat::expect_identical(
    target_dependency_hash_v2(base_cfg, 100000L, 4L),
    target_dependency_hash_v2(cfg_alt, 100000L, 4L)
  )
  testthat::expect_identical(
    pop_h_dependency_hash_v2(base_cfg, 100000L, 4L, 200L, 200^(-3 / 10)),
    pop_h_dependency_hash_v2(cfg_alt, 100000L, 4L, 200L, 200^(-3 / 10))
  )
})

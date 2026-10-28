pair_fit_root <- v2_require_root()
source(file.path(pair_fit_root, "scripts", "00_source_v2.R"), local = TRUE)
source_v2_module(pair_fit_root, "profile_v2", environment())
source_v2_module(pair_fit_root, "simulation_v2", environment())
withr::local_envvar(c(JPDMEQR_UNPENALISED_OPTIMIZER = "proximal_gradient",
                      JPDMEQR_PROXIMAL_SCALING = "none",
                      JPDMEQR_HESSIAN_REFRESH = "1", JPDMEQR_IDENTITY_DETAIL = "cluster",
                      JPDMEQR_DANTZIG_SOLVER = "clarabel_native"),
                    .local_envir = environment())

pair_fit_data <- function(q = 1L, p = 8L, n = 24L, seed = 20261023L) {
  generate_profile_qr_data_v2(
    n_clusters = n, p = p, s = 2L, tau = 0.25, q = q,
    m_values = 3:5, seed = seed
  )
}

pair_fit_real <- function(dat, targets = c("x00001", "x00004"), Lambda = 1,
                          control = list(), penalty_factor = NULL) {
  fit_profile_pair_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, tau = 0.25,
    h_inf = dat$n_clusters^(-3 / 10), lambda_gamma = Lambda,
    target_coordinates = targets, control = control,
    penalty_factor = penalty_factor
  )
}

testthat::test_that("paired fitting retains both penalty components and both loading passes", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  dat <- pair_fit_data()
  fit <- pair_fit_real(dat)
  stochastic <- lambda_0_n_v2(dat$n_clusters, ncol(dat$X))$lambda_0
  testthat::expect_true(fit$converged)
  testthat::expect_identical(fit$method_version, "profile_pair_20261022")
  testthat::expect_identical(fit$score_definition,
                            "profile_bandwidth_pair_2_half_minus_full")
  testthat::expect_equal(fit$lambda_stochastic, stochastic)
  testthat::expect_equal(fit$lambda_mean, fit$h_est)
  testthat::expect_equal(fit$lambda_base, stochastic + fit$h_est)
  testthat::expect_equal(fit$loadings$final,
                         pmax(fit$loadings$pass0, fit$loadings$pass1))
  testthat::expect_equal(fit$lambda_components$coordinate,
                         fit$lambda_base * fit$loadings$final)
  testthat::expect_lte(fit$final_kkt_normalized, 1e-3)
  testthat::expect_lte(fit$post_refit_gradient_max, 1e-7)
  testthat::expect_lte(fit$post_refit_nuisance_gradient_max, 1e-7)
  testthat::expect_gt(fit$post_refit_hessian_min_eigenvalue, 1e-8)
  testthat::expect_lt(fit$post_refit_hessian_condition_number, 1e10)
  testthat::expect_true(all(c("x00001", "x00004") %in% fit$refit_set_names))
  testthat::expect_equal(fit$beta, fit$beta_refit)
  testthat::expect_equal(fit$components_first_stage$h, fit$h_est)
  testthat::expect_equal(fit$components$components_half_h$h, fit$h_inf / 2)
  testthat::expect_equal(fit$components$components_h$h, fit$h_inf)
  testthat::expect_equal(fit$components$cluster_scores,
                         2 * fit$components$components_half_h$cluster_scores -
                           fit$components$components_h$cluster_scores)
})

testthat::test_that("paired inference uses paired fold curvature and its score covariance", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton",
                        JPDMEQR_DANTZIG_SOLVER = "clarabel_native"))
  dat <- pair_fit_data(q = 2L, p = 26L, n = 20L, seed = 20261024L)
  Lambda <- matrix(c(1, 0.2, 0.2, 0.8), 2L)
  fit <- pair_fit_real(dat, Lambda = Lambda)
  result <- infer_profile_pair_v2(fit, c("x00001", "x00004"), c(0.05, 0.15, 0.4))
  testthat::expect_gt(ncol(dat$X), dat$n_clusters)
  testthat::expect_true(all(result$table$feasible))
  testthat::expect_true(all(result$table$solver == "CLARABEL-native"))
  testthat::expect_true(all(result$table$solver_static_regularization_constant == 1e-6))
  testthat::expect_true(all(result$table$mu_source == "fold_cv"))
  testthat::expect_equal(fit$components$components_h$Lambda, Lambda)
  testthat::expect_equal(fit$components$components_half_h$Lambda, Lambda)
  testthat::expect_equal(Reduce(`+`, fit$components$hessian_fold_sums) /
                           dat$n_clusters, fit$components$hessian,
                         tolerance = 1e-11)
  for (j in seq_len(nrow(result$table))) {
    row <- result$table[j, ]
    omega <- result$directions[[j]]$omega
    G <- fit$components$cluster_scores
    projection <- as.numeric(G %*% omega)
    variance <- mean((projection - mean(projection))^2)
    testthat::expect_equal(row$se^2 * dat$n_clusters, variance,
                           tolerance = 1e-10)
    testthat::expect_equal(row$beta_tilde,
                           as.numeric(fit$beta_refit[row$index] -
                                        sum(omega * colMeans(G))), tolerance = 1e-10)
    testthat::expect_lte(row$dantzig_residual, row$mu + 1e-6)
    cv <- result$cv[[j]]
    eligible <- cv$candidates[cv$candidates$eligible, ]
    testthat::expect_equal(row$mu_cv_defect, min(eligible$defect_mean))
    testthat::expect_true(all(cv$folds$solver == "CLARABEL-native"))
    testthat::expect_true(all(cv$folds$static_regularization_constant == 1e-6))
    testthat::expect_equal(nrow(cv$folds), 2L * 3L)
  }
  # 同一 coefficient 上独立重算两个 profile，验证完整设计的计算位置。
  component <- profile_bandwidth_pair_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, fit$beta_refit,
    0.25, fit$h_inf, lambda_gamma = Lambda
  )
  testthat::expect_equal(component$hessian, fit$components$hessian,
                         tolerance = 1e-10)
  testthat::expect_equal(component$score, fit$components$score,
                         tolerance = 1e-10)
})

testthat::test_that("always included covariates retain zero penalty", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  dat <- pair_fit_data(p = 5L, n = 20L, seed = 20261025L)
  fit <- pair_fit_real(dat, targets = "x00004", penalty_factor = c(0, 1, 1, 1, 1))
  testthat::expect_equal(fit$lambda_components$coordinate[1], c(x00001 = 0))
  testthat::expect_true("x00001" %in% fit$refit_set_names)
  plain <- pair_fit_real(dat, targets = "x00004", penalty_factor = rep(0, 5L))
  testthat::expect_false(plain$first_stage_calibrated)
  testthat::expect_equal(plain$lambda_base, 0)
  testthat::expect_null(plain$beta_l1)
  testthat::expect_equal(plain$refit_set_names, colnames(dat$X))
  interval <- infer_profile_pair_v2(plain, "x00004", c(0.05, 0.15))
  testthat::expect_true(is.finite(interval$table$se))
})

testthat::test_that("paired failure conditions preserve the attempted stage", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  dat <- pair_fit_data(p = 8L, n = 12L, seed = 20261026L)
  failure <- tryCatch(pair_fit_real(dat, colnames(dat$X)),
                      profile_pair_error = identity)
  testthat::expect_s3_class(failure, "profile_pair_error")
  testthat::expect_identical(failure$failure_stage, "post_refit_dimension")
  testthat::expect_equal(failure$fit$refit_set_names, colnames(dat$X))
  testthat::expect_equal(failure$fit$lambda_mean, failure$fit$h_est)
  testthat::expect_error(
    pair_fit_real(dat, targets = "unknown_coordinate"), "Unknown coordinates"
  )
  testthat::expect_error(
    fit_profile_pair_v2(dat$y, dat$X, dat$Z, dat$cluster_id, 0.25, 0.4,
                        mean_score_multiplier = 0, target_coordinates = "x00001"),
    "positive scalars"
  )
})

testthat::test_that("paired entry points preserve numerical and score definition failures", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  dat <- pair_fit_data(p = 5L, n = 20L, seed = 20261027L)
  fit <- pair_fit_real(dat)
  invalid <- fit
  invalid$components <- fit$components_first_stage
  testthat::expect_error(infer_profile_pair_v2(invalid, "x00001", 0.15),
                         "converged profile_pair")
  invalid <- fit
  invalid$beta[1] <- invalid$beta[1] + 0.01
  testthat::expect_error(infer_profile_pair_v2(invalid, "x00001", 0.15),
                         "stored beta_refit")
  testthat::expect_error(infer_profile_pair_v2(fit, "x00005", 0.15),
                         "prespecified")
  testthat::expect_error(infer_profile_pair_v2(fit, 1.5, 0.15), "integer indices")
  degenerate <- tryCatch(infer_profile_pair_v2(fit, "x00001", c(1, 1.2)),
                         profile_pair_error = identity)
  testthat::expect_s3_class(degenerate, "profile_pair_error")
  testthat::expect_identical(degenerate$failure_stage, "precision_solver")
  testthat::expect_equal(degenerate$diagnostics$mu_grid_admissible, numeric())
  testthat::expect_equal(degenerate$diagnostics$mu_grid_requested, c(1, 1.2))
  testthat::expect_equal(degenerate$diagnostics$excluded_zero_direction$mu, c(1, 1.2))
  testthat::expect_true(all(degenerate$diagnostics$excluded_zero_direction$status ==
                             "excluded_zero_direction"))
  testthat::expect_equal(degenerate$fit$method_version, fit$method_version)
  mixed <- infer_profile_pair_v2(fit, "x00001", c(0.05, 0.15, 1, 1.2))
  admissible <- infer_profile_pair_v2(fit, "x00001", c(0.05, 0.15))
  testthat::expect_equal(mixed$table, admissible$table)
  testthat::expect_equal(mixed$mu_grid_requested, c(0.05, 0.15, 1, 1.2))
  testthat::expect_equal(mixed$mu_grid_admissible, c(0.05, 0.15))
  testthat::expect_equal(mixed$excluded_zero_direction$mu, c(1, 1.2))
  testthat::expect_true(all(mixed$cv[[1]]$folds$mu < 1))
  testthat::expect_equal(nrow(mixed$cv[[1]]$folds), 4L)
  testthat::expect_equal(mixed$cv[[1]]$excluded_zero_direction,
                         mixed$excluded_zero_direction)
  failed_refit <- tryCatch(
    pair_fit_real(dat, control = list(refit_control = list(max_iter = 1L))),
    profile_pair_error = identity
  )
  testthat::expect_s3_class(failed_refit, "profile_pair_error")
  testthat::expect_identical(failed_refit$failure_stage, "post_refit")
  testthat::expect_true(length(failed_refit$fit$refit_set) > 0L)
  testthat::expect_true(!is.null(failed_refit$diagnostics$attempted_fit))
  failed_profile <- tryCatch(
    profile_pair_run_stage_v2(
      "profile_geometry", fit,
      profile_bandwidth_pair_v2(
        dat$y, dat$X, dat$Z, dat$cluster_id, fit$beta_refit,
        0.25, fit$h_inf, lambda_gamma = fit$Lambda,
        nuisance_control = list(newton_maxit = 1L, maxit = 1L, grad_tol = 1e-30)
      )
    ), profile_pair_error = identity
  )
  testthat::expect_s3_class(failed_profile, "profile_pair_error")
  testthat::expect_identical(failed_profile$failure_stage, "profile_geometry")
  testthat::expect_s3_class(failed_profile$diagnostics$cause, "error")
  testthat::expect_equal(failed_profile$fit$beta_refit, fit$beta_refit)
  failed_solver <- tryCatch(
    infer_profile_pair_v2(fit, "x00001", c(0.05, 0.15),
                           solver_opts = list(control = list(max_iter = 1L))),
    profile_pair_error = identity
  )
  testthat::expect_s3_class(failed_solver, "profile_pair_error")
  testthat::expect_identical(failed_solver$failure_stage, "precision_solver")
  testthat::expect_identical(failed_solver$diagnostics$coordinate, 1L)
  testthat::expect_equal(failed_solver$diagnostics$mu, 0.05)
  testthat::expect_equal(failed_solver$diagnostics$fold, 1L)
  testthat::expect_true(is.data.frame(failed_solver$diagnostics$candidates))
  testthat::expect_true(is.data.frame(failed_solver$diagnostics$folds))
  testthat::expect_s3_class(failed_solver$diagnostics$cause, "profile_pair_error")
})

testthat::test_that("prospective penalty multipliers apply before coordinate loadings", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  dat <- pair_fit_data(p = 5L, n = 20L, seed = 20261028L)
  fit <- fit_profile_pair_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, 0.25, 20^(-3 / 10),
    mean_score_multiplier = 2, lambda_multiplier = 0.8,
    target_coordinates = c("x00001", "x00004")
  )
  testthat::expect_equal(fit$lambda_mean, 2 * fit$h_est)
  testthat::expect_equal(fit$lambda_beta,
                         0.8 * (fit$lambda_stochastic + 2 * fit$h_est))
  testthat::expect_equal(fit$lambda_components$coordinate,
                         fit$lambda_beta * fit$loadings$final)
  testthat::expect_identical(fit$lambda_components$mean_score_constant_status,
                            "prospective_model_class_condition")
})

testthat::test_that("paired CV retains the four fold rule at two hundred clusters", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  dat <- pair_fit_data(p = 4L, n = 200L, seed = 20261029L)
  component <- profile_bandwidth_pair_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, rep(0, 4L),
    0.25, 200^(-3 / 10)
  )
  cv <- profile_pair_mu_cv_v2(component, 1L, c(0.1, 0.3))
  testthat::expect_equal(cv$n_folds, 4L)
  testthat::expect_equal(nrow(cv$folds), 8L)
  testthat::expect_true(all(cv$folds$n_validation == 50L))
  testthat::expect_true(all(cv$folds$n_training == 150L))
  testthat::expect_true(all(cv$folds$solver == "CLARABEL-native"))
})

testthat::test_that("fixed solver controls certify rank deficient training folds", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  dat <- generate_profile_qr_data_v2(
    n_clusters = 30L, p = 40L, s = 3L, tau = 0.25, q = 2L,
    m_values = 3:5, seed = 698518065L
  )
  cluster <- match(as.character(dat$cluster_id), rownames(dat$random_effects))
  known <- rowSums(dat$Z * dat$random_effects[cluster, , drop = FALSE])
  fit <- fit_profile_pair_v2(
    dat$y - known, dat$X, matrix(0, length(dat$y), 1L), dat$cluster_id,
    0.25, 30^(-3 / 10),
    target_coordinates = c("x00001", "x00002", "x00004", "x00005")
  )
  grid <- c(0.02, 0.05, 0.1, 0.25, 0.5, 1, 2, 4) *
    (sqrt(log(40) / (30 * fit$h_inf)) + fit$h_inf^2)
  result <- infer_profile_pair_v2(fit, "x00002", grid)
  folds <- result$cv[[1]]$folds
  testthat::expect_true(any(!folds$feasible))
  testthat::expect_true(all(folds$static_regularization_constant == 1e-6))
  testthat::expect_true(all(is.na(folds$defect[!folds$feasible])))
  testthat::expect_true(all(result$table$feasible))
  testthat::expect_true(is.finite(result$table$se))
})

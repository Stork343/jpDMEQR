pair_accounting_root <- normalizePath(v2_require_root(), mustWork = TRUE)
root <- pair_accounting_root
source(file.path(root, "scripts", "00_source_v2.R"), local = TRUE)
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, environment())
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"), local = TRUE)

testthat::test_that("real failed paired inference retains requested coordinates and coverage denominators", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton",
                       JPDMEQR_DANTZIG_SOLVER = "clarabel_native",
                       JPDMEQR_ASSEMBLY = "sequential", JPDMEQR_CV_CORES = "1",
                       JPDMEQR_HESSIAN_REFRESH = "1", JPDMEQR_SHARED_FIRST_STAGE = "FALSE"))
  registry <- read_simulation_registry_v2(file.path(root, "config", "simulation_preflight_profile_pair_20261022.csv"))
  config <- config_row_to_list_v2(registry[registry$experiment_id == "PAIR01", , drop = FALSE])
  config$experiment_id <- "PAIR-ACCOUNTING-REAL-FAILURE"
  config$n_clusters <- 12L
  config$p <- 8L
  config$target_coordinate_count <- 8L
  config$methods <- "PROFILE-DQR"
  result <- run_one_replication_v2(root, config, 1L, base_seed = 20261028L,
                                   final = FALSE, run_id = "real_failure_accounting_test")
  requested <- result$method_diagnostics[["PROFILE-DQR"]]$requested_target_coordinates
  coordinates <- result$coordinates
  testthat::expect_identical(result$metrics$status, "failed")
  testthat::expect_gt(length(requested), floor(config$n_clusters / log(config$n_clusters)))
  testthat::expect_setequal(coordinates$coordinate, requested)
  testthat::expect_equal(nrow(coordinates), length(requested))
  testthat::expect_true(all(coordinates$status == "failed"))
  testthat::expect_true(all(is.na(coordinates$beta_hat) & is.na(coordinates$beta_tilde)))
  testthat::expect_true(all(is.na(coordinates$estimated_se) & is.na(coordinates$p_value)))
  testthat::expect_true(all(is.na(coordinates$covered) & !coordinates$feasible))
  testthat::expect_true(validate_paired_failed_coordinate_accounting_v2(result$metrics, coordinates))
  summary <- coordinate_coverage_summary_v2(coordinates)
  testthat::expect_equal(summary$coverage_denominator, length(requested))
  testthat::expect_equal(summary$coverage, 0)
  testthat::expect_equal(summary$evaluable, 0)
  testthat::expect_true(is.na(summary$coverage_successful))
  testthat::expect_equal(summary$failed_inference_coordinates, length(requested))
  testthat::expect_true(all(is.na(coordinates$covered)))
})

testthat::test_that("failed source intervals retain every requested coordinate", {
  registry <- read_simulation_registry_v2(file.path(root, "config", "simulation_pair_preflight_20261022.yml"))
  config <- config_row_to_list_v2(registry[registry$experiment_id == "PF01", , drop = FALSE])
  config$experiment_id <- "SOURCE-INTERVAL-UNDERIDENTIFIED"
  config$n_clusters <- 4L
  config$p <- 40L
  config$m_values <- 2:3
  config$methods <- "BIAS-ADJ-LQMM"
  result <- run_one_replication_v2(root, config, 1L, base_seed = 2026102203L,
                                   final = FALSE, run_id = "source_interval_accounting_test")
  requested <- result$method_diagnostics[["BIAS-ADJ-LQMM"]]$requested_target_coordinates
  testthat::expect_identical(result$metrics$status, "failed")
  testthat::expect_equal(result$metrics$n_clusters, 2L)
  testthat::expect_lt(result$metrics$total_rows, result$metrics$p)
  testthat::expect_setequal(result$coordinates$coordinate, requested)
  testthat::expect_true(all(is.na(result$coordinates$beta_tilde)))
  testthat::expect_true(all(is.na(result$coordinates$estimated_se)))
  testthat::expect_true(all(result$coordinates$status == "failed"))
  summary <- coordinate_coverage_summary_v2(result$coordinates)
  testthat::expect_equal(summary$coverage_denominator, length(requested))
  testthat::expect_equal(summary$evaluable, 0L)
  testthat::expect_equal(summary$coverage, 0)
})

testthat::test_that("real paired population means center matching score diagnostics", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton",
                       JPDMEQR_DANTZIG_SOLVER = "clarabel_native",
                       JPDMEQR_ASSEMBLY = "sequential", JPDMEQR_CV_CORES = "1"))
  config <- list(
    experiment_id = "PAIR-CENTERING-REAL-TEST", module = "paired_development",
    n_clusters = 30L, p = 1L, s = 1L, tau = 0.5, q = 1L,
    m_values = 3:4, m_rule = "uniform", m_multiset = NULL, design_type = "ar1",
    rho_x = 0.3, signal = 0.75, error_dist = "gaussian", heteroskedastic = FALSE,
    error_dependence = "independent", copula_rho = 0, random_effect_dist = "normal",
    sigma_b0 = 0.8, sigma_b1 = 0.5, rho_b = 0.4, x_b_corr = 0,
    informative_size = FALSE, nonlinear_re_strength = 0, fit_random_effects = "intercept",
    lambda_gamma = 1, h_multiplier = 1, target_mode = "structural"
  )
  target <- build_profile_pair_target_asset_v2(
    config, 60L, 2L, seed = 20261029L, development = TRUE
  )
  population <- build_profile_pair_direction_asset_v2(
    config, target, 60L, 2L, seed = 20271029L, development = TRUE
  )
  dat <- profile_pair_population_data_v2(config, config$n_clusters, 20281029L)
  h <- config$h_multiplier * config$n_clusters^(-3 / 10)
  centering <- paired_population_score_centering_v2(population, dat$beta0, h, config)
  testthat::expect_true(centering$matched)
  testthat::expect_identical(centering$mean, population$score_mean_population)
  testthat::expect_identical(paired_population_score_centering_v2(
    population, dat$beta0 + 1e-12, h, config)$reason, "population_score_beta_mismatch")
  testthat::expect_identical(paired_population_score_centering_v2(
    population, dat$beta0, h / 2, config)$reason, "population_score_bandwidth_mismatch")
  other_config <- config
  other_config$lambda_gamma <- 2
  testthat::expect_identical(paired_population_score_centering_v2(
    population, dat$beta0, h, other_config)$reason, "population_score_design_mismatch")
  fit <- fit_benchmark_profile_pair_pop_h_v2(
    dat, config$tau, names(dat$beta0),
    tuning = list(h_inf = h, h_est = 0.5, lambda_gamma = 1), seed = 20291029L,
    asset = population,
    control = list(final_assets = FALSE, penalty_factor = 0,
                   implementation_commit = population$implementation_commit),
    config = config
  )
  comp_target <- profile_bandwidth_pair_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, dat$beta0, config$tau, h, 1,
    nuisance_control = list(grad_tol = 1e-8), need_fold_sums = FALSE
  )
  fit$target_cluster_scores <- comp_target$cluster_scores
  fit$population_score_centering <- centering
  fit$delta_full <- fit$fit_object$beta - dat$beta0
  nm <- names(dat$beta0)[1]
  diagnostic <- pop_direction_diagnostics_v2(
    fit, population, nm, sum(abs(fit$delta_full)), config$n_clusters, config$tau,
    beta_star = dat$beta0[[nm]]
  )
  omega <- population$directions[[nm]]$omega
  projected <- as.numeric(comp_target$cluster_scores %*% omega)
  expected <- (sqrt(config$n_clusters) * (fit$beta_tilde[[nm]] - dat$beta0[[nm]]) +
    sum(projected - sum(population$score_mean_population * omega)) / sqrt(config$n_clusters)) /
    population$sigma0_pop[[nm]]
  testthat::expect_equal(diagnostic$bahadur_scaled, expected, tolerance = 1e-12)
  centered_scores <- sweep(comp_target$cluster_scores, 2, centering$mean)
  ladder <- variance_ladder_one_v2(dat, fit, dat$beta0, population, nm, config$tau, 1,
                                   Gs = centered_scores, h = h)
  testthat::expect_equal(ladder$L_k, -sum(centered_scores %*% omega) / sqrt(config$n_clusters),
                         tolerance = 1e-12)
  testthat::expect_equal(ladder$meat1, as.numeric(crossprod(omega, population$Sigma_population %*% omega)),
                         tolerance = 1e-12)
  wald <- coordinate_wald_statistics_v2(fit$beta_tilde[[nm]], fit$se[[nm]], dat$beta0[[nm]])
  testthat::expect_equal(wald$wald_z, (fit$beta_tilde[[nm]] - dat$beta0[[nm]]) / fit$se[[nm]])
  testthat::expect_equal(wald$p_value, 2 * stats::pnorm(-abs(fit$beta_tilde[[nm]] / fit$se[[nm]])))
  testthat::expect_gt(abs(wald$p_value - 2 * stats::pnorm(-abs(wald$wald_z))), 1e-6)
  fit$population_score_centering <- paired_population_score_centering_v2(population, dat$beta0, h / 2, config)
  testthat::expect_true(is.na(pop_direction_diagnostics_v2(
    fit, population, nm, sum(abs(fit$delta_full)), config$n_clusters, config$tau,
    beta_star = dat$beta0[[nm]])$bahadur_scaled))
  target_obj <- list(beta = dat$beta0, target_type = "structural", target_mc_se = 0)
  testthat::expect_warning(theory <- method_theory_diagnostics_v2(
    config, 1L, 20281029L, "PROFILE-DQR-POP-H", fit, target_obj, population, dat
  ), NA)
  testthat::expect_true(is.na(theory$max_bahadur_scaled))
  testthat::expect_true(is.na(theory$scaled_bahadur_remainder))
})

pair_benchmark_root <- v2_require_root()
source(file.path(pair_benchmark_root, "scripts", "00_source_v2.R"), local = TRUE)
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(pair_benchmark_root, module, environment())
}

pair_benchmark_config <- function(n = 40L, p = 6L, q = 1L, tau = 0.5) {
  list(
    experiment_id = "PAIR-BENCHMARK-ACCEPTANCE", n_clusters = n, p = p,
    s = min(2L, p), tau = tau, q = q, m_values = 3:5,
    m_rule = "uniform", m_multiset = NULL, design_type = "ar1", rho_x = 0.3,
    signal = 0.75, error_dist = "gaussian", heteroskedastic = FALSE,
    error_dependence = "independent", copula_rho = 0, random_effect_dist = "normal",
    sigma_b0 = 0.8, sigma_b1 = 0.5, rho_b = 0.4, x_b_corr = 0,
    informative_size = FALSE, nonlinear_re_strength = 0,
    fit_random_effects = if (q == 1L) "intercept" else "intercept+slope",
    lambda_gamma = 1, h_multiplier = 1,
    lambda_beta_multipliers = c(0.25, 0.5, 1, 2),
    dantzig_multipliers = c(0.02, 0.05, 0.10, 0.25, 0.50, 1, 2, 4),
    target_mode = "structural"
  )
}

pair_benchmark_tuning <- function(config) {
  tuning <- reference_tuning_values_v2(config)
  tuning$method_version <- "profile_pair_20261022"
  tuning$h_inf <- tuning$h
  tuning$mu_grid <- tuning$mu
  tuning$mean_score_multiplier <- 1
  tuning$lambda_multiplier <- 1
  tuning
}

pair_benchmark_schema <- function(answer, dat, targets) {
  required <- c("method_id", "status", "failure_stage", "failure_message",
                "beta_hat", "beta_tilde", "se", "ci_lower", "ci_upper",
                "selected", "prediction_function", "runtime_sec", "converged",
                "warning_messages", "implementation_version", "reference_identifier",
                "adapter_fidelity_status", "method_version", "score_definition")
  testthat::expect_true(all(required %in% names(answer)))
  testthat::expect_identical(answer$status, "ok")
  testthat::expect_true(answer$converged)
  testthat::expect_identical(names(answer$beta_hat), colnames(dat$X))
  testthat::expect_true(all(targets %in% names(answer$beta_tilde)))
  testthat::expect_true(all(is.finite(answer$beta_hat)))
  testthat::expect_true(all(is.finite(answer$beta_tilde)))
  testthat::expect_true(all(is.finite(answer$se) & answer$se > 0))
  testthat::expect_true(all(is.finite(answer$ci_lower)))
  testthat::expect_true(all(is.finite(answer$ci_upper)))
  testthat::expect_true(is.finite(answer$runtime_sec) && answer$runtime_sec >= 0)
  testthat::expect_identical(answer$method_version, "profile_pair_20261022")
  testthat::expect_identical(answer$score_definition,
                            "profile_bandwidth_pair_2_half_minus_full")
  prediction <- answer$prediction_function(dat$X)
  testthat::expect_equal(prediction, as.numeric(dat$X %*% answer$beta_hat),
                         tolerance = 1e-12)
  reversed <- dat$X[, rev(seq_len(ncol(dat$X))), drop = FALSE]
  testthat::expect_equal(answer$prediction_function(reversed), prediction,
                         tolerance = 1e-12)
}

pair_benchmark_score_fidelity <- function(answer, dat) {
  fit <- answer$fit_object
  component <- fit$components
  direct <- profile_bandwidth_pair_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, fit$beta_refit,
    fit$tau, fit$h_inf, lambda_gamma = fit$Lambda
  )
  testthat::expect_equal(component$score, direct$score, tolerance = 1e-10)
  testthat::expect_equal(component$hessian, direct$hessian, tolerance = 1e-10)
  for (j in seq_len(nrow(answer$inference_object$table))) {
    row <- answer$inference_object$table[j, ]
    omega <- answer$inference_object$directions[[j]]$omega
    projection <- as.numeric(component$cluster_scores %*% omega)
    testthat::expect_equal(row$beta_tilde,
                           unname(fit$beta_refit[row$coordinate] - sum(omega * direct$score)),
                           tolerance = 1e-10)
    testthat::expect_equal(row$se^2 * component$n_clusters,
                           mean((projection - mean(projection))^2), tolerance = 1e-10)
  }
}

testthat::test_that("PROFILE-DQR paired deterministic schema and fidelity", {
  config <- pair_benchmark_config()
  dat <- simulate_from_config_v2(config, 20261101L)
  tuning <- pair_benchmark_tuning(config)
  targets <- c("x00001", "x00004")
  first <- fit_benchmark_by_id_v2("PROFILE-DQR", dat, config$tau, targets,
                                  tuning, 20261102L)
  second <- fit_benchmark_by_id_v2("PROFILE-DQR", dat, config$tau, targets,
                                   tuning, 20261102L)
  pair_benchmark_schema(first, dat, targets)
  testthat::expect_identical(first$beta_hat, second$beta_hat)
  testthat::expect_identical(first$beta_tilde, second$beta_tilde)
  testthat::expect_equal(first$lambda_mean, first$fit_object$h_est)
  pair_benchmark_score_fidelity(first, dat)
})

testthat::test_that("PROFILE-DQR paired high dimensional slope limiting case", {
  config <- pair_benchmark_config(n = 20L, p = 26L, q = 2L, tau = 0.25)
  config$lambda_gamma <- matrix(c(1, 0.2, 0.2, 0.8), 2L)
  dat <- simulate_from_config_v2(config, 20261024L)
  tuning <- pair_benchmark_tuning(config)
  answer <- fit_benchmark_by_id_v2("PROFILE-DQR", dat, config$tau,
                                    c("x00001", "x00004"), tuning, 20261103L)
  pair_benchmark_schema(answer, dat, c("x00001", "x00004"))
  testthat::expect_gt(ncol(dat$X), dat$n_clusters)
  testthat::expect_equal(answer$fit_object$Lambda, config$lambda_gamma)
  testthat::expect_true(all(answer$inference_object$table$mu < 1))
  testthat::expect_equal(answer$inference_object$mu_grid_requested, tuning$mu_grid)
})

testthat::test_that("PROFILE-DQR-TRUE-SUPPORT paired deterministic schema limiting case and fidelity", {
  config <- pair_benchmark_config()
  dat <- simulate_from_config_v2(config, 20261104L)
  tuning <- pair_benchmark_tuning(config)
  targets <- c("x00001", "x00004")
  first <- fit_benchmark_by_id_v2("PROFILE-DQR-TRUE-SUPPORT", dat, config$tau,
                                  targets, tuning, 20261105L,
                                  context = list(active = dat$active))
  second <- fit_benchmark_by_id_v2("PROFILE-DQR-TRUE-SUPPORT", dat, config$tau,
                                   targets, tuning, 20261105L,
                                   context = list(active = dat$active))
  pair_benchmark_schema(first, dat, targets)
  testthat::expect_identical(first$beta_tilde, second$beta_tilde)
  keep <- sort(unique(c(dat$active, match(targets, colnames(dat$X)))))
  reduced <- subset_features_v2(dat, keep)
  direct <- fit_benchmark_by_id_v2("PROFILE-DQR", reduced, config$tau, targets,
                                   tuning, 20261105L,
                                   control = list(penalty_factor = rep(0, length(keep))))
  testthat::expect_equal(first$beta_hat[keep], direct$beta_hat, tolerance = 1e-10)
  testthat::expect_equal(first$beta_tilde, direct$beta_tilde, tolerance = 1e-10)
  testthat::expect_equal(first$fit_object$lambda_base, 0)
  pair_benchmark_score_fidelity(first, reduced)
})

testthat::test_that("PROFILE-DQR-TRUE-NUISANCE paired deterministic schema limiting case and fidelity", {
  config <- pair_benchmark_config(q = 2L)
  dat <- simulate_from_config_v2(config, 20261106L)
  tuning <- pair_benchmark_tuning(config)
  targets <- c("x00001", "x00004")
  first <- fit_benchmark_by_id_v2("PROFILE-DQR-TRUE-NUISANCE", dat, config$tau,
                                  targets, tuning, 20261107L)
  second <- fit_benchmark_by_id_v2("PROFILE-DQR-TRUE-NUISANCE", dat, config$tau,
                                   targets, tuning, 20261107L)
  pair_benchmark_schema(first, dat, targets)
  testthat::expect_identical(first$beta_tilde, second$beta_tilde)
  nuisance <- rowSums(dat$Z * dat$random_effects[
    match(dat$cluster_id, rownames(dat$random_effects)), , drop = FALSE])
  reduced <- dat
  reduced$y <- dat$y - nuisance
  reduced$Z <- matrix(0, nrow(dat$X), 1L)
  direct <- fit_benchmark_by_id_v2("PROFILE-DQR", reduced, config$tau, targets,
                                   tuning, 20261107L)
  testthat::expect_equal(first$beta_tilde, direct$beta_tilde, tolerance = 1e-10)
  testthat::expect_equal(first$oracle_nuisance_removed, nuisance, tolerance = 1e-12)
  testthat::expect_lt(max(abs(first$fit_object$gamma)), 1e-12)
  pair_benchmark_score_fidelity(first, reduced)
  omitted <- dat
  omitted$Z <- dat$Z[, 1, drop = FALSE]
  testthat::expect_error(fit_benchmark_by_id_v2(
    "PROFILE-DQR-TRUE-NUISANCE", omitted, config$tau, targets, tuning, 20261107L),
    "correctly specified")
})

testthat::test_that("POSTREFIT-EXACT-H paired deterministic schema limiting case and fidelity", {
  config <- pair_benchmark_config()
  dat <- simulate_from_config_v2(config, 20261108L)
  tuning <- pair_benchmark_tuning(config)
  targets <- c("x00001", "x00004")
  first <- fit_benchmark_by_id_v2("POSTREFIT-EXACT-H", dat, config$tau,
                                  targets, tuning, 20261109L)
  second <- fit_benchmark_by_id_v2("POSTREFIT-EXACT-H", dat, config$tau,
                                   targets, tuning, 20261109L)
  pair_benchmark_schema(first, dat, targets)
  testthat::expect_identical(first$beta_tilde, second$beta_tilde)
  fit <- first$fit_object
  selected <- fit$refit_set
  J <- fit$components$hessian[selected, selected, drop = FALSE]
  for (j in seq_along(targets)) {
    direction <- first$inference_object$directions[[j]]$omega
    canonical <- numeric(length(selected))
    canonical[match(targets[j], names(fit$beta)[selected])] <- 1
    testthat::expect_equal(as.numeric(J %*% direction[selected]), canonical,
                           tolerance = 1e-10)
    testthat::expect_equal(direction[-selected], rep(0, length(direction) - length(selected)))
  }
  pair_benchmark_score_fidelity(first, dat)
})

testthat::test_that("PROFILE-DQR-SPLIT paired deterministic schema limiting case and fidelity", {
  config <- pair_benchmark_config(n = 80L)
  dat <- simulate_from_config_v2(config, 20261110L)
  tuning <- pair_benchmark_tuning(config)
  targets <- c("x00001", "x00004")
  control <- list(screen_fraction = 0.5, screen_dim = 4L)
  first <- fit_benchmark_by_id_v2("PROFILE-DQR-SPLIT", dat, config$tau,
                                  targets, tuning, 20261111L, control = control)
  second <- fit_benchmark_by_id_v2("PROFILE-DQR-SPLIT", dat, config$tau,
                                   targets, tuning, 20261111L, control = control)
  pair_benchmark_schema(first, dat, targets)
  testthat::expect_identical(first$beta_tilde, second$beta_tilde)
  role <- first$screening
  testthat::expect_length(intersect(role$screen_clusters, role$fit_clusters), 0L)
  reduced <- subset_features_v2(subset_clusters_v2(dat, role$fit_clusters),
                                 role$selected_global)
  direct <- fit_benchmark_by_id_v2("PROFILE-DQR", reduced, config$tau, targets,
                                   tuning, 20261111L)
  testthat::expect_equal(first$beta_tilde, direct$beta_tilde, tolerance = 1e-10)
  testthat::expect_equal(first$fit_object$components$n_clusters,
                         length(role$fit_clusters))
  testthat::expect_error(fit_benchmark_by_id_v2(
    "PROFILE-DQR-SPLIT", dat, config$tau, targets, tuning, 20261111L,
    context = list(screening = role)), "inference clusters only")
  pair_benchmark_score_fidelity(first, reduced)
})

testthat::test_that("PROFILE-DQR-POP-H paired deterministic schema limiting case and fidelity", {
  config <- pair_benchmark_config(n = 60L, p = 2L)
  config$s <- 1L
  tuning <- pair_benchmark_tuning(config)
  dat <- simulate_from_config_v2(config, 20261112L)
  targets <- colnames(dat$X)
  target <- build_profile_pair_target_asset_v2(
    config, n_population = 60L, repeats = 2L, seed = 20261024L, development = TRUE)
  asset <- build_profile_pair_direction_asset_v2(
    config, target, n_population = 60L, repeats = 2L, seed = 20271024L,
    development = TRUE)
  control <- list(final_assets = FALSE, implementation_commit = "working_tree")
  context <- list(config = config, population_direction_asset = asset)
  first <- fit_benchmark_by_id_v2("PROFILE-DQR-POP-H", dat, config$tau,
                                  targets, tuning, 20261113L, context, control)
  second <- fit_benchmark_by_id_v2("PROFILE-DQR-POP-H", dat, config$tau,
                                   targets, tuning, 20261113L, context, control)
  pair_benchmark_schema(first, dat, targets)
  testthat::expect_identical(first$beta_tilde, second$beta_tilde)
  testthat::expect_length(intersect(target$seeds, asset$seeds), 0L)
  testthat::expect_equal(unname(asset$H_population %*% asset$Omega_population),
                         diag(2), tolerance = 1e-8)
  for (j in seq_along(targets)) {
    testthat::expect_equal(first$inference_object$directions[[j]]$omega,
                           as.numeric(asset$directions[[targets[j]]]$omega),
                           tolerance = 1e-12)
  }
  pair_benchmark_score_fidelity(first, dat)
  testthat::expect_error(fit_benchmark_by_id_v2(
    "PROFILE-DQR-POP-H", dat, config$tau, targets, tuning, 20261113L,
    list(population_direction_asset = asset), control), "context\\$config")
  changed <- config
  changed$lambda_gamma <- 2
  testthat::expect_error(fit_benchmark_by_id_v2(
    "PROFILE-DQR-POP-H", dat, config$tau, targets, tuning, 20261113L,
    modifyList(context, list(config = changed)), control), "identity")
  testthat::expect_error(fit_benchmark_by_id_v2(
    "PROFILE-DQR-POP-H", dat, config$tau, targets, tuning, 20261113L,
    context, modifyList(control, list(implementation_commit = "another_commit"))),
    "commit mismatch")
  testthat::expect_error(fit_benchmark_by_id_v2(
    "PROFILE-DQR-POP-H", dat, config$tau, targets, tuning, 20261113L,
    context, modifyList(control, list(final_assets = TRUE))), "commit identity")
})

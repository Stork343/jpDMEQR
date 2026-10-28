conditional_test_root <- v2_require_root()
source(file.path(conditional_test_root, "scripts", "00_source_v2.R"), local = TRUE)
for (module in c("profile_v2", "simulation_v2")) source_v2_module(conditional_test_root, module, environment())

conditional_test_config <- function(q = 1L) {
  table <- utils::read.csv(file.path(conditional_test_root, "config", "simulation_main.csv"), stringsAsFactors = FALSE)
  config <- config_row_to_list_v2(table[table$experiment_id == "B16", , drop = FALSE])
  config$target_mode <- "profile_mc"
  config$p <- 12L
  config$q <- q
  config$fit_random_effects <- if (q == 1L) "intercept" else "intercept+slope"
  config
}

testthat::test_that("conditional profile target preserves the complete signal integrand", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "compiled"))
  config <- conditional_test_config()
  dat <- simulate_from_config_v2(utils::modifyList(config, list(n_clusters = 80L)), seed = 202610231L)
  delta <- .09
  beta <- dat$beta0
  beta[1L] <- beta[1L] + delta
  control <- list(grad_tol = 1e-11)
  full <- profile_bandwidth_pair_v2(dat$y, dat$X, dat$Z, dat$cluster_id, beta, config$tau, .08, 1,
    nuisance_control = control)
  scalar <- profile_bandwidth_pair_v2(dat$y - drop(dat$X %*% dat$beta0), dat$X[, 1L, drop = FALSE],
    dat$Z, dat$cluster_id, delta, config$tau, .08, 1, nuisance_control = control)
  testthat::expect_equal(full$cluster_scores[, 1L], scalar$cluster_scores[, 1L], tolerance = 1e-9)
  testthat::expect_equal(full$hessian[1L, 1L], scalar$hessian[1L, 1L], tolerance = 1e-9)
  testthat::expect_equal(full$score[1L], scalar$score[1L], tolerance = 1e-9)
  testthat::expect_identical(dat$active, seq_len(8L))
  testthat::expect_identical(ar1_covariance_v2(1000L, config$rho_x)[1:11, 1:11], ar1_covariance_v2(11L, config$rho_x))
})

testthat::test_that("conditional profile target builds genuine deterministic full-coordinate roots", {
  config <- conditional_test_config()
  asset <- build_profile_pair_conditional_scale_target_v2(config, 400L, 2L, seed = 202610232L, development = TRUE)
  replay <- build_profile_pair_conditional_scale_target_v2(config, 400L, 2L, seed = 202610232L, development = TRUE)
  testthat::expect_true(validate_profile_pair_conditional_scale_target_v2(asset, config, final = FALSE)$valid)
  testthat::expect_identical(asset$estimates_by_repeat, replay$estimates_by_repeat)
  testthat::expect_identical(asset$asset_kind, "paired_conditional_scale_target")
  testthat::expect_identical(asset$population_engine, "conditional_gaussian_scale_v1")
  testthat::expect_identical(asset$finite_full_design_empirical_root, FALSE)
  testthat::expect_length(asset$beta_star_mc, 12L)
  testthat::expect_identical(asset$target_names, sprintf("x%05d", 1:12))
  testthat::expect_equal(asset$beta_star_mc[-1L], make_beta_v2(12L, 8L, .75)[-1L], tolerance = 0)
  testthat::expect_true(all(asset$target_mc_se_by_coordinate[-1L] == 0))
  testthat::expect_equal(asset$target_mc_se_by_coordinate, apply(asset$estimates_by_repeat, 2L, stats::sd) / sqrt(2), tolerance = 1e-14)
  testthat::expect_equal(asset$target_bandwidth_displacement_by_repeat,
    asset$refinement_estimates_by_repeat - asset$estimates_by_repeat, tolerance = 0)
  testthat::expect_equal(asset$bandwidths_refinement, c(asset$h_target / 4, asset$h_target / 2), tolerance = 0)
  testthat::expect_lt(asset$max_kkt_residual, 1e-7)
  testthat::expect_lt(asset$profile_identity_error, 1e-8)
  testthat::expect_lt(asset$max_nuisance_gradient, 1e-7)
  testthat::expect_identical(asset$repeat_assets[[1L]]$generated_complete_support, seq_len(8L))
  testthat::expect_equal(asset$repeat_assets[[1L]]$n_population, 400L)
  testthat::expect_error(validate_profile_pair_conditional_scale_target_v2(asset, config, final = TRUE))
  changed <- asset
  changed$target_mc_se_by_coordinate[1L] <- changed$target_mc_se_by_coordinate[1L] + .1
  testthat::expect_error(validate_profile_pair_conditional_scale_target_v2(changed, config, final = FALSE), "MCSE")
  changed <- asset
  changed$repeat_assets[[1L]]$target$conditional_geometry$full_score[2L] <- .1
  testthat::expect_error(validate_profile_pair_conditional_scale_target_v2(changed, config, final = FALSE), "diagnostics")
  changed <- asset
  changed$repeat_assets[[1L]]$target$conditional_geometry$hessian_coefficients[] <- 1e-15
  changed$repeat_assets[[1L]]$target$conditional_geometry$hessian_identity_coefficients[] <- 1e-15
  testthat::expect_error(validate_profile_pair_conditional_scale_target_v2(changed, config, final = FALSE), "diagnostics")
  changed <- asset
  changed$repeat_assets[[1L]]$analysis_score_at_target$max_nuisance_gradient <- 1
  testthat::expect_error(validate_profile_pair_conditional_scale_target_v2(changed, config, final = FALSE), "analysis-bandwidth")
  changed <- asset
  changed$source_paths <- changed$source_paths[-1L]
  testthat::expect_error(validate_profile_pair_conditional_scale_target_v2(changed, config, final = FALSE), "source identity")
})

testthat::test_that("conditional profile geometry covers random time slopes and non-diagonal Lambda", {
  config <- conditional_test_config(2L)
  config$lambda_gamma <- matrix(c(1, .2, .2, 1.3), 2L)
  asset <- build_profile_pair_conditional_scale_target_v2(config, 400L, 2L, seed = 202610233L, development = TRUE)
  testthat::expect_true(validate_profile_pair_conditional_scale_target_v2(asset, config, final = FALSE)$valid)
  testthat::expect_lt(asset$profile_identity_error, 1e-8)
  testthat::expect_true(all(vapply(asset$repeat_assets, function(record) {
    record$target$conditional_geometry$minimum_eigenvalue_lower_bound > 1e-8
  }, logical(1))))
  testthat::expect_identical(asset$config$q, 2L)
  invalid <- config
  invalid$x_b_corr <- .2
  testthat::expect_error(conditional_scale_target_conditions_v2(invalid), "conditions")
  invalid <- config
  invalid$heteroskedastic <- FALSE
  testthat::expect_error(conditional_scale_target_conditions_v2(invalid), "conditions")
})

testthat::test_that("own SQR marginal target agrees with an independent Normal-mixture equation", {
  config <- conditional_test_config()
  asset <- build_sqr_marginal_scale_target_v2(config, 400L, 2L, seed = 202610234L, development = TRUE)
  equation <- function(delta) {
    integrate_part <- function(interval) stats::integrate(function(x) {
      scale <- exp(.35 * pmax(pmin(x, 3), -3))
      x * (stats::pnorm((delta * x + scale * stats::qnorm(config$tau)) /
        sqrt(config$sigma_b0^2 + scale^2)) - config$tau) * stats::dnorm(x)
    }, interval[1L], interval[2L], rel.tol = 1e-10, abs.tol = 1e-12)$value
    sum(vapply(list(c(-Inf, -3), c(-3, 3), c(3, Inf)), integrate_part, numeric(1)))
  }
  independent <- stats::uniroot(equation, c(-2, 2), tol = 1e-10)$root
  testthat::expect_equal(asset$scalar_delta, independent, tolerance = 1e-7)
  testthat::expect_true(validate_sqr_marginal_scale_target_v2(asset, config, final = FALSE)$valid)
  testthat::expect_identical(asset$asset_kind, "sqr_marginal_qr_target")
  testthat::expect_identical(asset$target_scope, "own_pooled_visit_marginal_QR_population_projection")
  testthat::expect_identical(asset$beta_target, asset$beta_marginal)
  testthat::expect_true(all(asset$target_mc_se_by_coordinate == 0))
  testthat::expect_lt(asset$quadrature$score_residual, 1e-7)
  testthat::expect_lt(asset$quadrature_root_difference, 1e-7)
  testthat::expect_true(all(vapply(asset$repeat_assets, function(record) record$scalar_score_mcse > 0, logical(1))))
  testthat::expect_identical(asset$repeat_assets[[1L]]$generated_complete_support, seq_len(8L))
  testthat::expect_error(validate_sqr_marginal_scale_target_v2(asset, config, final = TRUE))
  changed <- asset
  changed$beta_target[1L] <- changed$beta_target[1L] + .1
  testthat::expect_error(validate_sqr_marginal_scale_target_v2(changed, config, final = FALSE), "definition")
  changed <- asset
  changed$repeat_assets[[1L]]$score_mcse[1L] <- 0
  testthat::expect_error(validate_sqr_marginal_scale_target_v2(changed, config, final = FALSE), "audit")
  changed <- asset
  changed$population_engine <- "conditional_gaussian_scale_v1"
  testthat::expect_error(validate_sqr_marginal_scale_target_v2(changed, config, final = FALSE), "definition")
})

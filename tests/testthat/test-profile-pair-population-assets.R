pair_asset_root <- v2_require_root()
source(file.path(pair_asset_root, "scripts", "00_source_v2.R"), local = TRUE)
source_v2_module(pair_asset_root, "profile_v2", environment())
source_v2_module(pair_asset_root, "simulation_v2", environment())

testthat::test_that("paired population assets build real roots and reject incompatible records", {
  config <- list(
    experiment_id = "PAIR-ASSET-TEST", n_clusters = 60L, p = 2L, s = 1L,
    tau = 0.5, q = 1L, m_values = 3:4, m_rule = "uniform", m_multiset = NULL,
    design_type = "ar1", rho_x = 0.3, signal = 0.75, error_dist = "gaussian",
    heteroskedastic = FALSE, error_dependence = "independent", copula_rho = 0,
    random_effect_dist = "normal", sigma_b0 = 0.8, sigma_b1 = 0.5, rho_b = 0.4,
    x_b_corr = 0, informative_size = FALSE, nonlinear_re_strength = 0,
    fit_random_effects = "intercept", lambda_gamma = 1, h_multiplier = 1,
    target_mode = "profile_mc"
  )
  target <- build_profile_pair_target_asset_v2(
    config, n_population = 60L, repeats = 2L, seed = 20261024L, development = TRUE
  )
  direction <- build_profile_pair_direction_asset_v2(
    config, target, n_population = 60L, repeats = 2L, seed = 20271024L, development = TRUE
  )
  testthat::expect_true(validate_profile_pair_target_asset_v2(target, config, final = FALSE)$valid)
  testthat::expect_true(validate_profile_pair_direction_asset_v2(direction, config, final = FALSE)$valid)
  testthat::expect_identical(target$score_definition, "profile_bandwidth_pair_2_half_minus_full")
  testthat::expect_identical(target$target_names, c("x00001", "x00002"))
  testthat::expect_lt(target$max_kkt_residual, 1e-7)
  testthat::expect_lt(target$profile_identity_error, 1e-8)
  testthat::expect_lt(direction$covariance_identity_error, 1e-8)
  testthat::expect_lt(max(abs(direction$H_population %*% direction$Omega_population - diag(2))), 1e-7)
  testthat::expect_identical(direction$sigma0_pop, direction$sigma_population)
  testthat::expect_identical(direction$target_source, "paired_population_target_mc")
  testthat::expect_identical(direction$evaluation_point, target$beta_star_mc)
  testthat::expect_identical(direction$beta_target, target$beta_star_mc)
  testthat::expect_length(intersect(target$seeds, direction$seeds), 0L)
  testthat::expect_equal(direction$Sigma_population,
                         4 * direction$Sigma_half + direction$Sigma_full -
                           2 * direction$cross_covariance_half_full -
                           2 * t(direction$cross_covariance_half_full), tolerance = 1e-12)
  destination <- tempfile("paired_asset_", tmpdir = file.path(pair_asset_root, "results", "scratch"))
  saveRDS(direction, destination)
  restored <- readRDS(destination)
  testthat::expect_true(validate_profile_pair_direction_asset_v2(restored, config, final = FALSE)$valid)
  commit <- system2("git", c("-C", shQuote(pair_asset_root), "rev-parse", "HEAD"), stdout = TRUE)
  testthat::expect_error(validate_profile_pair_target_asset_v2(target, config, final = TRUE), "commit identity")
  incompatible <- target
  incompatible$implementation_commit <- commit
  testthat::expect_error(validate_profile_pair_target_asset_v2(incompatible, config, final = TRUE), "Development")
  testthat::expect_error(build_profile_pair_target_asset_v2(
    config, 100000L, 4L, development = FALSE), "commit identity")
  testthat::expect_error(build_profile_pair_target_asset_v2(
    config, 100000L, 4L, development = FALSE, implementation_commit = ""), "commit identity")
  testthat::expect_error(build_profile_pair_direction_asset_v2(
    config, target, 100000L, 4L, development = FALSE), "commit identity")
  incompatible <- target
  incompatible$score_definition <- "ordinary_unsmoothed_quantile"
  testthat::expect_error(validate_profile_pair_target_asset_v2(incompatible, final = FALSE), "score_definition")
  incompatible <- direction
  incompatible$construction_version <- "pop-h-v2-analysis-bandwidth"
  testthat::expect_error(validate_profile_pair_direction_asset_v2(incompatible, final = FALSE), "construction_version")
  incompatible$construction_version <- "profile-pair-population-assets-20261022-v1"
  testthat::expect_error(validate_profile_pair_direction_asset_v2(incompatible, final = FALSE), "construction_version")
  incompatible <- direction
  incompatible$Sigma_population[1, 1] <- incompatible$Sigma_population[1, 1] + 0.1
  testthat::expect_error(validate_profile_pair_direction_asset_v2(incompatible, final = FALSE), "covariance")
  incompatible <- config
  incompatible$lambda_gamma <- 2
  testthat::expect_error(validate_profile_pair_direction_asset_v2(direction, incompatible, final = FALSE), "identity")
  incompatible <- target
  incompatible$config$lambda_gamma <- 2
  testthat::expect_error(validate_profile_pair_target_asset_v2(incompatible, final = FALSE), "working-design identity")
  incompatible <- direction
  incompatible$implementation_commit <- ""
  testthat::expect_error(validate_profile_pair_direction_asset_v2(incompatible, final = TRUE), "commit identity")
  testthat::expect_error(build_profile_pair_target_asset_v2(
    config, 60L, 2L, development = FALSE, implementation_commit = commit), "100000")
  testthat::expect_error(profile_pair_curvature_v2(diag(c(1, -1))), "curvature")
})

testthat::test_that("structural directions evaluate their own generated structural coefficient", {
  config <- list(
    experiment_id = "PAIR-ASSET-STRUCTURAL-TEST", n_clusters = 60L, p = 1L, s = 1L,
    tau = 0.5, q = 1L, m_values = 3:4, m_rule = "uniform", m_multiset = NULL,
    design_type = "ar1", rho_x = 0.3, signal = 0.75, error_dist = "gaussian",
    heteroskedastic = FALSE, error_dependence = "independent", copula_rho = 0,
    random_effect_dist = "normal", sigma_b0 = 0.8, sigma_b1 = 0.5, rho_b = 0.4,
    x_b_corr = 0, informative_size = FALSE, nonlinear_re_strength = 0,
    fit_random_effects = "intercept", lambda_gamma = 1, h_multiplier = 1,
    target_mode = "structural"
  )
  target <- build_profile_pair_target_asset_v2(
    config, n_population = 60L, repeats = 2L, seed = 20261025L, development = TRUE
  )
  direction <- build_profile_pair_direction_asset_v2(
    config, target, n_population = 60L, repeats = 2L, seed = 20271025L, development = TRUE
  )
  beta0 <- make_beta_v2(config$p, config$s, config$signal)
  testthat::expect_true(validate_profile_pair_direction_asset_v2(direction, config, final = FALSE)$valid)
  testthat::expect_identical(direction$target_source, "generated_structural_coefficient")
  testthat::expect_identical(direction$beta_target, beta0)
  testthat::expect_identical(direction$evaluation_point, beta0)
  testthat::expect_identical(target$target_source, "paired_population_target_mc")
  testthat::expect_identical(target$evaluation_point, target$beta_star_mc)
  testthat::expect_identical(direction$target_asset_identity$beta_star_mc, target$beta_star_mc)
  testthat::expect_true(all(vapply(direction$repeat_assets, function(x) {
    identical(x$generated_beta0, beta0) && identical(x$evaluation_point, beta0)
  }, logical(1))))
  incompatible <- direction
  incompatible$evaluation_point <- target$beta_star_mc
  testthat::expect_error(validate_profile_pair_direction_asset_v2(incompatible, final = FALSE), "evaluation point")
  incompatible <- direction
  incompatible$repeat_assets[[2]]$generated_beta0 <- target$beta_star_mc
  testthat::expect_error(validate_profile_pair_direction_asset_v2(incompatible, final = FALSE), "identity")
  testthat::expect_match(profile_pair_direction_asset_path_v2(pair_asset_root, config$experiment_id),
                         "construction_v3/directions")
})

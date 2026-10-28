root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, .GlobalEnv)
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"), local = TRUE)

testthat::test_that("paired fitting reuses only an identical numerical request", {
  dat <- generate_profile_qr_data_v2(
    n_clusters = 48L, p = 6L, s = 2L, tau = 0.5, q = 2L,
    m_values = 3:5, seed = 2026102201L
  )
  cache <- new.env(parent = emptyenv())
  tuning <- list(h_inf = 48^(-3 / 10), lambda_gamma = 1,
                 mean_score_multiplier = 1, lambda_multiplier = 1)
  control <- list(shared_first_stage = cache)
  targets <- colnames(dat$X)[c(1L, 4L)]
  first <- profile_pair_fit_from_data_v2(dat, 0.5, targets, tuning, control)
  reused <- profile_pair_fit_from_data_v2(dat, 0.5, targets, tuning, control)
  independent <- profile_pair_fit_from_data_v2(dat, 0.5, targets, tuning, list())
  testthat::expect_false(first$diagnostics$cache$reused)
  testthat::expect_true(reused$diagnostics$cache$reused)
  for (field in c("beta", "beta_l1", "loadings", "selected_support", "refit_set")) {
    testthat::expect_identical(reused[[field]], first[[field]])
    testthat::expect_identical(independent[[field]], first[[field]])
  }
  testthat::expect_identical(reused$components$score, independent$components$score)
  testthat::expect_identical(reused$components$hessian, independent$components$hessian)
  testthat::expect_equal(reused$first_stage_sec, 0)
  altered <- dat
  altered$y[1L] <- altered$y[1L] + 0.05
  changed <- profile_pair_fit_from_data_v2(altered, 0.5, targets, tuning, control)
  testthat::expect_false(changed$diagnostics$cache$reused)
  testthat::expect_identical(cache$paired_request$y, altered$y)
  changed_targets <- profile_pair_fit_from_data_v2(dat, 0.5, targets[1L], tuning, control)
  testthat::expect_false(changed_targets$diagnostics$cache$reused)
  changed_tuning <- utils::modifyList(tuning, list(lambda_multiplier = 1.1))
  changed_penalty <- profile_pair_fit_from_data_v2(dat, 0.5, targets, changed_tuning, control)
  testthat::expect_false(changed_penalty$diagnostics$cache$reused)
  testthat::expect_error(
    profile_pair_fit_from_data_v2(dat, 0.5, targets, tuning, list(shared_first_stage = list())),
    "environment"
  )
})

testthat::test_that("common low dimensional screening uses independent subject roles", {
  dat <- generate_profile_qr_data_v2(
    n_clusters = 40L, p = 30L, s = 3L, tau = 0.25, q = 1L,
    m_values = 3:6, seed = 2026102202L
  )
  targets <- colnames(dat$X)[c(1L, 12L)]
  first <- common_lqmm_subset_v2(dat, targets, 0.25, max_features = 8L, seed = 581L)
  repeated <- common_lqmm_subset_v2(dat, targets, 0.25, max_features = 8L, seed = 581L)
  testthat::expect_identical(first, repeated)
  testthat::expect_length(intersect(first$roles$screen_clusters, first$roles$fit_clusters), 0L)
  testthat::expect_setequal(c(first$roles$screen_clusters, first$roles$fit_clusters),
                            as.character(unique(dat$cluster_id)))
  testthat::expect_setequal(as.character(unique(first$data$cluster_id)), first$roles$fit_clusters)
  testthat::expect_equal(first$data$n_clusters, length(first$roles$fit_clusters))
  testthat::expect_true(all(targets %in% colnames(first$data$X)))
  testthat::expect_lte(ncol(first$data$X), 8L)
  altered <- dat
  fit_rows <- as.character(dat$cluster_id) %in% first$roles$fit_clusters
  altered$y[fit_rows] <- altered$y[fit_rows] + seq_len(sum(fit_rows)) / 10
  changed <- common_lqmm_subset_v2(altered, targets, 0.25, max_features = 8L, seed = 581L)
  testthat::expect_identical(changed$keep, first$keep)
  testthat::expect_identical(changed$roles, first$roles)
  oracle <- common_lqmm_subset_v2(dat, targets, 0.25, max_features = 8L, active = dat$active)
  testthat::expect_null(oracle$roles)
  testthat::expect_equal(oracle$data$n_clusters, 40L)
})

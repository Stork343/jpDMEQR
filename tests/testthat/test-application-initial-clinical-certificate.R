root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "application_v2", envir = .GlobalEnv)

testthat::test_that("真实训练数据的完整KKT计算与相同请求缓存保持一致", {
  path <- file.path(root, "results", "scratch", "initial_scaling_none_20261022_v1", "source_design.rds")
  testthat::skip_if_not(file.exists(path), "本地训练设计未提供。")
  observed <- readRDS(path)
  train <- observed$training
  design <- observed$design
  md <- observed$metadata
  cache <- new.env(parent = emptyenv())
  run <- function(h = NULL) fit_application_profile_initial_v2(
    md$log_c3[train], design$X[train, , drop = FALSE], design$Z[train, , drop = FALSE],
    md$subject_id[train], .5, h_est = h, penalty_factor = design$penalty_factor,
    control = list(clinical_acceleration = "clinical_kkt", clinical_cache = cache)
  )
  first <- run()
  repeated <- run()
  testthat::expect_identical(first$beta_l1, repeated$beta_l1)
  testthat::expect_identical(first$loadings, repeated$loadings)
  testthat::expect_length(cache$entries, 1L)
  testthat::expect_false(first$initial_optimizer_history$preliminary$cache_hit)
  testthat::expect_true(first$initial_optimizer_history$final$cache_hit)
  testthat::expect_true(repeated$initial_optimizer_history$preliminary$cache_hit)
  testthat::expect_equal(ncol(first$components$cluster_scores), ncol(design$X))
  testthat::expect_equal(dim(first$components$hessian), rep(ncol(design$X), 2L))
  testthat::expect_null(first$components$hessian_fold_sums)
  testthat::expect_lte(first$components$profile_identity_error, 1e-8)
  for (pass in first$initial_optimizer_history) {
    testthat::expect_true(pass$global_certificate_accepted)
    testthat::expect_named(pass$global_residual, colnames(design$X))
    testthat::expect_lte(max(pass$global_residual), 1e-7)
    testthat::expect_lte(pass$global_kkt_normalized, 1e-7)
    testthat::expect_gt(nrow(pass$clinical_history), 0L)
  }
  changed <- run(first$h_est * .75)
  testthat::expect_length(cache$entries, 2L)
  testthat::expect_false(changed$initial_optimizer_history$preliminary$cache_hit)
  record <- application_initial_fit_record_v2(first)
  testthat::expect_identical(record$initial_optimizer_history, first$initial_optimizer_history)
  testthat::expect_identical(record$initial_acceleration, "clinical_kkt")
})

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)

fold_fixture <- function(seed = 20260930L) {
  generate_profile_qr_data_v2(24L, 8L, 3L, 0.5, q = 2L, m_values = 3:5,
                              seed = seed)
}

testthat::test_that("omitting the fold sums leaves every other object bit-identical", {
  d <- fold_fixture()
  for (mode in c("sequential", "batch")) {
    full <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id,
                                  as.numeric(d$beta), 0.5, 0.4,
                                  need_hessian = TRUE, assembly = mode)
    lean <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id,
                                  as.numeric(d$beta), 0.5, 0.4,
                                  need_hessian = TRUE, assembly = mode,
                                  need_fold_sums = FALSE)

    testthat::expect_null(lean$hessian_fold_sums)
    testthat::expect_null(lean$hessian_fold_counts)
    testthat::expect_identical(lean$hessian, full$hessian)
    testthat::expect_identical(lean$gamma, full$gamma)
    testthat::expect_identical(lean$cluster_scores, full$cluster_scores)
    testthat::expect_identical(lean$unsmoothed_cluster_scores,
                               full$unsmoothed_cluster_scores)
    testthat::expect_identical(lean$objective, full$objective)
    testthat::expect_identical(lean$score, full$score)
    testthat::expect_identical(lean$residuals, full$residuals)
    testthat::expect_identical(lean$nuisance_diagnostics,
                               full$nuisance_diagnostics)
  }
})

testthat::test_that("omitting the fold sums keeps the identity audit intact", {
  d <- fold_fixture()
  full <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id,
                                as.numeric(d$beta), 0.5, 0.4,
                                need_hessian = TRUE, identity_check = TRUE,
                                assembly = "batch")
  lean <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id,
                                as.numeric(d$beta), 0.5, 0.4,
                                need_hessian = TRUE, identity_check = TRUE,
                                assembly = "batch", need_fold_sums = FALSE)

  testthat::expect_null(lean$hessian_fold_sums)
  testthat::expect_identical(lean$hessian_identity, full$hessian_identity)
  testthat::expect_identical(lean$profile_identity_error,
                             full$profile_identity_error)
  testthat::expect_identical(lean$cluster_identity_errors,
                             full$cluster_identity_errors)
})

testthat::test_that("a fit returns fold sums whose weights reproduce the Hessian", {
  d <- fold_fixture()
  fit <- fit_profile_lasso_v2(d$y, d$X, d$Z, d$cluster_id, 0.5, 0.4,
                              lambda_beta = 0.05, control = list(max_iter = 200L))
  n <- fit$components$n_clusters
  testthat::expect_identical(length(fit$components$hessian_fold_sums), 4L)
  testthat::expect_identical(sum(fit$components$hessian_fold_counts), n)
  testthat::expect_true(all(fit$components$hessian_fold_counts > 0L))
  assembled <- Reduce(`+`, fit$components$hessian_fold_sums) / n
  testthat::expect_lt(max(abs(assembled - fit$components$hessian)), 1e-14)
})
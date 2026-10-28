root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)

assembly_fixture <- function(seed = 20260922L) {
  generate_profile_qr_data_v2(24L, 8L, 3L, 0.5, q = 2L, m_values = 3:5,
                              seed = seed)
}

rel_diff <- function(a, b) max(abs(a - b)) / max(abs(b))

testthat::test_that("the batched assembly reproduces the sequential Hessian", {
  d <- assembly_fixture()
  n <- length(unique(d$cluster_id))
  seq_out <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id,
                                   as.numeric(d$beta), 0.5, 0.4,
                                   need_hessian = TRUE, assembly = "sequential")
  bat_out <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id,
                                   as.numeric(d$beta), 0.5, 0.4,
                                   need_hessian = TRUE, assembly = "batch")

  testthat::expect_lt(rel_diff(bat_out$hessian, seq_out$hessian), 1e-12)
  testthat::expect_lt(rel_diff(seq_out$hessian, seq_out$hessian), 1e-15)
  for (f in seq_len(4L)) {
    testthat::expect_lt(
      rel_diff(bat_out$hessian_fold_sums[[f]], seq_out$hessian_fold_sums[[f]]),
      1e-12
    )
  }
  testthat::expect_identical(bat_out$hessian_fold_counts,
                             seq_out$hessian_fold_counts)
  testthat::expect_equal(sum(bat_out$hessian_fold_counts), n)

  # 装配模式不进入簇求解，逐簇对象按位相同
  testthat::expect_identical(bat_out$gamma, seq_out$gamma)
  testthat::expect_identical(bat_out$cluster_scores, seq_out$cluster_scores)
  testthat::expect_identical(bat_out$objective, seq_out$objective)
  testthat::expect_identical(bat_out$residuals, seq_out$residuals)
  testthat::expect_identical(bat_out$nuisance_diagnostics,
                             seq_out$nuisance_diagnostics)
})

testthat::test_that("the batched assembly reproduces the identity audit", {
  withr::local_envvar(c(JPDMEQR_IDENTITY_DETAIL = "cluster"))
  d <- assembly_fixture()
  seq_out <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id,
                                   as.numeric(d$beta), 0.5, 0.4,
                                   need_hessian = TRUE, identity_check = TRUE,
                                   assembly = "sequential")
  bat_out <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id,
                                   as.numeric(d$beta), 0.5, 0.4,
                                   need_hessian = TRUE, identity_check = TRUE,
                                   assembly = "batch")

  testthat::expect_lt(rel_diff(bat_out$hessian_identity,
                               seq_out$hessian_identity), 1e-12)
  testthat::expect_lt(bat_out$profile_identity_error, 1e-9)
  testthat::expect_identical(bat_out$cluster_identity_errors,
                             seq_out$cluster_identity_errors)
})

testthat::test_that("the assembly mode is selected by the environment", {
  d <- assembly_fixture()
  env_out <- withr::with_envvar(
    c(JPDMEQR_ASSEMBLY = "batch"),
    profile_components_v2(d$y, d$X, d$Z, d$cluster_id, as.numeric(d$beta), 0.5,
                          0.4, need_hessian = FALSE)
  )
  testthat::expect_identical(env_out$cluster_scores,
                             profile_components_v2(d$y, d$X, d$Z, d$cluster_id,
                                                   as.numeric(d$beta), 0.5, 0.4,
                                                   need_hessian = FALSE)$cluster_scores)
})

testthat::test_that("an unknown assembly mode raises an error", {
  d <- assembly_fixture()
  testthat::expect_error(
    profile_components_v2(d$y, d$X, d$Z, d$cluster_id, as.numeric(d$beta), 0.5,
                          0.4, need_hessian = TRUE, assembly = "stream"),
    "Unknown assembly mode: stream"
  )
})

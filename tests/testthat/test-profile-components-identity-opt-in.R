root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)

profile_identity_fixture <- function(seed = 20260921L) {
  dat <- generate_profile_qr_data_v2(14L, 6L, 2L, 0.5, q = 1L,
                                     m_values = 3:5, seed = seed)
  list(dat = dat, h = 0.4,
       beta = dat$beta0 + c(0.06, -0.03, 0.02, -0.01, 0, 0))
}

testthat::test_that("identity Hessian is opt-in and leaves inferential objects untouched", {
  fx <- profile_identity_fixture()
  d <- fx$dat
  base <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id, fx$beta, 0.5, fx$h,
                                need_hessian = TRUE)
  audit <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id, fx$beta, 0.5, fx$h,
                                 need_hessian = TRUE, identity_check = TRUE)

  testthat::expect_false(base$identity_checked)
  testthat::expect_true(is.null(base$hessian_identity))
  testthat::expect_true(is.na(base$profile_identity_error))
  testthat::expect_true(all(is.na(base$cluster_identity_errors)))
  testthat::expect_false(is.null(base$hessian))
  testthat::expect_equal(base$hessian, t(base$hessian), tolerance = 0)

  testthat::expect_true(audit$identity_checked)
  testthat::expect_lt(audit$profile_identity_error, 1e-9)
  testthat::expect_lt(max(abs(audit$hessian - audit$hessian_identity)), 1e-9)

  # 审计开关只增加两个等价的 Hessian 表达式，不改变任何推断对象，逐位相同
  testthat::expect_identical(audit$hessian, base$hessian)
  testthat::expect_identical(audit$score, base$score)
  testthat::expect_identical(audit$score0, base$score0)
  testthat::expect_identical(audit$objective, base$objective)
  testthat::expect_identical(audit$gamma, base$gamma)
  testthat::expect_identical(audit$cluster_scores, base$cluster_scores)
  testthat::expect_identical(audit$unsmoothed_cluster_scores, base$unsmoothed_cluster_scores)
  testthat::expect_identical(audit$hessian_fold_sums, base$hessian_fold_sums)
  testthat::expect_identical(audit$hessian_fold_counts, base$hessian_fold_counts)
  testthat::expect_identical(audit$nuisance_diagnostics, base$nuisance_diagnostics)
})

testthat::test_that("cluster bookkeeping keeps its column types and fold assignment", {
  fx <- profile_identity_fixture()
  d <- fx$dat
  comp <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id, fx$beta, 0.5, fx$h,
                                need_hessian = TRUE)
  diag <- comp$nuisance_diagnostics

  testthat::expect_identical(
    names(diag),
    c("cluster", "m_i", "converged", "gradient_max", "optim_convergence",
      "function_evaluations", "gradient_evaluations", "newton_iterations")
  )
  testthat::expect_identical(diag$cluster, unique(as.character(d$cluster_id)))
  testthat::expect_identical(rownames(comp$gamma), diag$cluster)
  testthat::expect_type(diag$converged, "logical")
  testthat::expect_type(diag$m_i, "integer")
  testthat::expect_type(diag$optim_convergence, "integer")
  testthat::expect_type(diag$function_evaluations, "integer")
  testthat::expect_type(diag$gradient_evaluations, "integer")
  testthat::expect_type(diag$gradient_max, "double")
  testthat::expect_identical(comp$nuisance_converged, all(diag$converged))
  testthat::expect_identical(comp$max_nuisance_gradient, max(diag$gradient_max))
  testthat::expect_identical(sum(comp$hessian_fold_counts), nrow(diag))
  testthat::expect_true(all(comp$hessian_fold_counts > 0L))
  total <- Reduce(`+`, comp$hessian_fold_sums)
  testthat::expect_lt(max(abs(total / nrow(diag) - comp$hessian)), 1e-14)
})

testthat::test_that("fitted components carry the audit identity", {
  fx <- profile_identity_fixture()
  d <- fx$dat
  fit <- fit_profile_lasso_v2(d$y, d$X, d$Z, d$cluster_id, 0.5, fx$h,
                              lambda_beta = 0.05, control = list(max_iter = 200L))
  testthat::expect_true(fit$components$identity_checked)
  testthat::expect_true(is.finite(fit$components$profile_identity_error))
  testthat::expect_lt(fit$components$profile_identity_error, 1e-8)
  testthat::expect_true(is.matrix(fit$components$hessian_identity))
  testthat::expect_identical(sum(fit$components$hessian_fold_counts),
                             fit$components$n_clusters)
})

testthat::test_that("parallel cluster branch reproduces the serial accumulation", {
  testthat::skip_on_os("windows")
  fx <- profile_identity_fixture()
  d <- fx$dat
  old <- getOption("jpDMEQR.cluster_cores")
  on.exit(options(jpDMEQR.cluster_cores = old), add = TRUE)

  options(jpDMEQR.cluster_cores = 1L)
  serial <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id, fx$beta, 0.5, fx$h,
                                  need_hessian = TRUE, identity_check = TRUE)
  options(jpDMEQR.cluster_cores = 3L)
  par <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id, fx$beta, 0.5, fx$h,
                               need_hessian = TRUE, identity_check = TRUE)

  testthat::expect_identical(par$hessian, serial$hessian)
  testthat::expect_identical(par$hessian_fold_sums, serial$hessian_fold_sums)
  testthat::expect_identical(par$hessian_fold_counts, serial$hessian_fold_counts)
  testthat::expect_identical(par$score, serial$score)
  testthat::expect_identical(par$score0, serial$score0)
  testthat::expect_identical(par$objective, serial$objective)
  testthat::expect_identical(par$cluster_identity_errors, serial$cluster_identity_errors)
  testthat::expect_identical(par$nuisance_diagnostics, serial$nuisance_diagnostics)
})

testthat::test_that("bookkeeping columns survive the warm-start nuisance path", {
  fx <- profile_identity_fixture()
  d <- fx$dat
  reference <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id, fx$beta, 0.5, fx$h,
                                     need_hessian = TRUE)
  old <- Sys.getenv("JPDMEQR_CLUSTER_SOLVER", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("JPDMEQR_CLUSTER_SOLVER") else
      Sys.setenv(JPDMEQR_CLUSTER_SOLVER = old)
  }, add = TRUE)
  Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "newton")
  fast <- profile_components_v2(d$y, d$X, d$Z, d$cluster_id, fx$beta, 0.5, fx$h,
                                need_hessian = TRUE, identity_check = TRUE)

  testthat::expect_true(fast$nuisance_converged)
  testthat::expect_identical(
    vapply(fast$nuisance_diagnostics, function(x) class(x)[1], character(1)),
    vapply(reference$nuisance_diagnostics, function(x) class(x)[1], character(1))
  )
  testthat::expect_identical(fast$nuisance_diagnostics$cluster,
                             reference$nuisance_diagnostics$cluster)
  testthat::expect_identical(nrow(fast$nuisance_diagnostics), nrow(reference$nuisance_diagnostics))
  testthat::expect_identical(sum(fast$hessian_fold_counts), nrow(fast$nuisance_diagnostics))
  testthat::expect_lt(fast$profile_identity_error, 1e-9)
  testthat::expect_lt(max(abs(fast$hessian - reference$hessian)), 1e-6)
})
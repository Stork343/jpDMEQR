root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled", JPDMEQR_UNPENALISED_OPTIMIZER = "newton_armijo",
           JPDMEQR_PROXIMAL_SCALING = "cluster_rms", JPDMEQR_ASSEMBLY = "batch",
           JPDMEQR_HESSIAN_REFRESH = "20", JPDMEQR_MATRIX_BACKEND = "Eigen_single_thread",
           JPDMEQR_IDENTITY_DETAIL = "aggregate", JPDMEQR_DANTZIG_SOLVER = "dual_warm",
           JPDMEQR_SHARED_FIRST_STAGE = "TRUE")
options(jpDMEQR.matrix_backend = "Eigen_single_thread", jpDMEQR.identity_detail = "aggregate")
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
v2_require_root <- function() root
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a new certificate numerical run ID.")
output <- file.path(root, "results", "scratch", args[1])
if (dir.exists(output)) stop("The certificate numerical directory already exists.")
dir.create(output, recursive = TRUE)
source_path <- file.path(root, "results", "scratch", "initial_scaling_none_20261022_v1", "source_design.rds")
observed <- readRDS(source_path)
md <- observed$metadata
train <- observed$training
design <- observed$design
request <- list(source_design = source_path, scope = "observed_training_numerical_certificate_and_fixed_penalty_stationarity",
                actual_source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                evaluated_initial_function = fit_application_profile_initial_v2,
                solver_configuration = solver_configuration_v2(),
                research_low_penalty_multiplier = .01,
                production_prediction_penalty_grid_unchanged = TRUE)
saveRDS(request, file.path(output, "validation_request.rds"))
unit <- as.data.frame(testthat::test_file(
  file.path(root, "tests", "testthat", "test-application-initial-clinical-certificate.R"), reporter = "summary"
))
write.csv(unit[c("test", "nb", "failed", "skipped", "error", "warning", "passed")],
           file.path(output, "unit_test_summary.csv"), row.names = FALSE)
if (any(unit$failed > 0 | unit$error | unit$skipped | unit$warning > 0)) stop("The real clinical certificate unit check failed.")
reporter <- testthat::ListReporter$new()
summary <- list()
withCallingHandlers({
for (q in 1:2) for (tau in c(.25, .5, .75)) {
  Z <- if (q == 2L) cbind(intercept = 1, centered_time = md$cumulative_time - design$clinical$time_center) else design$Z
  fit <- fit_application_profile_initial_v2(
    md$log_c3[train], design$X[train, , drop = FALSE], Z[train, , drop = FALSE], md$subject_id[train], tau,
    penalty_factor = design$penalty_factor, control = list(clinical_acceleration = "clinical_kkt")
  )
  reference <- fit_profile_lasso_v2(
    md$log_c3[train], design$X[train, , drop = FALSE], Z[train, , drop = FALSE], md$subject_id[train], tau,
    fit$h_est, fit$lambda_beta, lambda_gamma = fit$lambda_gamma, penalty_factor = fit$penalty_factor,
    beta_init = fit$beta_l1,
    control = list(max_iter = 2000L, kkt_normalized_tol = 1e-7, beta_tol = 1e-9,
                    nuisance_control = list(grad_tol = 1e-11))
  )
  profile_pair_check_fit_v2(reference, "fixed_penalty_stationarity_check", 1e-7, 1e-9)
  result <- data.frame(q = q, tau = tau, beta_difference = max(abs(fit$beta_l1 - reference$beta)),
                       score_difference = max(abs(fit$components$score - reference$components$score)),
                       hessian_difference = max(abs(fit$components$hessian - reference$components$hessian)),
                       hessian_relative_difference = max(abs(fit$components$hessian - reference$components$hessian)) /
                         max(1, max(abs(fit$components$hessian))),
                       certified_identity_error = fit$components$profile_identity_error,
                       reference_identity_error = reference$components$profile_identity_error,
                       objective_difference = abs(fit$composite_objective - reference$composite_objective),
                       prediction_difference = max(abs(design$X %*% (fit$beta_l1 - reference$beta))),
                       reference_kkt = reference$kkt_normalized,
                       certified_kkt = fit$kkt_normalized,
                       scope = "fixed_penalty_stationarity_at_globally_certified_point")
  saveRDS(list(certified = fit, fixed_penalty_reference = reference, summary = result),
           file.path(output, sprintf("fixed_penalty_q%d_tau%02d.rds", q, 100 * tau)))
  summary[[length(summary) + 1L]] <- result
  testthat::with_reporter(reporter, {
    testthat::test_that(sprintf("真实q%d tau%.2f完整固定penalty接受点检查", q, tau), {
      testthat::expect_lte(result$beta_difference, 1e-6)
      testthat::expect_lte(result$score_difference, 1e-7)
      testthat::expect_lte(result$hessian_relative_difference, 1e-10)
      testthat::expect_lte(result$certified_identity_error, 1e-8)
      testthat::expect_lte(result$reference_identity_error, 1e-8)
      testthat::expect_lte(result$objective_difference, 1e-8)
      testthat::expect_lte(result$prediction_difference, 1e-6)
      testthat::expect_lte(result$reference_kkt, 1e-7)
      testthat::expect_lte(result$certified_kkt, 1e-7)
      testthat::expect_true(all(vapply(fit$initial_optimizer_history,
                                       function(pass) pass$global_certificate_accepted, logical(1))))
    })
  })
}
write.csv(do.call(rbind, summary), file.path(output, "fixed_penalty_stationarity_summary.csv"), row.names = FALSE)
low <- fit_application_profile_initial_v2(
  md$log_c3[train], design$X[train, , drop = FALSE], design$Z[train, , drop = FALSE], md$subject_id[train], .5,
  lambda_multiplier = .01, penalty_factor = design$penalty_factor,
  control = list(clinical_acceleration = "clinical_kkt")
)
saveRDS(low, file.path(output, "actual_low_penalty_continuation.rds"))
testthat::with_reporter(reporter, {
  testthat::test_that("真实低penalty的完整proximal继续分支保留原始载荷和接受检查", {
    testthat::expect_true(low$converged)
    testthat::expect_lte(low$kkt_normalized, 1e-3)
    testthat::expect_gt(sum(abs(low$beta_l1[design$penalty_factor > 0]) > 1e-10), 0L)
    testthat::expect_identical(low$loadings$final, pmax(low$loadings$pass0, low$loadings$pass1))
    testthat::expect_equal(unname(low$coordinate_penalty), low$lambda_beta * low$penalty_factor, tolerance = 0)
    for (pass in low$initial_optimizer_history) {
      testthat::expect_false(pass$global_certificate_accepted)
      testthat::expect_identical(pass$optimizer, "proximal_gradient")
      testthat::expect_equal(unname(pass$proximal_initial_beta[design$penalty_factor > 0]), rep(0, sum(design$penalty_factor > 0)))
      testthat::expect_equal(pass$proximal_initial_beta[design$penalty_factor == 0], pass$clinical_beta, tolerance = 0)
      testthat::expect_gt(nrow(pass$proximal_history), 0L)
      testthat::expect_lte(tail(pass$proximal_history$composite_objective, 1),
                            pass$proximal_history$composite_objective[1] + 1e-12)
      testthat::expect_identical(pass$coordinate_penalty,
                                  setNames(low$lambda_beta * design$penalty_factor *
                                             if (pass$pass == "preliminary") low$loadings$pass0 else low$loadings$final,
                                           colnames(design$X)))
    }
  })
})
}, error = function(condition) saveRDS(condition, file.path(output, "failure_condition.rds")))
checks <- as.data.frame(reporter$get_results())
write.csv(checks[c("test", "nb", "failed", "skipped", "error", "warning", "passed")],
           file.path(output, "test_summary.csv"), row.names = FALSE)
if (any(checks$failed > 0 | checks$error | checks$skipped | checks$warning > 0)) stop("The actual initial certificate checks failed.")
jsonlite::write_json(list(status = "completed", unit_assertions = sum(unit$passed),
                          stationarity_and_continuation_assertions = sum(checks$passed),
                          actual_source_commit = request$actual_source_commit,
                          solver_configuration = request$solver_configuration),
                     file.path(output, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

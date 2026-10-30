.libPaths(c(normalizePath("results/scratch/dependencies/library", mustWork = TRUE), .libPaths()))
if (utils::packageVersion("highs") < "1.15.1-1") stop("Run the declared dependency installation before comparison fitting.")
source("research/submission_completion_20261027/load.R")
source("research/submission_completion_20261027/clime_primal.R")
source("research/submission_completion_20261027/clime_parametric.R")
source("research/submission_completion_20261027/clime_certified.R")
source("research/submission_completion_20261027/clime_parametric_certified.R")
source("research/submission_completion_20261027/clime_qp_candidate.R")
source("research/submission_completion_20261027/sqr_cv_pruning.R")
source("research/submission_completion_20261027/qgee_lapack.R")
source("research/submission_completion_20261027/qgee_algebra.R")
source("research/submission_completion_20261027/quantile_column_generation.R")
source("research/submission_completion_20261027/quantile_highs.R")
source("research/submission_completion_20261027/load_compiled_clime.R")
submission_numerical_specification <- yaml::read_yaml("config/submission_comparator_numerics_20261030.yml")
submission_sqr_specification <- submission_numerical_specification[["SQR-DEBIASED-IID"]]
if (!is.list(submission_sqr_specification) ||
    !identical(submission_sqr_specification$precision_optimizer, "compiled_certified") ||
    !isTRUE(submission_sqr_specification$complete_gamma_grid) ||
    !identical(submission_sqr_specification$native_path_maximum_pivots, "4p") ||
    !identical(as.integer(submission_sqr_specification$lambda_subject_folds), 5L) ||
    !identical(as.integer(submission_sqr_specification$gamma_subject_folds), 5L))
  stop("The SQR numerical specification does not match the loaded implementation.")
compile_comparator_geometry <- function() {
  lock <- filelock::lock("results/scratch/submission_completion_20261027/qgee_compile.lock", timeout = 60000)
  if (is.null(lock)) stop("Unable to acquire the QGEE compilation lock.")
  on.exit(filelock::unlock(lock), add = TRUE)
  Rcpp::sourceCpp("research/submission_completion_20261027/qgee_algebra.cpp",
    cacheDir = "results/scratch/submission_completion_20261027/qgee_cache", rebuild = FALSE,
    showOutput = FALSE, env = .GlobalEnv)
}
compile_comparator_geometry()
rm(compile_comparator_geometry)
submission_comparator_reference <- fit_benchmark_frozen_tuning_v2
fit_benchmark_frozen_tuning_v2 <- function(method_id, train, tau, target_coords, tuning, seed,
                                          context = list(), control = list()) {
  if (method_id == "SQR-DEBIASED-IID") {
    control$solver_opts <- list(algorithm = "compiled_certified", native_guard = TRUE, column_loss_bounds = TRUE)
    control$exact_cv_pruning <- TRUE
  }
  if (method_id == "POOLED-QR-LASSO") control$pooled_optimizer <- "dual_lp"
  if (method_id == "DOUBLE-PEN-QLMM") control$raw_sum_penalty_scales <- TRUE
  if (method_id == "QGEE-SCAD") {
    control$batch_geometry <- TRUE
    control$standardize <- FALSE
  }
  submission_comparator_reference(method_id, train, tau, target_coords, tuning, seed, context, control)
}

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
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[1] %in% c("parity", "strict", "grid")) stop("Supply parity, strict or grid and a new run ID.")
output <- file.path(root, "results", "scratch", args[2])
if (dir.exists(output)) stop("The initial estimator numerical directory already exists.")
dir.create(output, recursive = TRUE)
source_path <- file.path(root, "results", "scratch", "initial_scaling_none_20261022_v1", "source_design.rds")
observed <- readRDS(source_path)
metadata <- observed$metadata
training <- observed$training
base_design <- observed$design
request <- list(scope = "observed_training_fold_numerical_validation", source_design = source_path,
                actual_source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                evaluated_initial_function = fit_application_profile_initial_v2,
                training_subjects = sort(unique(as.character(metadata$subject_id[training]))),
                solver_configuration = solver_configuration_v2(), modes = c("none", "clinical_kkt"))
saveRDS(request, file.path(output, "validation_request.rds"))
designs <- list(`100` = base_design)
if (args[1] == "grid") {
  rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
  members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
  signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
  annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
  cohort <- readRDS(file.path(root, "data", "derived", "GSE65391", "audit_20261022_bootstrap_v1", "analysis_cohort.rds"))
  expression <- readRDS(cohort$expression_file)[, rownames(metadata), drop = FALSE]
  full <- prepare_gse65391_prediction_design_v2(
    metadata, expression, annotation, signatures, members, rule, training,
    include_sledai = TRUE, filter_size = 5000L, screening_dimension = 400L, screening_tau = .5
  )
  clinical <- which(full$penalty_factor == 0)
  genes <- which(full$penalty_factor > 0)
  for (dimension in c(100L, 200L, 400L)) {
    value <- full
    keep <- c(clinical, genes[seq_len(dimension)])
    value$X <- full$X[, keep, drop = FALSE]
    value$penalty_factor <- full$penalty_factor[keep]
    value$highdim_preprocessing$screening_dimension <- dimension
    value$highdim_preprocessing$selected_genes <- colnames(value$X)[value$penalty_factor > 0]
    designs[[as.character(dimension)]] <- value
  }
  if (!identical(designs[["100"]]$X, base_design$X)) stop("The real source training design changed.")
  saveRDS(designs, file.path(output, "three_observed_training_designs.rds"))
}
grid <- if (args[1] %in% c("parity", "strict")) {
  expand.grid(q = 1:2, tau = c(.25, .5, .75), h_multiplier = 1, lambda_gamma = 1,
               lambda_multiplier = 1, screen_dimension = 100L)
} else {
  execution <- yaml::read_yaml(file.path(root, "config", "application_execution_20261022.yml"))
  expand.grid(q = 1L, tau = .5,
               h_multiplier = execution$prediction$h_multiplier_grid,
               lambda_gamma = execution$prediction$lambda_gamma_grid,
               lambda_multiplier = execution$prediction$lambda_multiplier_grid,
               screen_dimension = execution$prediction$screening_dimensions)
}
grid$candidate_id <- seq_len(nrow(grid))
write.csv(grid, file.path(output, "scheduled_candidates.csv"), row.names = FALSE)
summary <- comparisons <- list()
cache <- new.env(parent = emptyenv())
reporter <- testthat::ListReporter$new()
withCallingHandlers({
for (index in seq_len(nrow(grid))) {
  candidate <- grid[index, ]
  design <- designs[[as.character(candidate$screen_dimension)]]
  Z <- if (candidate$q == 2L) {
    cbind(intercept = 1, centered_time = metadata$cumulative_time - design$clinical$time_center)
  } else design$Z
  h <- candidate$h_multiplier * (log(sum(design$penalty_factor > 0)) /
                                  length(unique(metadata$subject_id[training])))^(1 / 4)
  fits <- list()
  for (mode in c("none", "clinical_kkt")) {
    started <- proc.time()[[3]]
    fit <- fit_application_profile_initial_v2(
      metadata$log_c3[training], design$X[training, , drop = FALSE], Z[training, , drop = FALSE],
      metadata$subject_id[training], candidate$tau, h_est = h,
      lambda_gamma = candidate$lambda_gamma, lambda_multiplier = candidate$lambda_multiplier,
      penalty_factor = design$penalty_factor,
      control = list(clinical_acceleration = mode, clinical_cache = cache,
                     fit_control = if (args[1] == "strict") {
                       list(max_iter = 20000L, kkt_normalized_tol = 1e-7, beta_tol = 1e-9)
                     } else list())
    )
    elapsed <- proc.time()[[3]] - started
    fits[[mode]] <- fit
    saveRDS(fit, file.path(output, sprintf("candidate%03d_%s.rds", index, mode)))
    summary[[length(summary) + 1L]] <- data.frame(
      candidate, mode = mode, runtime_sec = elapsed, iterations = fit$iterations,
      kkt_normalized = fit$kkt_normalized, kkt_absolute = fit$kkt_residual,
      objective = fit$composite_objective, selected_genes = sum(abs(fit$beta_l1[design$penalty_factor > 0]) > 1e-10),
      preliminary_certificate = fit$initial_optimizer_history$preliminary$global_certificate_accepted,
      final_certificate = fit$initial_optimizer_history$final$global_certificate_accepted,
      preliminary_cache_hit = fit$initial_optimizer_history$preliminary$cache_hit,
      final_cache_hit = fit$initial_optimizer_history$final$cache_hit
    )
    write.csv(do.call(rbind, summary), file.path(output, "actual_fit_summary.csv"), row.names = FALSE)
  }
  original <- fits$none
  accelerated <- fits$clinical_kkt
  difference <- data.frame(
    candidate, beta_difference = max(abs(original$beta_l1 - accelerated$beta_l1)),
    objective_difference = abs(original$composite_objective - accelerated$composite_objective),
    prediction_difference = max(abs(design$X %*% (original$beta_l1 - accelerated$beta_l1))),
    loading_difference = max(abs(original$loadings$final - accelerated$loadings$final)),
    selection_identical = identical(abs(original$beta_l1[design$penalty_factor > 0]) > 1e-10,
                                    abs(accelerated$beta_l1[design$penalty_factor > 0]) > 1e-10)
  )
  comparisons[[length(comparisons) + 1L]] <- difference
  write.csv(do.call(rbind, comparisons), file.path(output, "actual_comparisons.csv"), row.names = FALSE)
  testthat::with_reporter(reporter, {
    testthat::test_that(sprintf("真实训练设计candidate %d完整KKT计算与原路径记录", index), {
      testthat::expect_true(original$converged && accelerated$converged)
      testthat::expect_lte(original$kkt_normalized, 1e-3)
      testthat::expect_lte(accelerated$kkt_normalized, 1e-3)
      if (args[1] != "grid") {
        testthat::expect_lte(difference$beta_difference, 1e-6)
        testthat::expect_lte(difference$objective_difference, 1e-8)
        testthat::expect_lte(difference$prediction_difference, 1e-6)
      }
      testthat::expect_true(difference$selection_identical)
      testthat::expect_identical(original$loadings$pass0, accelerated$loadings$pass0)
      for (pass in accelerated$initial_optimizer_history) {
        if (pass$global_certificate_accepted) {
          testthat::expect_lte(pass$global_kkt_normalized, pass$global_certificate_threshold)
          testthat::expect_lte(max(pass$global_residual), pass$global_certificate_threshold)
          testthat::expect_lte(pass$full_profile_identity_error, 1e-8)
        }
      }
    })
  })
  cat(sprintf("completed actual candidate %d/%d\n", index, nrow(grid)))
}
}, error = function(condition) {
  saveRDS(list(condition = condition, candidate = if (exists("candidate")) candidate else NULL),
          file.path(output, "failure_condition.rds"))
})
checks <- as.data.frame(reporter$get_results())
write.csv(checks[c("test", "nb", "failed", "skipped", "error", "warning", "passed")],
           file.path(output, "test_summary.csv"), row.names = FALSE)
if (any(checks$failed > 0 | checks$error | checks$skipped | checks$warning > 0)) stop("Initial clinical numerical comparison failed.")
jsonlite::write_json(list(status = "completed", candidates = nrow(grid), actual_fits = nrow(grid) * 2L,
                          scope = request$scope, actual_source_commit = request$actual_source_commit,
                          comparison_scope = if (args[1] == "grid") {
                            "full_frozen_grid_global_KKT_and_actual_timing_with_zero_start_differences_recorded"
                          } else "strict_zero_start_numerical_comparison",
                          solver_configuration = request$solver_configuration),
                     file.path(output, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

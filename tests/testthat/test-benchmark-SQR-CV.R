root <- v2_require_root()
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"), local = TRUE)
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, environment())
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"), local = TRUE)

sqr_cv_actual_example <- function() {
  registry <- read_simulation_registry_v2(file.path(root, "tests", "fixtures", "simulation_streaming_pilot_20261022.yml"))
  config <- config_row_to_list_v2(registry[1L, , drop = FALSE])
  config$n_clusters <- 50L
  config$p <- 6L
  config$s <- 2L
  config$experiment_id <- "SQR-CV-ACTUAL-DGP"
  config$methods <- "SQR-DEBIASED-IID"
  dat <- simulate_from_config_v2(config, 20261023L)
  list(data = dat, config = config, tuning = build_tuning_for_analysis_v2(config, dat),
       targets = colnames(dat$X)[c(1L, 4L)])
}

testthat::test_that("SQR-DEBIASED-IID source CV uses the complete grids within independent subjects", {
  example <- sqr_cv_actual_example()
  actual <- fit_benchmark_frozen_tuning_v2(
    "SQR-DEBIASED-IID", example$data, .5, example$targets, example$tuning,
    20261024L, context = list(config = example$config)
  )
  testthat::expect_identical(actual$status, "ok")
  testthat::expect_identical(actual$method_version, "source_20261023_CV")
  diagnostic <- actual$tuning_diagnostics$gamma_cv
  testthat::expect_identical(diagnostic$configured_multipliers, example$config$dantzig_multipliers)
  testthat::expect_equal(nrow(diagnostic$cells), 5L * length(example$config$dantzig_multipliers))
  testthat::expect_equal(nrow(diagnostic$candidates), length(example$config$dantzig_multipliers))
  testthat::expect_equal(sum(diagnostic$candidates$selected), 1L)
  testthat::expect_true(all(diagnostic$candidates$successful_folds[diagnostic$candidates$eligible] == 5L))
  testthat::expect_true(any(diagnostic$cells$status == "excluded_zero_direction"))
  nesting <- vapply(seq_len(5L), function(fold) {
    outer <- diagnostic$roles[[fold]]
    inner <- diagnostic$inner_lambda_cv[[fold]]
    length(intersect(outer$training_subjects, outer$validation_subjects)) == 0L &&
      setequal(inner$training_subjects, outer$training_subjects) &&
      all(vapply(inner$roles, function(role) {
        all(c(role$training_subjects, role$validation_subjects) %in% outer$training_subjects) &&
          !any(c(role$training_subjects, role$validation_subjects) %in% outer$validation_subjects)
      }, logical(1))) &&
      identical(inner$configured_multipliers, example$config$lambda_beta_multipliers)
  }, logical(1))
  testthat::expect_true(all(nesting))
  testthat::expect_length(actual$tuning_values$mu_grid, 1L)
  testthat::expect_true(all(actual$inference_object$table$mu == diagnostic$selected_gamma))
  eligible <- diagnostic$candidates[diagnostic$candidates$eligible, ]
  selected <- eligible[order(eligible$loss, eligible$multiplier), ][1L, ]
  testthat::expect_equal(diagnostic$selected_multiplier, selected$multiplier)
  for (fold in seq_len(5L)) {
    rows <- diagnostic$cells[diagnostic$cells$fold == fold, ]
    source <- diagnostic$fold_initial_estimators[[fold]]
    testthat::expect_true(all(rows$lambda_source[rows$status == "ok"] == source$lambda))
  }
  good <- which(diagnostic$cells$status == "ok")[1L]
  cell <- diagnostic$cells[good, ]
  role <- diagnostic$roles[[cell$fold]]
  training <- subset_clusters_v2(example$data, role$training_subjects)
  validation <- subset_clusters_v2(example$data, role$validation_subjects)
  initial <- diagnostic$fold_initial_estimators[[cell$fold]]
  H_train <- sqr_iid_hessian_fn_v2(training$y, training$X, initial$beta, .5, initial$h)
  H_validation <- sqr_iid_hessian_fn_v2(validation$y, validation$X, initial$beta, .5, initial$h)
  precision <- solve_clime_matrix_iid_v2(H_train, cell$gamma)
  independent_loss <- sum((H_validation %*% precision$omega_symmetric - diag(ncol(training$X)))^2) /
    ncol(training$X)
  testthat::expect_equal(cell$loss, independent_loss, tolerance = 1e-10)
  testthat::expect_equal(nrow(diagnostic$precision_records[[good]]$result$actual_LP_attempts), ncol(training$X))
  testthat::expect_equal(actual$tuning_values$lambda_beta_multiplier,
                         actual$tuning_diagnostics$selected_multiplier)
})

testthat::test_that("SQR-DEBIASED-IID source CV is deterministic and independent of generated truth", {
  example <- sqr_cv_actual_example()
  select <- function(data) {
    select_sqr_gamma_cv_v2(
      data, .5, example$targets, example$config, example$tuning, 20261025L,
      lambda_selector = function(training, seed) {
        select_frozen_method_tuning_v2("SQR-DEBIASED-IID", training, .5, example$targets,
                                        example$config, build_tuning_for_analysis_v2(example$config, training), seed)
      }, fold_tuning_builder = function(training) build_tuning_for_analysis_v2(example$config, training)
    )
  }
  a <- select(example$data)
  altered <- example$data
  altered$beta0 <- NULL
  altered$active <- NULL
  altered$random_effects <- NULL
  b <- select(altered)
  testthat::expect_equal(a$gamma, b$gamma, tolerance = 0)
  testthat::expect_equal(a$multiplier, b$multiplier, tolerance = 0)
  testthat::expect_identical(a$diagnostics$cells, b$diagnostics$cells)
  testthat::expect_identical(a$diagnostics$candidates, b$diagnostics$candidates)
  testthat::expect_equal(a$diagnostics$fold_initial_estimators, b$diagnostics$fold_initial_estimators,
                         tolerance = 0)
})

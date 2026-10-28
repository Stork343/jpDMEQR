tuning_root <- normalizePath(v2_require_root(), mustWork = TRUE)
root <- tuning_root
source(file.path(root, "scripts", "00_source_v2.R"), local = TRUE)
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, environment())
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"), local = TRUE)
withr::local_envvar(c(JPDMEQR_UNPENALISED_OPTIMIZER = "proximal_gradient",
                      JPDMEQR_PROXIMAL_SCALING = "none",
                      JPDMEQR_HESSIAN_REFRESH = "1", JPDMEQR_IDENTITY_DETAIL = "cluster",
                      JPDMEQR_DANTZIG_SOLVER = "clarabel_native"),
                    .local_envir = environment())
frozen_tuning_test_results <- getOption("jpDMEQR.frozen_tuning_test_results",
                                       new.env(parent = emptyenv()))

frozen_tuning_test_data <- function(n = 40L, p = 4L, seed = 20261103L) {
  registry <- read_simulation_registry_v2(file.path(root, "config", "simulation_pair_main_20261022.yml"))
  cfg <- config_row_to_list_v2(registry[registry$experiment_id == "E07", , drop = FALSE])
  cfg$n_clusters <- n
  cfg$p <- p
  cfg$s <- 1L
  cfg$target_coordinate_count <- 1L
  dat <- simulate_from_config_v2(cfg, seed)
  list(config = cfg, data = dat, tuning = build_tuning_for_analysis_v2(cfg, dat))
}

testthat::test_that("E07 executes all paired multipliers and the selected complete training fit", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton",
                        JPDMEQR_DANTZIG_SOLVER = "clarabel_native",
                        JPDMEQR_ASSEMBLY = "sequential", JPDMEQR_CV_CORES = "1",
                        JPDMEQR_HESSIAN_REFRESH = "1", JPDMEQR_SHARED_FIRST_STAGE = "FALSE"))
  ex <- frozen_tuning_test_data()
  fit <- fit_benchmark_frozen_tuning_v2(
    "PROFILE-DQR", ex$data, ex$config$tau, colnames(ex$data$X)[1],
    ex$tuning, 20261104L, context = list(config = ex$config, active = ex$data$active)
  )
  d <- fit$tuning_diagnostics
  frozen_tuning_test_results$paired <- fit
  testthat::expect_identical(fit$status, "ok")
  testthat::expect_equal(nrow(d$candidates), 8L)
  testthat::expect_equal(nrow(d$fold_losses), 40L)
  testthat::expect_equal(d$candidates$multiplier, sort(unique(ex$config$lambda_beta_multipliers)))
  testthat::expect_true(any(d$candidates$eligible))
  testthat::expect_equal(sum(d$candidates$selected), 1L)
  testthat::expect_equal(fit$fit_object$lambda_multiplier, d$selected_multiplier)
  testthat::expect_identical(fit$fit_object$score_definition,
                             "profile_bandwidth_pair_2_half_minus_full")
  testthat::expect_equal(fit$fit_object$lambda_components$coordinate,
                         d$selected_multiplier * fit$fit_object$lambda_base *
                           fit$fit_object$penalty_factor)
  testthat::expect_setequal(d$full_training$subject_ids, unique(as.character(ex$data$cluster_id)))
  testthat::expect_equal(d$full_training$multiplier, d$selected_multiplier)
  testthat::expect_true(all(vapply(d$roles, function(role) {
    !length(intersect(role$training_subjects, role$validation_subjects)) &&
      setequal(c(role$training_subjects, role$validation_subjects), d$training_subjects)
  }, logical(1))))
  first <- which(d$fold_losses$status == "ok")[1]
  row <- d$fold_losses[first, ]
  role <- d$roles[[row$fold]]
  changed <- ex$data
  changed$y[changed$cluster_id %in% role$validation_subjects] <-
    changed$y[changed$cluster_id %in% role$validation_subjects] + 7
  training <- subset_clusters_v2(changed, role$training_subjects)
  ft <- build_tuning_for_analysis_v2(ex$config, training)
  repeated <- fit_frozen_tuning_training_v2(
    "PROFILE-DQR", training, ex$config$tau, colnames(ex$data$X)[1], ft,
    row$multiplier, row$seed, context = list(active = ex$data$active)
  )
  testthat::expect_equal(repeated$beta, d$training_fits[[first]]$fit$beta, tolerance = 0)
  testthat::expect_equal(repeated$fit_summary$lambda_components,
                         d$training_fits[[first]]$fit$lambda_components, tolerance = 0)
  cfg_fixed <- ex$config
  cfg_fixed$experiment_id <- "E01"
  fixed <- select_frozen_method_tuning_v2(
    "PROFILE-DQR", ex$data, ex$config$tau, colnames(ex$data$X)[1], cfg_fixed,
    ex$tuning, 20261104L
  )
  testthat::expect_equal(fixed$diagnostics$selected_multiplier, 1)
  testthat::expect_null(fixed$diagnostics$fold_losses)
})

testthat::test_that("pooled and SQR use all frozen multipliers with training-only formulas", {
  library_path <- file.path(root, "results", "scratch", "profile_pair_package_20261022", "library")
  if (dir.exists(library_path)) withr::local_libpaths(c(library_path, .libPaths()))
  ex <- frozen_tuning_test_data(n = 35L, seed = 20261105L)
  ex$config$experiment_id <- "P01"
  ex$config$lambda_beta_multipliers <- c(0.25, 0.5, 1, 2)
  for (method in c("POOLED-QR-LASSO", "SQR-DEBIASED-IID")) {
    first <- select_frozen_method_tuning_v2(
      method, ex$data, ex$config$tau, colnames(ex$data$X)[1], ex$config, ex$tuning, 20261106L
    )
    repeated <- select_frozen_method_tuning_v2(
      method, ex$data, ex$config$tau, colnames(ex$data$X)[1], ex$config, ex$tuning, 20261106L
    )
    d <- first$diagnostics
    frozen_tuning_test_results[[method]] <- first
    testthat::expect_equal(nrow(d$candidates), 4L)
    testthat::expect_equal(nrow(d$fold_losses), 20L)
    testthat::expect_equal(d$folds, repeated$diagnostics$folds, tolerance = 0)
    testthat::expect_equal(d$fold_losses, repeated$diagnostics$fold_losses, tolerance = 0)
    testthat::expect_equal(d$candidates, repeated$diagnostics$candidates, tolerance = 0)
    testthat::expect_equal(first$tuning, repeated$tuning, tolerance = 0)
    full <- fit_benchmark_frozen_tuning_v2(
      method, ex$data, ex$config$tau, colnames(ex$data$X)[1], ex$tuning,
      20261106L, context = list(config = ex$config)
    )
    frozen_tuning_test_results[[paste0(method, "_full")]] <- full
    testthat::expect_identical(full$status, "ok")
    testthat::expect_equal(full$tuning_diagnostics$selected_multiplier, d$selected_multiplier)
    if (method == "SQR-DEBIASED-IID") {
      testthat::expect_equal(full$fit_object$lambda, d$selected_multiplier *
        sqrt(ex$config$tau * (1 - ex$config$tau) * log(ncol(ex$data$X)) / length(ex$data$y)))
    } else {
      testthat::expect_equal(unique(full$fit_object$penalty),
                             d$selected_multiplier * ex$tuning$lambda_beta)
    }
    rows <- d$fold_losses
    testthat::expect_true(all(d$candidates$eligible))
    testthat::expect_equal(rows$raw_lambda_anchor,
                           sqrt(log(ncol(ex$data$X)) / rows$n_training))
    if (method == "SQR-DEBIASED-IID") {
      testthat::expect_equal(rows$h_source,
        pmax(sqrt(ex$config$tau * (1 - ex$config$tau)) *
               sqrt(log(ncol(ex$data$X))) * rows$N_training^(-3 / 10), 0.05))
      testthat::expect_equal(rows$lambda_source, rows$multiplier *
        sqrt(ex$config$tau * (1 - ex$config$tau) * log(ncol(ex$data$X)) / rows$N_training))
    } else {
      testthat::expect_equal(rows$lambda_source, rows$multiplier * rows$raw_lambda_anchor)
    }
    for (candidate in seq_len(4L)) {
      r <- rows[rows$candidate == candidate, ]
      testthat::expect_equal(d$candidates$loss[candidate], sum(r$loss_sum) / sum(r$n_validation))
    }
    changed <- ex$data
    validation_ids <- d$roles[[1]]$validation_subjects
    changed$y[changed$cluster_id %in% validation_ids] <-
      changed$y[changed$cluster_id %in% validation_ids] + 7
    counterfactual <- select_frozen_method_tuning_v2(
      method, changed, ex$config$tau, colnames(changed$X)[1], ex$config, ex$tuning, 20261106L
    )
    fold1 <- which(rows$fold == 1L)
    for (j in fold1) {
      testthat::expect_equal(d$training_fits[[j]]$fit$beta,
                             counterfactual$diagnostics$training_fits[[j]]$fit$beta, tolerance = 0)
    }
    testthat::expect_gt(max(abs(rows$loss_mean[fold1] -
      counterfactual$diagnostics$fold_losses$loss_mean[fold1])), 0.1)
  }
})

testthat::test_that("runner retains selected CV records and the complete Doublepen multiplier grid", {
  library_path <- file.path(root, "results", "scratch", "profile_pair_package_20261022", "library")
  if (dir.exists(library_path)) withr::local_libpaths(c(library_path, .libPaths()))
  ex <- frozen_tuning_test_data(n = 35L, seed = 20261111L)
  cfg <- ex$config
  cfg$experiment_id <- "TUNING-RUNNER-REAL"
  cfg$methods <- c("POOLED-QR-LASSO", "SQR-DEBIASED-IID", "DOUBLE-PEN-QLMM")
  cfg$lambda_beta_multipliers <- c(0.25, 0.5, 1, 2)
  result <- run_one_replication_v2(root, cfg, 1L, base_seed = 20261112L,
                                   final = FALSE, run_id = "frozen_tuning_actual_dgp")
  frozen_tuning_test_results$runner <- result
  testthat::expect_equal(nrow(result$metrics), 3L)
  for (method in cfg$methods[1:2]) {
    d <- result$method_diagnostics[[method]]$tuning
    testthat::expect_equal(nrow(d$candidates), 4L)
    testthat::expect_equal(nrow(d$fold_losses), 20L)
    testthat::expect_true(d$selected_multiplier %in% cfg$lambda_beta_multipliers)
    testthat::expect_equal(sum(d$candidates$selected), 1L)
    testthat::expect_identical(result$metrics$status[result$metrics$method_id == method], "ok")
    testthat::expect_equal(result$metrics$lambda_beta[result$metrics$method_id == method],
                           d$full_training$lambda_beta)
    if (method == "SQR-DEBIASED-IID") {
      row <- result$metrics[result$metrics$method_id == method, ]
      testthat::expect_equal(row$lambda_beta, d$selected_multiplier *
        sqrt(cfg$tau * (1 - cfg$tau) * log(cfg$p) / row$total_rows))
      testthat::expect_equal(row$h,
        max(sqrt(cfg$tau * (1 - cfg$tau)) * sqrt(log(cfg$p)) * row$total_rows^(-3 / 10), 0.05))
    }
  }
  double <- result$method_diagnostics[["DOUBLE-PEN-QLMM"]]$tuning
  testthat::expect_equal(double$configured_multipliers, cfg$lambda_beta_multipliers)
  testthat::expect_equal(double$lambda_alpha_multipliers, c(0.5, 1))
  if (!is.null(double$source_block_tuning)) {
    N <- result$metrics$total_rows[result$metrics$method_id == "DOUBLE-PEN-QLMM"]
    testthat::expect_equal(double$source_block_tuning$lambda_beta_grid,
                           cfg$lambda_beta_multipliers * sqrt(log(cfg$p) / N))
    testthat::expect_true(double$selected_multiplier %in% cfg$lambda_beta_multipliers)
    testthat::expect_equal(result$metrics$lambda_beta[result$metrics$method_id == "DOUBLE-PEN-QLMM"],
                           double$source_block_tuning$lambda_beta)
  }
})

testthat::test_that("E07 oracle support retains all multiplier records with zero L1 penalty", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton", JPDMEQR_ASSEMBLY = "sequential"))
  ex <- frozen_tuning_test_data(n = 30L, seed = 20261109L)
  selected <- select_frozen_method_tuning_v2(
    "PROFILE-DQR-TRUE-SUPPORT", ex$data, ex$config$tau, colnames(ex$data$X)[1],
    ex$config, ex$tuning, 20261110L, context = list(active = ex$data$active)
  )
  d <- selected$diagnostics
  frozen_tuning_test_results$oracle_support <- selected
  testthat::expect_equal(nrow(d$fold_losses), 40L)
  testthat::expect_true(all(d$candidates$eligible))
  testthat::expect_equal(length(unique(d$candidates$loss)), 1L)
  testthat::expect_equal(d$selected_multiplier, min(unique(ex$config$lambda_beta_multipliers)))
  testthat::expect_true(all(vapply(d$training_fits, function(value) {
    all(value$fit$lambda_components$coordinate == 0)
  }, logical(1))))
})

testthat::test_that("nonconverged real candidates retain every fold and terminate tuning", {
  ex <- frozen_tuning_test_data(n = 35L, seed = 20261107L)
  ex$config$lambda_beta_multipliers <- c(0.25, 0.5, 1, 2)
  failure <- tryCatch(select_frozen_method_tuning_v2(
    "SQR-DEBIASED-IID", ex$data, ex$config$tau, colnames(ex$data$X)[1],
    ex$config, ex$tuning, 20261108L, control = list(fit_control = list(max_iter = 1L))
  ), error = identity)
  testthat::expect_s3_class(failure, "frozen_tuning_error")
  testthat::expect_identical(failure$failure_stage, "tuning")
  d <- failure$diagnostics$tuning
  frozen_tuning_test_results$failed_candidates <- d
  testthat::expect_equal(nrow(d$fold_losses), 20L)
  testthat::expect_false(any(d$candidates$eligible))
  testthat::expect_true(all(is.na(d$candidates$loss)))
  testthat::expect_true(all(d$fold_losses$status == "failed"))
  testthat::expect_true(all(nzchar(d$fold_losses$failure_message)))
  testthat::expect_true(all(vapply(d$training_fits, function(value) {
    !is.null(value$fit$beta) && !value$fit$converged
  }, logical(1))))
})

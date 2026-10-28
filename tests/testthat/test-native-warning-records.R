warning_root <- normalizePath(v2_require_root(), mustWork = TRUE)
root <- warning_root
source(file.path(root, "scripts", "00_source_v2.R"), local = TRUE)
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, environment())
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"), local = TRUE)

testthat::test_that("real source warnings survive runner diagnostics and CSV with unchanged inference", {
  registry <- read_simulation_registry_v2(file.path(root, "config", "simulation_pair_pilot_20261022.yml"))
  cfg <- config_row_to_list_v2(registry[registry$experiment_id == "P01", , drop = FALSE])
  cfg$experiment_id <- "NATIVE-WARNING-BIAS-REAL"
  cfg$module <- "source_warning_development"
  cfg$n_clusters <- 40L
  cfg$p <- 5L
  cfg$s <- 2L
  cfg$m_values <- 3:5
  cfg$sigma_b0 <- 1.2
  cfg$copula_rho <- 0.4
  cfg$methods <- "BIAS-ADJ-LQMM"
  cfg$target_coordinate_count <- 4L
  cfg$replications <- 1L
  cfg$screen_method <- "none"
  seed_offset <- seed_from_id_v2(1L, cfg$experiment_id, 1L, "training") - 1
  base_seed <- as.integer((202 - seed_offset) %% 2147483000)
  seed <- seed_from_id_v2(base_seed, cfg$experiment_id, 1L, "training")
  testthat::expect_equal(seed, 202L)
  dat0 <- simulate_from_config_v2(cfg, seed)
  known <- generate_profile_qr_data_v2(
    n_clusters = 40L, p = 5L, s = 2L, tau = 0.5, q = 1L,
    m_values = 3:5, m_rule = "uniform", sigma_b0 = 1.2,
    error_dist = "gaussian", seed = 202L
  )
  testthat::expect_identical(dat0[c("y", "X", "Z", "cluster_id")],
                             known[c("y", "X", "Z", "cluster_id")])
  runner_conditions <- list()
  result <- withCallingHandlers(run_one_replication_v2(
    root, cfg, 1L, base_seed = base_seed, final = FALSE,
    run_id = "actual_source_warning_recording"
  ), warning = function(condition) {
    runner_conditions[[length(runner_conditions) + 1L]] <<- condition
  })
  prepared <- prepare_analysis_data_v2(dat0, cfg, seed)
  dat <- prepared$data
  targets <- resolve_target_coordinates_v2(cfg, dat)
  low <- common_lqmm_subset_v2(
    dat, targets, cfg$tau, max_features = 20L,
    seed = seed_from_id_v2(base_seed, cfg$experiment_id, 1L, "common-lowdim-roles"),
    screen_fraction = cfg$screen_fraction
  )
  targets <- intersect(targets, colnames(low$data$X))
  direct_conditions <- list()
  direct <- withCallingHandlers(fit_benchmark_bias_adj_lqmm_v2(
    low$data, cfg$tau, targets, build_tuning_for_analysis_v2(cfg, dat),
    seed = seed_from_id_v2(base_seed, cfg$experiment_id, 1L, "method-BIAS-ADJ-LQMM"),
    control = list(B = 25L, nK = 7L, ci_level = 0.95)
  ), warning = function(condition) {
    direct_conditions[[length(direct_conditions) + 1L]] <<- condition
  })
  diag <- result$method_diagnostics[["BIAS-ADJ-LQMM"]]
  records <- diag$native_warning_records
  messages <- vapply(records, function(record) record$message, character(1))
  direct_messages <- vapply(direct_conditions, conditionMessage, character(1))
  testthat::expect_gt(length(records), 0L)
  testthat::expect_identical(messages, direct_messages)
  testthat::expect_identical(messages, vapply(runner_conditions, conditionMessage, character(1)))
  testthat::expect_true(any(grepl("non-positive fis|variance-covariance|nonunique", messages)))
  testthat::expect_identical(vapply(records, function(record) record$sequence, integer(1)),
                             seq_along(records))
  for (j in seq_along(records)) {
    testthat::expect_identical(records[[j]]$condition_class, class(runner_conditions[[j]]))
    testthat::expect_identical(records[[j]]$call,
      paste(deparse(conditionCall(runner_conditions[[j]])), collapse = " "))
    testthat::expect_identical(records[[j]]$stage, "benchmark_fitting_and_tuning")
  }
  testthat::expect_identical(result$metrics$status, "warning")
  testthat::expect_identical(diag$status, "warning")
  testthat::expect_identical(direct$status, "ok")
  expected_messages <- c(direct$warning_messages, direct_messages)
  testthat::expect_equal(result$metrics$warning_count, length(expected_messages))
  testthat::expect_identical(result$metrics$warning_messages, paste(expected_messages, collapse = " | "))
  testthat::expect_equal(diag$beta_hat, direct$beta_hat, tolerance = 0)
  testthat::expect_identical(diag$training_design$common_screen_roles, low$roles)
  coordinates <- result$coordinates[match(targets, result$coordinates$coordinate), ]
  testthat::expect_equal(coordinates$beta_hat, as.numeric(direct$beta_hat[targets]), tolerance = 0)
  for (column in c("beta_tilde", "ci_lower", "ci_upper")) {
    testthat::expect_equal(coordinates[[column]], as.numeric(direct[[column]][targets]), tolerance = 0)
  }
  testthat::expect_equal(coordinates$estimated_se, as.numeric(direct$se[targets]), tolerance = 0)
  testthat::expect_true(all(is.finite(coordinates$estimated_se) & coordinates$estimated_se > 0))
  output <- getOption("jpDMEQR.unit_artifact_dir")
  if (is.null(output)) {
    output <- file.path(root, "results", "scratch", "native_warning_records_20261022",
                        paste0("test_", Sys.getpid()))
    dir.create(output, recursive = TRUE)
  }
  utils::write.csv(result$metrics, file.path(output, "replication_metrics.csv"), row.names = FALSE)
  utils::write.csv(result$coordinates, file.path(output, "coordinate_metrics.csv"), row.names = FALSE)
  csv <- utils::read.csv(file.path(output, "replication_metrics.csv"), stringsAsFactors = FALSE)
  testthat::expect_equal(csv$warning_count, length(expected_messages))
  testthat::expect_identical(csv$warning_messages, paste(expected_messages, collapse = " | "))
  testthat::expect_identical(csv$status, "warning")
  testthat::expect_identical(names(csv), names(result$metrics))
  saveRDS(list(config = cfg, base_seed = base_seed, raw_data = dat0,
               result = result, direct = direct, runner_conditions = runner_conditions,
               direct_conditions = direct_conditions, training_roles = low$roles),
           file.path(output, "actual_warning_comparison.rds"))
})

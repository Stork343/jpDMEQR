root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
sys.source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"),
           envir = environment())

testthat::test_that("interval helpers build symmetric bounds and 90% companions", {
  z95 <- normal_interval_multiplier_v2(0.95)
  z90 <- normal_interval_multiplier_v2(0.90)
  testthat::expect_equal(z95, stats::qnorm(0.975))
  testthat::expect_equal(z90, stats::qnorm(0.95))
  testthat::expect_true(z95 > z90)

  bounds <- interval_bounds_v2(1, 0.5, 0.95)
  testthat::expect_equal(unname(bounds[["lower"]]), 1 - z95 * 0.5)
  testthat::expect_equal(unname(bounds[["upper"]]), 1 + z95 * 0.5)

  negative <- interval_bounds_v2(-2, 1, 0.95)
  testthat::expect_equal(unname(negative[["lower"]]), -2 - z95)

  unusable <- interval_bounds_v2(NA_real_, 0.25, 0.95)
  testthat::expect_true(is.na(unusable[["lower"]]) && is.na(unusable[["upper"]]))
  negative_se <- interval_bounds_v2(1, -1, 0.95)
  testthat::expect_true(is.na(negative_se[["lower"]]))

  testthat::expect_true(coverage_indicator_v2(
    unname(bounds[["lower"]]), unname(bounds[["upper"]]), 1
  ))
  testthat::expect_false(coverage_indicator_v2(
    unname(bounds[["lower"]]), unname(bounds[["upper"]]), 5
  ))
  testthat::expect_identical(coverage_indicator_v2(NA_real_, 1, 0), NA)

  narrow <- interval_bounds_v2(1, 0.5, 0.90)
  testthat::expect_true(narrow[["upper"]] - narrow[["lower"]] <
                          bounds[["upper"]] - bounds[["lower"]])
})

testthat::test_that("replicate coverage aggregates at the replication level", {
  covered <- c(1, 1, 0, 1, 0, 1, 1, 1)
  replicate <- c(1, 1, 1, 1, 2, 2, 2, 2)
  summary <- replicate_coverage_summary_v2(covered, replicate)
  testthat::expect_identical(unname(summary[["replications"]]), 2)
  testthat::expect_equal(unname(summary[["rate"]]), 0.75)
  testthat::expect_equal(unname(summary[["mcse"]]), 0)

  mixed <- replicate_coverage_summary_v2(c(1, 1, 0, 0), c(1, 1, 2, 2))
  testthat::expect_equal(unname(mixed[["rate"]]), 0.5)
  testthat::expect_equal(unname(mixed[["mcse"]]), 0.5)

  empty <- replicate_coverage_summary_v2(c(NA_real_, NA_real_), c(1, 2))
  testthat::expect_equal(unname(empty[["replications"]]), 0)
  testthat::expect_true(is.na(empty[["rate"]]))
})

testthat::test_that("null p-value bins and the uniform statistic behave", {
  uniform_p <- seq(0.001, 0.999, length.out = 1000)
  bins <- null_pvalue_bins_v2(uniform_p)
  testthat::expect_equal(sum(bins$count), 1000)
  testthat::expect_true(all(abs(bins$expected_fraction - 0.1) < 1e-8))
  testthat::expect_equal(sum(bins$fraction), 1)
  testthat::expect_true(max(abs(bins$fraction - 0.1)) < 0.01)

  ks <- ks_uniform_v2(uniform_p)
  testthat::expect_true(ks[["statistic"]] < 0.01)
  testthat::expect_equal(ks[["n"]], 1000)
  testthat::expect_equal(ks[["mcse_bound"]], 0.5 / sqrt(1000))

  skewed <- rep(0.001, 500)
  testthat::expect_warning(skewed_ks <- ks_uniform_v2(skewed), "ties should not be present")
  testthat::expect_true(skewed_ks[["statistic"]] > 0.5)
})

testthat::test_that("studentised quantiles carry the normal reference", {
  z <- stats::rnorm(2000)
  quantiles <- studentised_quantiles_v2(z, probs = c(0.05, 0.5, 0.95))
  testthat::expect_identical(names(quantiles), c("probability", "empirical", "normal"))
  testthat::expect_equal(quantiles$normal, stats::qnorm(c(0.05, 0.5, 0.95)))
  testthat::expect_true(all(diff(quantiles$empirical) > 0))
  testthat::expect_true(abs(quantiles$empirical[2]) < 0.2)

  empty <- studentised_quantiles_v2(c(NA_real_, NaN), probs = c(0.5))
  testthat::expect_true(is.na(empty$empirical[1]))
})

testthat::test_that("log-log slopes recover a known rate", {
  n <- c(100, 200, 400, 800)
  value <- 5 * n^(-0.5)
  testthat::expect_warning(slope <- loglog_slope_v2(n, value), "essentially perfect fit")
  testthat::expect_equal(unname(slope[["slope"]]), -0.5, tolerance = 1e-8)
  testthat::expect_equal(unname(slope[["points"]]), 4)

  too_few <- loglog_slope_v2(c(100, 200), c(0.1, 0.05))
  testthat::expect_true(is.na(too_few[["slope"]]))
  testthat::expect_identical(unname(too_few[["points"]]), 2)

  nonpositive <- loglog_slope_v2(n, c(0.1, 0, 0.01, 0.005))
  testthat::expect_identical(unname(nonpositive[["points"]]), 3)
})

testthat::test_that("selection stability counts frequencies exactly", {
  selected <- c("1|2|3", "1|2", "1|9", "1|2|3")
  stability <- selection_stability_v2(selected, active_names = c("1", "2", "3"),
                                      null_names = c("9"))
  testthat::expect_equal(unname(stability[["active_selection_frequency_mean"]]), 0.75)
  testthat::expect_equal(unname(stability[["active_selection_frequency_max"]]), 1)
  testthat::expect_equal(unname(stability[["stable_active_fraction"]]), 1)
  testthat::expect_equal(unname(stability[["null_selection_frequency_mean"]]), 0.25)
  testthat::expect_equal(unname(stability[["unstable_null_count"]]), 0)

  none_selected <- selection_stability_v2(rep("", 3), c("1"), c("2"))
  testthat::expect_equal(unname(none_selected[["active_selection_frequency_mean"]]), 0)
  testthat::expect_equal(unname(none_selected[["stable_active_fraction"]]), 0)
  testthat::expect_equal(unname(none_selected[["unstable_null_count"]]), 0)
})

testthat::test_that("cluster loss quantiles and memory snapshots are well formed", {
  y <- c(1, 2, 3, 4)
  qhat <- rep(2, 4)
  cluster_id <- c("a", "a", "b", "b")
  quantiles <- cluster_loss_quantiles_v2(y, qhat, cluster_id, tau = 0.5)
  testthat::expect_identical(names(quantiles), c("q10", "q50", "q90"))
  testthat::expect_true(quantiles[["q10"]] <= quantiles[["q50"]])
  testthat::expect_true(quantiles[["q50"]] <= quantiles[["q90"]])
  testthat::expect_equal(unname(quantiles[["q50"]]), 0.5)

  unusable <- cluster_loss_quantiles_v2(y, rep(NA_real_, 4), cluster_id, tau = 0.5)
  testthat::expect_true(all(is.na(unusable)))

  snapshot <- memory_snapshot_v2()
  testthat::expect_true(all(c("rss_mb", "peak_wset_mb") %in% names(snapshot)))
  testthat::expect_true(is.numeric(snapshot[["rss_mb"]]))
  testthat::expect_true(is.numeric(snapshot[["peak_wset_mb"]]))
})

testthat::test_that("safe aggregation helpers keep NA instead of zero", {
  all_na <- c(NA_real_, NA_real_)
  testthat::expect_true(is.na(mean_or_na_v2(all_na)))
  testthat::expect_true(is.na(median_or_na_v2(all_na)))
  testthat::expect_true(is.na(max_or_na_v2(all_na)))
  testthat::expect_true(is.na(mean_or_na_v2(numeric())))
  testthat::expect_equal(mean_or_na_v2(c(1, 2, NA)), 1.5)
  testthat::expect_equal(max_or_na_v2(c(1, 2, NA)), 2)
})

testthat::test_that("the named gate manifest requires the final flag", {
  dir <- file.path(tempdir(), "jpDMEQR_named_gate_final")
  if (dir.exists(dir)) unlink(dir, recursive = TRUE)
  dir.create(file.path(dir, "results", "preflight"), recursive = TRUE)
  commit <- current_commit_v2(root)
  writeLines(commit, file.path(dir, "implementation_commit.txt"))
  path <- micro_preflight_manifest_path_v2(dir)

  write_json_v2(list(pass = TRUE, commit_sha = commit, final = FALSE), path)
  without_final <- validate_named_gate_manifest_v2(dir, path, "micro-preflight",
                                                   require_final = TRUE)
  testthat::expect_false(without_final$pass)
  testthat::expect_match(without_final$reason, "final flag")
  testthat::expect_true(validate_named_gate_manifest_v2(dir, path, "micro-preflight")$pass)

  write_json_v2(list(pass = TRUE, commit_sha = commit, final = TRUE), path)
  with_final <- validate_named_gate_manifest_v2(dir, path, "micro-preflight",
                                                require_final = TRUE)
  testthat::expect_true(with_final$pass)
  testthat::expect_identical(with_final$reason, "accepted")

  write_json_v2(list(pass = TRUE, commit_sha = "another-commit", final = TRUE), path)
  stale <- validate_named_gate_manifest_v2(dir, path, "micro-preflight",
                                           require_final = TRUE)
  testthat::expect_false(stale$pass)
})

testthat::test_that("the output contract documents every required column", {
  contract <- readLines(file.path(root, "docs", "RESULTS_CONTRACT.md"), warn = FALSE)
  contract <- paste(contract, collapse = "\n")
  required <- c(
    "experiment_id", "method_id", "replicate", "status", "converged",
    "n_clusters", "n_clusters_original", "p", "p_original",
    "l1_error", "l2_error", "pinball_loss", "pinball_loss_intercept",
    "pinball_skill_pooled", "pinball_skill_intercept", "quantile_calibration",
    "cluster_loss_q10", "cluster_loss_q50", "cluster_loss_q90",
    "selected_size", "tpr", "fdp", "selected_names", "objective",
    "runtime_sec", "first_stage_sec", "refit_sec", "inference_sec",
    "peak_memory_mb", "process_peak_wset_mb", "outer_iterations",
    "nuisance_iterations", "dantzig_feasible_fraction", "warning_count"
  )
  missing <- required[!vapply(required, grepl, logical(1), x = contract, fixed = TRUE)]
  testthat::expect_identical(missing, character())

  coordinate_required <- c(
    "experiment_id", "method_id", "replicate", "coordinate", "coordinate_type",
    "beta_tilde", "target_value", "estimated_se", "p_value",
    "ci_lower", "ci_upper", "covered",
    "ci_level_secondary", "ci_lower_secondary", "ci_upper_secondary",
    "covered_secondary", "interval_length_secondary"
  )
  coordinate_missing <- coordinate_required[
    !vapply(coordinate_required, grepl, logical(1), x = contract, fixed = TRUE)
  ]
  testthat::expect_identical(coordinate_missing, character())

  summaries <- c(
    "method_summary.csv", "coverage_summary.csv", "coverage_by_scenario.csv",
    "inference_summary.csv", "null_pvalue_uniformity.csv",
    "studentised_quantiles.csv", "selection_stability.csv",
    "rate_slopes.csv", "computation_summary.csv",
    "paired_method_differences.csv", "registry_accounting.csv"
  )
  summary_missing <- summaries[
    !vapply(summaries, grepl, logical(1), x = contract, fixed = TRUE)
  ]
  testthat::expect_identical(summary_missing, character())
})

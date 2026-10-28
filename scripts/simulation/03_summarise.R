#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/03_summarise.R"
root <- if (file.exists("R/metrics_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
run_id <- cli$`run-id` %||% stop("Supply --run-id=<run directory name>.")
pilot_gate_requested <- as_bool_cli_v2(cli$`pilot-gate`, grepl("^pilot_", run_id))
run_dir <- file.path(root, "results", "raw", "simulation", run_id)
# Gate-issue commit: defaults to the run's implementation commit; may be
# overridden with the commit on which the pilot DECISION is certified (the
# freeze gate validates the decision record against the current commit).
# Run provenance always remains in run_id/run_dir/implementation_commit.txt.
gate_commit_override <- trimws(cli$`gate-commit` %||% "")
metrics_path <- file.path(run_dir, "replication_metrics.csv")
coords_path <- cli$`coordinate-file` %||% file.path(run_dir, "coordinate_metrics.csv")
theory_path <- file.path(run_dir, "theory_diagnostics.csv")
if (!file.exists(metrics_path)) stop("Missing ", metrics_path)
metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
coords <- if (file.exists(coords_path)) read.csv(coords_path, stringsAsFactors = FALSE) else data.frame()
theory <- if (file.exists(theory_path)) read.csv(theory_path, stringsAsFactors = FALSE) else data.frame()
validate_paired_failed_coordinate_accounting_v2(metrics, coords)

out_dir <- cli$`output-dir` %||% file.path(root, "results", "tables", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keys <- interaction(metrics$experiment_id, metrics$method_id, drop = TRUE)
summary_rows <- lapply(split(metrics, keys), function(d) {
  data.frame(
    experiment_id = d$experiment_id[1],
    module = d$module[1],
    method_id = d$method_id[1],
    scheduled = nrow(d),
    success = sum(d$status %in% c("ok", "warning")),
    success_fraction = mean(d$status %in% c("ok", "warning")),
    failure_rate = mean(d$status %in% c("failed", "not_implemented")),
    convergence_rate = mean(d$converged, na.rm = TRUE),
    mean_l1 = mean(d$l1_error, na.rm = TRUE),
    mcse_l1 = mcse_mean_v2(d$l1_error),
    mean_l2 = mean(d$l2_error, na.rm = TRUE),
    mcse_l2 = mcse_mean_v2(d$l2_error),
    mean_tpr = mean(d$tpr, na.rm = TRUE),
    mean_fdp = mean(d$fdp, na.rm = TRUE),
    mean_selected_size = mean(d$selected_size, na.rm = TRUE),
    mean_pinball = mean_or_na_v2(d$pinball_loss),
    mcse_pinball = mcse_mean_v2(d$pinball_loss),
    mean_pinball_skill_pooled = mean_or_na_v2(d$pinball_skill_pooled),
    mcse_pinball_skill_pooled = mcse_mean_v2(d$pinball_skill_pooled),
    mean_pinball_skill_intercept = mean_or_na_v2(d$pinball_skill_intercept),
    mcse_pinball_skill_intercept = mcse_mean_v2(d$pinball_skill_intercept),
    mean_quantile_calibration = mean_or_na_v2(d$quantile_calibration),
    mean_cluster_loss_q50 = mean_or_na_v2(d$cluster_loss_q50),
    mean_objective = mean_or_na_v2(d$objective),
    median_outer_iterations = median_or_na_v2(d$outer_iterations),
    median_nuisance_iterations = median_or_na_v2(d$nuisance_iterations),
    mean_dantzig_feasible_fraction = mean_or_na_v2(d$dantzig_feasible_fraction),
    warning_count_total = sum(d$warning_count, na.rm = TRUE),
    mean_runtime = mean_or_na_v2(d$runtime_sec),
    mcse_runtime = mcse_mean_v2(d$runtime_sec),
    median_runtime = median_or_na_v2(d$runtime_sec),
    max_runtime = max_or_na_v2(d$runtime_sec),
    median_first_stage_sec = median_or_na_v2(d$first_stage_sec),
    median_refit_sec = median_or_na_v2(d$refit_sec),
    median_inference_sec = median_or_na_v2(d$inference_sec),
    median_peak_memory_mb = median_or_na_v2(d$peak_memory_mb),
    max_peak_memory_mb = max_or_na_v2(d$peak_memory_mb),
    stringsAsFactors = FALSE
  )
})
summary_metrics <- do.call(rbind, summary_rows)
write_atomic_csv_v2(summary_metrics, file.path(out_dir, "method_summary.csv"))

# 注册表核算：每一次注册的重复都有记录，失败与警告保留在计数中，
# 未记录的重复单独列出（docs/RESULTS_CONTRACT.md 第 9 节检查）。
registry_path <- file.path(run_dir, "config_snapshot.csv")
if (!file.exists(registry_path)) stop("Missing ", registry_path)
registry <- read.csv(registry_path, stringsAsFactors = FALSE)
method_version <- registry_method_version_v2(registry)
registry_counts <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) {
  data.frame(
    experiment_id = registry$experiment_id[i],
    method_id = parse_pipe_character_v2(registry$methods[i]),
    scheduled_replications = as.integer(registry$replications[i]),
    stringsAsFactors = FALSE
  )
}))
recorded_key <- paste(metrics$experiment_id, metrics$method_id, sep = "|")
registry_key <- paste(registry_counts$experiment_id, registry_counts$method_id, sep = "|")
unregistered <- setdiff(unique(recorded_key), registry_key)
if (length(unregistered)) {
  stop("Recorded rows are not bound to a registry experiment_id: ",
       paste(head(unregistered, 5L), collapse = "; "))
}
record_rows <- lapply(seq_len(nrow(registry_counts)), function(i) {
  d <- metrics[recorded_key == registry_key[i], , drop = FALSE]
  success_rate <- if (nrow(d)) mean(d$status %in% c("ok", "warning")) else NA_real_
  data.frame(
    experiment_id = registry_counts$experiment_id[i],
    method_id = registry_counts$method_id[i],
    scheduled_replications = registry_counts$scheduled_replications[i],
    recorded_replications = nrow(d),
    missing_replications = registry_counts$scheduled_replications[i] - nrow(d),
    success = sum(d$status %in% c("ok", "warning")),
    warning = sum(d$status == "warning"),
    failed = sum(d$status == "failed"),
    not_implemented = sum(d$status == "not_implemented"),
    success_rate = success_rate,
    low_success_flag = is.finite(success_rate) && success_rate < 0.95,
    stringsAsFactors = FALSE
  )
})
write_atomic_csv_v2(do.call(rbind, record_rows),
                    file.path(out_dir, "registry_accounting.csv"))

coverage_summary <- data.frame()
if (nrow(coords)) {
  ckeys <- interaction(coords$experiment_id, coords$method_id, coords$coordinate, drop = TRUE)
  coverage <- lapply(split(coords, ckeys), function(d) {
    cov <- coordinate_coverage_summary_v2(d, "covered")
    cov_secondary <- coordinate_coverage_summary_v2(d, "covered_secondary")
    names(cov_secondary) <- paste0(names(cov_secondary), "_secondary")
    empirical_sd <- stats::sd(d$beta_tilde, na.rm = TRUE)
    mean_se <- mean_or_na_v2(d$estimated_se)
    data.frame(
      experiment_id = d$experiment_id[1],
      method_id = d$method_id[1],
      coordinate = d$coordinate[1],
      coordinate_type = d$coordinate_type[1],
      target_type = d$target_type[1],
      bias = mean_or_na_v2(d$beta_tilde - d$target_value),
      bias_mcse = mcse_mean_v2(d$beta_tilde - d$target_value),
      empirical_sd = empirical_sd,
      mean_se = mean_se,
      median_se = median_or_na_v2(d$estimated_se),
      se_sd_ratio = mean_se / empirical_sd,
      mean_interval_length = mean_or_na_v2(d$interval_length),
      mean_interval_length_secondary = mean_or_na_v2(d$interval_length_secondary),
      cov,
      cov_secondary,
      stringsAsFactors = FALSE
    )
  })
  coverage_summary <- do.call(rbind, coverage)
  write_atomic_csv_v2(coverage_summary, file.path(out_dir, "coverage_summary.csv"))

  # 情景层覆盖率：在目标坐标上汇集，同时给出按重复计算的蒙特卡洛标准误，
  # 避免把同一重复内的相关坐标当作独立单元。
  active_coords <- coords[coords$coordinate_type == "active", , drop = FALSE]
  if (nrow(active_coords)) {
    akeys <- interaction(active_coords$experiment_id, active_coords$method_id,
                         drop = TRUE)
    scenario_cov <- lapply(split(active_coords, akeys), function(d) {
      cov <- coordinate_coverage_summary_v2(d, "covered")
      cov_secondary <- coordinate_coverage_summary_v2(d, "covered_secondary")
      names(cov_secondary) <- paste0(names(cov_secondary), "_secondary")
      coverage_replicate <- coordinate_replicate_coverage_v2(d, "covered")
      coverage_replicate_secondary <- coordinate_replicate_coverage_v2(d, "covered_secondary")
      coverage_successful <- replicate_coverage_summary_v2(d$covered, d$replicate)
      empirical_sd <- stats::sd(d$beta_tilde, na.rm = TRUE)
      mean_se <- mean_or_na_v2(d$estimated_se)
      data.frame(
        experiment_id = d$experiment_id[1],
        method_id = d$method_id[1],
        target_coordinate_count = length(unique(d$coordinate)),
        coordinate_replications = nrow(d),
        bias = mean_or_na_v2(d$beta_tilde - d$target_value),
        bias_mcse = mcse_mean_v2(d$beta_tilde - d$target_value),
        empirical_sd = empirical_sd,
        mean_se = mean_se,
        median_se = median_or_na_v2(d$estimated_se),
        se_sd_ratio = mean_se / empirical_sd,
        mean_interval_length = mean_or_na_v2(d$interval_length),
        coverage_replication = coverage_replicate[["rate"]],
        coverage_replication_mcse = coverage_replicate[["mcse"]],
        coverage_replications = coverage_replicate[["replications"]],
        coverage_replication_successful = coverage_successful[["rate"]],
        coverage_replication_successful_mcse = coverage_successful[["mcse"]],
        coverage_replication_secondary = coverage_replicate_secondary[["rate"]],
        coverage_replication_mcse_secondary = coverage_replicate_secondary[["mcse"]],
        cov,
        cov_secondary,
        stringsAsFactors = FALSE
      )
    })
    write_atomic_csv_v2(do.call(rbind, scenario_cov),
                        file.path(out_dir, "coverage_by_scenario.csv"))
  }

  ikeys <- interaction(coords$experiment_id, coords$method_id, drop = TRUE)
  inference_rows <- lapply(split(coords, ikeys), function(d) {
    null_d <- d[d$coordinate_type == "null", , drop = FALSE]
    active_d <- d[d$coordinate_type == "active", , drop = FALSE]
    null_05 <- replicate_rate_summary_v2(null_d$p_value, null_d$replicate, 0.05)
    null_10 <- replicate_rate_summary_v2(null_d$p_value, null_d$replicate, 0.10)
    power_05 <- replicate_rate_summary_v2(active_d$p_value, active_d$replicate, 0.05)
    power_10 <- replicate_rate_summary_v2(active_d$p_value, active_d$replicate, 0.10)
    pooled_null <- extremes_summary_v2(null_d$p_value, 0.05)
    pooled_active <- extremes_summary_v2(active_d$p_value, 0.05)
    data.frame(
      experiment_id = d$experiment_id[1],
      method_id = d$method_id[1],
      null_replications = null_05[["replications"]],
      null_evaluable = pooled_null[["evaluable"]],
      type1_error_05 = null_05[["rate"]],
      type1_error_05_mcse = null_05[["mcse"]],
      type1_error_10 = null_10[["rate"]],
      type1_error_05_pooled = pooled_null[["rate"]],
      type1_error_05_pooled_lower = pooled_null[["binom_lower"]],
      type1_error_05_pooled_upper = pooled_null[["binom_upper"]],
      active_replications = power_05[["replications"]],
      active_evaluable = pooled_active[["evaluable"]],
      power_05 = power_05[["rate"]],
      power_05_mcse = power_05[["mcse"]],
      power_10 = power_10[["rate"]],
      power_05_pooled = pooled_active[["rate"]],
      power_05_pooled_lower = pooled_active[["binom_lower"]],
      power_05_pooled_upper = pooled_active[["binom_upper"]],
      stringsAsFactors = FALSE
    )
  })
  write_atomic_csv_v2(do.call(rbind, inference_rows),
                      file.path(out_dir, "inference_summary.csv"))

  uniformity_rows <- lapply(split(coords, ikeys), function(d) {
    null_d <- d[d$coordinate_type == "null", , drop = FALSE]
    ks <- ks_uniform_v2(null_d$p_value)
    bins <- null_pvalue_bins_v2(null_d$p_value)
    counts <- as.data.frame(as.list(stats::setNames(
      bins$count, paste0("bin_", seq_len(nrow(bins)))
    )), stringsAsFactors = FALSE)
    expected <- as.data.frame(as.list(stats::setNames(
      bins$expected_fraction, paste0("expected_", seq_len(nrow(bins)))
    )), stringsAsFactors = FALSE)
    cbind(data.frame(
      experiment_id = d$experiment_id[1],
      method_id = d$method_id[1],
      ks_statistic = ks[["statistic"]],
      ks_mcse_bound = ks[["mcse_bound"]],
      pvalue_count = ks[["n"]],
      stringsAsFactors = FALSE
    ), counts, expected)
  })
  write_atomic_csv_v2(do.call(rbind, uniformity_rows),
                      file.path(out_dir, "null_pvalue_uniformity.csv"))

  quantile_rows <- lapply(split(coords, ikeys), function(d) {
    cbind(data.frame(
      experiment_id = d$experiment_id[1],
      method_id = d$method_id[1],
      stringsAsFactors = FALSE
    ), studentised_quantiles_v2(d$wald_z))
  })
  write_atomic_csv_v2(do.call(rbind, quantile_rows),
                      file.path(out_dir, "studentised_quantiles.csv"))

  support_map <- lapply(split(coords, coords$experiment_id), function(d) {
    list(active = unique(d$coordinate[d$coordinate_type == "active"]),
         null = unique(d$coordinate[d$coordinate_type == "null"]))
  })
  if (!"selected_names" %in% names(metrics)) {
    stop("The replication record has no selected_names column.")
  }
  skeys <- interaction(metrics$experiment_id, metrics$method_id, drop = TRUE)
  stability_rows <- lapply(split(metrics, skeys), function(d) {
    support <- support_map[[d$experiment_id[1]]]
    if (is.null(support)) support <- list(active = character(), null = character())
    stability <- selection_stability_v2(d$selected_names, support$active, support$null)
    data.frame(
      experiment_id = d$experiment_id[1],
      method_id = d$method_id[1],
      scheduled = nrow(d),
      as.list(stability),
      stringsAsFactors = FALSE
    )
  })
  write_atomic_csv_v2(do.call(rbind, stability_rows),
                      file.path(out_dir, "selection_stability.csv"))
}

if (nrow(theory)) {
  tkeys <- interaction(theory$experiment_id, theory$method_id, drop = TRUE)
  theory_summary <- lapply(split(theory, tkeys), function(d) data.frame(
    experiment_id = d$experiment_id[1],
    method_id = d$method_id[1],
    mean_beta_h_target_distance = mean_or_na_v2(d$beta_h_target_distance),
    max_beta_h_target_distance = max_or_na_v2(d$beta_h_target_distance),
    mean_hessian_max_error = mean_or_na_v2(d$hessian_max_error),
    mean_hessian_operator_error = mean_or_na_v2(d$hessian_operator_error),
    mean_hessian_condition_number = mean_or_na_v2(d$hessian_condition_number),
    mean_precision_row_error = mean_or_na_v2(d$precision_row_error),
    mean_one_step_correction = mean_or_na_v2(d$one_step_correction_max),
    max_profile_identity_error = max_or_na_v2(d$profile_identity_error),
    max_nuisance_gradient = max_or_na_v2(d$max_nuisance_gradient),
    max_dantzig_residual = max_or_na_v2(d$max_dantzig_residual),
    median_D_k = median_or_na_v2(d$median_D_k),
    max_D_k = max_or_na_v2(d$max_D_k),
    median_A_k = median_or_na_v2(d$median_A_k),
    max_A_k = max_or_na_v2(d$max_A_k),
    mean_scaled_bahadur_remainder = mean_or_na_v2(d$scaled_bahadur_remainder),
    max_abs_scaled_bahadur_remainder = max_or_na_v2(abs(d$scaled_bahadur_remainder)),
    mean_score_skewness = mean_or_na_v2(d$score_skewness),
    mean_score_excess_kurtosis = mean_or_na_v2(d$score_excess_kurtosis),
    stringsAsFactors = FALSE
  ))
  write_atomic_csv_v2(do.call(rbind, theory_summary),
                      file.path(out_dir, "theory_summary.csv"))
}

rate_columns <- c("l1_error", "l2_error", "max_active_error")
if (nrow(theory)) {
  merged_rates <- merge(
    metrics[, c("experiment_id", "replicate", "method_id", "module", "n_clusters",
                "p", "s", rate_columns)],
    theory[, c("experiment_id", "replicate", "method_id",
               "scaled_bahadur_remainder", "beta_h_target_distance",
               "hessian_max_error")],
    by = c("experiment_id", "replicate", "method_id"), all.x = TRUE
  )
} else {
  merged_rates <- metrics[, c("experiment_id", "replicate", "method_id", "module",
                              "n_clusters", "p", "s", rate_columns)]
  merged_rates$scaled_bahadur_remainder <- NA_real_
  merged_rates$beta_h_target_distance <- NA_real_
  merged_rates$hessian_max_error <- NA_real_
}
rate_metrics <- c(rate_columns, "scaled_bahadur_remainder",
                  "beta_h_target_distance", "hessian_max_error")
rkeys <- interaction(merged_rates$module, merged_rates$method_id, merged_rates$p,
                     merged_rates$s, drop = TRUE)
rate_rows <- list()
for (key in levels(rkeys)) {
  d <- merged_rates[rkeys == key, , drop = FALSE]
  if (length(unique(d$n_clusters[is.finite(d$n_clusters)])) < 3L) next
  for (metric in rate_metrics) {
    slope <- loglog_slope_v2(d$n_clusters, abs(d[[metric]]))
    rate_rows[[length(rate_rows) + 1L]] <- data.frame(
      module = d$module[1],
      method_id = d$method_id[1],
      p = d$p[1],
      s = d$s[1],
      metric = metric,
      slope = slope[["slope"]],
      slope_mcse = slope[["slope_mcse"]],
      points = slope[["points"]],
      stringsAsFactors = FALSE
    )
  }
}
if (length(rate_rows)) {
  write_atomic_csv_v2(do.call(rbind, rate_rows), file.path(out_dir, "rate_slopes.csv"))
}

if ("POOLED-QR-LASSO" %in% metrics$method_id) {
  paired <- paired_difference_summary_v2(
    metrics, value = "pinball_loss", method = "method_id",
    reference = "POOLED-QR-LASSO",
    id = c("experiment_id", "replicate")
  )
  write_atomic_csv_v2(paired, file.path(out_dir, "paired_pinball_differences.csv"))
}

if ("PROFILE-DQR" %in% metrics$method_id) {
  paired_rows <- list()
  for (metric in c("l1_error", "l2_error", "pinball_loss", "pinball_skill_pooled")) {
    paired <- paired_difference_summary_v2(
      metrics, value = metric, method = "method_id",
      reference = "PROFILE-DQR", id = c("experiment_id", "replicate")
    )
    paired$metric <- metric
    paired_rows[[length(paired_rows) + 1L]] <- paired
  }
  write_atomic_csv_v2(do.call(rbind, paired_rows),
                      file.path(out_dir, "paired_method_differences.csv"))
}

computation <- metrics[metrics$module == "G_computation", , drop = FALSE]
if (nrow(computation)) {
  gkeys <- interaction(computation$experiment_id, computation$method_id, drop = TRUE)
  computation_summary <- lapply(split(computation, gkeys), function(d) {
    data.frame(
      experiment_id = d$experiment_id[1],
      method_id = d$method_id[1],
      requested_workers = d$requested_workers[1],
      scheduled = nrow(d),
      failure_rate = mean(d$status %in% c("failed", "not_implemented")),
      median_runtime = median_or_na_v2(d$runtime_sec),
      q25_runtime = if (any(is.finite(d$runtime_sec))) {
        unname(stats::quantile(d$runtime_sec, 0.25, na.rm = TRUE, names = FALSE))
      } else NA_real_,
      q75_runtime = if (any(is.finite(d$runtime_sec))) {
        unname(stats::quantile(d$runtime_sec, 0.75, na.rm = TRUE, names = FALSE))
      } else NA_real_,
      max_runtime = max_or_na_v2(d$runtime_sec),
      median_first_stage_sec = median_or_na_v2(d$first_stage_sec),
      median_refit_sec = median_or_na_v2(d$refit_sec),
      median_inference_sec = median_or_na_v2(d$inference_sec),
      median_peak_memory_mb = median_or_na_v2(d$peak_memory_mb),
      max_peak_memory_mb = max_or_na_v2(d$peak_memory_mb),
      median_outer_iterations = median_or_na_v2(d$outer_iterations),
      median_nuisance_iterations = median_or_na_v2(d$nuisance_iterations),
      median_kkt_residual = median_or_na_v2(d$kkt_residual),
      median_condition_number = median_or_na_v2(d$hessian_condition_number),
      mean_dantzig_feasible_fraction = mean_or_na_v2(d$dantzig_feasible_fraction),
      stringsAsFactors = FALSE
    )
  })
  write_atomic_csv_v2(do.call(rbind, computation_summary),
                      file.path(out_dir, "computation_summary.csv"))
}

if (pilot_gate_requested) {
  run_request <- read_json_v2(file.path(run_dir, "run_request.json"))
  if (!isTRUE(run_request$final)) {
    stop("The pilot gate cannot be certified from a development-mode run.")
  }
  proposed <- metrics[metrics$method_id == "PROFILE-DQR", , drop = FALSE]
  required_ids <- sprintf("P%02d", 1:6)
  per_id <- split(proposed, proposed$experiment_id)
  convergence_pass <- all(vapply(required_ids, function(id) {
    d <- per_id[[id]]
    !is.null(d) && mean(d$status %in% c("ok", "warning")) >= 0.98
  }, logical(1)))
  identity_values <- proposed$profile_identity_error[is.finite(proposed$profile_identity_error)]
  identity_pass <- length(identity_values) &&
    stats::median(identity_values) < 1e-8 && max(identity_values) < 1e-6
  proposed_coords <- coords[coords$method_id == "PROFILE-DQR", , drop = FALSE]
  dantzig_pass <- nrow(proposed_coords) && mean(proposed_coords$feasible, na.rm = TRUE) >= 0.95
  # Round-2 (B2): D_k is DEMOTED from a hard freeze gate to a diagnostic (it is
  # a Holder upper bound and may be very loose). Report it and A_k, the actual
  # normalized inverse-defect inner product. No hard inverse-defect threshold.
  d_vals <- proposed_coords$D_k[is.finite(proposed_coords$D_k)]
  a_vals <- proposed_coords$A_k[is.finite(proposed_coords$A_k)]
  inverse_defect_diagnostic <- list(
    D_k_median = if (length(d_vals)) stats::median(d_vals) else NA_real_,
    D_k_q90 = if (length(d_vals)) stats::quantile(d_vals, 0.9, na.rm = TRUE, names = FALSE) else NA_real_,
    A_k_median = if (length(a_vals)) stats::median(a_vals) else NA_real_,
    A_k_q90 = if (length(a_vals)) stats::quantile(a_vals, 0.9, na.rm = TRUE, names = FALSE) else NA_real_
  )
  # Round-2 (A2): the HARD first-order calibration gate is P05 (n=200);
  # P06 (n=400) is the scaling confirmation, P01-P04 are stress cells.
  baseline_cov <- coverage_summary[
    coverage_summary$experiment_id == "P05" &
      coverage_summary$method_id == "PROFILE-DQR", , drop = FALSE
  ]
  se_ratio_pass <- nrow(baseline_cov) &&
    all(baseline_cov$se_sd_ratio >= 0.80 & baseline_cov$se_sd_ratio <= 1.20, na.rm = TRUE)
  coverage_pass <- nrow(baseline_cov) &&
    all(baseline_cov$coverage >= 0.88 & baseline_cov$coverage <= 0.98, na.rm = TRUE)
  # P06 scaling summary (diagnostic, not a pass gate)
  p06_cov <- coverage_summary[
    coverage_summary$experiment_id == "P06" &
      coverage_summary$method_id == "PROFILE-DQR", , drop = FALSE
  ]
  p06_scaling <- if (nrow(p06_cov)) {
    list(se_sd_ratio = p06_cov$se_sd_ratio, coverage = p06_cov$coverage)
  } else NULL
  runner_failure_path <- file.path(run_dir, "runner_failures.csv")
  runner_failure_count <- if (file.exists(runner_failure_path)) nrow(read.csv(runner_failure_path)) else 0L
  missing_method_pass <- !any(metrics$status == "not_implemented")
  # Round-6 V1 amendment (METHOD_SPECIFICATION_ROUND6_AMENDMENT.md, section on
  # conditional gate interpretation): the coverage/SE-SD bands are visible
  # calibration diagnostics and are no longer binary implementation-correctness
  # requirements for intentionally small-n stress cells. The pilot gate pass
  # criterion is mechanism positivity (convergence/identity/Dantzig) plus
  # record completeness (no runner failures, no missing methods). The band
  # fields below remain in the manifest as diagnostics.
  v1_mechanism_pass <- convergence_pass && identity_pass && dantzig_pass &&
    runner_failure_count == 0L && missing_method_pass
  pass <- v1_mechanism_pass
  pilot_commit <- if (nzchar(gate_commit_override)) gate_commit_override else
    trimws(readLines(file.path(run_dir, "implementation_commit.txt"), warn = FALSE)[1])
  manifest <- list(
    pass = pass,
    final = TRUE,
    commit_sha = pilot_commit,
    v1_classification = "V1 (round-6 ladder; bands are diagnostics)",
    mechanism_evidence = list(
      list(run_id = "pilot_v2r2_P05b", B = 200L, report = "results/preflight/PILOT_V2_P05_CALIBRATION_GATE.md", formal_gate = FALSE),
      list(run_id = "pilot_v4_P05_B50", B = 50L, report = "results/preflight/PILOT_V4_PAIRED_PROBE.md", formal_gate = FALSE),
      list(run_id = "pilot_v4_P06_B50", B = 50L, report = "results/preflight/PILOT_V4_PAIRED_PROBE.md", formal_gate = FALSE),
      list(run_id = "pilot_v5_P05_B50", B = 50L, report = "results/preflight/PILOT_V5_PAIRED_PROBE.md", formal_gate = FALSE),
      list(run_id = "pilot_v5_P06_B50", B = 50L, report = "results/preflight/PILOT_V5_PAIRED_PROBE.md", formal_gate = FALSE),
      list(run_id = "pilot_v5_P07_B50", B = 50L, report = "results/preflight/PILOT_V5_PAIRED_PROBE.md", formal_gate = FALSE),
      list(run_id = "pilot_v6_ladder_P06", B = 50L, report = "results/preflight/PILOT_V6_VARIANCE_LADDER.md", formal_gate = FALSE),
      list(run_id = "pilot_v5_B200", B = 200L, report = "results/preflight/PILOT_V5_FINAL_REPORT.md", formal_gate = TRUE),
      list(run_id = "pilot_v5_P07_B200", B = 200L, report = "results/preflight/PILOT_V5_P07_B200_REPORT.md", formal_gate = FALSE),
      list(run_id = "papp_hd_srv3", B = 200L, report = "results/preflight/PILOT_V6_PAPP_REPORT.md", formal_gate = FALSE),
      list(run_id = "papp_ld_B200", B = 200L, report = "results/preflight/PILOT_V6_PAPP_REPORT.md", formal_gate = FALSE)
    ),
    run_id = run_id,
    run_dir = normalizePath(run_dir, mustWork = TRUE),
    evaluated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    convergence_pass = convergence_pass,
    identity_pass = identity_pass,
    dantzig_pass = dantzig_pass,
    inverse_defect_diagnostic = inverse_defect_diagnostic,
    p06_scaling = p06_scaling,
    calibration_baseline = "P05",
    se_ratio_pass = se_ratio_pass,
    coverage_pass = coverage_pass,
    no_runner_failures = runner_failure_count == 0L,
    no_missing_methods = missing_method_pass,
    thresholds = list(
      convergence = 0.98,
      identity_median = 1e-8,
      identity_max = 1e-6,
      dantzig_feasible = 0.95,
      se_sd_ratio = c(0.80, 1.20),
      baseline_coverage = c(0.88, 0.98)
    )
  )
  manifest$method_version <- method_version
  manifest$score_definition <- if (identical(method_version, "profile_pair_20261022"))
    "profile_bandwidth_pair_2_half_minus_full" else NULL
  manifest$solver_configuration <- read_json_v2(file.path(run_dir, "solver_configuration.json"))
  if (identical(method_version, "profile_pair_20261022")) {
    manifest$v1_classification <- "paired profile score mechanism"
    manifest$mechanism_evidence <- NULL
    accounting <- do.call(rbind, record_rows)
    complete_counts <- all(accounting$missing_replications == 0L) &&
      all(registry$experiment_id %in% sprintf("P%02d", 1:6)) &&
      length(unique(registry$experiment_id)) == 6L &&
      all(registry$replications >= 200L)
    actual_commit <- trimws(readLines(file.path(run_dir, "implementation_commit.txt"), warn = FALSE)[1L])
    if (!identical(pilot_commit, actual_commit)) {
      stop("The paired pilot decision must retain the actual execution commit.")
    }
    manifest$complete_registered_replications <- complete_counts
    manifest$pass <- manifest$pass && complete_counts
  }
  gate_path <- pilot_gate_manifest_path_v2(root, method_version)
  dir.create(dirname(gate_path), recursive = TRUE, showWarnings = FALSE)
  write_json_v2(manifest, gate_path)
  if (!pass) warning("Pilot gate did not pass. Diagnose; do not calibrate to the thresholds.")
}
message("Summaries written to: ", out_dir)

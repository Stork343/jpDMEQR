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
metrics_path <- file.path(run_dir, "replication_metrics.csv")
coords_path <- file.path(run_dir, "coordinate_metrics.csv")
theory_path <- file.path(run_dir, "theory_diagnostics.csv")
if (!file.exists(metrics_path)) stop("Missing ", metrics_path)
metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
coords <- if (file.exists(coords_path)) read.csv(coords_path, stringsAsFactors = FALSE) else data.frame()
theory <- if (file.exists(theory_path)) read.csv(theory_path, stringsAsFactors = FALSE) else data.frame()

out_dir <- file.path(root, "results", "tables", run_id)
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
    mean_pinball = mean(d$pinball_loss, na.rm = TRUE),
    mcse_pinball = mcse_mean_v2(d$pinball_loss),
    mean_runtime = mean(d$runtime_sec, na.rm = TRUE),
    mcse_runtime = mcse_mean_v2(d$runtime_sec),
    stringsAsFactors = FALSE
  )
})
summary_metrics <- do.call(rbind, summary_rows)
write_atomic_csv_v2(summary_metrics, file.path(out_dir, "method_summary.csv"))

coverage_summary <- data.frame()
if (nrow(coords)) {
  ckeys <- interaction(coords$experiment_id, coords$method_id, coords$coordinate, drop = TRUE)
  coverage <- lapply(split(coords, ckeys), function(d) {
    cov <- summarise_coverage_v2(d$covered, scheduled = nrow(d))
    empirical_sd <- stats::sd(d$beta_tilde, na.rm = TRUE)
    mean_se <- mean(d$estimated_se, na.rm = TRUE)
    data.frame(
      experiment_id = d$experiment_id[1],
      method_id = d$method_id[1],
      coordinate = d$coordinate[1],
      coordinate_type = d$coordinate_type[1],
      target_type = d$target_type[1],
      bias = mean(d$beta_tilde - d$target_value, na.rm = TRUE),
      bias_mcse = mcse_mean_v2(d$beta_tilde - d$target_value),
      empirical_sd = empirical_sd,
      mean_se = mean_se,
      se_sd_ratio = mean_se / empirical_sd,
      mean_interval_length = mean(d$interval_length, na.rm = TRUE),
      cov,
      stringsAsFactors = FALSE
    )
  })
  coverage_summary <- do.call(rbind, coverage)
  write_atomic_csv_v2(coverage_summary, file.path(out_dir, "coverage_summary.csv"))
}

if (nrow(theory)) {
  tkeys <- interaction(theory$experiment_id, theory$method_id, drop = TRUE)
  theory_summary <- lapply(split(theory, tkeys), function(d) data.frame(
    experiment_id = d$experiment_id[1],
    method_id = d$method_id[1],
    mean_hessian_max_error = mean(d$hessian_max_error, na.rm = TRUE),
    mean_precision_row_error = mean(d$precision_row_error, na.rm = TRUE),
    mean_one_step_correction = mean(d$one_step_correction_max, na.rm = TRUE),
    max_profile_identity_error = max(d$profile_identity_error, na.rm = TRUE),
    max_nuisance_gradient = max(d$max_nuisance_gradient, na.rm = TRUE),
    max_dantzig_residual = max(d$max_dantzig_residual, na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
  write_atomic_csv_v2(do.call(rbind, theory_summary),
                      file.path(out_dir, "theory_summary.csv"))
}

if ("POOLED-QR-LASSO" %in% metrics$method_id) {
  paired <- paired_difference_summary_v2(
    metrics, value = "pinball_loss", method = "method_id",
    reference = "POOLED-QR-LASSO",
    id = c("experiment_id", "replicate")
  )
  write_atomic_csv_v2(paired, file.path(out_dir, "paired_pinball_differences.csv"))
}

if (pilot_gate_requested) {
  proposed <- metrics[metrics$method_id == "PROFILE-DQR", , drop = FALSE]
  required_ids <- sprintf("P%02d", 1:4)
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
  # Inverse-defect / row-accuracy gate (theory decision section 3.5): where a
  # POP-H asset exists, require median(D_k) < 0.5 and Q0.9(D_k) < 1 for the
  # baseline; a violation means the precision approximation is not yet in the
  # first-order regime, not that the interval should be widened.
  d_vals <- proposed_coords$D_k[is.finite(proposed_coords$D_k)]
  inverse_defect_pass <- length(d_vals) > 0 &&
    stats::median(d_vals) < 0.5 &&
    stats::quantile(d_vals, 0.9, na.rm = TRUE, names = FALSE) < 1
  baseline_cov <- coverage_summary[
    coverage_summary$experiment_id == "P01" &
      coverage_summary$method_id == "PROFILE-DQR", , drop = FALSE
  ]
  se_ratio_pass <- nrow(baseline_cov) &&
    all(baseline_cov$se_sd_ratio >= 0.80 & baseline_cov$se_sd_ratio <= 1.20, na.rm = TRUE)
  coverage_pass <- nrow(baseline_cov) &&
    all(baseline_cov$coverage >= 0.88 & baseline_cov$coverage <= 0.98, na.rm = TRUE)
  runner_failure_path <- file.path(run_dir, "runner_failures.csv")
  runner_failure_count <- if (file.exists(runner_failure_path)) nrow(read.csv(runner_failure_path)) else 0L
  missing_method_pass <- !any(metrics$status == "not_implemented")
  pass <- convergence_pass && identity_pass && dantzig_pass && inverse_defect_pass &&
    se_ratio_pass && coverage_pass && runner_failure_count == 0L && missing_method_pass
  manifest <- list(
    pass = pass,
    commit_sha = trimws(readLines(file.path(run_dir, "implementation_commit.txt"), warn = FALSE)[1]),
    run_id = run_id,
    run_dir = normalizePath(run_dir, mustWork = TRUE),
    evaluated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    convergence_pass = convergence_pass,
    identity_pass = identity_pass,
    dantzig_pass = dantzig_pass,
    inverse_defect_pass = inverse_defect_pass,
    se_ratio_pass = se_ratio_pass,
    coverage_pass = coverage_pass,
    no_runner_failures = runner_failure_count == 0L,
    no_missing_methods = missing_method_pass,
    thresholds = list(
      convergence = 0.98,
      identity_median = 1e-8,
      identity_max = 1e-6,
      dantzig_feasible = 0.95,
      inverse_defect_median = 0.5,
      inverse_defect_q90 = 1,
      se_sd_ratio = c(0.80, 1.20),
      baseline_coverage = c(0.88, 0.98)
    )
  )
  write_json_v2(manifest, pilot_gate_manifest_path_v2(root))
  if (!pass) warning("Pilot gate did not pass. Diagnose; do not calibrate to the thresholds.")
}
message("Summaries written to: ", out_dir)

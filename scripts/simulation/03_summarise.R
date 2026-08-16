#!/usr/bin/env Rscript
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "scripts/simulation/03_summarise.R"
root <- if (file.exists("R/metrics_v2.R")) "." else normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
run_id <- cli$`run-id` %||% stop("Supply --run-id=<run directory name>.")
run_dir <- file.path(root, "results", "raw", "simulation", run_id)
metrics_path <- file.path(run_dir, "replication_metrics.csv")
coords_path <- file.path(run_dir, "coordinate_metrics.csv")
if (!file.exists(metrics_path)) stop("Missing ", metrics_path)
metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
coords <- if (file.exists(coords_path)) read.csv(coords_path, stringsAsFactors = FALSE) else data.frame()

out_dir <- file.path(root, "results", "tables", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keys <- interaction(metrics$experiment_id, metrics$method_id, drop = TRUE)
summary_rows <- lapply(split(metrics, keys), function(d) {
  data.frame(
    experiment_id = d$experiment_id[1],
    method_id = d$method_id[1],
    scheduled = nrow(d),
    success = sum(d$status %in% c("ok", "warning")),
    failure_rate = mean(d$status %in% c("failed", "not_implemented")),
    mean_l2 = mean(d$l2_error, na.rm = TRUE),
    mcse_l2 = mcse_mean_v2(d$l2_error),
    mean_pinball = mean(d$pinball_loss, na.rm = TRUE),
    mcse_pinball = mcse_mean_v2(d$pinball_loss),
    mean_runtime = mean(d$runtime_sec, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
summary_metrics <- do.call(rbind, summary_rows)
write.csv(summary_metrics, file.path(out_dir, "method_summary.csv"), row.names = FALSE)

if (nrow(coords)) {
  ckeys <- interaction(coords$experiment_id, coords$method_id, coords$coordinate, drop = TRUE)
  coverage <- lapply(split(coords, ckeys), function(d) {
    cov <- summarise_coverage_v2(d$covered, scheduled = nrow(d))
    data.frame(experiment_id = d$experiment_id[1], method_id = d$method_id[1],
               coordinate = d$coordinate[1], cov, stringsAsFactors = FALSE)
  })
  write.csv(do.call(rbind, coverage), file.path(out_dir, "coverage_summary.csv"), row.names = FALSE)
}
message("Summaries written to: ", out_dir)

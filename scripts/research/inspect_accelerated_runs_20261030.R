library(data.table)
root <- "results/raw/simulation/submission_20261028_main_sqr"
transitions <- list.files(root, "^numerical_transition_[0-9]+[.]rds$", full.names = TRUE)
if (!length(transitions)) stop("The numerical transition record is unavailable.")
transition <- readRDS(transitions[which.max(file.info(transitions)$mtime)])
dispatchers <- fread("results/scratch/sqr_acceleration_20261030/active_dispatchers.csv",
  colClasses = list(character = "started"))
started <- as.POSIXct(dispatchers[queue == "sqr", started], tz = "Asia/Shanghai")
if (length(started) != 1L || is.na(started)) stop("The active SQR start time is unavailable.")
observed <- Sys.time()
paths <- list.files(file.path(root, "replications"), "^result[.]rds$", recursive = TRUE, full.names = TRUE)
modified <- file.info(paths)$mtime
keep <- modified <= observed
paths <- paths[keep]
modified <- modified[keep]
recent <- which(modified > transition$created)
records <- rbindlist(lapply(recent, function(i) {
  result <- readRDS(paths[i])
  row <- as.data.table(result$metrics)
  if (nrow(row) != 1L || row$method_id != "SQR-DEBIASED-IID") stop("An accelerated result has a different method identity.")
  row <- row[, .(experiment_id, replicate, n_clusters, p_original, status, converged,
    runtime_sec, failure_stage, failure_message, implementation_commit)]
  execution <- result$method_diagnostics[[1L]]$tuning$gamma_cv$execution
  row[, `:=`(completed = as.character(modified[i]), source = paths[i],
    numerical_specification = result$numerical_specification$version,
    execution = if (is.null(execution)) NA_character_ else execution)]
  row
}), fill = TRUE)
if (!nrow(records)) stop("The accelerated queue has not yet produced a completed result.")
if (any(records$numerical_specification != "submission_comparator_numerics_20261030") ||
    any(records$converged & (is.na(records$execution) | records$execution != "certified_column_and_fold_loss_bounds")))
  stop("A resumed result does not use the accepted numerical implementation.")
request <- readRDS(file.path(root, "request.rds"))
elapsed <- as.numeric(difftime(observed, started, units = "hours"))
if (!is.finite(elapsed) || elapsed <= 0) stop("The observation time must follow the dispatcher start time.")
summary <- records[, .(completed = .N, converged = sum(converged),
  median_fitting_seconds = median(runtime_sec), maximum_fitting_seconds = max(runtime_sec)),
  by = .(n_clusters, p_original)]
output <- "results/preflight/sqr_acceleration_20261030"
fwrite(records, file.path(output, "resumed_records.csv"))
fwrite(summary, file.path(output, "resumed_runtime_summary.csv"))
report <- list(observed = as.character(observed), started = as.character(started),
  completed_before_transition = length(paths) - length(recent), completed_after_transition = nrow(records),
  converged_after_transition = sum(records$converged), elapsed_hours = elapsed,
  completed_per_hour = nrow(records) / elapsed, registered_datasets = nrow(request$jobs),
  remaining_datasets = nrow(request$jobs) - length(paths), source = root)
jsonlite::write_json(report, file.path(output, "resumed_progress.json"), auto_unbox = TRUE, pretty = TRUE)
print(report)
print(summary)
if (any(!records$converged)) print(records[!converged, .(experiment_id, replicate, failure_stage, failure_message)])

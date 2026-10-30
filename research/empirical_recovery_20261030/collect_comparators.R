library(data.table)
source("research/empirical_recovery_20261030/queue.R")
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1L] else "audit"
if (!mode %in% c("audit", "final")) stop("Select audit or final collection.")
cohort <- readRDS("data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds")
metadata <- cohort$eligible_metadata
request <- readRDS(file.path(empirical_comparator_root, "request.rds"))
jobs <- as.data.table(request$jobs)
status <- fread(file.path(empirical_queue_root, "status.csv"))
expected <- file.path(empirical_comparator_root, sprintf("job_%04d", jobs$job_id))
position <- match(expected, status$output)
if (anyNA(position) || nrow(jobs) != 1350L || anyDuplicated(position))
  stop("The comparison registry and dispatch records do not agree.")
jobs[, `:=`(state = status$state[position], exit_status = status$exit_status[position], output = expected)]
if (mode == "final" && any(!jobs$state %in% c("completed", "failed")))
  stop("The registered prediction comparisons are still running.")
identity <- c("method", "model", "tau", "processing", "repetition", "fold")
accounting <- predictions <- list()
for (i in which(jobs$state %in% c("completed", "failed"))) {
  job <- jobs[i]
  result_path <- file.path(job$output, "result.rds")
  failure_path <- file.path(job$output, "failure.rds")
  successful <- identical(job$state, "completed") && job$exit_status == 0L
  if (successful != file.exists(result_path) || successful == file.exists(failure_path))
    stop("A comparison output does not match its process record: ", job$output)
  value <- readRDS(if (successful) result_path else failure_path)
  if (!isTRUE(all.equal(as.data.frame(job[, ..identity]), value$job[, identity], check.attributes = FALSE)))
    stop("A comparison output differs from its registered specification.")
  row <- job[, ..identity]
  row[, `:=`(successful = successful, failure_stage = "", failure_message = "",
    source_file = if (successful) result_path else failure_path)]
  if (successful) {
    data <- as.data.table(value$predictions)
    heldout <- cohort$subject_folds$subject_id[cohort$subject_folds$repetition == job$repetition &
      cohort$subject_folds$fold == job$fold]
    expected_samples <- rownames(metadata)[metadata$subject_id %in% heldout]
    index <- match(data$sample_id, rownames(metadata))
    if (!isTRUE(value$fit$converged) || anyDuplicated(data$sample_id) || anyNA(index) ||
        !setequal(data$sample_id, expected_samples) ||
        any(as.character(data$subject_id) != as.character(metadata$subject_id[index])))
      stop("The comparison predictions do not cover the declared heldout patients.")
    for (field in identity) if (any(data[[field]] != job[[field]])) stop("A prediction identity differs.")
    residual <- metadata$log_c3[index] - data$predicted_quantile
    loss <- residual * (job$tau - (residual < 0))
    if (any(!is.finite(data$predicted_quantile)) ||
        max(abs(data$observed_log_c3 - metadata$log_c3[index])) > 1e-12 ||
        max(abs(data$loss - loss)) > 1e-12) stop("The saved comparison losses do not reproduce the observed outcomes.")
    row[, `:=`(heldout_subjects = uniqueN(data$subject_id), heldout_visits = nrow(data),
      selected_dimension = value$selected$screen_dimension, selected_multiplier = value$selected$lambda_multiplier)]
    predictions[[length(predictions) + 1L]] <- data
  } else row[, `:=`(failure_stage = value$stage, failure_message = conditionMessage(value$condition),
    heldout_subjects = 0L, heldout_visits = 0L, selected_dimension = NA_integer_, selected_multiplier = NA_real_)]
  accounting[[length(accounting) + 1L]] <- row
}
accounting <- rbindlist(accounting)
predictions <- rbindlist(predictions)
if (!nrow(accounting) || anyDuplicated(accounting[, ..identity]) ||
    anyDuplicated(predictions[, c(identity, "sample_id"), with = FALSE])) stop("Comparison records are empty or duplicated.")
counts <- accounting[, .(attempted = .N, successful = sum(successful), failed = sum(!successful)),
  by = .(processing, method, model, tau)]
if (mode == "final" && (nrow(accounting) != 1350L || nrow(counts) != 54L || any(counts$attempted != 25L)))
  stop("The full comparison study does not account for every registered fit.")
output <- if (mode == "audit") "results/scratch/empirical_recovery_20261030/comparison_audit" else
  "results/tables/empirical_recovery_20261030/application"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(accounting, file.path(output, "job_accounting.csv"))
fwrite(counts, file.path(output, "setting_accounting.csv"))
if (mode == "final") {
  fwrite(predictions, file.path(output, "heldout_predictions.csv"))
  saveRDS(list(accounting = accounting, predictions = predictions,
    subject_folds = cohort$subject_folds, collected = Sys.time()), file.path(output, "records.rds"))
}
cat("Verified", nrow(accounting), "comparison outcomes and", nrow(predictions), "observed visit predictions\n")

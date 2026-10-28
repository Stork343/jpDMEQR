library(data.table)
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1L] else "final"
if (!mode %in% c("audit", "final")) stop("Unknown application collection mode.")
roots <- c(mixed = "results/raw/application_comparison/submission_20261027_v2",
  qgee = "results/raw/application_comparison/submission_20261028_qgee")
cohort <- readRDS("data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds")
metadata <- cohort$eligible_metadata
folds <- as.data.table(cohort$subject_folds)
identity_fields <- c("method", "model", "tau", "processing", "repetition", "fold")
accounting <- predictions <- vector("list", 1350L)
count <- 0L
for (family in names(roots)) {
  root <- roots[[family]]
  request <- readRDS(file.path(root, "request.rds"))
  if (mode == "final" && !file.exists(file.path(root, "completed.rds"))) stop("An application queue is incomplete: ", root)
  jobs <- as.data.table(request$jobs)
  if (!"job_id" %in% names(jobs)) jobs[, job_id := seq_len(.N)]
  jobs <- if (family == "qgee") jobs[method == "QGEE-SCAD"] else jobs[method != "QGEE-SCAD"]
  status <- fread(file.path(root, "status.csv"))
  if (!"job_id" %in% names(status)) status[, job_id := seq_len(.N)]
  jobs <- merge(jobs[, c(identity_fields, "job_id"), with = FALSE],
    status[, .(job_id, state, exit_status)], by = "job_id", all.x = TRUE, sort = FALSE)
  jobs[, record_root := root]
  if (family == "qgee") {
    recovery_root <- "results/raw/application_comparison/submission_20261028_qgee_recovery"
    if (file.exists(file.path(recovery_root, "completed.rds"))) {
      recovered <- fread(file.path(recovery_root, "status.csv"))
      for (id in recovered$job_id) {
        value <- recovered[job_id == id]
        jobs[job_id == id, `:=`(state = value$state, exit_status = value$exit_status, record_root = recovery_root)]
      }
    }
  }
  if (anyNA(jobs$state) || anyDuplicated(jobs$job_id)) stop("Application dispatch accounting differs.")
  if (mode == "final" && any(jobs$state != "completed")) stop("Unfinished application jobs remain.")
  for (i in which(jobs$state == "completed")) {
    job <- jobs[i]
    path <- file.path(job$record_root, sprintf("job_%04d", job$job_id))
    result_file <- file.path(path, "result.rds")
    failure_file <- file.path(path, "failure.rds")
    if (mode == "audit" && job$exit_status < 0L && !file.exists(failure_file)) next
    success <- job$exit_status == 0L
    if (success != file.exists(result_file) || success == file.exists(failure_file))
      stop("Application result and dispatch status differ: ", path)
    result <- readRDS(if (success) result_file else failure_file)
    if (!isTRUE(all.equal(as.data.frame(job[, ..identity_fields]),
      as.data.frame(result$job[, identity_fields]), check.attributes = FALSE))) stop("Application job identity differs: ", path)
    count <- count + 1L
    row <- job[, ..identity_fields]
    row[, `:=`(successful = success, failure_stage = "", failure_message = "", source_file = if (success) result_file else failure_file)]
    if (success) {
      if (!isTRUE(result$fit$converged)) stop("A successful application record has an unconverged fit.")
      data <- as.data.table(result$predictions)
      expected_subjects <- folds[repetition == job$repetition & fold == job$fold, subject_id]
      expected_samples <- rownames(metadata)[metadata$subject_id %in% expected_subjects]
      index <- match(data$sample_id, rownames(metadata))
      if (anyDuplicated(data$sample_id) || !setequal(data$sample_id, expected_samples) || anyNA(index) ||
          any(as.character(data$subject_id) != as.character(metadata$subject_id[index]))) stop("Heldout application subjects differ.")
      for (field in identity_fields) if (any(data[[field]] != job[[field]])) stop("Prediction identity differs.")
      if (any(!is.finite(data$predicted_quantile)) ||
          max(abs(data$observed_log_c3 - metadata$log_c3[index])) > 1e-12) stop("Prediction outcome values differ.")
      residual <- data$observed_log_c3 - data$predicted_quantile
      computed_loss <- residual * (job$tau - (residual < 0))
      if (any(!is.finite(data$loss)) || max(abs(data$loss - computed_loss)) > 1e-12) stop("Saved pinball losses differ.")
      row[, `:=`(heldout_subjects = uniqueN(data$subject_id), heldout_visits = nrow(data),
        selected_dimension = result$selected$screen_dimension, selected_multiplier = result$selected$lambda_multiplier)]
      predictions[[count]] <- data
    } else {
      row[, `:=`(failure_stage = result$stage, failure_message = conditionMessage(result$condition),
        heldout_subjects = 0L, heldout_visits = 0L, selected_dimension = NA_integer_, selected_multiplier = NA_real_)]
    }
    accounting[[count]] <- row
  }
}
accounting <- rbindlist(accounting, fill = TRUE)
predictions <- rbindlist(predictions, fill = TRUE)
if (anyDuplicated(accounting[, ..identity_fields]) ||
    anyDuplicated(predictions[, c(identity_fields, "sample_id"), with = FALSE])) stop("Duplicated application results.")
if (mode == "final" && nrow(accounting) != 1350L) stop("The application comparison requires 1350 job outcomes.")
counts <- accounting[, .(attempted = .N, successful = sum(successful), failed = sum(!successful)),
  by = .(processing, method, model, tau)]
if (mode == "final" && (nrow(counts) != 54L || any(counts$attempted != 25L))) stop("Incomplete application setting counts.")
output <- if (mode == "final") "results/tables/submission_20261028/application" else
  "results/scratch/submission_completion_20261027/application_collection_audit"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(accounting, file.path(output, "job_accounting.csv"))
fwrite(counts, file.path(output, "setting_accounting.csv"))
if (mode == "final") {
  fwrite(predictions, file.path(output, "heldout_predictions.csv"))
  saveRDS(list(accounting = accounting, predictions = predictions, roots = roots,
    subject_folds = folds, collected = Sys.time()), file.path(output, "records.rds"))
}
cat("Verified", nrow(accounting), "application job outcomes and", nrow(predictions), "heldout visit predictions\n")

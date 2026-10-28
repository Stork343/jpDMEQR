source("research/submission_completion_20261027/application.R")
source("research/submission_completion_20261027/application_state.R")
args <- commandArgs(trailingOnly = TRUE)
root <- args[1L]
job_id <- as.integer(args[2L])
request <- readRDS(file.path(root, "request.rds"))
job <- request$jobs[request$jobs$job_id == job_id, c("method", "model", "tau", "processing", "repetition", "fold")]
if (nrow(job) != 1L || anyNA(job) || job$method != "QGEE-SCAD") stop("Invalid QGEE application job.")
output <- file.path(root, sprintf("job_%04d", job_id))
if (dir.exists(output)) stop("The QGEE application job already exists.")
dir.create(output, recursive = TRUE)
saveRDS(job, file.path(output, "request.rds"))
stage <- "preprocessing"
globalCallingHandlers(error = function(condition) saveRDS(list(job = job, stage = stage,
  condition = condition, finished = Sys.time()), file.path(output, "failure.rds")))
source_output <- file.path(request$source_root, sprintf("job_%04d", job_id))
source_finished <- submission_application_completed(request$source_root, job_id)
cohort <- readRDS("data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds")
metadata <- cohort$eligible_metadata
heldout <- cohort$subject_folds$subject_id[cohort$subject_folds$repetition == job$repetition & cohort$subject_folds$fold == job$fold]
inner_roles <- read.csv("data/derived/GSE65391/execution_20261022_v1/inner_subject_folds.csv")
inner_roles <- inner_roles[inner_roles$outer_repetition == job$repetition & inner_roles$outer_fold == job$fold, ]
grid <- data.frame(lambda_multiplier = 1, screen_dimension = c(100L, 200L, 400L), candidate_id = 1:3)
losses <- matrix(NA_real_, 5L, 3L)
valid <- matrix(FALSE, 5L, 3L)
native_warnings <- list()
valid_saved_fit <- function(fit, design) {
  !is.null(fit$beta_prediction) && identical(names(fit$beta_prediction), colnames(design$X)) &&
    all(is.finite(fit$beta_prediction)) && !is.null(fit$lambda_path) && nrow(fit$lambda_path) == 30L &&
    (identical(fit$penalty_failure_policy, "source_HBIC_individual_penalties") || all(is.finite(fit$lambda_path$hbic)))
}
fit_one <- function(design, train, seed) {
  withCallingHandlers(tryCatch(
    submission_qgee(metadata, design, train, job$tau, seed, allow_failed_penalties = TRUE),
    error = function(condition) list(converged = FALSE, status = "failed", condition = condition,
      failure_message = conditionMessage(condition))), warning = function(condition) {
    native_warnings[[length(native_warnings) + 1L]] <<- list(stage = stage,
      message = conditionMessage(condition), call = deparse(conditionCall(condition)))
  })
}
started <- proc.time()[[3L]]
reused <- 0L
for (g in 1:5) {
  base <- submission_application_base(job$processing, job$repetition, job$fold, g, metadata)
  training <- seq_len(nrow(metadata)) %in% base$training_rows
  testing <- metadata$subject_id %in% inner_roles$subject_id[inner_roles$inner_fold == g]
  if (any(training & testing) || any(metadata$subject_id[training] %in% heldout)) stop("QGEE subject roles overlap.")
  for (j in 1:3) {
    stage <- paste0("inner_fold_", g, "_candidate_", j)
    design <- submission_application_design(base, metadata, job$model, job$tau, grid$screen_dimension[j])
    source_file <- file.path(source_output, sprintf("inner_%d_candidate_%02d.rds", g, j))
    old <- if (source_finished && file.exists(source_file)) readRDS(source_file) else NULL
    reuse <- !is.null(old) && valid_saved_fit(old$fit, design)
    if (reuse) {
      if (!isTRUE(all.equal(old$clinical_transform, design$clinical_transform, tolerance = 1e-12)) ||
          !setequal(old$training_subjects, unique(metadata$subject_id[training]))) stop("Cached QGEE training definitions differ.")
      fitted <- old$fit
      reused <- reused + 1L
    } else fitted <- fit_one(design, training,
      20260817L + job$repetition * 100000L + job$fold * 10000L + g * 100L)
    valid[g, j] <- isTRUE(fitted$converged)
    subject_loss <- setNames(rep(NA_real_, length(unique(metadata$subject_id[testing]))), unique(metadata$subject_id[testing]))
    if (valid[g, j]) {
      prediction <- drop(design$X[testing, , drop = FALSE] %*% fitted$beta_prediction)
      residual <- metadata$log_c3[testing] - prediction
      subject_loss <- tapply(residual * (job$tau - (residual < 0)), metadata$subject_id[testing], mean)
      if (reuse && max(abs(subject_loss - old$subject_loss[match(names(subject_loss), names(old$subject_loss))])) > 1e-12)
        stop("Cached QGEE validation losses differ.")
      losses[g, j] <- sum(subject_loss)
    }
    saveRDS(list(fit = fitted, clinical_transform = design$clinical_transform,
      training_subjects = unique(metadata$subject_id[training]), validation_subjects = names(subject_loss),
      subject_loss = subject_loss, reused_source = if (reuse) source_file else NULL),
      file.path(output, sprintf("inner_%d_candidate_%02d.rds", g, j)))
  }
}
stage <- "candidate_selection"
grid$loss <- colSums(losses) / length(unique(metadata$subject_id[!metadata$subject_id %in% heldout]))
grid$eligible <- colSums(valid) == 5L
write.csv(grid, file.path(output, "candidate_scores.csv"), row.names = FALSE)
if (!any(grid$eligible)) stop("No QGEE screen candidate converged in all five inner folds.")
eligible <- which(grid$eligible)
selected <- eligible[order(grid$loss[eligible], grid$candidate_id[eligible])][1L]
base <- submission_application_base(job$processing, job$repetition, job$fold, 0L, metadata)
design <- submission_application_design(base, metadata, job$model, job$tau, grid$screen_dimension[selected])
training <- seq_len(nrow(metadata)) %in% base$training_rows
testing <- metadata$subject_id %in% heldout
if (any(training & testing) || !all(training | testing)) stop("Incomplete outer QGEE subject roles.")
stage <- "outer_fit"
old_path <- file.path(source_output, "result.rds")
old <- if (source_finished && file.exists(old_path)) readRDS(old_path) else NULL
reuse_outer <- !is.null(old) && old$selected$screen_dimension == grid$screen_dimension[selected] && valid_saved_fit(old$fit, design)
if (reuse_outer) {
  if (!isTRUE(all.equal(old$clinical_transform, design$clinical_transform, tolerance = 1e-12))) stop("Cached outer QGEE units differ.")
  fitted <- old$fit
  reused <- reused + 1L
} else fitted <- fit_one(design, training, 20260817L + job$repetition * 100000L + job$fold * 10000L)
saveRDS(fitted, file.path(output, "outer_fit.rds"))
if (!isTRUE(fitted$converged)) stop("The selected outer QGEE fit did not converge.")
prediction <- drop(design$X[testing, , drop = FALSE] %*% fitted$beta_prediction)
if (any(!is.finite(prediction))) stop("Nonfinite heldout QGEE prediction.")
if (reuse_outer && max(abs(prediction - old$predictions$predicted_quantile[
    match(rownames(metadata)[testing], old$predictions$sample_id)])) > 1e-12) stop("Cached outer QGEE predictions differ.")
residual <- metadata$log_c3[testing] - prediction
predictions <- data.frame(sample_id = rownames(metadata)[testing], subject_id = metadata$subject_id[testing],
  observed_log_c3 = metadata$log_c3[testing], predicted_quantile = prediction,
  loss = residual * (job$tau - (residual < 0)), repetition = job$repetition, fold = job$fold,
  model = job$model, tau = job$tau, method = job$method, processing = job$processing,
  filter_size = 5000L, random_slope = FALSE)
saveRDS(list(job = job, predictions = predictions, fit = fitted, selected = grid[selected, ],
  clinical_transform = design$clinical_transform, native_warnings = native_warnings,
  reused_candidate_fits = reused, new_calculation_seconds = proc.time()[[3L]] - started,
  penalty_failure_policy = "source_HBIC_individual_penalties", completed = Sys.time()), file.path(output, "result.rds"))
write.csv(predictions, file.path(output, "predictions.csv"), row.names = FALSE)
cat("Completed source QGEE job", job_id, "\n")

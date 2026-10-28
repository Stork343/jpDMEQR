source("research/score_validity_20261026/load_joint.R")
source_v2_module(".", "application_v2", .GlobalEnv)
source("research/submission_completion_20261027/application.R")
args <- commandArgs(trailingOnly = TRUE)
root <- args[1L]
job_id <- as.integer(args[2L])
request <- readRDS(file.path(root, "request.rds"))
job <- request$jobs[job_id, c("method", "model", "tau", "processing", "repetition", "fold")]
if (nrow(job) != 1L || anyNA(job)) stop("Invalid comparator job.")
output <- file.path(root, sprintf("job_%04d", job_id))
if (dir.exists(output)) stop("The comparator job already exists.")
dir.create(output, recursive = TRUE)
saveRDS(job, file.path(output, "request.rds"))
stage <- "preprocessing"
globalCallingHandlers(error = function(condition) saveRDS(list(job = job, stage = stage,
  condition = condition, finished = Sys.time()), file.path(output, "failure.rds")))
cohort <- readRDS("data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds")
metadata <- cohort$eligible_metadata
heldout <- cohort$subject_folds$subject_id[cohort$subject_folds$repetition == job$repetition & cohort$subject_folds$fold == job$fold]
inner_roles <- read.csv("data/derived/GSE65391/execution_20261022_v1/inner_subject_folds.csv")
inner_roles <- inner_roles[inner_roles$outer_repetition == job$repetition & inner_roles$outer_fold == job$fold, ]
grid <- if (job$method == "POOLED-QR-LASSO") {
  expand.grid(lambda_multiplier = c(.25, .5, 1, 2), screen_dimension = c(100L, 200L, 400L))
} else if (job$method == "QGEE-SCAD") {
  data.frame(lambda_multiplier = 1, screen_dimension = c(100L, 200L, 400L))
} else data.frame(lambda_multiplier = 1, screen_dimension = 0L)
grid$candidate_id <- seq_len(nrow(grid))
losses <- matrix(NA_real_, 5L, nrow(grid))
valid <- matrix(FALSE, 5L, nrow(grid))
native_warnings <- list()
fit_one <- function(design, train, multiplier, seed) {
  withCallingHandlers({
    if (job$method == "QGEE-SCAD") {
      candidate_path <- file.path(output, stage)
      dir.create(candidate_path)
      saveRDS(list(metadata = metadata, design = design, training = train, tau = job$tau, seed = seed),
               file.path(candidate_path, "input.rds"))
      process <- processx::run(file.path(R.home("bin"), "Rscript"),
        c("research/submission_completion_20261027/run_qgee_candidate.R", candidate_path),
        error_on_status = FALSE, stdout = file.path(candidate_path, "solver.log"), stderr = "2>&1")
      if (process$status != 0L) {
        return(list(converged = FALSE, status = "failed", process_status = process$status,
                    failure_file = file.path(candidate_path, "failure.rds")))
      }
      readRDS(file.path(candidate_path, "result.rds"))
    } else if (job$method == "LQMM") {
      submission_lqmm(metadata, design, train, job$tau, seed)
    } else
      fit_application_comparator_v2(job$method, metadata, design, train, job$tau, multiplier, seed)
  }, warning = function(condition) {
    native_warnings[[length(native_warnings) + 1L]] <<- list(stage = stage,
      message = conditionMessage(condition), call = deparse(conditionCall(condition)))
  })
}
started <- proc.time()[[3L]]
for (g in 1:5) {
  base <- submission_application_base(job$processing, job$repetition, job$fold, g, metadata)
  training <- seq_len(nrow(metadata)) %in% base$training_rows
  testing <- metadata$subject_id %in% inner_roles$subject_id[inner_roles$inner_fold == g]
  if (any(training & testing) || any(metadata$subject_id[training] %in% heldout)) stop("Comparator subject roles overlap.")
  for (dimension in unique(grid$screen_dimension)) {
    design <- submission_application_design(base, metadata, job$model, job$tau, dimension)
    for (j in which(grid$screen_dimension == dimension)) {
      stage <- paste0("inner_fold_", g, "_candidate_", j)
      fitted <- fit_one(design, training, grid$lambda_multiplier[j],
        20260817L + job$repetition * 100000L + job$fold * 10000L + g * 100L)
      valid[g, j] <- isTRUE(fitted$converged)
      subject_loss <- setNames(rep(NA_real_, length(unique(metadata$subject_id[testing]))), unique(metadata$subject_id[testing]))
      if (valid[g, j]) {
        prediction <- drop(design$X[testing, , drop = FALSE] %*% fitted$beta_prediction)
        residual <- metadata$log_c3[testing] - prediction
        loss <- residual * (job$tau - (residual < 0))
        subject_loss <- tapply(loss, metadata$subject_id[testing], mean)
        losses[g, j] <- sum(subject_loss)
      }
      saveRDS(list(fit = fitted, clinical_transform = design$clinical_transform,
                   training_subjects = unique(metadata$subject_id[training]),
                   validation_subjects = names(subject_loss), subject_loss = subject_loss),
               file.path(output, sprintf("inner_%d_candidate_%02d.rds", g, j)))
    }
  }
}
stage <- "candidate_selection"
grid$loss <- colSums(losses) / length(unique(metadata$subject_id[!metadata$subject_id %in% heldout]))
grid$eligible <- colSums(valid) == 5L
write.csv(grid, file.path(output, "candidate_scores.csv"), row.names = FALSE)
if (!any(grid$eligible)) stop("No comparator candidate converged in every inner subject fold.")
eligible <- which(grid$eligible)
selected <- eligible[order(grid$loss[eligible], grid$candidate_id[eligible])][1L]
base <- submission_application_base(job$processing, job$repetition, job$fold, 0L, metadata)
design <- submission_application_design(base, metadata, job$model, job$tau, grid$screen_dimension[selected])
training <- seq_len(nrow(metadata)) %in% base$training_rows
testing <- metadata$subject_id %in% heldout
if (any(training & testing) || !all(training | testing)) stop("Incomplete outer subject roles.")
stage <- "outer_fit"
fitted <- fit_one(design, training, grid$lambda_multiplier[selected],
                   20260817L + job$repetition * 100000L + job$fold * 10000L)
saveRDS(fitted, file.path(output, "outer_fit.rds"))
if (!isTRUE(fitted$converged)) stop("The selected outer comparator fit did not converge.")
prediction <- drop(design$X[testing, , drop = FALSE] %*% fitted$beta_prediction)
if (any(!is.finite(prediction))) stop("Nonfinite heldout prediction.")
residual <- metadata$log_c3[testing] - prediction
predictions <- data.frame(sample_id = rownames(metadata)[testing], subject_id = metadata$subject_id[testing],
  observed_log_c3 = metadata$log_c3[testing], predicted_quantile = prediction,
  loss = residual * (job$tau - (residual < 0)), repetition = job$repetition, fold = job$fold,
  model = job$model, tau = job$tau, method = job$method, processing = job$processing,
  filter_size = 5000L, random_slope = FALSE)
saveRDS(list(job = job, predictions = predictions, fit = fitted, selected = grid[selected, ],
             clinical_transform = design$clinical_transform, native_warnings = native_warnings,
             runtime = proc.time()[[3L]] - started, completed = Sys.time()), file.path(output, "result.rds"))
write.csv(predictions, file.path(output, "predictions.csv"), row.names = FALSE)
cat("Completed comparator job", job_id, job$method, "\n")

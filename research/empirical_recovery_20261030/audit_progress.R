library(data.table)
source("research/empirical_recovery_20261030/queue.R")
cohort <- readRDS("data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds")
metadata <- cohort$eligible_metadata
jobs <- empirical_registry()$jobs
checks <- list()
for (j in which(jobs$phase != "comparator_checks")) {
  job <- jobs[j, ]
  paths <- list.files(job$output, "^result[.]rds$",
    recursive = !job$phase %in% c("comparators", "processing_modules"), full.names = TRUE)
  if (!length(paths)) next
  request <- readRDS(file.path(job$output, "request.rds"))
  for (path in paths) {
    value <- readRDS(path)
    difference <- 0
    if (job$phase %in% c("prediction", "slope_prediction", "comparators")) {
      predictions <- value$predictions
      identity <- unique(predictions[, c("method", "model", "tau", "repetition", "fold", "processing")])
      if (nrow(identity) != 1L || anyDuplicated(predictions$sample_id)) stop("Prediction result identities are incomplete.")
      rows <- match(predictions$sample_id, rownames(metadata))
      heldout <- cohort$subject_folds$subject_id[cohort$subject_folds$repetition == identity$repetition &
        cohort$subject_folds$fold == identity$fold]
      expected <- which(metadata$subject_id %in% heldout)
      if (anyNA(rows) || !setequal(rows, expected) ||
          !identical(as.character(predictions$subject_id), as.character(metadata$subject_id[rows])) ||
          any(!is.finite(predictions$predicted_quantile))) stop("Prediction observations differ from the registered subject fold.")
      residual <- metadata$log_c3[rows] - predictions$predicted_quantile
      loss <- residual * (identity$tau - (residual < 0))
      difference <- max(abs(c(predictions$observed_log_c3 - metadata$log_c3[rows], predictions$loss - loss)))
      if (difference > 1e-12) stop("Saved prediction losses differ from the observed C3 measurements.")
      if (job$phase != "comparators") {
        if (length(intersect(request$selection_subjects, request$heldout_subjects)) ||
            !setequal(request$heldout_subjects, heldout)) stop("The prediction training and heldout subjects overlap.")
        candidates <- value$fit$candidates
        selected <- candidates[order(candidates$loss, candidates$size, candidates$lambda)[1L], ]
        if (!isTRUE(all.equal(value$fit$selected, selected)) ||
            max(value$fit$diagnostics$max_gap, value$fit$final_gap) > 1.01e-8)
          stop("The selected prediction candidate or convex certificate differs.")
      } else if (!isTRUE(value$fit$converged)) stop("A completed comparator result has a nonconverged fit.")
      units <- nrow(predictions)
    } else if (job$phase == "exploratory") {
      if (length(intersect(request$selection_subjects, request$inference_subjects)) ||
          !setequal(c(request$selection_subjects, request$inference_subjects), unique(metadata$subject_id)) ||
          length(value$selected) != 3L || nrow(value$table) != 3L ||
          !setequal(value$selected, value$table$coordinate) ||
          value$fit$diagnostics$n_clusters != length(request$inference_subjects))
        stop("The exploratory inference does not match its independent subject design.")
      influence <- value$fit$influence
      centered <- sweep(influence, 2L, colMeans(influence), "-")
      variance <- crossprod(centered) / length(request$inference_subjects)^2
      difference <- max(abs(variance - value$fit$covariance))
      if (difference > 1e-12 || any(!is.finite(value$table$se)) || any(value$table$se <= 0))
        stop("The exploratory covariance differs from its subject score contributions.")
      units <- nrow(value$table)
    } else {
      estimates <- value$estimates
      if (nrow(estimates) != 129L || ncol(estimates) != 5L ||
          !setequal(value$omitted_subjects, unique(metadata$subject_id)) ||
          anyDuplicated(value$omitted_subjects) || nrow(value$diagnostics) != 130L)
        stop("The complete module deletion workflow is incomplete.")
      centered <- sweep(estimates, 2L, colMeans(estimates), "-")
      se <- sqrt(128 / 129 * colSums(centered^2))
      difference <- max(abs(se - value$table$workflow_se))
      critical <- qnorm(1 - .05 / 60)
      difference <- max(difference, abs(value$table$workflow_simultaneous_lower -
        (value$table$beta_hat - critical * se)), abs(value$table$workflow_simultaneous_upper -
        (value$table$beta_hat + critical * se)))
      if (difference > 1e-12) stop("The module covariance or simultaneous intervals differ.")
      units <- nrow(value$table)
    }
    checks[[length(checks) + 1L]] <- data.table(job_id = job$job_id, phase = job$phase,
      processing = job$processing, numerical_difference = difference, records = units, source = path)
  }
}
checks <- rbindlist(checks)
if (!nrow(checks)) stop("No complete empirical outputs are available for validation.")
summary <- checks[, .(completed_fits = .N, records = sum(records), maximum_difference = max(numerical_difference)),
  by = .(phase, processing)]
output <- "results/preflight/empirical_recovery_20261030"
fwrite(checks, file.path(output, "completed_result_checks.csv"))
fwrite(summary, file.path(output, "progress_summary.csv"))
jsonlite::write_json(list(passed = TRUE, inspected = as.character(Sys.time()),
  prediction_fits = sum(checks$phase %in% c("prediction", "slope_prediction")),
  prediction_required = 900L, exploratory_fits = sum(checks$phase == "exploratory"),
  exploratory_required = 360L, processing_workflows = sum(checks$phase == "processing_modules"),
  processing_required = 12L, comparator_fits = sum(checks$phase == "comparators"), comparator_required = 1350L,
  maximum_numerical_difference = max(checks$numerical_difference)),
  file.path(output, "progress_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
print(summary)

root <- normalizePath(".", mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[3] %in% c("progress", "complete")) {
  stop("Supply the frozen workload directory, a new summary ID and progress or complete.")
}
workload <- normalizePath(args[1], mustWork = TRUE)
output <- file.path(root, "results", "tables", "application", args[2])
if (dir.exists(output)) stop("The immutable application summary directory already exists.")
registries <- lapply(c("module", "validation", "exploratory"), function(phase) {
  read.csv(file.path(workload, paste0(phase, "_jobs.csv")), stringsAsFactors = FALSE)
})
names(registries) <- c("module", "validation", "exploratory")
if (!identical(vapply(registries, nrow, integer(1)), c(module = 40L, validation = 1050L, exploratory = 120L))) {
  stop("The complete frozen application matrix is required.")
}
bind_tables <- function(tables) {
  if (!length(tables)) return(data.frame())
  fields <- unique(unlist(lapply(tables, names), use.names = FALSE))
  do.call(rbind, lapply(tables, function(table) {
    for (field in setdiff(fields, names(table))) table[[field]] <- NA
    table[fields]
  }))
}
accounting <- modules <- clinical <- metrics <- predictions <- list()
for (phase in names(registries)) {
  jobs <- registries[[phase]]
  for (i in seq_len(nrow(jobs))) {
    job <- jobs[i, , drop = FALSE]
    marker <- file.path(root, job$expected_manifest)
    directory <- dirname(marker)
    status <- "not_started"
    failure_stage <- failure_message <- NA_character_
    actual_source_commit <- numeric_certificate_commit <- NA_character_
    if (file.exists(marker)) {
      manifest <- jsonlite::read_json(marker, simplifyVector = TRUE)
      status <- "ok"
      implementation <- if (phase == "exploratory") manifest$implementation else manifest
      actual_source_commit <- implementation$actual_source_commit
      numeric_certificate_commit <- implementation$numeric_certificate_commit
      if (is.null(actual_source_commit) || is.null(numeric_certificate_commit)) {
        stop("A completed application job lacks its actual numerical certificate identity.")
      }
      if (phase == "module") {
        table <- read.csv(file.path(directory, "confirmatory_module_inference.csv"), stringsAsFactors = FALSE)
        covariates <- read.csv(file.path(directory, "clinical_coordinate_inference.csv"), stringsAsFactors = FALSE)
        if (nrow(table) != 5L || nrow(covariates) != (if (job$model == "A") 11L else 12L) ||
            any(table$n_clusters != 129L | table$jackknife_completed != 129L) ||
            any(table$score_definition != "profile_bandwidth_pair_2_half_minus_full")) {
          stop("A complete observed module inference table is invalid.")
        }
        for (field in c("job_id", "run_id", "specification")) {
          table[[field]] <- job[[field]]
          covariates[[field]] <- job[[field]]
        }
        table$actual_source_commit <- covariates$actual_source_commit <- actual_source_commit
        table$numeric_certificate_commit <- covariates$numeric_certificate_commit <- numeric_certificate_commit
        modules[[length(modules) + 1L]] <- table
        clinical[[length(clinical) + 1L]] <- covariates
      } else if (phase == "validation") {
        table <- read.csv(file.path(directory, "fold_metrics.csv"), stringsAsFactors = FALSE)
        heldout <- read.csv(file.path(directory, "predictions.csv"), stringsAsFactors = FALSE)
        if (nrow(table) != 1L || anyDuplicated(heldout$sample_id) ||
            any(!is.finite(heldout$observed_log_c3) | !is.finite(heldout$predicted_quantile) | !is.finite(heldout$loss))) {
          stop("A complete held-out prediction table is invalid.")
        }
        for (field in c("job_id", "run_id", "method", "model", "tau", "filter_size", "random_slope", "repetition", "fold")) {
          table[[field]] <- heldout[[field]] <- job[[field]]
        }
        table$actual_source_commit <- actual_source_commit
        table$numeric_certificate_commit <- numeric_certificate_commit
        metrics[[length(metrics) + 1L]] <- table
        predictions[[length(predictions) + 1L]] <- heldout
      }
    } else if (file.exists(file.path(directory, "failed_manifest.json"))) {
      manifest <- jsonlite::read_json(file.path(directory, "failed_manifest.json"), simplifyVector = TRUE)
      status <- "failed"
      failure_stage <- manifest$failure_stage
      failure_message <- manifest$message
      actual_source_commit <- manifest$implementation$actual_source_commit
      numeric_certificate_commit <- manifest$implementation$numeric_certificate_commit
    } else if (file.exists(file.path(directory, "run_request.rds")) || file.exists(file.path(directory, "request.rds"))) {
      status <- "started_incomplete"
    }
    accounting[[length(accounting) + 1L]] <- data.frame(
      phase = phase, job_id = job$job_id, run_id = job$run_id, status = status,
      failure_stage = failure_stage, failure_message = failure_message,
      actual_source_commit = actual_source_commit, numeric_certificate_commit = numeric_certificate_commit,
      expected_manifest = job$expected_manifest, stringsAsFactors = FALSE
    )
  }
}
accounting <- bind_tables(accounting)
if (args[3] == "complete" && any(!accounting$status %in% c("ok", "failed"))) {
  stop("Every scheduled real application job requires an actual result before complete summarisation.")
}
dir.create(output, recursive = TRUE)
write.csv(accounting, file.path(output, "application_registry_accounting.csv"), row.names = FALSE)
write.csv(accounting[accounting$status == "failed", ], file.path(output, "application_failure_log.csv"), row.names = FALSE)
modules <- bind_tables(modules)
clinical <- bind_tables(clinical)
if (nrow(modules)) {
  write.csv(modules, file.path(output, "confirmatory_module_inference_all_specifications.csv"), row.names = FALSE)
  write.csv(modules[modules$specification == "reference", ], file.path(output, "confirmatory_module_inference.csv"), row.names = FALSE)
  write.csv(modules[modules$specification != "reference", ], file.path(output, "module_sensitivity_inference.csv"), row.names = FALSE)
  write.csv(clinical, file.path(output, "clinical_coordinate_inference.csv"), row.names = FALSE)
}
metrics <- bind_tables(metrics)
predictions <- bind_tables(predictions)
if (nrow(metrics)) write.csv(metrics, file.path(output, "fold_level_metrics.csv"), row.names = FALSE)
if (nrow(predictions)) {
  write.csv(predictions, file.path(output, "fold_level_predictions.csv"), row.names = FALSE)
  predictions$calibration <- as.numeric(predictions$observed_log_c3 <= predictions$predicted_quantile)
  grouping <- c("method", "model", "tau", "filter_size", "random_slope")
  subjects <- aggregate(predictions[c("loss", "calibration")],
                         by = predictions[c(grouping, "repetition", "subject_id")], FUN = mean)
  write.csv(subjects, file.path(output, "heldout_subject_losses.csv"), row.names = FALSE)
  groups <- unique(subjects[grouping])
  match_group <- function(table, group) Reduce(`&`, lapply(grouping, function(field) table[[field]] == group[[field]]))
  repeated <- overall <- paired <- list()
  baseline <- subjects[subjects$method == "POOLED-QR-LASSO" & subjects$filter_size == 5000L & !subjects$random_slope, ]
  for (i in seq_len(nrow(groups))) {
    group <- groups[i, , drop = FALSE]
    values <- subjects[match_group(subjects, group), ]
    fold <- metrics[match_group(metrics, group), ]
    person <- aggregate(values[c("loss", "calibration")], by = values["subject_id"], FUN = mean)
    for (r in seq_len(5L)) {
      rows <- values[values$repetition == r, ]
      runs <- fold[fold$repetition == r, ]
      repeated[[length(repeated) + 1L]] <- data.frame(
        group, repetition = r, registered_subjects = 129L, scored_subjects = nrow(rows), completed_folds = nrow(runs),
        pinball_loss = if (nrow(rows)) mean(rows$loss) else NA_real_,
        quantile_calibration = if (nrow(rows)) mean(rows$calibration) else NA_real_,
        intercept_pinball_loss = if (nrow(runs)) stats::weighted.mean(runs$intercept_pinball_loss, runs$n_test_subjects) else NA_real_,
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
    base_loss <- stats::weighted.mean(fold$intercept_pinball_loss, fold$n_test_subjects)
    overall[[length(overall) + 1L]] <- data.frame(
      group, registered_subjects = 129L, scored_subjects = nrow(person), completed_outer_jobs = nrow(fold),
      scheduled_outer_jobs = 25L, pinball_loss = mean(person$loss),
      quantile_calibration = mean(person$calibration), intercept_pinball_loss = base_loss,
      intercept_skill = 1 - mean(person$loss) / base_loss,
      mean_selected_genes = mean(fold$selected_genes), median_selected_genes = stats::median(fold$selected_genes),
      mean_runtime_sec = mean(fold$runtime_sec), known_converged_jobs = sum(fold$converged %in% TRUE),
      known_nonconverged_jobs = sum(fold$converged %in% FALSE), convergence_unreported_jobs = sum(is.na(fold$converged)),
      uncertainty_scope = "descriptive_heldout_subjects_and_repeated_subject_folds", stringsAsFactors = FALSE, row.names = NULL
    )
    reference <- baseline[baseline$model == group$model & baseline$tau == group$tau, ]
    if (nrow(reference) && !(group$method == "POOLED-QR-LASSO" && group$filter_size == 5000L && !group$random_slope)) {
      comparison <- merge(values[c("repetition", "subject_id", "loss")],
                            reference[c("repetition", "subject_id", "loss")],
                            by = c("repetition", "subject_id"), suffixes = c("_method", "_pooled"))
      comparison$difference <- comparison$loss_method - comparison$loss_pooled
      if (!nrow(comparison)) next
      patient <- aggregate(difference ~ subject_id, comparison, mean)
      paired[[length(paired) + 1L]] <- data.frame(
        group, registered_subjects = 129L, paired_subjects = nrow(patient),
        mean_difference = mean(patient$difference), median_difference = stats::median(patient$difference),
        difference_q25 = stats::quantile(patient$difference, .25, names = FALSE),
        difference_q75 = stats::quantile(patient$difference, .75, names = FALSE),
        subjects_with_lower_loss = sum(patient$difference < 0),
        reference = "POOLED-QR-LASSO_filter5000_random_intercept", stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }
  write.csv(bind_tables(repeated), file.path(output, "repetition_level_prediction_summary.csv"), row.names = FALSE)
  write.csv(bind_tables(overall), file.path(output, "application_performance_summary.csv"), row.names = FALSE)
  write.csv(bind_tables(paired), file.path(output, "paired_subject_pinball_differences.csv"), row.names = FALSE)
}
execution <- yaml::read_yaml(file.path(workload, "application_execution.yml"))
audit <- file.path(root, execution$cohort_directory)
audit_files <- c(cohort_flow = "cohort_flow.csv", outcome_audit = "outcome_audit_summary.csv",
                 outcome_audit_gates = "outcome_audit_gates.csv", visit_distribution = "visit_distribution.csv",
                 missingness_by_variable = "missingness_by_variable.csv", batch_qc = "batch_outcome_summary.csv",
                 local_curvature = "local_curvature.csv", subject_folds = "subject_folds.csv")
for (name in names(audit_files)) {
  source <- file.path(audit, audit_files[[name]])
  if (!file.exists(source)) stop("An actual frozen cohort audit table is missing: ", source)
  file.copy(source, file.path(output, paste0(name, ".csv")), overwrite = FALSE)
}
exploratory_run <- unique(registries$exploratory$run_id)
if (length(exploratory_run) != 1L) stop("The shared exploratory run identity changed.")
exploratory <- file.path(root, "results", "raw", "application", exploratory_run)
for (name in c("exploratory_role_metrics.csv", "exploratory_role_intervals.csv", "exploratory_selection_frequency.csv",
               "exploratory_coefficient_median_iqr.csv")) {
  source <- file.path(exploratory, name)
  if (file.exists(source)) file.copy(source, file.path(output, name), overwrite = FALSE)
}
jsonlite::write_json(list(
  stage = args[3], workload_directory = workload, registered_jobs = 1210L,
  successful_jobs = sum(accounting$status == "ok"), failed_jobs = sum(accounting$status == "failed"),
  incomplete_jobs = sum(!accounting$status %in% c("ok", "failed")),
  observed_numeric_certificates = unique(accounting$numeric_certificate_commit[!is.na(accounting$numeric_certificate_commit)]),
  summary_source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "application_summary_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

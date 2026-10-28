library(dplyr)

output <- "results/tables/joint_profile_manuscript_20261026"
generated <- "manuscript/generated/joint_profile"
root <- "results/raw/application_prediction/joint_profile_20261026"
routes <- c("reference", "official_R1_R2_normexp_training_quantile", "processed_median_across_mapped_probes", "reference_conic_slope")
records <- list()
fits <- list()
for (route in routes) {
  paths <- list.files(file.path(root, route), "^result[.]rds$", recursive = TRUE, full.names = TRUE)
  for (path in paths) {
    result <- readRDS(path)
    if (route == "reference" && result$job$random_slope) next
    records[[length(records) + 1L]] <- result$predictions
    job <- result$job
    job$processing <- result$predictions$processing[1L]
    job$repetition <- result$predictions$repetition[1L]
    job$fold <- result$predictions$fold[1L]
    job$runtime <- result$fit$runtime
    job$selected_size <- result$fit$selected$size
    job$selected_lambda <- result$fit$selected$lambda
    job$final_gap <- result$fit$final_gap
    job$source_file <- path
    fits[[length(fits) + 1L]] <- job
  }
}
if (length(records) != 900L) stop(paste("Prediction fits completed", length(records), "of 900"))
predictions <- bind_rows(records)
fits <- bind_rows(fits)
keys <- c("processing", "model", "tau", "filter_size", "random_slope")
if (anyDuplicated(predictions[c(keys, "repetition", "sample_id")]) ||
    anyDuplicated(fits[c(keys, "repetition", "fold")])) stop("Duplicate final prediction records.")
predictions$calibrated <- predictions$observed_log_c3 <= predictions$predicted_quantile
subject <- predictions %>% group_by(across(all_of(c(keys, "repetition", "subject_id")))) %>%
  summarise(loss = mean(loss), calibration = mean(calibrated), visits = n(), .groups = "drop")
subject_vectors <- subject %>% group_by(across(all_of(c(keys, "subject_id")))) %>%
  summarise(loss = mean(loss), calibration = mean(calibration), repetitions = n(), .groups = "drop")
if (any(subject_vectors$repetitions != 5L)) stop("A subject lacks a prediction repetition.")
summary <- subject_vectors %>% group_by(across(all_of(keys))) %>%
  summarise(subjects = n(), mean_loss = mean(loss), calibration = mean(calibration), .groups = "drop")
if (any(summary$subjects != 129L) || nrow(summary) != 36L) stop("Prediction cohort sizes differ.")
timings <- fits %>% group_by(across(all_of(keys))) %>% summarise(fits = n(), runtime_median = median(runtime),
  runtime_maximum = max(runtime), selected_size_median = median(selected_size), maximum_gap = max(final_gap), .groups = "drop")
summary <- summary %>% left_join(timings, by = keys)

set.seed(20261026L)
subjects <- sort(unique(subject$subject_id))
weights <- rmultinom(10000L, length(subjects), rep(1 / length(subjects), length(subjects)))
saveRDS(list(subject_id = subjects, weights = weights, seed = 20261026L), file.path(output, "prediction_subject_weights.rds"))
reference <- subject_vectors %>% filter(processing == "reference", filter_size == 5000L, !random_slope) %>%
  select(model, tau, subject_id, reference_loss = loss)
paired <- subject_vectors %>% left_join(reference, by = c("model", "tau", "subject_id")) %>%
  mutate(difference = loss - reference_loss)
paired_summary <- paired %>% group_by(across(all_of(keys))) %>% group_modify(function(data, key) {
  values <- data$difference[match(subjects, data$subject_id)]
  if (anyNA(values)) stop("A paired prediction subject is missing.")
  draws <- drop(crossprod(values, weights)) / length(subjects)
  data.frame(difference = mean(values), lower = quantile(draws, .025), upper = quantile(draws, .975), subjects = length(values))
}) %>% ungroup()

old_root <- "results/tables/application/paired_loss_readback_committed_d83f62b_v1"
baselines <- read.csv(file.path(old_root, "all_available_heldout_predictions.csv")) %>%
  filter(method %in% c("POOLED-QR-LASSO", "QGEE-SCAD", "LQMM"),
         processing == "processed_highest_training_mad_reference", filter_size == 5000L, !random_slope)
available <- predictions %>% filter(processing == "reference", filter_size == 5000L, !random_slope) %>%
  select(sample_id, subject_id, repetition, fold, model, tau, proposed_loss = loss, observed_log_c3)
comparison <- baselines %>% inner_join(available, by = c("sample_id", "subject_id", "repetition", "fold", "model", "tau"), suffix = c("", "_new"))
if (nrow(comparison) != nrow(baselines) || any(abs(comparison$observed_log_c3 - comparison$observed_log_c3_new) > 1e-12)) stop("Saved comparator observations differ.")
paired_baseline <- comparison %>% group_by(method, model, tau, subject_id, repetition) %>%
  summarise(difference = mean(proposed_loss - loss), .groups = "drop") %>%
  group_by(method, model, tau, subject_id) %>% summarise(repetitions = n(), difference = mean(difference), .groups = "drop")
baseline_weights <- readRDS(file.path(old_root, "subject_resampling_weights.rds"))[[1L]]
baseline_summary <- paired_baseline %>% group_by(method, model, tau) %>% group_modify(function(data, key) {
  values <- data$difference[match(baseline_weights$subject_id, data$subject_id)]
  if (anyNA(values) || nrow(data) != length(values)) stop("Saved comparator resampling subjects differ.")
  draws <- drop(crossprod(values, baseline_weights$weights)) / length(values)
  data.frame(subjects = nrow(data), repetitions = min(data$repetitions), difference = mean(values),
             lower = quantile(draws, .025), upper = quantile(draws, .975))
}) %>% ungroup()
baseline_status <- read.csv(file.path(old_root, "all_2250_validation_job_accounting.csv")) %>%
  filter(method %in% c("POOLED-QR-LASSO", "QGEE-SCAD", "LQMM"),
         processing == "processed_highest_training_mad_reference")
accounting <- baseline_status %>% count(method, status, reported_fit_status, reported_converged, name = "fits")
returned_status <- baseline_status %>% filter(status == "ok") %>%
  select(method, model, tau, reported_fit_status, reported_converged)
baseline_summary <- baseline_summary %>% left_join(returned_status, by = c("method", "model", "tau"))
if (nrow(baseline_summary) != 18L) stop("Comparator fit records differ.")
write.csv(baseline_summary, file.path(output, "returned_comparator_records.csv"), row.names = FALSE)
baseline_summary$fit_status <- ifelse(is.na(baseline_summary$reported_converged), "Returned",
                                     ifelse(baseline_summary$reported_converged, "Converged", "Unconverged"))
unconverged <- !is.na(baseline_summary$reported_converged) & !baseline_summary$reported_converged
baseline_summary[unconverged, c("difference", "lower", "upper")] <- NA_real_
write.csv(predictions, file.path(output, "prediction_records.csv"), row.names = FALSE)
write.csv(fits, file.path(output, "prediction_fits.csv"), row.names = FALSE)
write.csv(subject_vectors, file.path(output, "prediction_subject_vectors.csv"), row.names = FALSE)
write.csv(summary, file.path(output, "prediction_summary.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(output, "prediction_sensitivity_differences.csv"), row.names = FALSE)
write.csv(baseline_summary, file.path(output, "available_comparator_differences.csv"), row.names = FALSE)
write.csv(accounting, file.path(output, "comparator_availability.csv"), row.names = FALSE)

fmt <- function(x) ifelse(is.na(x), "--", formatC(ifelse(abs(x) < .0005, 0, x), digits = 3L, format = "f"))
write_table <- function(table, caption, label, file, long = FALSE) {
  columns <- paste0("l", paste(rep("r", ncol(table) - 1L), collapse = ""))
  rows <- apply(table, 1L, function(x) paste0(paste(x, collapse = " & "), " \\\\"))
  header <- paste0(paste(names(table), collapse = " & "), " \\\\")
  if (long) {
    content <- c(paste0("\\begin{longtable}{", columns, "}"), paste0("\\caption{", caption, "}\\label{", label, "}\\\\"),
                 "\\toprule", header, "\\midrule\\endfirsthead", "\\toprule", header, "\\midrule\\endhead", rows,
                 "\\bottomrule", "\\end{longtable}")
  } else {
    content <- c("\\begin{table}[t]", paste0("\\caption{", caption, "}\\label{", label, "}"), "\\centering\\small",
                 paste0("\\begin{tabular}{", columns, "}"), "\\toprule", header, "\\midrule", rows,
                 "\\bottomrule", "\\end{tabular}", "\\end{table}")
  }
  writeLines(content, file.path(generated, paste0(file, ".tex")))
}
selected <- summary %>% filter(processing == "reference", filter_size == 5000L, !random_slope)
write_table(data.frame(Model = selected$model, Quantile = fmt(selected$tau), Loss = fmt(selected$mean_loss),
                       Calibration = fmt(selected$calibration), Genes = selected$selected_size_median,
                       Seconds = fmt(selected$runtime_median)),
            "Complete subject validation with five thousand candidate genes", "tab:joint-prediction", "prediction")
all <- left_join(summary, paired_summary, by = keys, suffix = c("", "_paired"))
all$route <- ifelse(all$processing == "reference", "Reference", ifelse(grepl("official", all$processing), "Raw arrays", "Median probes"))
write_table(data.frame(Processing = all$route, Model = all$model, Quantile = fmt(all$tau), Genes = all$filter_size,
                      Slope = ifelse(all$random_slope, "Yes", "No"), Loss = fmt(all$mean_loss),
                      Difference = fmt(all$difference), Range = paste0("[", fmt(all$lower), ", ", fmt(all$upper), "]")),
            "Prediction sensitivity differences against the reference random intercept learner", "tab:joint-prediction-all", "prediction_sensitivity", long = TRUE)
write_table(data.frame(Method = ifelse(baseline_summary$method == "POOLED-QR-LASSO", "Pooled Lasso", baseline_summary$method),
                      Model = baseline_summary$model, Quantile = fmt(baseline_summary$tau),
                      Subjects = baseline_summary$subjects, Status = baseline_summary$fit_status,
                      Difference = fmt(baseline_summary$difference),
                      Range = ifelse(unconverged, "--", paste0("[", fmt(baseline_summary$lower), ", ", fmt(baseline_summary$upper), "]"))),
            "Available paired comparator predictions on their recorded heldout subjects", "tab:joint-comparators", "available_comparators", long = TRUE)
macros <- c(paste0("\\newcommand{\\JointPredictionFits}{", nrow(fits), "}"),
            paste0("\\newcommand{\\JointPredictionSubjects}{", length(subjects), "}"),
            paste0("\\newcommand{\\JointComparatorSubjects}{", unique(baseline_summary$subjects), "}"),
            paste0("\\newcommand{\\JointComparatorUnconverged}{", sum(unconverged), "}"))
writeLines(macros, file.path(generated, "prediction_values.tex"))
jsonlite::write_json(list(method = "joint_profile_prediction_20261026", generated = as.character(Sys.time()),
  final_fits = nrow(fits), heldout_visit_records = nrow(predictions), settings = nrow(summary),
  subjects_per_setting = unique(summary$subjects), repetitions_per_subject = unique(subject_vectors$repetitions),
  random_slope_solver = "clarabel_selected_qp", random_slope_source = file.path(root, "reference_conic_slope"),
  numerical_attempt_failure_records = length(list.files(root, "^failure[.]rds$", recursive = TRUE)),
  comparator_accounting = accounting, comparator_nonconverged = sum(unconverged), source_tables = output,
  resampling_seed = 20261026L, resampling_replicates = ncol(weights), fitted_models = "held_fixed"),
  file.path(generated, "prediction_manifest.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
print(selected)
print(accounting)
cat("Completed prediction summaries for", nrow(fits), "fits and", nrow(predictions), "heldout visit records\n")

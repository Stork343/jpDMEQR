library(data.table)
output <- "results/tables/empirical_recovery_20261030/application"
collected <- readRDS(file.path(output, "records.rds"))
accounting <- as.data.table(collected$accounting)
predictions <- as.data.table(collected$predictions)
if (nrow(accounting) != 1350L) stop("The application comparison is incomplete.")
source_root <- "results/raw/application_prediction/joint_profile_20261026"
routes <- c("reference", "official_R1_R2_normexp_training_quantile", "processed_median_across_mapped_probes", "reference_conic_slope")
own <- list()
own_jobs <- list()
for (route in routes) {
  paths <- list.files(file.path(source_root, route), "^result[.]rds$", recursive = TRUE, full.names = TRUE)
  for (path in paths) {
    result <- readRDS(path)
    if (route == "reference" && result$job$random_slope) next
    own[[length(own) + 1L]] <- as.data.table(result$predictions)
    own_jobs[[length(own_jobs) + 1L]] <- data.table(source_file = path,
      processing = result$predictions$processing[1L], model = result$job$model, tau = result$job$tau,
      repetition = result$predictions$repetition[1L], fold = result$predictions$fold[1L],
      filter_size = result$job$filter_size, random_slope = result$job$random_slope)
  }
}
own <- rbindlist(own)
own_jobs <- rbindlist(own_jobs)
if (nrow(own_jobs) != 900L) stop("The joint prediction study requires 900 successful fits.")
own_keys <- c("processing", "model", "tau", "filter_size", "random_slope", "repetition", "fold")
if (anyDuplicated(own_jobs[, ..own_keys]) || anyDuplicated(own[, c(own_keys, "sample_id"), with = FALSE]))
  stop("Duplicated joint prediction records.")
own[, calibrated := observed_log_c3 <= predicted_quantile]
own_subjects <- own[, .(visits = .N, loss = mean(loss), calibration = mean(calibrated)),
  by = .(processing, model, tau, filter_size, random_slope, repetition, subject_id)]
complete_own <- own_subjects[, .(repetitions = .N), by = .(processing, model, tau, filter_size, random_slope, subject_id)]
if (nrow(complete_own) != 36L * 129L || any(complete_own$repetitions != 5L)) stop("Joint prediction subjects are incomplete.")
own_reference <- own[filter_size == 5000L & !random_slope,
  .(processing, model, tau, repetition, fold, sample_id, subject_id,
    proposed_loss = loss, proposed_observed = observed_log_c3,
    proposed_calibration = calibrated)]
join_keys <- c("processing", "model", "tau", "repetition", "fold", "sample_id", "subject_id")
paired <- merge(predictions, own_reference, by = join_keys, all.x = TRUE, sort = FALSE)
if (nrow(paired) != nrow(predictions) || anyNA(paired$proposed_loss) ||
    any(abs(paired$proposed_observed - paired$observed_log_c3) > 1e-12)) stop("Paired prediction observations differ.")
paired[, `:=`(difference = proposed_loss - loss, calibrated = observed_log_c3 <= predicted_quantile)]
subjects <- paired[, .(difference = mean(difference), loss = mean(loss), calibration = mean(calibrated),
  proposed_loss = mean(proposed_loss), proposed_calibration = mean(proposed_calibration), visits = .N),
  by = .(processing, method, model, tau, subject_id, repetition)]
subject_vectors <- subjects[, .(repetitions = .N, visits = sum(visits), difference = mean(difference),
  loss = mean(loss), calibration = mean(calibration), proposed_loss = mean(proposed_loss),
  proposed_calibration = mean(proposed_calibration)), by = .(processing, method, model, tau, subject_id)]
if (any(subject_vectors$repetitions > 5L)) stop("A comparator subject has duplicated repetitions.")
settings <- unique(accounting[, .(processing, method, model, tau)])
setorder(settings, processing, method, model, tau)
records <- vector("list", 2L * nrow(settings))
weights_records <- vector("list", length(records))
for (i in seq_len(nrow(settings))) {
  setting <- settings[i]
  key_fields <- names(setting)
  available <- subject_vectors[setting, on = key_fields, nomatch = 0L]
  fitted <- accounting[setting, on = key_fields]
  for (j in 1:2) {
    data <- if (j == 1L) available else available[repetitions == 5L]
    setorder(data, subject_id)
    index <- 2L * (i - 1L) + j
    row <- copy(setting)
    row[, `:=`(scope = if (j == 1L) "available_repetitions" else "five_repetitions",
      attempted_fits = nrow(fitted), successful_fits = sum(fitted$successful), subjects = nrow(data),
      subject_repetitions = sum(data$repetitions),
      minimum_repetitions = if (nrow(data)) min(data$repetitions) else 0L,
      maximum_repetitions = if (nrow(data)) max(data$repetitions) else 0L)]
    values <- c(loss = NA_real_, calibration = NA_real_, proposed_loss = NA_real_, proposed_calibration = NA_real_,
      difference = NA_real_, lower = NA_real_, upper = NA_real_)
    if (nrow(data)) {
      values[c("loss", "calibration", "proposed_loss", "proposed_calibration", "difference")] <-
        vapply(data[, .(loss, calibration, proposed_loss, proposed_calibration, difference)], mean, numeric(1))
      seed <- 202610280L + index
      set.seed(seed)
      weights <- rmultinom(10000L, nrow(data), rep(1 / nrow(data), nrow(data)))
      draws <- drop(crossprod(data$difference, weights)) / nrow(data)
      values[c("lower", "upper")] <- quantile(draws, c(.025, .975), names = FALSE)
      weights_records[[index]] <- list(setting = setting, scope = row$scope, subject_id = data$subject_id,
        seed = seed, weights = weights)
    }
    records[[index]] <- cbind(row, as.data.table(as.list(values)))
  }
}
summary <- rbindlist(records)
fwrite(subjects, file.path(output, "paired_subject_repetitions.csv"))
fwrite(subject_vectors, file.path(output, "paired_subject_vectors.csv"))
fwrite(summary, file.path(output, "paired_comparison_summary.csv"))
fwrite(own_jobs, file.path(output, "joint_prediction_jobs.csv"))
saveRDS(weights_records, file.path(output, "subject_resampling.rds"), compress = "xz")
jsonlite::write_json(list(registered_comparator_jobs = nrow(accounting), successful_comparator_jobs = sum(accounting$successful),
  registered_joint_jobs = nrow(own_jobs), comparator_settings = nrow(settings),
  resampling_replicates = 10000L, resampling_unit = "subject", fitted_models = "held_fixed",
  difference_definition = "joint_profile_loss_minus_comparator_loss", generated = as.character(Sys.time())),
  file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
print(summary[processing == "reference" & scope == "available_repetitions"])

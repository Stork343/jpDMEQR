library(yaml)
library(jsonlite)
library(data.table)
library(matrixStats)
source("research/clinical_search_20261105/panel.R")
source("research/empirical_selection_20261101/prediction_core.R")
config <- read_yaml("config/clinical_quantile_20261105.yml")
records <- checks <- list()
for (path in unlist(config$source_studies)) {
  study <- read_yaml(path)
  root <- file.path(config$output_root, study$study_id)
  status <- read_json(file.path(root, "status.json"), simplifyVector = TRUE)
  if (status$status != "completed") next
  panel <- cp_panel(study$input_directory)
  for (fold in seq_len(study$outer_folds)) {
    rows <- fread(file.path(study$output_root, paste0("outer_", fold), "evaluation.csv"))
    stopifnot(all(rows$role == "discovery_evaluation"))
    predictions <- matrix(NA_real_, nrow(rows), length(config$quantiles))
    for (j in seq_along(config$quantiles)) {
      tau <- unlist(config$quantiles)[j]
      directory <- file.path(root, paste0("outer_", fold), sprintf("tau%03d", as.integer(100 * tau)))
      fit <- readRDS(file.path(directory, "fit.rds"))
      stopifnot(!length(intersect(fit$recipe$training_strains, rows$strain_id)))
      features <- cp_features(rows, panel, fit$recipe)
      reconstructed <- drop(features$X[, names(fit$beta), drop = FALSE] %*% fit$beta) *
        fit$recipe$y_scale + fit$recipe$y_center
      stored <- fread(file.path(directory, "prediction.csv"))
      stopifnot(identical(stored$sample_id, rows$sample_id))
      error <- max(abs(reconstructed - stored$prediction))
      stopifnot(error < 1e-8)
      predictions[, j] <- reconstructed
      checks[[length(checks) + 1L]] <- data.table(study_id = study$study_id, fold, tau,
        reconstruction_error = error, patients = uniqueN(rows$strain_id))
    }
    sorted <- t(apply(predictions, 1L, sort))
    sorted[sorted < 0] <- 0
    for (j in seq_along(config$quantiles)) {
      tau <- unlist(config$quantiles)[j]
      for (rule in c("raw", "nonnegative_sorted_quantiles")) {
        value <- if (rule == "raw") predictions[, j] else sorted[, j]
        records[[length(records) + 1L]] <- data.table(study_id = study$study_id, fold, tau, rule,
          strain_id = rows$strain_id, loss = ep_loss(rows$y, value, tau))
      }
    }
  }
}
patient <- rbindlist(records)[, .(loss = mean(loss)), by = .(study_id, tau, rule, strain_id)]
summary <- patient[, .(loss = mean(loss), patients = .N), by = .(study_id, tau, rule)]
output <- "results/preflight/clinical_quantile_20261105"
saved <- fread(file.path(output, "discovery_losses.csv"))
comparison <- merge(summary, saved, by = c("study_id", "tau", "rule"), suffixes = c("_R", "_Python"))
stopifnot(nrow(comparison) == nrow(summary), all(comparison$patients_R == comparison$patients_Python))
comparison[, loss_difference := abs(loss_R - loss_Python)]
stopifnot(max(comparison$loss_difference) < 1e-10)
fwrite(rbindlist(checks), file.path(output, "prediction_reconstruction_checks.csv"))
fwrite(comparison, file.path(output, "independent_loss_checks.csv"))
write_json(list(status = "passed", completed_studies = unique(summary$study_id),
  maximum_reconstruction_error = max(rbindlist(checks)$reconstruction_error),
  maximum_loss_difference = max(comparison$loss_difference)), file.path(output, "prediction_acceptance.json"),
  auto_unbox = TRUE, pretty = TRUE)
print(comparison)

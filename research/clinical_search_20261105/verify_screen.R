library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/clinical_search_20261105/panel.R")
config <- read_yaml(commandArgs(trailingOnly = TRUE)[1L])
root <- config$output_root
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
panel <- cp_panel(config$input_directory)
jobs <- fread(file.path(root, "progress.csv"))
stopifnot(all(jobs$status == "completed"), all(jobs$exit_code[!is.na(jobs$exit_code)] == 0L))
records <- checks <- list()
for (j in seq_len(nrow(jobs))) {
  job <- jobs[j]
  directory <- job$directory
  fit <- readRDS(file.path(directory, "fit.rds"))
  rows <- fread(file.path(root, job$trait, "evaluation.csv"))
  stopifnot(!length(intersect(rows$strain_id, fit$recipe$training_strains)))
  features <- if (fit$offset_model) cp_offset_features(rows, panel, fit$recipe) else cp_features(rows, panel, fit$recipe)
  if (job$method == "CLINICAL_QR") {
    columns <- rownames(fit$beta)
    stopifnot(length(columns) == nrow(fit$beta), all(columns %in% colnames(features$X)))
    prediction <- features$X[, columns, drop = FALSE] %*% fit$beta * fit$recipe$y_scale + fit$recipe$y_center
  } else {
    base <- drop(features$X %*% fit$beta) * fit$recipe$y_scale + fit$recipe$y_center
    if (fit$offset_model) base <- base + features$offset
    prediction <- vapply(fit$calibration, function(value) base + value[rows$condition], numeric(nrow(rows)))
  }
  error <- max(abs(prediction - fit$prediction))
  stopifnot(is.finite(error), error < 1e-10)
  expected <- rbindlist(lapply(seq_along(fit$quantiles), function(k)
    data.table(sample_id = rows$sample_id, strain_id = rows$strain_id, method = job$method,
      tau = fit$quantiles[k], prediction = prediction[, k], y = rows$y)))
  saved <- fread(file.path(directory, "prediction.csv"))
  paired <- merge(expected, saved, by = c("sample_id", "strain_id", "method", "tau"), suffixes = c("_new", "_saved"))
  serialized_error <- max(abs(paired$prediction_new - paired$prediction_saved))
  stopifnot(nrow(paired) == nrow(expected), serialized_error < 1e-10)
  records[[j]] <- expected
  checks[[j]] <- data.table(run_id = basename(directory), outer_fold = job$trait,
    prediction_error = error, serialized_error = serialized_error)
}
prediction <- rbindlist(records)
stopifnot(!anyDuplicated(prediction[, .(sample_id, method, tau)]))
raw <- copy(prediction)
raw[, scoring := "raw"]
primary <- copy(prediction)
setorder(primary, method, sample_id, tau)
primary[, prediction := sort(pmax(prediction, 0)), by = .(method, sample_id)]
primary[, scoring := "primary"]
all <- rbindlist(list(raw, primary))
all[, residual := y - prediction]
all[, loss := abs(residual) * fifelse(residual >= 0, tau, 1 - tau)]
patient <- all[, .(loss = mean(loss)), by = .(scoring, method, tau, strain_id)]
summary <- patient[, .(loss = mean(loss), se = sd(loss) / sqrt(.N)), by = .(scoring, method, tau)]
saved <- fread(file.path(config$report_directory, "discovery_loss.csv"))
paired <- merge(summary, saved, by = c("scoring", "method", "tau"), suffixes = c("_new", "_saved"))
error <- max(abs(paired$loss_new - paired$loss_saved), abs(paired$se_new - paired$se_saved))
stopifnot(error < 1e-10, nrow(paired) == nrow(summary))
fwrite(rbindlist(checks), file.path(config$report_directory, "prediction_reconstruction_checks.csv"))
write_json(list(status = "passed", models_checked = nrow(jobs),
  maximum_prediction_error = max(rbindlist(checks)$prediction_error),
  independent_loss_error = error, formal_validation_used = FALSE),
  file.path(config$report_directory, "screen_audit.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("Checked", nrow(jobs), "models; loss error", error, "\n")

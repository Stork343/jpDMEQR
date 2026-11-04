library(data.table)
library(jsonlite)
library(yaml)
library(BEDMatrix)
source("research/yeast_1011_application_20261104/features.R")
source("research/yeast_stress_application_20261104/features.R")
source("research/empirical_selection_20261101/prediction_core.R")

config <- read_yaml("config/discovery_yeast_stress_20261104.yml")
queue <- fread("results/raw/application/yeast_stress_20261104/exact_comparisons_20261104/progress.csv")
stopifnot(nrow(queue) > 0L, all(queue$status == "completed"))
panel <- y1011_panel(config$genomic_directory)
screen <- fread("results/preflight/empirical_selection_20261101/yeast_stress_prediction/comparisons.csv")
checks <- comparisons <- list()
set.seed(20261105L)
for (j in seq_len(nrow(queue))) {
  assay_name <- queue$trait[j]
  directory <- queue$directory[j]
  evaluation <- fread(file.path(config$input_directory, assay_name, "discovery.csv"))[exploration_role == "evaluation"]
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  stopifnot(!request$formal_validation_used, all(evaluation$role == "discovery"))
  fit <- readRDS(file.path(directory, "fit.rds"))
  stopifnot(!length(intersect(fit$recipe$training_strains, evaluation$strain_id)))
  features <- stress_features(evaluation, panel, fit$recipe)
  prediction <- drop(ep_predict(features, fit$beta)) + fit$calibration[[1L]][evaluation$condition]
  saved <- fread(file.path(directory, "predictions.csv"))
  stopifnot(identical(saved$sample_id, evaluation$sample_id), all(saved$tau == 0.25))
  error <- max(abs(saved$y - evaluation$y), abs(saved$prediction - prediction),
               abs(saved$loss - ep_loss(evaluation$y, prediction, 0.25)))
  stopifnot(is.finite(error), error < 1e-9)
  checks[[j]] <- data.table(assay = assay_name, method = "POOLED", maximum_error = error)
  pooled <- saved[, .(pooled = mean(loss)), by = strain_id]
  joint <- fread(file.path(config$output_root, assay_name, "JOINT_prediction_tau025", "strain_loss.csv"))
  paired <- merge(pooled, joint[, .(strain_id, joint = loss)], by = "strain_id")
  stopifnot(nrow(paired) == 148L, !anyNA(paired))
  weights <- rmultinom(10000L, nrow(paired), rep(1 / nrow(paired), nrow(paired))) / nrow(paired)
  difference <- paired$pooled - paired$joint
  bootstrap <- drop(crossprod(weights, difference)) / drop(crossprod(weights, paired$pooled))
  gain <- mean(difference) / mean(paired$pooled)
  mean_gain <- screen[assay == assay_name]$relative_gain
  stopifnot(length(mean_gain) == 1L)
  comparisons[[j]] <- data.table(assay = assay_name, tau = 0.25, strains = nrow(paired),
    joint_loss = mean(paired$joint), comparator_loss = mean(paired$pooled),
    relative_gain = gain, lower = quantile(bootstrap, 0.025), upper = quantile(bootstrap, 0.975),
    forest_comparisons_required = gain > 0 && mean_gain > 0, manuscript_adoption = FALSE)
  rm(features, fit)
  gc()
}
output <- "results/preflight/empirical_selection_20261101/yeast_stress_exact"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(comparisons), file.path(output, "comparisons.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
write_json(list(stage = "discovery", completed_runs = nrow(queue), formal_validation_used = FALSE,
                manuscript_adoption = FALSE), file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(comparisons))

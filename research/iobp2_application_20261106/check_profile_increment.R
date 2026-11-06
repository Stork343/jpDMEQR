library(data.table)
library(glmnet)
library(jsonlite)
library(yaml)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
study <- read_yaml(args[1L])
root <- study$output_root
implementation <- file.path(root, "implementation")
source(file.path(implementation, "panel.R"))
source(file.path(implementation, "prediction_core.R"))
panel <- cp_panel(study$input_directory)
panel$X <- panel$X[, startsWith(colnames(panel$X), "C_"), drop = FALSE]
directory <- file.path(root, "profile_increment")
dir.create(directory, showWarnings = FALSE)
progress <- function(stage, details) cat(stage, unlist(details), "\n")
records <- list()
for (fold in seq_len(study$outer_folds)) {
  outer <- file.path(root, paste0("outer_", fold))
  rows <- fread(file.path(outer, "training.csv"))
  evaluation <- fread(file.path(outer, "evaluation.csv"))
  stopifnot(all(evaluation$role == "discovery_evaluation"),
    !length(intersect(rows$strain_id, evaluation$strain_id)))
  fit_dir <- file.path(directory, paste0("outer_", fold))
  dir.create(fit_dir, showWarnings = FALSE)
  path <- file.path(fit_dir, "fit.rds")
  if (!file.exists(path)) {
    fit <- ep_fit_prediction(rows, evaluation, panel, "ENET", NA_real_, study, progress, fit_dir,
      cp_recipe, cp_features, "strain_id", "condition", ep_predict, "oof")
    saveRDS(fit, path)
  }
  clinical <- readRDS(path)
  full <- readRDS(file.path(outer, "models", "ENET_all", "fit.rds"))
  stopifnot(identical(sort(unique(rows$strain_id)), clinical$recipe$training_strains),
    identical(clinical$quantiles, full$quantiles))
  scores <- function(fit) {
    pred <- t(apply(fit$prediction, 1L, function(x) sort(pmax(0, x))))
    index <- match(study$decision_quantile, fit$quantiles)
    ep_loss(evaluation$y, pred[, index], study$decision_quantile)
  }
  records[[fold]] <- data.table(strain_id = evaluation$strain_id, fold, clinical_loss = scores(clinical),
    curve_loss = scores(full))
}
patient <- rbindlist(records)[, .(clinical_loss = mean(clinical_loss), curve_loss = mean(curve_loss)), by = .(strain_id, fold)]
stopifnot(!anyDuplicated(patient$strain_id))
result <- patient[, .(patients = .N, clinical_loss = mean(clinical_loss), curve_loss = mean(curve_loss))]
result[, relative_gain_with_curve := 1 - curve_loss / clinical_loss]
stored <- fread(file.path(study$report_directory, "discovery_loss.csv"))[scoring == "primary" &
  method == "ENET" & tau == study$decision_quantile]
stopifnot(nrow(stored) == 1L, abs(result$curve_loss - stored$loss) < 1e-10)
fwrite(patient, file.path(directory, "patient_loss.csv"))
fwrite(result, file.path(study$report_directory, "exploratory_profile_increment.csv"))
print(result)

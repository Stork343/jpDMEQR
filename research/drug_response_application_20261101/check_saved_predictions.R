library(data.table)
library(jsonlite)
source("research/drug_response_application_20261101/data.R")
source("research/drug_response_application_20261101/fbio_adapter.R")
source("research/drug_response_application_20261101/prediction.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
data <- readRDS("data/derived/drug_response_20261101/discovery_cases.rds")
expression <- readRDS("data/derived/drug_response_20261101/expression.rds")
records <- list()
for (directory in args) {
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  fit <- readRDS(file.path(directory, "fit.rds"))
  evaluation <- data[DRUG_NAME == request$drug & analysis_curve & exploration_role == "evaluation"]
  stopifnot(!length(intersect(fit$recipe$training_patients, evaluation$patient_id)))
  if (request$method == "JOINT") {
    features <- drug_features(evaluation, expression, fit$recipe)
    prediction <- matrix(drop(features$X %*% fit$beta) +
      drug_apply_calibration(fit$calibration, evaluation$CONC), ncol = 1L)
    quantiles <- request$tau
  } else {
    features <- drug_features(evaluation, expression, fit$recipe, expand = FALSE)
    expected_mean <- drug_fbio_predict(fit$model, features$genes, features$dose, features$clinical)
    prediction <- vapply(fit$calibration, function(calibration) expected_mean +
      drug_apply_calibration(calibration, evaluation$CONC), numeric(nrow(evaluation)))
    quantiles <- unlist(request$config$quantiles)
  }
  saved <- fread(file.path(directory, "predictions.csv"))
  errors <- vapply(seq_along(quantiles), function(j) {
    current <- saved[tau == quantiles[j]]
    stopifnot(identical(current$patient_id, evaluation$patient_id),
              max(abs(current$y - evaluation$normalized_intensity)) < 1e-12)
    max(abs(current$prediction - prediction[, j]),
        abs(current$loss - drug_loss(evaluation$normalized_intensity, prediction[, j], quantiles[j])))
  }, numeric(1))
  stopifnot(max(errors) < 1e-10)
  record <- list(status = "prediction_artifacts_verified", directory = directory,
                   maximum_error = max(errors), patients = uniqueN(evaluation$patient_id),
                   training_overlap = 0L, formal_validation_used = FALSE)
  write_json(record, file.path(directory, "artifact_checks.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  records[[length(records) + 1L]] <- as.data.table(record)
  rm(features, fit, prediction)
  gc()
}
print(rbindlist(records))

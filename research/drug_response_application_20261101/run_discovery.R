library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/drug_response_application_20261101/data.R")
source("research/drug_response_application_20261101/fbio_adapter.R")
source("research/drug_response_application_20261101/prediction.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3L)
drug <- args[1L]
method <- args[2L]
run_id <- args[3L]
config <- read_yaml("config/discovery_gdsc_20261101.yml")
stopifnot(drug %in% config$drugs, method %in% c("JOINT", "FBIO"), grepl("^[a-zA-Z0-9_]+$", run_id))
tau <- if (method == "JOINT") as.numeric(args[4L]) else NA_real_
if (method == "JOINT") stopifnot(tau %in% config$quantiles)
directory <- file.path("results/raw/application/drug_response_20261101/discovery", drug,
                       paste0(method, "_", run_id))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, detail = list()) {
  record <- c(list(stage = stage, drug = drug, method = method,
                   time = format(Sys.time(), tz = "UTC", usetz = TRUE)), detail)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(detail), collapse = " "), "\n")
  flush.console()
}
progress("loading")
data <- readRDS("data/derived/drug_response_20261101/discovery_cases.rds")
data <- data[DRUG_NAME == drug & analysis_curve == TRUE]
stopifnot(all(data$role == "discovery"))
training <- data[exploration_role == "training"]
evaluation <- data[exploration_role == "evaluation"]
stopifnot(!length(intersect(training$patient_id, evaluation$patient_id)))
expression <- readRDS("data/derived/drug_response_20261101/expression.rds")
request <- list(drug = drug, method = method, run_id = run_id, tau = tau, config = config,
                 code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                 training_patients = uniqueN(training$patient_id),
                 evaluation_patients = uniqueN(evaluation$patient_id),
                 formal_validation_used = FALSE, manuscript_adoption = FALSE)
write_json(request, file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
started <- proc.time()[["elapsed"]]
if (method == "JOINT") {
  source("research/score_validity_20261026/load_joint.R")
  source("research/joint_profile_manuscript_20261026/conic_search.R")
  result <- drug_joint_prediction(training, evaluation, expression, config, tau, progress, directory)
  quantiles <- tau
  predicted <- matrix(result$prediction, ncol = 1L)
} else {
  result <- drug_fbio_prediction(training, evaluation, expression, config, progress, directory)
  quantiles <- unlist(config$quantiles)
  predicted <- result$prediction
}
elapsed <- proc.time()[["elapsed"]] - started
stopifnot(nrow(predicted) == nrow(evaluation), all(is.finite(predicted)))
saveRDS(result, file.path(directory, "fit.rds"))
losses <- rbindlist(lapply(seq_along(quantiles), function(j) {
  data.table(patient_id = evaluation$patient_id, concentration = evaluation$CONC,
             y = evaluation$normalized_intensity, prediction = predicted[, j], tau = quantiles[j],
             loss = drug_loss(evaluation$normalized_intensity, predicted[, j], quantiles[j]))
}))
fwrite(losses, file.path(directory, "predictions.csv"))
patient_loss <- losses[, .(loss = mean(loss)), by = .(patient_id, tau)]
fwrite(patient_loss, file.path(directory, "patient_loss.csv"))
summary <- patient_loss[, .(patients = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = tau]
summary[, `:=`(drug = drug, method = method, runtime_sec = elapsed)]
fwrite(summary, file.path(directory, "summary.csv"))
write_json(list(status = "completed", runtime_sec = elapsed,
                 formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
progress("completed", list(runtime_sec = elapsed))
print(summary)

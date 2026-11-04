library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L, args[1L] == "ALT")
method <- args[2L]
run_id <- args[3L]
implementation <- normalizePath(args[5L], mustWork = TRUE)
stopifnot(method %in% c("JOINT", "ENET", "POOLED"), grepl("^[A-Za-z0-9_]+$", run_id))
config <- read_yaml(file.path(implementation, "config.yml"))
tau <- if (method == "ENET") NA_real_ else as.numeric(args[4L])
if (method != "ENET") stopifnot(tau %in% config$quantiles)
source(file.path(implementation, "yeast_prediction.R"))
source(file.path(implementation, "yeast_comparators.R"))
source(file.path(implementation, "prediction.R"))
if (method == "ENET") source(file.path(implementation, "prediction_core.R"))
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
if (method == "POOLED") source("research/soynam_application_20261031/load_pooled_lp.R")
rows <- readRDS(config$input_growth)
stopifnot(all(rows$role == "discovery"), all(is.finite(rows$y)))
training <- rows[exploration_role == "training"]
evaluation <- rows[exploration_role == "evaluation"]
stopifnot(uniqueN(training$compound_id) >= 100L, uniqueN(evaluation$compound_id) >= 30L,
          !length(intersect(training$compound_id, evaluation$compound_id)))
panel <- readRDS(config$input_panel)
annotation <- fread(config$input_annotation)
panel$allowed_probes <- intersect(annotation[TYPE == config$probe_type, ID], colnames(panel$X))
stopifnot(length(panel$allowed_probes) >= 1000L)
directory <- file.path(config$output_root, "ALT", paste0(method, "_", run_id))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, method = method, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(method = method, tau = tau, config = config,
                implementation_directory = implementation,
                code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                training_compounds = uniqueN(training$compound_id), evaluation_compounds = uniqueN(evaluation$compound_id),
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
started <- proc.time()[["elapsed"]]
fit <- dm_fit_prediction(training, evaluation, panel, method, tau, config, progress, directory)
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(directory, "fit.rds"))
predictions <- rbindlist(lapply(seq_along(fit$quantiles), function(j)
  data.table(sample_id = evaluation$sample_id, compound_id = evaluation$compound_id,
    study_id = evaluation$study_id, y = evaluation$y, tau = fit$quantiles[j],
    prediction = fit$prediction[, j], loss = yeast_loss(evaluation$y, fit$prediction[, j], fit$quantiles[j]))))
fwrite(predictions, file.path(directory, "predictions.csv"))
losses <- predictions[, .(loss = mean(loss)), by = .(compound_id, tau)]
fwrite(losses, file.path(directory, "compound_loss.csv"))
summary <- losses[, .(compounds = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = tau]
summary[, `:=`(method = method, runtime_sec = elapsed)]
fwrite(summary, file.path(directory, "summary.csv"))
write_json(list(status = "completed", runtime_sec = elapsed, formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
progress("completed", list(runtime_sec = elapsed))
print(summary)

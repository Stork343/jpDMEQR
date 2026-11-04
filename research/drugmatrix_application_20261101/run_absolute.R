library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
method <- args[1L]
run_id <- args[2L]
implementation <- normalizePath(args[3L], mustWork = TRUE)
stopifnot(method %in% c("JOINT", "POOLED"), grepl("^[A-Za-z0-9_]+$", run_id))
config <- read_yaml(file.path(implementation, "config.yml"))
source(file.path(implementation, "prediction_core.R"))
source(file.path(implementation, "prediction.R"))
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
if (method == "POOLED") source("research/soynam_application_20261031/load_pooled_lp.R")
rows <- readRDS(config$input_growth)
stopifnot(all(rows$role == "discovery"), all(rows$alt_u_l > 0))
rows[, y := alt_u_l]
training <- rows[exploration_role == "training"]
evaluation <- rows[exploration_role == "evaluation"]
stopifnot(!length(intersect(training$compound_id, evaluation$compound_id)))
panel <- readRDS(config$input_panel)
annotation <- fread(config$input_annotation)
panel$allowed_probes <- intersect(annotation[TYPE == "DISCOVERY", ID], colnames(panel$X))
raw_recipe <- function(rows, panel, config) {
  recipe <- dm_recipe(rows, panel)
  weights <- 1 / as.numeric(table(rows$compound_id)[rows$compound_id])
  quantiles <- vapply(c(0.25, 0.5, 0.75), function(tau)
    as.numeric(quantreg::rq.fit(matrix(weights, ncol = 1L), rows$y * weights, tau = tau, method = "br")$coefficients), numeric(1))
  recipe$y_center <- quantiles[2L]
  recipe$y_scale <- (quantiles[3L] - quantiles[1L]) / (qnorm(0.75) - qnorm(0.25))
  stopifnot(is.finite(recipe$y_scale), recipe$y_scale > 0)
  recipe$normalization <- config$normalization
  recipe
}
directory <- file.path(config$output_root, paste0(method, "_", run_id))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, method = method, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(method = method, tau = config$decision_quantile, config = config,
                implementation_directory = implementation, code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                training_compounds = uniqueN(training$compound_id), evaluation_compounds = uniqueN(evaluation$compound_id),
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE)
started <- proc.time()[["elapsed"]]
fit <- ep_fit_prediction(training, evaluation, panel, method, config$decision_quantile, config,
                         progress, directory, raw_recipe, dm_features, "compound_id")
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(directory, "fit.rds"))
prediction <- drop(fit$prediction)
saved <- data.table(sample_id = evaluation$sample_id, compound_id = evaluation$compound_id,
  study_id = evaluation$study_id, y = evaluation$y, tau = config$decision_quantile,
  prediction = prediction, loss = ep_loss(evaluation$y, prediction, config$decision_quantile))
fwrite(saved, file.path(directory, "predictions.csv"))
losses <- saved[, .(loss = mean(loss)), by = compound_id]
fwrite(losses, file.path(directory, "compound_loss.csv"))
summary <- losses[, .(compounds = .N, loss = mean(loss), se = sd(loss) / sqrt(.N))]
summary[, `:=`(method = method, tau = config$decision_quantile, runtime_sec = elapsed)]
fwrite(summary, file.path(directory, "summary.csv"))
write_json(list(status = "completed", runtime_sec = elapsed, formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
progress("completed", list(runtime_sec = elapsed))
print(summary)

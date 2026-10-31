library(data.table)
library(jsonlite)
library(yaml)
library(glmnet)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
trait_name <- args[1L]
method <- args[2L]
run_id <- args[3L]
implementation_directory <- normalizePath(args[5L], mustWork = TRUE)
source(file.path(implementation_directory, "prediction.R"))
stopifnot(method %in% c("JOINT", "ENET"), grepl("^[A-Za-z0-9_]+$", run_id))
config <- read_yaml(file.path(implementation_directory, "config.yml"))
tau <- if (method == "JOINT") as.numeric(args[4L]) else NA_real_
if (method == "JOINT") stopifnot(tau %in% config$quantiles)
data <- readRDS("data/derived/yeast_20261101/growth_discovery.rds")
stopifnot(trait_name %in% data$trait, all(data$role == "discovery"))
data <- data[trait == trait_name]
stopifnot(uniqueN(data$trait) == 1L, identical(unique(data$trait), trait_name))
training <- data[exploration_role == "training"]
evaluation <- data[exploration_role == "evaluation"]
stopifnot(uniqueN(training$strain_id) >= 100L, uniqueN(evaluation$strain_id) >= 50L,
          !length(intersect(training$strain_id, evaluation$strain_id)))
panel <- readRDS("data/derived/yeast_20261101/panel.rds")
directory <- file.path("results/raw/application/yeast_20261101/discovery", trait_name, paste0(method, "_", run_id))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, trait = trait_name, method = method,
                   time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(trait = trait_name, method = method, tau = tau, config = config,
                 implementation_directory = implementation_directory,
                 training_strains = uniqueN(training$strain_id), evaluation_strains = uniqueN(evaluation$strain_id),
                 formal_validation_used = FALSE, manuscript_adoption = FALSE,
                 code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
           file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
started <- proc.time()[["elapsed"]]
if (method == "JOINT") {
  source("research/score_validity_20261026/load_joint.R")
  source("research/joint_profile_manuscript_20261026/conic_search.R")
  result <- yeast_profile_prediction(training, evaluation, panel, config, tau, progress, directory)
  quantiles <- tau
  predicted <- matrix(result$prediction, ncol = 1L)
} else {
  result <- yeast_enet_prediction(training, evaluation, panel, config, progress, directory)
  quantiles <- unlist(config$quantiles)
  predicted <- result$prediction
}
elapsed <- proc.time()[["elapsed"]] - started
stopifnot(all(is.finite(predicted)), nrow(predicted) == nrow(evaluation))
saveRDS(result, file.path(directory, "fit.rds"))
predictions <- rbindlist(lapply(seq_along(quantiles), function(j) {
  data.table(strain_id = evaluation$strain_id, replicate = evaluation$replicate, y = evaluation$y,
             tau = quantiles[j], prediction = predicted[, j], loss = yeast_loss(evaluation$y, predicted[, j], quantiles[j]))
}))
fwrite(predictions, file.path(directory, "predictions.csv"))
losses <- predictions[, .(loss = mean(loss)), by = .(strain_id, tau)]
fwrite(losses, file.path(directory, "strain_loss.csv"))
summary <- losses[, .(strains = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = tau]
summary[, `:=`(trait = trait_name, method = method, runtime_sec = elapsed)]
fwrite(summary, file.path(directory, "summary.csv"))
write_json(list(status = "completed", runtime_sec = elapsed, manuscript_adoption = FALSE),
           file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
progress("completed", list(runtime_sec = elapsed))
print(summary)

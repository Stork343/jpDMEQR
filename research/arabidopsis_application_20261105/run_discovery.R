library(data.table)
library(jsonlite)
library(yaml)
library(glmnet)
library(BEDMatrix)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L, args[1L] == "Flowering")
method <- args[2L]
run_id <- args[3L]
implementation <- normalizePath(args[5L], mustWork = TRUE)
config <- read_yaml(file.path(implementation, "config.yml"))
stopifnot(method %in% c("JOINT", "ENET", "POOLED"), grepl("^[A-Za-z0-9_]+$", run_id))
tau <- if (method == "ENET") NA_real_ else as.numeric(args[4L])
if (method != "ENET") stopifnot(tau %in% config$quantiles)
source(file.path(implementation, "prediction_core.R"))
source(file.path(implementation, "features.R"))
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
if (method == "POOLED") source("research/soynam_application_20261031/load_pooled_lp.R")
check <- read_json("results/preflight/empirical_selection_20261101/arabidopsis/input_checks.json", simplifyVector = TRUE)
stopifnot(check$status == "passed")
rows <- fread(file.path(config$input_directory, "discovery.csv"), colClasses = c(accession_id = "character"))
stopifnot(all(rows$role == "discovery"))
training <- rows[exploration_role == "training"]
evaluation <- rows[exploration_role == "evaluation"]
stopifnot(uniqueN(training$accession_id) == 509L, uniqueN(evaluation$accession_id) == 170L,
          !length(intersect(training$accession_id, evaluation$accession_id)))
panel <- arab_panel(config$input_directory)
directory <- file.path(config$output_root, "Flowering", paste0(method, "_", run_id))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, method = method, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(method = method, tau = tau, config = config, implementation_directory = implementation,
                code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
started <- proc.time()[["elapsed"]]
fit <- ep_fit_prediction(training, evaluation, panel, method, tau, config, progress, directory,
                         arab_recipe, arab_features, "accession_id", "condition")
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(directory, "fit.rds"))
predictions <- rbindlist(lapply(seq_along(fit$quantiles), function(j)
  data.table(sample_id = evaluation$sample_id, accession_id = evaluation$accession_id,
    condition = evaluation$condition, y = evaluation$y, tau = fit$quantiles[j],
    prediction = fit$prediction[, j], loss = ep_loss(evaluation$y, fit$prediction[, j], fit$quantiles[j]))))
fwrite(predictions, file.path(directory, "predictions.csv"))
losses <- predictions[, .(loss = mean(loss)), by = .(accession_id, tau)]
fwrite(losses, file.path(directory, "accession_loss.csv"))
summary <- losses[, .(accessions = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = tau]
summary[, `:=`(method = method, runtime_sec = elapsed)]
fwrite(summary, file.path(directory, "summary.csv"))
write_json(list(status = "completed", runtime_sec = elapsed, formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
capture.output(sessionInfo(), file = file.path(directory, "session.txt"))
progress("completed", list(runtime_sec = elapsed))
print(summary)

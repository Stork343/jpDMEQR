library(data.table)
library(jsonlite)
library(yaml)
library(grf)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
assay_name <- args[1L]
method <- args[2L]
run_id <- args[3L]
implementation <- normalizePath(args[5L], mustWork = TRUE)
config <- read_yaml(file.path(implementation, "config.yml"))
stopifnot(assay_name %in% config$assays, method %in% c("QRF_RAW", "QRF_OFFSET"), run_id == "native")
source(file.path(implementation, "prediction_core.R"))
source(file.path(implementation, "ridge_direction.R"))
source(file.path(implementation, "contrasts.R"))
source(file.path(implementation, "features.R"))
source(file.path(implementation, "forest.R"))
checks <- read_json("results/preflight/empirical_selection_20261101/ridge_direction/forest_checks.json", simplifyVector = TRUE)
stopifnot(checks$status == "passed")
screen <- fread(config$screening_decision)
stopifnot(nrow(screen[assay == assay_name]) == 1L, screen[assay == assay_name]$further_comparisons_required)
panel <- proteome_panel(config$input_directory)
recipes <- readRDS(file.path(config$source_recipes, assay_name, "recipes.rds"))
rows <- fread(file.path(config$input_directory, assay_name, "discovery.csv"))
stopifnot(all(rows$role == "discovery"))
training <- rows[exploration_role == "training"]
evaluation <- rows[exploration_role == "evaluation"]
stopifnot(!length(intersect(training$strain_id, evaluation$strain_id)))
directory <- file.path(config$output_root, assay_name, paste0(method, "_", run_id))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, method = method, assay = assay_name,
                   time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(method = method, config = config, phase = "method_development",
  implementation_directory = implementation, code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  formal_validation_used = FALSE, manuscript_adoption = FALSE),
  file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE)
started <- proc.time()[["elapsed"]]
fit <- rs_forest_prediction(training, evaluation, panel, recipes, method == "QRF_OFFSET", progress, directory)
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(directory, "fit.rds"))
predictions <- rbindlist(lapply(seq_along(fit$quantiles), function(j)
  data.table(sample_id = evaluation$sample_id, strain_id = evaluation$strain_id, condition = evaluation$condition,
    y = evaluation$y, tau = fit$quantiles[j], prediction = fit$prediction[, j],
    loss = ep_loss(evaluation$y, fit$prediction[, j], fit$quantiles[j]))))
fwrite(predictions, file.path(directory, "predictions.csv"))
losses <- predictions[, .(loss = mean(loss)), by = .(strain_id, tau)]
fwrite(losses, file.path(directory, "strain_loss.csv"))
summary <- losses[, .(strains = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = tau]
summary[, `:=`(method = method, assay = assay_name, runtime_sec = elapsed, phase = "method_development")]
fwrite(summary, file.path(directory, "summary.csv"))
capture.output(sessionInfo(), file = file.path(directory, "session.txt"))
write_json(list(status = "completed", runtime_sec = elapsed, formal_validation_used = FALSE,
                manuscript_adoption = FALSE), file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
progress("completed", list(runtime_sec = elapsed))
print(summary)

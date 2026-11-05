library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
library(grf)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L, args[1L] == "hemoglobin", args[2L] %in% c("QRF_RAW", "QRF_OFFSET"))
implementation <- normalizePath(args[5L], mustWork = TRUE)
config <- read_yaml(file.path(implementation, "config.yml"))
for (file in c("prediction_core.R", "ridge_direction.R", "features.R", "forest_base.R", "forest_train.R")) source(file.path(implementation, file))
rows <- fread(file.path(config$input_directory, "training.csv"))
panel <- ch_panel(config$input_directory)
recipes <- readRDS(file.path(config$output_root, "recipes.rds"))
directory <- file.path(config$output_root, "models", paste0(args[2L], "_all"))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(method = args[2L], stage = "training", config = config,
  code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  implementation_directory = implementation, formal_validation_used = FALSE),
  file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE)
started <- proc.time()[["elapsed"]]
fit <- ch_train_forest(rows, panel, recipes, config, args[2L] == "QRF_OFFSET", progress, directory)
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(directory, "fit.rds"))
capture.output(sessionInfo(), file = file.path(directory, "session.txt"))
capture.output(warnings(), file = file.path(directory, "warnings.txt"))
write_json(list(status = "completed", stage = "training", runtime_sec = elapsed,
  formal_validation_used = FALSE, manuscript_adoption = FALSE), file.path(directory, "status.json"),
  auto_unbox = TRUE, pretty = TRUE)
progress("training_completed", list(runtime_sec = elapsed))

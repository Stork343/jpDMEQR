library(data.table)
library(jsonlite)
library(yaml)
config_path <- "config/application_charite_hb_20261105.yml"
config <- read_yaml(config_path)
check <- read_json(file.path(config$report_directory, "forest_checks.json"), simplifyVector = TRUE)
stopifnot(check$status == "passed", !check$formal_validation_used)
directory <- file.path(config$output_root, "forest_implementation")
if (dir.exists(directory)) stop("实现目录已经存在")
dir.create(directory)
paths <- c(config_path, "research/charite_application_20261105/features.R",
  "research/empirical_selection_20261101/prediction_core.R", "research/empirical_selection_20261101/ridge_direction.R",
  "research/ridge_direction_20261105/forest.R", "research/charite_application_20261105/forest_train.R",
  "research/charite_application_20261105/run_forest.R")
stopifnot(all(file.copy(paths, file.path(directory, c("config.yml", "features.R", "prediction_core.R",
  "ridge_direction.R", "forest_base.R", "forest_train.R", "run_discovery.R")))))
queue <- data.table(trait = "hemoglobin", method = c("QRF_RAW", "QRF_OFFSET"), tau = NA_real_,
  run_id = "all", job_id = 1:2, status = "pending")
queue[, directory := file.path(config$output_root, "models", paste0(method, "_all"))]
dir.create(file.path(config$output_root, "forest_queue"))
fwrite(queue, file.path(config$output_root, "forest_queue", "jobs.csv"))
write_json(list(code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), validation_used = FALSE),
  file.path(directory, "source.json"), auto_unbox = TRUE, pretty = TRUE)

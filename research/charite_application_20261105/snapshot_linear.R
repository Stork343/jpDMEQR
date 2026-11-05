library(jsonlite)
library(yaml)
config_path <- "config/application_charite_hb_20261105.yml"
config <- read_yaml(config_path)
check <- read_json(file.path(config$report_directory, "training_checks.json"), simplifyVector = TRUE)
stopifnot(check$status == "passed", !check$formal_validation_used)
directory <- file.path(config$output_root, "linear_implementation")
if (dir.exists(directory)) stop("实现目录已经存在")
dir.create(directory)
paths <- c(config_path, "research/charite_application_20261105/features.R",
  "research/charite_application_20261105/train_only.R", "research/empirical_selection_20261101/prediction_core.R",
  "research/empirical_selection_20261101/ridge_direction.R", "research/charite_application_20261105/run_linear.R")
stopifnot(all(file.copy(paths, file.path(directory, c("config.yml", "features.R", "train_only.R",
  "prediction_core.R", "ridge_direction.R", "run_discovery.R")))))
write_json(list(code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  validation_used = FALSE), file.path(directory, "source.json"), auto_unbox = TRUE, pretty = TRUE)

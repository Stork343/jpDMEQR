library(data.table)
library(jsonlite)
library(yaml)
config_path <- "config/development_ibdmdb_forest_20261105.yml"
config <- read_yaml(config_path)
comparison <- fread("results/preflight/empirical_selection_20261101/ibdmdb_median/comparisons.csv")[tau == config$decision_quantile]
check <- read_json("results/preflight/empirical_selection_20261101/ibdmdb_cohort/forest_checks.json", simplifyVector = TRUE)
stopifnot(nrow(comparison) == 4L, all(comparison$relative_gain > 0), check$status == "passed")
root <- config$output_root
if (dir.exists(root)) stop("运行目录已经存在")
dir.create(root, recursive = TRUE)
implementation <- file.path(root, "implementation")
dir.create(implementation)
paths <- c(config_path, "research/empirical_selection_20261101/prediction_core.R",
  "research/empirical_selection_20261101/ridge_direction.R", "research/ibdmdb_application_20261105/features.R",
  "research/ridge_direction_20261105/forest.R", "research/ibdmdb_application_20261105/run_forest.R")
stopifnot(all(file.copy(paths, file.path(implementation,
  c("config.yml", "prediction_core.R", "ridge_direction.R", "features.R", "forest.R", "run_discovery.R")))))
queue <- CJ(trait = sprintf("fold%02d", 1:5), method = unlist(config$methods), sorted = FALSE)
queue[, `:=`(tau = NA_real_, run_id = "native", job_id = seq_len(.N))]
queue[, directory := file.path(root, trait, paste0(method, "_", run_id))]
queue[, status := "pending"]
fwrite(queue, file.path(root, "jobs.csv"))

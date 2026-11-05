library(data.table)
library(jsonlite)
library(yaml)
config_path <- "config/development_ibdmdb_median_20261105.yml"
config <- read_yaml(config_path)
comparison <- fread("results/preflight/empirical_selection_20261101/ibdmdb_prediction/comparisons.csv")
criterion <- comparison[tau == config$decision_quantile & comparator %in% c("ENET", "PRESS_RIDGE")]
stopifnot(nrow(criterion) == 2L, all(criterion$relative_gain > 0))
root <- config$output_root
if (dir.exists(root)) stop("运行目录已经存在")
dir.create(root, recursive = TRUE)
jobs <- list()
for (fold in 1:5) {
  fold_name <- sprintf("fold%02d", fold)
  directory <- file.path(root, fold_name)
  dir.create(directory)
  paths <- file.path(config$source_root, fold_name, c("recipes.rds", "training_roles.csv"))
  stopifnot(all(file.copy(paths, directory)))
  jobs[[fold]] <- data.table(trait = fold_name, method = c("POOLED_OFFSET", "JOINT_DIRECT"), tau = 0.5, run_id = "tau050")
}
implementation <- file.path(root, "implementation")
dir.create(implementation)
paths <- file.path(config$source_root, "implementation", c("prediction_core.R", "ridge_direction.R", "features.R", "run_discovery.R"))
stopifnot(all(file.copy(paths, implementation)), file.copy(config_path, file.path(implementation, "config.yml")))
queue <- rbindlist(jobs)
queue[, job_id := seq_len(.N)]
queue[, directory := file.path(root, trait, paste0(method, "_", run_id))]
queue[, status := "pending"]
fwrite(queue, file.path(root, "jobs.csv"))
write_json(list(status = "prepared", selected_quantile = config$decision_quantile,
  selection_basis = config$selection_basis, formal_validation_used = FALSE, manuscript_adoption = FALSE),
  file.path(root, "status.json"), auto_unbox = TRUE, pretty = TRUE)

library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/clinical_search_20261105/panel.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
config <- read_yaml(args[1L])
root <- config$output_root
if (dir.exists(root)) stop("运行目录已经存在")
dir.create(root, recursive = TRUE)
rows <- fread(file.path(config$input_directory, "discovery.csv"))
stopifnot(all(rows$role == "discovery"))
panel <- cp_panel(config$input_directory)
jobs <- summaries <- list()
for (outer in seq_len(config$outer_folds)) {
  label <- paste0("outer_", outer)
  directory <- file.path(root, label)
  dir.create(directory)
  fitting <- copy(rows[outer_fold != outer])
  evaluation <- copy(rows[outer_fold == outer])
  ids <- sort(unique(fitting$strain_id))
  set.seed(config$role_seed + outer)
  shuffled <- sample(ids)
  assignment <- setNames(rep(seq_len(config$inner_folds), length.out = length(ids)), shuffled)
  fitting[, `:=`(inner_fold = assignment[strain_id], role = "training")]
  evaluation[, `:=`(inner_fold = 0L, role = "discovery_evaluation")]
  stopifnot(!length(intersect(fitting$strain_id, evaluation$strain_id)))
  fwrite(fitting, file.path(directory, "training.csv"))
  fwrite(evaluation, file.path(directory, "evaluation.csv"))
  recipes <- list()
  for (key in c("training", paste0("inner_", seq_len(config$inner_folds)))) {
    current <- if (key == "training") fitting else fitting[inner_fold != as.integer(sub("inner_", "", key))]
    started <- proc.time()[["elapsed"]]
    direction <- ep_direction_recipe(current, panel, config, cp_recipe, cp_features)
    recipe <- ep_offset_recipe(current, panel, direction, cp_features)
    recipes[[key]] <- recipe
    features <- cp_offset_features(current, panel, recipe)
    j <- recipe$penalized
    error <- max(abs(colSums(features$X[, j, drop = FALSE] * features$weight) / features$n),
      abs(colSums(features$X[, j, drop = FALSE]^2 * features$weight) / features$n - 1))
    stopifnot(error < 1e-10)
    summaries[[length(summaries) + 1L]] <- data.table(outer, panel = key, patients = features$n,
      records = nrow(features$X), variables = ncol(features$X), standardization_error = error,
      preparation_sec = proc.time()[["elapsed"]] - started)
  }
  saveRDS(recipes, file.path(directory, "recipes.rds"))
  queue <- rbindlist(list(CJ(method = c("JOINT_OFFSET", "JOINT_DIRECT", "POOLED", "POOLED_OFFSET"),
    tau = unlist(config$quantiles), sorted = FALSE),
    data.table(method = c("ENET", "PRESS_RIDGE", "CLINICAL_QR"), tau = NA_real_)))
  queue[, run_id := ifelse(is.na(tau), "all", sprintf("tau%03d", as.integer(100 * tau)))]
  queue[, `:=`(trait = label, status = "pending")]
  queue[, directory := file.path(directory, "models", paste0(method, "_", run_id))]
  jobs[[outer]] <- queue
  cat(label, "prepared\n")
  flush.console()
}
jobs <- rbindlist(jobs)
jobs[, job_id := seq_len(.N)]
fwrite(jobs, file.path(root, "jobs.csv"))
fwrite(rbindlist(summaries), file.path(root, "preparation_summary.csv"))
implementation <- file.path(root, "implementation")
dir.create(implementation)
paths <- c(args[1L], "research/clinical_search_20261105/panel.R",
  "research/clinical_search_20261105/run_discovery.R", "research/empirical_selection_20261101/prediction_core.R",
  "research/empirical_selection_20261101/ridge_direction.R")
stopifnot(all(file.copy(paths, file.path(implementation,
  c("config.yml", "panel.R", "run_discovery.R", "prediction_core.R", "ridge_direction.R")))))
write_json(list(status = "prepared", study_id = config$study_id,
  code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  formal_validation_used = FALSE, manuscript_adoption = FALSE), file.path(root, "preparation.json"),
  auto_unbox = TRUE, pretty = TRUE)

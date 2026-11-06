library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/clinical_search_20261105/panel.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1L, 2L))
stage <- if (length(args) == 2L) args[2L] else "complete"
stopifnot(stage %in% c("necessary", "complete", "offsets"))
config <- read_yaml(args[1L])
root <- config$output_root
if (stage == "offsets") {
  prepared <- read_json(file.path(root, "preparation.json"))
  stopifnot(identical(prepared$stage, "necessary"), identical(prepared$study_id, config$study_id))
} else {
  if (dir.exists(root)) stop("运行目录已经存在")
  dir.create(root, recursive = TRUE)
}
rows <- fread(file.path(config$input_directory, "discovery.csv"))
stopifnot(all(rows$role == "discovery"))
panel <- cp_panel(config$input_directory)
jobs <- summaries <- list()
for (outer in seq_len(config$outer_folds)) {
  label <- paste0("outer_", outer)
  directory <- file.path(root, label)
  if (stage != "offsets") dir.create(directory)
  fitting <- copy(rows[outer_fold != outer])
  evaluation <- copy(rows[outer_fold == outer])
  ids <- sort(unique(fitting$strain_id))
  set.seed(config$role_seed + outer)
  shuffled <- sample(ids)
  assignment <- setNames(rep(seq_len(config$inner_folds), length.out = length(ids)), shuffled)
  fitting[, `:=`(inner_fold = assignment[strain_id], role = "training")]
  evaluation[, `:=`(inner_fold = 0L, role = "discovery_evaluation")]
  stopifnot(!length(intersect(fitting$strain_id, evaluation$strain_id)))
  if (stage == "offsets") {
    for (pair in list(list("training.csv", fitting), list("evaluation.csv", evaluation))) {
      saved <- fread(file.path(directory, pair[[1L]]))
      stopifnot(isTRUE(all.equal(saved, pair[[2L]], tolerance = 1e-13, check.attributes = FALSE)))
    }
    direct <- readRDS(file.path(directory, "recipes.rds"))
  } else {
    fwrite(fitting, file.path(directory, "training.csv"))
    fwrite(evaluation, file.path(directory, "evaluation.csv"))
  }
  recipes <- list()
  for (key in c("training", paste0("inner_", seq_len(config$inner_folds)))) {
    current <- if (key == "training") fitting else fitting[inner_fold != as.integer(sub("inner_", "", key))]
    started <- proc.time()[["elapsed"]]
    if (stage == "necessary") {
      base <- cp_recipe(current, panel, config)
      recipe <- list(direction_recipe = list(base_recipe = base), training_strains = base$training_strains)
      features <- cp_features(current, panel, base)
      j <- which(features$penalty_factor == 1)
    } else {
      direction <- ep_direction_recipe(current, panel, config, cp_recipe, cp_features)
      recipe <- ep_offset_recipe(current, panel, direction, cp_features)
      if (stage == "offsets") {
        stopifnot(identical(recipe$direction_recipe$base_recipe, direct[[key]]$direction_recipe$base_recipe))
      }
      features <- cp_offset_features(current, panel, recipe)
      j <- recipe$penalized
    }
    recipes[[key]] <- recipe
    error <- max(abs(colSums(features$X[, j, drop = FALSE] * features$weight) / features$n),
      abs(colSums(features$X[, j, drop = FALSE]^2 * features$weight) / features$n - 1))
    stopifnot(error < 1e-10)
    summaries[[length(summaries) + 1L]] <- data.table(outer, panel = key, patients = features$n,
      records = nrow(features$X), variables = ncol(features$X), standardization_error = error,
      preparation_sec = proc.time()[["elapsed"]] - started)
  }
  if (stage == "offsets") {
    destination <- file.path(directory, "offset_recipes.rds")
    if (file.exists(destination)) {
      stopifnot(identical(readRDS(destination), recipes))
    } else {
      saveRDS(recipes, paste0(destination, ".partial"))
      stopifnot(file.rename(paste0(destination, ".partial"), destination))
    }
    next
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
if (stage == "offsets") {
  fwrite(rbindlist(summaries), file.path(root, "offset_preparation_summary.csv"))
  quit(save = "no")
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
write_json(list(status = "prepared", stage = stage, study_id = config$study_id,
  code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  formal_validation_used = FALSE, manuscript_adoption = FALSE), file.path(root, "preparation.json"),
  auto_unbox = TRUE, pretty = TRUE)

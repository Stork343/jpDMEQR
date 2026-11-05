library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/charite_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
config <- read_yaml("config/application_charite_hb_20261105.yml")
root <- config$output_root
if (dir.exists(root)) stop("运行目录已经存在")
dir.create(root, recursive = TRUE)
rows <- fread(file.path(config$input_directory, "training.csv"))
stopifnot(all(rows$role == "training"))
panel <- ch_panel(config$input_directory)
recipes <- list()
summary <- list()
for (key in c("training", "inner_1", "inner_2", "inner_3")) {
  current <- if (key == "training") rows else rows[inner_fold != as.integer(sub("inner_", "", key))]
  started <- proc.time()[["elapsed"]]
  direction <- ep_direction_recipe(current, panel, config, ch_recipe, ch_features)
  recipe <- ep_offset_recipe(current, panel, direction, ch_features)
  recipes[[key]] <- recipe
  features <- ch_offset_features(current, panel, recipe)
  j <- recipe$penalized
  error <- max(abs(colSums(features$X[, j, drop = FALSE] * features$weight) / features$n),
    abs(colSums(features$X[, j, drop = FALSE]^2 * features$weight) / features$n - 1))
  stopifnot(error < 1e-10)
  summary[[key]] <- data.table(panel = key, patients = features$n, records = nrow(features$X),
    variables = ncol(features$X), standardization_error = error,
    background_sec = proc.time()[["elapsed"]] - started)
  for (mode in c("raw", "offset")) {
    directory <- file.path(root, "exports", mode, key)
    dir.create(directory, recursive = TRUE)
    sets <- list(fitting = current)
    if (key != "training") sets$validation <- rows[inner_fold == as.integer(sub("inner_", "", key))]
    for (name in names(sets)) {
      selected_rows <- sets[[name]]
      value <- if (mode == "offset") ch_offset_features(selected_rows, panel, recipe) else
        ch_features(selected_rows, panel, direction$base_recipe)
      offset <- if (mode == "offset") value$offset else rep(0, nrow(selected_rows))
      fwrite(as.data.table(value$X[, -1L, drop = FALSE]), file.path(directory, paste0(name, "_X.csv")))
      fwrite(data.table(sample_id = selected_rows$sample_id, strain_id = selected_rows$strain_id,
        y = value$y, raw_y = value$raw_y, weight = value$weight, offset = offset),
        file.path(directory, paste0(name, "_meta.csv")))
      write_json(list(y_scale = value$recipe$y_scale, y_center = value$recipe$y_center,
        columns = colnames(value$X)[-1L]), file.path(directory, "scaling.json"), auto_unbox = TRUE, pretty = TRUE)
    }
  }
  cat(key, "prepared\n")
  flush.console()
}
saveRDS(recipes, file.path(root, "recipes.rds"))
fwrite(rbindlist(summary), file.path(root, "recipe_summary.csv"))
queue <- rbindlist(list(CJ(method = c("JOINT_OFFSET", "JOINT_DIRECT", "POOLED", "POOLED_OFFSET"),
  tau = unlist(config$quantiles), sorted = FALSE), data.table(method = c("ENET", "PRESS_RIDGE", "CLINICAL_QR"), tau = NA_real_)))
queue[, run_id := ifelse(is.na(tau), "all", sprintf("tau%03d", as.integer(100 * tau)))]
queue[, `:=`(trait = "hemoglobin", job_id = seq_len(.N), status = "pending")]
queue[, directory := file.path(root, "models", paste0(method, "_", run_id))]
fwrite(queue, file.path(root, "jobs.csv"))
write_json(list(status = "prepared", training_patients = uniqueN(rows$strain_id), formal_validation_used = FALSE),
  file.path(root, "preparation.json"), auto_unbox = TRUE, pretty = TRUE)

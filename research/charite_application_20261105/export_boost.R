library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/charite_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
config <- read_yaml("config/application_charite_hb_20261105.yml")
rows <- fread(file.path(config$input_directory, "training.csv"))
panel <- ch_panel(config$input_directory)
recipes <- readRDS(file.path(config$output_root, "recipes.rds"))
root <- file.path(config$output_root, "exports_full_precision")
if (dir.exists(root)) stop("导出目录已经存在")
checks <- list()
for (key in names(recipes)) {
  sets <- if (key == "training") list(fitting = rows) else {
    fold <- as.integer(sub("inner_", "", key))
    list(fitting = rows[inner_fold != fold], validation = rows[inner_fold == fold])
  }
  for (mode in c("raw", "offset")) {
    directory <- file.path(root, mode, key)
    dir.create(directory, recursive = TRUE)
    for (part in names(sets)) {
      current <- sets[[part]]
      value <- if (mode == "offset") ch_offset_features(current, panel, recipes[[key]]) else
        ch_features(current, panel, recipes[[key]]$direction_recipe$base_recipe)
      offset <- if (mode == "offset") value$offset else rep(0, nrow(current))
      reconstructed <- value$y * value$recipe$y_scale + value$recipe$y_center + offset
      error <- max(abs(reconstructed - current$y))
      stopifnot(error < 1e-10)
      fwrite(as.data.table(value$X[, -1L, drop = FALSE]), file.path(directory, paste0(part, "_X.csv")))
      fwrite(data.table(sample_id = current$sample_id, strain_id = current$strain_id, y = value$y,
        raw_y = value$raw_y, weight = value$weight, offset = offset), file.path(directory, paste0(part, "_meta.csv")))
      write_json(list(y_scale = value$recipe$y_scale, y_center = value$recipe$y_center, columns = colnames(value$X)[-1L]),
        file.path(directory, "scaling.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
      checks[[length(checks) + 1L]] <- data.table(panel = key, mode, part, maximum_error = error)
    }
  }
}
fwrite(rbindlist(checks), file.path(config$report_directory, "response_scale_checks.csv"))

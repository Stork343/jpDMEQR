ch_train_forest <- function(rows, panel, recipes, config, offset_model, progress, directory) {
  stopifnot(all(rows$role == "training"), config$forest$trees == 2000L)
  grid <- expand.grid(regression_splitting = unlist(config$forest$regression_splitting),
    mtry_multiplier = unlist(config$forest$mtry_multiplier), min_node_size = unlist(config$forest$min_node_size),
    honesty = unlist(config$forest$honesty))
  settings <- rbindlist(lapply(seq_len(nrow(grid)), function(g)
    data.table(forest_index = g, calibration = unlist(config$forest$calibration))))
  quantiles <- unlist(config$quantiles)
  stopifnot(identical(quantiles, c(0.25, 0.5, 0.75)))
  loss <- matrix(0, nrow(settings), length(quantiles))
  features <- function(current, recipe) {
    value <- if (offset_model) ch_offset_features(current, panel, recipe) else
      ch_features(current, panel, recipe$direction_recipe$base_recipe)
    value$sample_id <- current$sample_id
    value
  }
  for (fold in sort(unique(rows$inner_fold))) {
    fitting_rows <- rows[inner_fold != fold]
    held_rows <- rows[inner_fold == fold]
    recipe <- recipes[[paste0("inner_", fold)]]
    fitting <- features(fitting_rows, recipe)
    held <- features(held_rows, recipe)
    for (g in seq_len(nrow(grid))) {
      progress("forest_inner", list(fold = fold, grid = g))
      model <- rs_forest_fit(fitting, grid[g, ], 20261105L + fold, offset_model)
      prediction <- rs_forest_predict(model, held)
      calibration_prediction <- list(training = rs_forest_predict(model, fitting),
        oob = rs_forest_predict(model, fitting, oob = TRUE))
      calibrations <- list()
      for (s in which(settings$forest_index == g)) {
        mode <- settings$calibration[s]
        for (k in seq_along(quantiles)) {
          calibration <- if (mode == "none") c(hemoglobin = 0) else
            ep_calibrate(fitting_rows, calibration_prediction[[mode]][, k], quantiles[k], "strain_id", "condition")
          adjusted <- prediction[, k] + calibration[held_rows$condition]
          loss[s, k] <- loss[s, k] + sum(ep_loss(held$raw_y, adjusted, quantiles[k]) * held$weight)
          calibrations[[paste(s, k, sep = "_")]] <- calibration
        }
      }
      saveRDS(list(model = model, recipe = recipe, prediction = prediction, calibration = calibrations,
        sample_ids = held_rows$sample_id), file.path(directory, sprintf("inner_%d_forest_%d.rds", fold, g)))
      rm(model)
      gc()
    }
  }
  selected <- apply(loss, 2L, which.min)
  recipe <- recipes$training
  fitting <- features(rows, recipe)
  models <- list()
  calibration <- vector("list", length(quantiles))
  for (g in unique(settings$forest_index[selected])) {
    progress("forest_training", list(grid = g))
    model <- rs_forest_fit(fitting, grid[g, ], 20261115L, offset_model)
    prediction <- list(training = rs_forest_predict(model, fitting), oob = rs_forest_predict(model, fitting, oob = TRUE))
    for (k in which(settings$forest_index[selected] == g)) {
      mode <- settings$calibration[selected[k]]
      calibration[[k]] <- if (mode == "none") c(hemoglobin = 0) else
        ep_calibrate(rows, prediction[[mode]][, k], quantiles[k], "strain_id", "condition")
    }
    models[[as.character(g)]] <- model
  }
  list(models = models, recipe = recipe, quantiles = quantiles, calibration = calibration,
    selected = selected, settings = settings, grid = grid, loss = loss, offset_model = offset_model,
    model_type = "forest")
}

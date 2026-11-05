rs_forest_fit <- function(features, settings, seed, offset_model, compute_oob = TRUE) {
  X <- features$X[, -1L, drop = FALSE]
  stopifnot(length(features$sample_id) == nrow(X), !anyDuplicated(features$sample_id))
  model <- grf::quantile_forest(X, features$y, num.trees = 2000L,
    quantiles = c(0.25, 0.5, 0.75), regression.splitting = settings$regression_splitting,
    clusters = factor(features$cluster_id), equalize.cluster.weights = TRUE,
    mtry = min(ncol(X), as.integer(ceiling(settings$mtry_multiplier * (sqrt(ncol(X)) + 20)))),
    min.node.size = as.integer(settings$min_node_size), honesty = settings$honesty,
    sample.fraction = 0.5, compute.oob.predictions = compute_oob, num.threads = 1L, seed = seed)
  list(model = model, columns = colnames(X), settings = settings, offset_model = offset_model,
       sample_ids = features$sample_id, grf_version = as.character(utils::packageVersion("grf")))
}

rs_forest_predict <- function(fit, features, oob = FALSE) {
  X <- features$X[, -1L, drop = FALSE]
  stopifnot(identical(colnames(X), fit$columns))
  if (oob) stopifnot(identical(features$sample_id, fit$sample_ids))
  prediction <- predict(fit$model, if (oob) NULL else X, quantiles = c(0.25, 0.5, 0.75), num.threads = 1L)$predictions
  prediction <- prediction * features$recipe$y_scale + features$recipe$y_center
  if (fit$offset_model) {
    stopifnot(length(features$offset) == nrow(X))
    prediction <- prediction + features$offset
  }
  stopifnot(all(is.finite(prediction)))
  prediction
}

rs_forest_prediction <- function(rows, evaluation, panel, recipes, offset_model, progress, directory,
                                 base_feature_function = proteome_features) {
  grid <- expand.grid(regression_splitting = c(FALSE, TRUE), mtry_multiplier = c(1, 5),
                       min_node_size = c(5L, 10L), honesty = c(TRUE, FALSE))
  settings <- data.table::rbindlist(lapply(seq_len(nrow(grid)), function(j)
    data.table::data.table(forest_index = j, calibration = c("none", "training", "oob"))))
  quantiles <- c(0.25, 0.5, 0.75)
  loss <- matrix(0, nrow(settings), length(quantiles))
  feature_function <- function(current, recipe) {
    result <- if (offset_model) ep_offset_features(current, panel, recipe, base_feature_function) else
      base_feature_function(current, panel, recipe$direction_recipe$base_recipe)
    result$sample_id <- current$sample_id
    result
  }
  for (fold in sort(unique(rows$inner_fold))) {
    fitting_rows <- rows[inner_fold != fold]
    validation_rows <- rows[inner_fold == fold]
    recipe <- recipes[[paste0("inner_", fold)]]
    fitting <- feature_function(fitting_rows, recipe)
    validation <- feature_function(validation_rows, recipe)
    for (g in seq_len(nrow(grid))) {
      progress("forest_inner", list(fold = fold, grid = g, clusters = fitting$n))
      model <- rs_forest_fit(fitting, grid[g, ], 20261105L + fold, offset_model)
      predictions <- rs_forest_predict(model, validation)
      calibration_predictions <- list(training = rs_forest_predict(model, fitting),
                                      oob = rs_forest_predict(model, fitting, oob = TRUE))
      calibrations <- list()
      for (j in which(settings$forest_index == g)) {
        mode <- settings$calibration[j]
        for (k in seq_along(quantiles)) {
          calibration <- if (mode == "none") setNames(rep(0, length(unique(rows$condition))), sort(unique(rows$condition))) else
            ep_calibrate(fitting_rows, calibration_predictions[[mode]][, k], quantiles[k], "strain_id", "condition")
          adjusted <- predictions[, k] + calibration[validation_rows$condition]
          loss[j, k] <- loss[j, k] + sum(ep_loss(validation$raw_y, adjusted, quantiles[k]) * validation$weight)
          calibrations[[paste(j, k, sep = "_")]] <- calibration
        }
      }
      saveRDS(list(model = model, recipe = recipe, prediction = predictions, calibration = calibrations),
        file.path(directory, sprintf("inner_%d_forest_%d.rds", fold, g)))
      rm(model, calibration_predictions)
      gc()
    }
  }
  selected <- apply(loss, 2L, which.min)
  recipe <- recipes$training
  fitting <- feature_function(rows, recipe)
  test <- feature_function(evaluation, recipe)
  predictions <- matrix(NA_real_, nrow(evaluation), length(quantiles))
  calibrations <- vector("list", length(quantiles))
  models <- list()
  for (g in unique(settings$forest_index[selected])) {
    progress("forest_training", list(grid = g))
    model <- rs_forest_fit(fitting, grid[g, ], 20261115L, offset_model)
    predicted <- rs_forest_predict(model, test)
    calibration_predictions <- list(training = rs_forest_predict(model, fitting),
                                    oob = rs_forest_predict(model, fitting, oob = TRUE))
    for (k in which(settings$forest_index[selected] == g)) {
      mode <- settings$calibration[selected[k]]
      calibrations[[k]] <- if (mode == "none") setNames(rep(0, length(unique(rows$condition))), sort(unique(rows$condition))) else
        ep_calibrate(rows, calibration_predictions[[mode]][, k], quantiles[k], "strain_id", "condition")
      predictions[, k] <- predicted[, k] + calibrations[[k]][evaluation$condition]
    }
    models[[as.character(g)]] <- model
  }
  list(models = models, recipe = recipe, prediction = predictions, quantiles = quantiles,
       calibration = calibrations, selected = selected, settings = settings, grid = grid, loss = loss,
       offset_model = offset_model)
}

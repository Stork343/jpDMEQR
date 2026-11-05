ch_train_prediction <- function(rows, panel, method, tau, config, progress, directory,
                                 recipe_function, feature_function, prediction_function,
                                 calibration_mode) {
  stopifnot(all(rows$role == "training"), calibration_mode %in% c("oof", "training"))
  if (calibration_mode == "training") stopifnot(method != "ENET")
  folds <- sort(unique(rows$inner_fold))
  recipes <- lapply(folds, function(fold) recipe_function(rows[inner_fold != fold], panel, config))
  counts <- vapply(folds, function(fold) uniqueN(rows[inner_fold != fold]$strain_id), integer(1))
  maximum <- min(vapply(seq_along(folds), function(j)
    min(recipes[[j]]$search_dimension, floor(counts[j] / log(counts[j]))), numeric(1)))
  loss <- oof <- NULL
  for (f in seq_along(folds)) {
    fitting_rows <- rows[inner_fold != folds[f]]
    held_rows <- rows[inner_fold == folds[f]]
    training <- feature_function(fitting_rows, panel, recipes[[f]])
    validation <- feature_function(held_rows, panel, recipes[[f]])
    progress("inner_fit", list(fold = folds[f], patients = training$n, variables = ncol(training$X)))
    fit <- ep_path(training, method, tau, config, maximum)
    predicted <- as.matrix(prediction_function(validation, fit$beta))
    if (calibration_mode == "training") {
      fitted_prediction <- as.matrix(prediction_function(training, fit$beta))
      for (j in seq_len(ncol(predicted))) {
        calibration <- ep_calibrate(fitting_rows, fitted_prediction[, j], tau, "strain_id", "condition")
        predicted[, j] <- predicted[, j] + calibration[held_rows$condition]
      }
    }
    if (f == 1L) {
      grid <- fit$grid
      oof <- matrix(NA_real_, nrow(rows), ncol(predicted))
      loss <- numeric(ncol(predicted))
    }
    stopifnot(identical(as.list(grid), as.list(fit$grid)))
    values <- if (method == "ENET") (validation$raw_y - predicted)^2 else ep_loss(validation$raw_y, predicted, tau)
    loss <- loss + colSums(values * validation$weight)
    oof[which(rows$inner_fold == folds[f]), ] <- predicted
    saveRDS(list(path = fit, recipe = recipes[[f]], prediction = predicted), file.path(directory, paste0("inner_", folds[f], ".rds")))
  }
  stopifnot(all(is.finite(oof)))
  choice <- which.min(loss)
  recipe <- recipe_function(rows, panel, config)
  training <- feature_function(rows, panel, recipe)
  progress("training_fit", list(selected_index = choice))
  final <- ep_path(training, method, tau, config, maximum, as.list(grid[choice]))
  stopifnot(ncol(final$beta) == 1L)
  beta <- setNames(as.numeric(final$beta[, 1L]), colnames(training$X))
  quantiles <- if (method == "ENET") unlist(config$quantiles) else tau
  calibration_prediction <- if (calibration_mode == "oof") oof[, choice] else drop(prediction_function(training, beta))
  calibration <- lapply(quantiles, function(value) ep_calibrate(rows, calibration_prediction, value, "strain_id", "condition"))
  list(beta = beta, recipe = recipe, quantiles = quantiles, calibration = calibration, choice = choice,
    grid = grid, loss = loss, oof = oof[, choice], final_diagnostics = final$diagnostics,
    calibration_mode = calibration_mode)
}

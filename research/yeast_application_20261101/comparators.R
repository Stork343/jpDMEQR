yeast_pooled_path <- function(features, tau, lambda = NULL) {
  row_weight <- nrow(features$X) * features$weight / features$n
  fitted <- soy_pooled_lp_cpp(features$X * row_weight, features$y * row_weight,
    features$penalty_factor, tau, 40L, 0.01, if (is.null(lambda)) numeric() else lambda, TRUE)
  rownames(fitted$beta) <- colnames(features$X)
  fitted
}

yeast_pooled_prediction <- function(rows, evaluation, panel, config, tau, progress, directory) {
  loss <- numeric(40L)
  oof <- matrix(NA_real_, nrow(rows), 40L)
  diagnostics <- list()
  for (fold in sort(unique(rows$inner_fold))) {
    recipe <- yeast_recipe(rows[inner_fold != fold], panel)
    training <- yeast_features(rows[inner_fold != fold], panel, recipe)
    validation <- yeast_features(rows[inner_fold == fold], panel, recipe)
    progress("pooled_inner", list(fold = fold, strains = training$n, variables = ncol(training$X)))
    fit <- yeast_pooled_path(training, tau)
    prediction <- as.matrix(yeast_raw_prediction(validation, fit$beta))
    loss <- loss + colSums(yeast_loss(validation$raw_y, prediction, tau) * validation$weight)
    oof[which(rows$inner_fold == fold), ] <- prediction
    diagnostics[[fold]] <- fit[c("primal_dual_gap", "dual_violation", "iterations")]
    saveRDS(list(fit = fit, recipe = recipe, prediction = prediction),
            file.path(directory, paste0("inner_", fold, ".rds")))
  }
  choice <- which.min(loss)
  stopifnot(all(is.finite(oof)))
  recipe <- yeast_recipe(rows, panel)
  calibration <- yeast_calibrate(rows, oof[, choice], tau, recipe$replicate_levels)
  training <- yeast_features(rows, panel, recipe)
  progress("pooled_training", list(penalty_index = choice))
  fit <- yeast_pooled_path(training, tau)
  beta <- fit$beta[, choice]
  test <- yeast_features(evaluation, panel, recipe)
  prediction <- drop(yeast_raw_prediction(test, beta)) + calibration[as.character(evaluation$replicate)]
  list(beta = beta, recipe = recipe, prediction = prediction, calibration = calibration,
       oof = oof[, choice], losses = loss, selected_penalty_index = choice,
       diagnostics = diagnostics, final_diagnostics = fit[c("primal_dual_gap", "dual_violation", "iterations")])
}

yeast_forest_fit <- function(features, settings, seed) {
  X <- features$X[, -1L, drop = FALSE]
  model <- grf::quantile_forest(X, features$y, num.trees = 2000L,
    quantiles = c(0.25, 0.5, 0.75), regression.splitting = settings$regression_splitting,
    clusters = factor(features$cluster_id), equalize.cluster.weights = TRUE,
    mtry = min(ncol(X), as.integer(ceiling(settings$mtry_multiplier * (sqrt(ncol(X)) + 20)))),
    min.node.size = as.integer(settings$min_node_size), honesty = settings$honesty,
    sample.fraction = 0.5, compute.oob.predictions = FALSE, num.threads = 1L, seed = seed)
  list(model = model, columns = colnames(X), settings = settings)
}

yeast_forest_predict <- function(fit, features) {
  X <- features$X[, -1L, drop = FALSE]
  stopifnot(identical(colnames(X), fit$columns))
  predicted <- predict(fit$model, X, quantiles = c(0.25, 0.5, 0.75), num.threads = 1L)$predictions
  stopifnot(all(is.finite(predicted)))
  predicted * features$recipe$y_scale + features$recipe$y_center
}

yeast_forest_prediction <- function(rows, evaluation, panel, config, progress, directory) {
  grid <- expand.grid(regression_splitting = c(FALSE, TRUE),
                       mtry_multiplier = c(1, 5), min_node_size = c(5L, 10L), honesty = c(TRUE, FALSE))
  quantiles <- unlist(config$quantiles)
  loss <- matrix(0, nrow(grid), length(quantiles))
  oof <- array(NA_real_, c(nrow(rows), nrow(grid), length(quantiles)))
  for (fold in sort(unique(rows$inner_fold))) {
    recipe <- yeast_recipe(rows[inner_fold != fold], panel)
    training <- yeast_features(rows[inner_fold != fold], panel, recipe)
    validation <- yeast_features(rows[inner_fold == fold], panel, recipe)
    for (g in seq_len(nrow(grid))) {
      progress("forest_inner", list(fold = fold, grid = g, strains = training$n))
      fit <- yeast_forest_fit(training, grid[g, ], 20261201L + fold)
      predicted <- yeast_forest_predict(fit, validation)
      for (j in seq_along(quantiles)) {
        loss[g, j] <- loss[g, j] + sum(yeast_loss(validation$raw_y, predicted[, j], quantiles[j]) * validation$weight)
        oof[which(rows$inner_fold == fold), g, j] <- predicted[, j]
      }
      saveRDS(list(fit = fit, recipe = recipe, prediction = predicted),
              file.path(directory, sprintf("inner_%d_grid_%d.rds", fold, g)))
      rm(fit, predicted)
      gc()
    }
  }
  stopifnot(all(is.finite(oof)))
  selected <- apply(loss, 2L, which.min)
  recipe <- yeast_recipe(rows, panel)
  training <- yeast_features(rows, panel, recipe)
  test <- yeast_features(evaluation, panel, recipe)
  predictions <- matrix(NA_real_, nrow(evaluation), length(quantiles))
  calibrations <- vector("list", length(quantiles))
  models <- list()
  for (g in sort(unique(selected))) {
    progress("forest_training", list(grid = g))
    fit <- yeast_forest_fit(training, grid[g, ], 20261299L)
    predicted <- yeast_forest_predict(fit, test)
    for (j in which(selected == g)) {
      calibrations[[j]] <- yeast_calibrate(rows, oof[, g, j], quantiles[j], recipe$replicate_levels)
      predictions[, j] <- predicted[, j] + calibrations[[j]][as.character(evaluation$replicate)]
    }
    models[[as.character(g)]] <- fit
  }
  list(models = models, recipe = recipe, prediction = predictions, calibration = calibrations,
       losses = loss, selected = selected, grid = grid, oof = oof)
}

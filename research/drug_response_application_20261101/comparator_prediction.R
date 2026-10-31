drug_pooled_prediction <- function(rows, evaluation, expression, config, tau, progress, directory) {
  folds <- sort(unique(rows$inner_fold))
  nlambda <- 40L
  ratios <- exp(seq(0, log(0.01), length.out = nlambda))
  loss <- numeric(nlambda)
  oof <- matrix(NA_real_, nrow(rows), nlambda)
  diagnostics <- list()
  for (fold in folds) {
    training_rows <- rows[inner_fold != fold]
    validation_rows <- rows[inner_fold == fold]
    recipe <- drug_recipe(training_rows, expression, config)
    training <- drug_features(training_rows, expression, recipe)
    validation <- drug_features(validation_rows, expression, recipe)
    progress("pooled_inner_fit", list(fold = fold, patients = training$n, variables = ncol(training$X)))
    fit <- drug_pooled_path(training, tau, nlambda = nlambda, lambda_min = 0.01, progress = TRUE)
    stopifnot(max(abs(fit$lambda / fit$lambda[1L] - ratios)) < 1e-10)
    predicted <- validation$X %*% fit$beta
    loss <- loss + colSums(drug_loss(validation$y, predicted, tau) * validation$weight)
    oof[which(rows$inner_fold == fold), ] <- predicted
    diagnostics[[fold]] <- fit$diagnostics
    saveRDS(list(recipe = recipe, fit = fit, prediction = predicted),
            file.path(directory, paste0("inner_", fold, ".rds")))
    rm(training, validation, fit, predicted)
    gc()
  }
  selected <- which.min(loss)
  stopifnot(all(is.finite(oof)))
  calibration <- drug_calibration(rows, oof[, selected], tau)
  recipe <- drug_recipe(rows, expression, config)
  training <- drug_features(rows, expression, recipe)
  progress("pooled_training_fit", list(selected_ratio = ratios[selected], patients = training$n))
  fit <- drug_pooled_path(training, tau, nlambda = nlambda, lambda_min = 0.01, progress = TRUE)
  beta <- fit$beta[, selected]
  final_diagnostics <- fit$diagnostics
  rm(training, fit)
  gc()
  test <- drug_features(evaluation, expression, recipe)
  prediction <- drop(test$X %*% beta)
  list(beta = beta, recipe = recipe, prediction = prediction +
         drug_apply_calibration(calibration, evaluation$CONC),
       uncalibrated_prediction = prediction, calibration = calibration,
       oof = oof[, selected], loss = loss, ratios = ratios, selected = selected,
       diagnostics = diagnostics, final_diagnostics = final_diagnostics)
}

drug_qrf_design <- function(features) {
  cbind(drug_dose = features$dose, features$clinical[, -1L, drop = FALSE], features$genes)
}

drug_qrf_fit <- function(features, mtry_multiplier, min_node_size, seed, trees = 1000L) {
  X <- drug_qrf_design(features)
  patients <- sort(unique(features$cluster_id))
  group <- match(features$cluster_id, patients)
  set.seed(seed)
  inbag <- lapply(seq_len(trees), function(j)
    tabulate(sample.int(length(patients), length(patients), replace = TRUE), nbins = length(patients))[group])
  fit <- ranger::ranger(x = X, y = features$y,
    num.trees = trees, mtry = min(ncol(X), as.integer(ceiling(mtry_multiplier * sqrt(ncol(X))))),
    min.node.size = as.integer(min_node_size), quantreg = TRUE, inbag = inbag,
    always.split.variables = c("drug_dose", colnames(features$clinical)[-1L]),
    oob.error = FALSE, num.threads = 1L, seed = seed, na.action = "na.fail", verbose = FALSE)
  list(fit = fit, columns = colnames(X), mtry_multiplier = mtry_multiplier,
       min_node_size = min_node_size, trees = trees, bootstrap_unit = "patient_id")
}

drug_qrf_predict <- function(model, features, quantiles) {
  X <- drug_qrf_design(features)
  stopifnot(identical(colnames(X), model$columns))
  result <- predict(model$fit, data = X, type = "quantiles", quantiles = quantiles, num.threads = 1L)$predictions
  stopifnot(all(is.finite(result)))
  result
}

drug_qrf_prediction <- function(rows, evaluation, expression, config, progress, directory) {
  grid <- expand.grid(mtry_multiplier = c(1, 5), min_node_size = c(5L, 10L, 20L))
  quantiles <- unlist(config$quantiles)
  loss <- matrix(0, nrow(grid), length(quantiles))
  oof <- array(NA_real_, c(nrow(rows), nrow(grid), length(quantiles)))
  for (fold in sort(unique(rows$inner_fold))) {
    training_rows <- rows[inner_fold != fold]
    validation_rows <- rows[inner_fold == fold]
    recipe <- drug_recipe(training_rows, expression, config)
    training <- drug_features(training_rows, expression, recipe, expand = FALSE)
    validation <- drug_features(validation_rows, expression, recipe, expand = FALSE)
    weight <- 1 / as.numeric(table(validation$cluster_id)[validation$cluster_id])
    for (g in seq_len(nrow(grid))) {
      progress("qrf_inner_fit", list(fold = fold, grid = g,
               mtry_multiplier = grid$mtry_multiplier[g], min_node_size = grid$min_node_size[g]))
      model <- drug_qrf_fit(training, grid$mtry_multiplier[g], grid$min_node_size[g], 20261200L + fold)
      predicted <- drug_qrf_predict(model, validation, quantiles)
      for (j in seq_along(quantiles)) {
        loss[g, j] <- loss[g, j] + sum(weight * drug_loss(validation$y, predicted[, j], quantiles[j]))
        oof[which(rows$inner_fold == fold), g, j] <- predicted[, j]
      }
      saveRDS(list(recipe = recipe, model = model, prediction = predicted),
              file.path(directory, sprintf("inner_%d_grid_%d.rds", fold, g)))
      rm(model, predicted)
      gc()
    }
    rm(training, validation)
    gc()
  }
  stopifnot(all(is.finite(oof)))
  selected <- apply(loss, 2L, which.min)
  recipe <- drug_recipe(rows, expression, config)
  training <- drug_features(rows, expression, recipe, expand = FALSE)
  test <- drug_features(evaluation, expression, recipe, expand = FALSE)
  predictions <- matrix(NA_real_, nrow(evaluation), length(quantiles))
  calibrations <- vector("list", length(quantiles))
  models <- list()
  for (g in sort(unique(selected))) {
    progress("qrf_training_fit", list(grid = g))
    model <- drug_qrf_fit(training, grid$mtry_multiplier[g], grid$min_node_size[g], 20261299L)
    index <- which(selected == g)
    predicted <- drug_qrf_predict(model, test, quantiles)
    for (j in index) {
      calibrations[[j]] <- drug_calibration(rows, oof[, g, j], quantiles[j])
      predictions[, j] <- predicted[, j] + drug_apply_calibration(calibrations[[j]], evaluation$CONC)
    }
    models[[as.character(g)]] <- model
  }
  list(models = models, recipe = recipe, prediction = predictions, calibration = calibrations,
       loss = loss, selected = selected, grid = grid, oof = oof)
}

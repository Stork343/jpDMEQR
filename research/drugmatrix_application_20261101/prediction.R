dm_recipe <- function(rows, panel) {
  ids <- rows$sample_id
  weight <- 1 / as.numeric(table(rows$compound_id)[rows$compound_id])
  weight <- weight / sum(weight)
  raw <- panel$X[ids, panel$allowed_probes, drop = FALSE]
  availability <- colSums(is.finite(raw) * weight)
  variables <- which(availability >= 0.8)
  stopifnot(length(variables) > 0L)
  raw <- raw[, variables, drop = FALSE]
  medians <- matrixStats::colMedians(raw, na.rm = TRUE)
  for (j in seq_len(ncol(raw))) raw[!is.finite(raw[, j]), j] <- medians[j]
  center <- colSums(raw * weight)
  scale <- sqrt(colSums(sweep(raw, 2L, center, "-")^2 * weight))
  keep <- is.finite(scale) & scale > 1e-8
  low <- panel$C[ids, , drop = FALSE]
  decomposition <- qr(low)
  low_columns <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  y_center <- sum(rows$y * weight)
  y_scale <- sqrt(sum((rows$y - y_center)^2 * weight))
  stopifnot(y_scale > 0, is.finite(y_scale), low_columns[1L] == 1L)
  list(variables = colnames(raw)[keep], medians = medians[keep],
       center = center[keep], scale = scale[keep], low_columns = low_columns,
       y_center = y_center, y_scale = y_scale, training_compounds = sort(unique(rows$compound_id)))
}

dm_features <- function(rows, panel, recipe) {
  raw <- panel$X[rows$sample_id, recipe$variables, drop = FALSE]
  for (j in seq_len(ncol(raw))) raw[!is.finite(raw[, j]), j] <- recipe$medians[j]
  high <- sweep(sweep(raw, 2L, recipe$center, "-"), 2L, recipe$scale, "/")
  low <- panel$C[rows$sample_id, recipe$low_columns, drop = FALSE]
  X <- cbind(low, high)
  stopifnot(all(is.finite(X)), all(X[, 1L] == 1), !anyDuplicated(colnames(X)))
  list(X = X, Z = matrix(1, nrow(X), 1L), y = (rows$y - recipe$y_center) / recipe$y_scale,
       raw_y = rows$y, cluster_id = rows$compound_id, n = length(unique(rows$compound_id)),
       weight = 1 / as.numeric(table(rows$compound_id)[rows$compound_id]),
       penalty_factor = c(rep(0, ncol(low)), rep(1, ncol(high))), recipe = recipe)
}

dm_path <- function(features, method, tau, config, maximum) {
  if (method == "JOINT") {
    fits <- lapply(config$working_penalties, function(lambda) {
      prepared <- pc_prepare_initial(features$y, features$X, features$Z, features$cluster_id,
        tau, matrix(lambda, 1L, 1L), features$penalty_factor, fit_always_included = FALSE)
      pc_orthogonal_path(prepared, tau, 0:maximum, solver = jp_joint_qp, gap_tolerance = config$gap_tolerance)
    })
    grid <- rbindlist(lapply(config$working_penalties, function(lambda)
      data.table(Lambda = lambda, size = 0:maximum)))
    list(beta = do.call(cbind, lapply(fits, `[[`, "coefficients")), grid = grid,
         diagnostics = lapply(fits, `[[`, "diagnostics"))
  } else if (method == "ENET") {
    glmnet::glmnet.control(fdev = 0, devmax = 1)
    fits <- lapply(config$elastic_net_alpha, function(alpha) yeast_enet_path(features, alpha, config))
    grid <- rbindlist(lapply(config$elastic_net_alpha, function(alpha)
      data.table(alpha = alpha, penalty_index = seq_len(config$elastic_net_path_length))))
    list(beta = as.matrix(do.call(cbind, lapply(fits, function(fit) fit$beta))), grid = grid,
         diagnostics = lapply(fits, function(fit) list(jerr = fit$jerr, lambda = fit$lambda)))
  } else {
    stopifnot(method == "POOLED")
    fit <- yeast_pooled_path(features, tau)
    list(beta = fit$beta, grid = data.table(penalty_index = seq_len(ncol(fit$beta))),
         diagnostics = fit[c("primal_dual_gap", "dual_violation", "iterations")])
  }
}

dm_predict <- function(features, beta) {
  features$X %*% beta * features$recipe$y_scale + features$recipe$y_center
}

dm_fit_prediction <- function(rows, evaluation, panel, method, tau, config, progress, directory) {
  folds <- sort(unique(rows$inner_fold))
  inner <- lapply(folds, function(fold) {
    recipe <- dm_recipe(rows[inner_fold != fold], panel)
    list(train = dm_features(rows[inner_fold != fold], panel, recipe),
         validation = dm_features(rows[inner_fold == fold], panel, recipe))
  })
  maximum <- min(vapply(inner, function(fold)
    max(0, floor(fold$train$n / log(fold$train$n)) - sum(fold$train$penalty_factor == 0)), numeric(1)))
  loss <- NULL
  oof <- NULL
  for (f in seq_along(folds)) {
    progress("inner_fit", list(fold = folds[f], compounds = inner[[f]]$train$n,
             variables = ncol(inner[[f]]$train$X), additions = maximum))
    fit <- dm_path(inner[[f]]$train, method, tau, config, maximum)
    predicted <- as.matrix(dm_predict(inner[[f]]$validation, fit$beta))
    if (f == 1L) {
      grid <- fit$grid
      oof <- matrix(NA_real_, nrow(rows), ncol(predicted))
      loss <- numeric(ncol(predicted))
    }
    stopifnot(identical(as.list(grid), as.list(fit$grid)), ncol(predicted) == ncol(oof))
    validation <- inner[[f]]$validation
    values <- if (method == "ENET") (validation$raw_y - predicted)^2 else yeast_loss(validation$raw_y, predicted, tau)
    loss <- loss + colSums(values * validation$weight)
    oof[which(rows$inner_fold == folds[f]), ] <- predicted
    saveRDS(list(path = fit, recipe = inner[[f]]$train$recipe, prediction = predicted),
            file.path(directory, paste0("inner_", folds[f], ".rds")))
  }
  stopifnot(all(is.finite(oof)))
  choice <- which.min(loss)
  recipe <- dm_recipe(rows, panel)
  training <- dm_features(rows, panel, recipe)
  test <- dm_features(evaluation, panel, recipe)
  progress("training_fit", list(compounds = training$n, selected_index = choice))
  final <- dm_path(training, method, tau, config, maximum)
  stopifnot(identical(as.list(final$grid), as.list(grid)))
  beta <- final$beta[, choice]
  names(beta) <- colnames(training$X)
  quantiles <- if (method == "ENET") unlist(config$quantiles) else tau
  calibration <- vapply(quantiles, function(value) {
    fit <- quantreg::rq.fit(matrix(training$weight, ncol = 1L),
      (rows$y - oof[, choice]) * training$weight, tau = value, method = "br")
    as.numeric(fit$coefficients)
  }, numeric(1))
  predicted <- drop(dm_predict(test, beta))
  prediction <- vapply(calibration, function(value) predicted + value, numeric(nrow(evaluation)))
  stopifnot(all(is.finite(prediction)))
  list(beta = beta, recipe = recipe, prediction = prediction, quantiles = quantiles,
       calibration = calibration, choice = choice, grid = grid, loss = loss,
       oof = oof[, choice], final_diagnostics = final$diagnostics)
}

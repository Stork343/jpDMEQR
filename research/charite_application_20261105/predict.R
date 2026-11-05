ch_design <- function(rows, panel, recipe) {
  values <- ch_values(rows, panel, recipe)[, recipe$columns, drop = FALSE]
  values <- sweep(sweep(values, 2L, recipe$center, "-"), 2L, recipe$scale, "/")
  X <- cbind("(Intercept)" = 1, values)
  stopifnot(all(is.finite(X)), !anyDuplicated(colnames(X)))
  X
}

ch_prediction_features <- function(rows, panel, recipe, offset_model) {
  if (!offset_model) return(list(X = ch_design(rows, panel, recipe), recipe = recipe, sample_id = rows$sample_id))
  direction <- recipe$direction_recipe
  base <- direction$base_recipe
  X <- ch_design(rows, panel, base)
  score <- drop(X %*% direction$background$beta)
  D <- (matrix(score, ncol = 1L) - X[, direction$low_columns, drop = FALSE] %*% direction$projection)
  D <- sweep(D, 2L, direction$direction_scale, "/")
  colnames(D) <- "RIDGE_mean"
  X <- cbind(X, D)
  X[, recipe$penalized] <- sweep(sweep(X[, recipe$penalized, drop = FALSE], 2L,
    recipe$center, "-"), 2L, recipe$scale, "/")
  list(X = X, recipe = recipe, offset = score * base$y_scale + base$y_center, sample_id = rows$sample_id)
}

ch_predict <- function(fit, rows, panel) {
  if (fit$model_type == "forest") {
    recipe <- if (fit$offset_model) fit$recipe else fit$recipe$direction_recipe$base_recipe
    features <- ch_prediction_features(rows, panel, recipe, fit$offset_model)
    prediction <- matrix(NA_real_, nrow(rows), length(fit$quantiles))
    for (g in unique(fit$settings$forest_index[fit$selected])) {
      native <- rs_forest_predict(fit$models[[as.character(g)]], features)
      indices <- which(fit$settings$forest_index[fit$selected] == g)
      prediction[, indices] <- native[, indices, drop = FALSE]
    }
  } else if (fit$model_type == "clinical") {
    X <- ch_design(rows, panel, fit$recipe)[, fit$clinical_columns, drop = FALSE]
    prediction <- X %*% fit$beta * fit$recipe$y_scale + fit$recipe$y_center
  } else {
    features <- ch_prediction_features(rows, panel, fit$recipe, fit$offset_model)
    base <- drop(if (fit$offset_model) ep_predict_offset(features, fit$beta) else ep_predict(features, fit$beta))
    prediction <- matrix(base, nrow(rows), length(fit$quantiles))
  }
  for (k in seq_along(fit$quantiles)) prediction[, k] <- prediction[, k] + fit$calibration[[k]][rows$condition]
  stopifnot(all(is.finite(prediction)))
  prediction
}

ch_expand_offset <- function(beta, recipe, calibration) {
  direction <- recipe$direction_recipe
  base <- direction$base_recipe
  background <- direction$background$beta
  dimension <- length(background)
  stopifnot(length(beta) == dimension + 1L, length(direction$direction_scale) == 1L)
  scale <- rep(1, length(beta))
  scale[recipe$penalized] <- recipe$scale
  correction <- beta * recipe$y_scale / scale
  projection <- numeric(dimension)
  projection[direction$low_columns] <- drop(direction$projection)
  expanded <- background * base$y_scale + correction[seq_len(dimension)] +
    (background - projection) * correction[dimension + 1L] / direction$direction_scale[1L]
  names(expanded) <- names(background)
  expanded[1L] <- expanded[1L] + base$y_center - sum(correction[recipe$penalized] * recipe$center) + calibration
  expanded
}

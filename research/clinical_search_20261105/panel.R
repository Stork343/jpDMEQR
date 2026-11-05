cp_panel <- function(directory) {
  frame <- data.table::fread(file.path(directory, "features.csv.gz"))
  X <- as.matrix(frame[, -1L, with = FALSE])
  rownames(X) <- frame$sample_id
  stopifnot(!anyDuplicated(rownames(X)), !anyDuplicated(colnames(X)),
    all(is.finite(X[!is.na(X)])))
  list(X = X)
}

cp_values <- function(rows, panel, recipe) {
  X <- panel$X[rows$sample_id, recipe$features, drop = FALSE]
  missing <- is.na(X)
  for (j in seq_along(recipe$features)) X[missing[, j], j] <- recipe$impute[j]
  M <- missing[, recipe$missing_features, drop = FALSE] * 1
  colnames(M) <- paste0("M_", recipe$missing_features)
  cbind(X, M)
}

cp_recipe <- function(rows, panel, config) {
  ids <- sort(unique(rows$strain_id))
  weight <- 1 / as.numeric(table(rows$strain_id)[rows$strain_id])
  X <- panel$X[rows$sample_id, , drop = FALSE]
  observed <- colSums(!is.na(X) * weight) / length(ids)
  keep <- observed >= config$minimum_observed_fraction |
    (startsWith(colnames(X), "C_") & observed > 0)
  features <- colnames(X)[keep]
  impute <- vapply(features, function(column) matrixStats::weightedMedian(X[, column], weight, na.rm = TRUE), numeric(1))
  missing_features <- features[colSums(is.na(X[, features, drop = FALSE])) > 0L]
  recipe <- list(features = features, impute = impute, missing_features = missing_features, observed = observed[features])
  values <- cp_values(rows, panel, recipe)
  center <- colSums(values * weight) / length(ids)
  scale <- sqrt(colSums(sweep(values, 2L, center, "-")^2 * weight) / length(ids))
  columns <- names(scale)[is.finite(scale) & scale > 1e-8]
  y_center <- sum(rows$y * weight) / length(ids)
  y_scale <- sqrt(sum((rows$y - y_center)^2 * weight) / length(ids))
  stopifnot(y_scale > 0, length(columns) > length(ids), length(unique(rows$condition)) == 1L)
  c(recipe, list(columns = columns, center = center[columns], scale = scale[columns],
    y_center = y_center, y_scale = y_scale, training_strains = ids, search_dimension = length(columns),
    contrasts = matrix(numeric(), 1L, 0L, dimnames = list(unique(rows$condition), character()))))
}

cp_features <- function(rows, panel, recipe) {
  values <- cp_values(rows, panel, recipe)[, recipe$columns, drop = FALSE]
  values <- sweep(sweep(values, 2L, recipe$center, "-"), 2L, recipe$scale, "/")
  X <- cbind("(Intercept)" = 1, values)
  stopifnot(all(is.finite(X)), !anyDuplicated(colnames(X)))
  list(X = X, Z = matrix(1, nrow(rows), 1L), y = (rows$y - recipe$y_center) / recipe$y_scale,
    raw_y = rows$y, cluster_id = rows$strain_id, n = data.table::uniqueN(rows$strain_id),
    weight = 1 / as.numeric(table(rows$strain_id)[rows$strain_id]),
    penalty_factor = c(0, rep(1, ncol(values))), recipe = recipe, sample_id = rows$sample_id)
}

cp_cached_recipe <- function(rows, panel, config, offset) {
  ids <- sort(unique(rows$strain_id))
  keys <- names(panel$recipes)[vapply(panel$recipes, function(recipe) identical(recipe$training_strains, ids), logical(1))]
  stopifnot(length(keys) == 1L)
  recipe <- panel$recipes[[keys]]
  if (offset) recipe else recipe$direction_recipe$base_recipe
}

cp_offset_features <- function(rows, panel, recipe) ep_offset_features(rows, panel, recipe, cp_features)

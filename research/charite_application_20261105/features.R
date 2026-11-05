ch_panel <- function(root) {
  protein <- data.table::fread(file.path(root, "protein_abundance.csv.gz"))
  values <- as.matrix(protein[, -1L, with = FALSE])
  rownames(values) <- protein$sample_id
  stopifnot(all(values[!is.na(values)] > 0))
  list(protein = log2(values), annotation = data.table::fread(file.path(root, "protein_catalog.csv")))
}

ch_clinical <- function(rows, bmi_median) {
  bmi <- rows$bmi
  missing <- is.na(bmi)
  bmi[missing] <- bmi_median
  cbind(C_age = rows$age, C_bmi = bmi, C_charlson = rows$charlson,
    C_log_day = log1p(rows$admission_day), C_who4 = as.numeric(rows$who == 4L),
    C_who5 = as.numeric(rows$who == 5L), C_who6 = as.numeric(rows$who == 6L),
    C_who7 = as.numeric(rows$who == 7L), C_bmi_missing = as.numeric(missing))
}

ch_values <- function(rows, panel, recipe) {
  P <- panel$protein[rows$sample_id, recipe$genes, drop = FALSE]
  missing <- is.na(P)
  for (j in seq_along(recipe$genes)) P[missing[, j], j] <- recipe$impute[j]
  colnames(P) <- paste0("P_", recipe$genes)
  M <- missing[, recipe$missing_genes, drop = FALSE] * 1
  colnames(M) <- paste0("M_", recipe$missing_genes)
  cbind(ch_clinical(rows, recipe$bmi_median), P, M)
}

ch_recipe <- function(rows, panel, config) {
  ids <- sort(unique(rows$strain_id))
  weight <- 1 / as.numeric(table(rows$strain_id)[rows$strain_id])
  P <- panel$protein[rows$sample_id, , drop = FALSE]
  observed <- colSums(!is.na(P) * weight) / length(ids)
  genes <- colnames(P)[observed >= config$protein_observed_training_fraction]
  impute <- vapply(genes, function(gene) matrixStats::weightedMedian(P[, gene], weight, na.rm = TRUE), numeric(1))
  missing_genes <- genes[colSums(is.na(P[, genes, drop = FALSE])) > 0L]
  recipe <- list(genes = genes, impute = impute, missing_genes = missing_genes,
    bmi_median = matrixStats::weightedMedian(rows$bmi, weight, na.rm = TRUE), observed = observed[genes])
  values <- ch_values(rows, panel, recipe)
  center <- colSums(values * weight) / length(ids)
  scale <- sqrt(colSums(sweep(values, 2L, center, "-")^2 * weight) / length(ids))
  columns <- names(scale)[is.finite(scale) & scale > 1e-8]
  y_center <- sum(rows$y * weight) / length(ids)
  y_scale <- sqrt(sum((rows$y - y_center)^2 * weight) / length(ids))
  stopifnot(y_scale > 0, length(columns) > length(ids))
  c(recipe, list(columns = columns, center = center[columns], scale = scale[columns],
    y_center = y_center, y_scale = y_scale, training_strains = ids, search_dimension = length(columns),
    contrasts = matrix(numeric(), 1L, 0L, dimnames = list("hemoglobin", character()))))
}

ch_features <- function(rows, panel, recipe) {
  values <- ch_values(rows, panel, recipe)[, recipe$columns, drop = FALSE]
  values <- sweep(sweep(values, 2L, recipe$center, "-"), 2L, recipe$scale, "/")
  X <- cbind("(Intercept)" = 1, values)
  stopifnot(all(is.finite(X)), !anyDuplicated(colnames(X)))
  list(X = X, Z = matrix(1, nrow(rows), 1L), y = (rows$y - recipe$y_center) / recipe$y_scale,
    raw_y = rows$y, cluster_id = rows$strain_id, n = data.table::uniqueN(rows$strain_id),
    weight = 1 / as.numeric(table(rows$strain_id)[rows$strain_id]),
    penalty_factor = c(0, rep(1, ncol(values))), recipe = recipe, sample_id = rows$sample_id)
}

ch_cached_recipe <- function(rows, panel, config, offset) {
  ids <- sort(unique(rows$strain_id))
  keys <- names(panel$recipes)[vapply(panel$recipes, function(recipe) identical(recipe$training_strains, ids), logical(1))]
  stopifnot(length(keys) == 1L)
  recipe <- panel$recipes[[keys]]
  if (offset) recipe else recipe$direction_recipe$base_recipe
}

ch_offset_features <- function(rows, panel, recipe) ep_offset_features(rows, panel, recipe, ch_features)

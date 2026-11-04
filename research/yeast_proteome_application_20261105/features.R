proteome_panel <- function(root) {
  protein <- data.table::fread(file.path(root, "protein_abundance.csv.gz"))
  genes <- protein[[1L]]
  values <- t(as.matrix(protein[, -1L, with = FALSE]))
  colnames(values) <- genes
  stopifnot(all(values[!is.na(values)] > 0))
  values <- log2(values)
  signal <- apply(values, 1L, median, na.rm = TRUE)
  list(protein = values - signal, signal = signal,
       population = readRDS(file.path(root, "population_controls.rds")))
}

proteome_low <- function(rows, panel, key, controls, contrasts) {
  values <- cbind(log_harvest_od = log2(rows$harvest_od), median_signal = panel$signal[rows$strain_id])
  values <- sweep(sweep(values, 2L, controls$center, "-"), 2L, controls$scale, "/")
  base <- cbind(ploidy = rows$ploidy - 2, ploidy_squared = (rows$ploidy - 2)^2,
    aneuploid = rows$aneuploid, panel$population[[key]]$scores[rows$strain_id, , drop = FALSE], values)
  condition <- contrasts[rows$condition, , drop = FALSE]
  interaction <- lapply(seq_len(ncol(condition)), function(j) {
    result <- base * condition[, j]
    colnames(result) <- paste0(colnames(base), "_", colnames(condition)[j])
    result
  })
  do.call(cbind, c(list("(Intercept)" = rep(1, nrow(rows)), condition, base), interaction))
}

proteome_recipe <- function(rows, panel, config) {
  ids <- sort(unique(rows$strain_id))
  keys <- names(panel$population)[vapply(panel$population,
    function(value) identical(value$eligible_training_strains, ids), logical(1))]
  stopifnot(length(keys) == 1L)
  conditions <- unlist(config$conditions)
  stopifnot(all(table(rows$strain_id, rows$condition)[, conditions, drop = FALSE] == 1L))
  P <- panel$protein[ids, , drop = FALSE]
  keep <- colMeans(!is.na(P)) >= config$protein_observed_fraction_minimum
  P <- P[, keep, drop = FALSE]
  impute <- apply(P, 2L, median, na.rm = TRUE)
  missing <- which(is.na(P), arr.ind = TRUE)
  if (nrow(missing)) P[missing] <- impute[missing[, 2L]]
  center <- colMeans(P)
  scale <- sqrt(colMeans(sweep(P, 2L, center, "-")^2))
  keep <- is.finite(scale) & scale > 1e-8
  controls <- cbind(log_harvest_od = log2(rows$harvest_od), median_signal = panel$signal[rows$strain_id])
  control_center <- colMeans(controls)
  control_scale <- sqrt(colMeans(sweep(controls, 2L, control_center, "-")^2))
  stopifnot(all(control_scale > 1e-8))
  controls <- list(center = control_center, scale = control_scale)
  contrasts <- stress_contrasts(conditions)
  low <- proteome_low(rows, panel, keys, controls, contrasts)
  decomposition <- qr(low)
  low_columns <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  y_center <- mean(rows$y)
  y_scale <- sqrt(mean((rows$y - y_center)^2))
  stopifnot(y_scale > 0, low_columns[1L] == 1L)
  list(panel_key = keys, genes = colnames(P)[keep], impute = impute[keep], center = center[keep],
       scale = scale[keep], controls = controls, contrasts = contrasts, low_columns = low_columns,
       search_dimension = sum(keep) * length(conditions), y_center = y_center, y_scale = y_scale,
       training_strains = ids)
}

proteome_features <- function(rows, panel, recipe) {
  P <- panel$protein[rows$strain_id, recipe$genes, drop = FALSE]
  missing <- which(is.na(P), arr.ind = TRUE)
  if (nrow(missing)) P[missing] <- recipe$impute[missing[, 2L]]
  P <- sweep(sweep(P, 2L, recipe$center, "-"), 2L, recipe$scale, "/")
  colnames(P) <- paste0("P_", recipe$genes)
  low <- proteome_low(rows, panel, recipe$panel_key, recipe$controls, recipe$contrasts)[, recipe$low_columns, drop = FALSE]
  condition <- recipe$contrasts[rows$condition, , drop = FALSE]
  interaction <- lapply(seq_len(ncol(condition)), function(j) {
    result <- P * condition[, j]
    colnames(result) <- paste0(colnames(P), "_", colnames(condition)[j])
    result
  })
  X <- do.call(cbind, c(list(low, P), interaction))
  stopifnot(all(is.finite(X)), all(X[, 1L] == 1), !anyDuplicated(colnames(X)))
  list(X = X, Z = matrix(1, nrow(rows), 1L), y = (rows$y - recipe$y_center) / recipe$y_scale,
       raw_y = rows$y, cluster_id = rows$strain_id, n = data.table::uniqueN(rows$strain_id),
       weight = 1 / as.numeric(table(rows$strain_id)[rows$strain_id]),
       penalty_factor = c(rep(0, ncol(low)), rep(1, ncol(X) - ncol(low))), recipe = recipe)
}

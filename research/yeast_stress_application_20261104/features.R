stress_contrasts <- function(conditions) {
  result <- stats::contr.helmert(length(conditions))
  result <- sweep(result, 2L, sqrt(colMeans(result^2)), "/")
  dimnames(result) <- list(conditions, paste0("condition", seq_len(ncol(result))))
  result
}

stress_low <- function(rows, scores, contrasts) {
  condition <- contrasts[rows$condition, , drop = FALSE]
  base <- cbind(ploidy = rows$ploidy - 2, ploidy_squared = (rows$ploidy - 2)^2,
                aneuploid = rows$aneuploid, scores[rows$strain_id, , drop = FALSE])
  interactions <- lapply(seq_len(ncol(condition)), function(j) {
    result <- base * condition[, j]
    colnames(result) <- paste0(colnames(base), "_", colnames(condition)[j])
    result
  })
  do.call(cbind, c(list("(Intercept)" = rep(1, nrow(rows)), condition, base), interactions))
}

stress_recipe <- function(rows, panel, config) {
  ids <- sort(unique(rows$strain_id))
  conditions <- unlist(config$conditions)
  stopifnot(all(table(rows$strain_id, rows$condition)[, conditions, drop = FALSE] == 1L))
  recipe <- y1011_genomic_recipe(ids, panel, config)
  scores <- y1011_scores(ids, panel, recipe)
  contrasts <- stress_contrasts(conditions)
  low <- stress_low(rows, scores, contrasts)
  decomposition <- qr(low)
  low_columns <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  y_center <- mean(rows$y)
  y_scale <- sqrt(mean((rows$y - y_center)^2))
  stopifnot(is.finite(y_scale), y_scale > 0, low_columns[1L] == 1L)
  c(recipe, list(contrasts = contrasts, search_dimension = length(conditions) * recipe$base_dimension,
                 low_columns = low_columns, y_center = y_center, y_scale = y_scale))
}

stress_features <- function(rows, panel, recipe) {
  ids <- sort(unique(rows$strain_id))
  scores <- y1011_scores(ids, panel, recipe)
  low <- stress_low(rows, scores, recipe$contrasts)[, recipe$low_columns, drop = FALSE]
  condition <- recipe$contrasts[rows$condition, , drop = FALSE]
  pg <- length(recipe$markers)
  pc <- length(recipe$copy_names)
  p <- pg + pc
  d <- ncol(low)
  X <- matrix(NA_real_, nrow(rows), d + p * (1L + ncol(condition)))
  X[, seq_len(d)] <- low
  for (start in seq.int(1L, pg, by = 4096L)) {
    j <- start:min(start + 4095L, pg)
    G <- y1011_genotypes(rows$strain_id, panel, recipe$panel_key, recipe$markers[j], recipe$g_center[j], recipe$g_scale[j])
    X[, d + j] <- G
    for (k in seq_len(ncol(condition))) X[, d + k * p + j] <- G * condition[, k]
  }
  C <- panel$copy[rows$strain_id, recipe$copy_names, drop = FALSE]
  missing <- which(is.na(C), arr.ind = TRUE)
  if (nrow(missing)) C[missing] <- recipe$c_median[missing[, 2L]]
  C <- sweep(sweep(C, 2L, recipe$c_center, "-"), 2L, recipe$c_scale, "/")
  X[, d + pg + seq_len(pc)] <- C
  for (k in seq_len(ncol(condition))) X[, d + k * p + pg + seq_len(pc)] <- C * condition[, k]
  names <- c(paste0("G_", recipe$markers), paste0("C_", recipe$copy_names))
  colnames(X) <- c(colnames(low), names, unlist(lapply(colnames(condition), function(value) paste0(names, "_", value))))
  stopifnot(all(is.finite(X)), all(X[, 1L] == 1), !anyDuplicated(colnames(X)))
  list(X = X, Z = matrix(1, nrow(rows), 1L), y = (rows$y - recipe$y_center) / recipe$y_scale,
       raw_y = rows$y, cluster_id = rows$strain_id, n = length(ids),
       weight = 1 / as.numeric(table(rows$strain_id)[rows$strain_id]),
       penalty_factor = c(rep(0, d), rep(1, ncol(X) - d)), recipe = recipe)
}

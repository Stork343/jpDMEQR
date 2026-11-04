arab_panel <- function(root, load_recipes = TRUE, panel_subdirectory = "ld20") {
  groups <- jsonlite::read_json(file.path(root, "training_sets.json"), simplifyVector = TRUE)
  beds <- lapply(names(groups), function(key)
    BEDMatrix::BEDMatrix(file.path(root, panel_subdirectory, paste0(key, "_panel")), simple_names = TRUE))
  names(beds) <- names(groups)
  recipes <- if (load_recipes) readRDS(file.path(root, "training_recipes.rds")) else NULL
  list(groups = groups, beds = beds, recipes = recipes)
}

arab_key <- function(rows, panel) {
  ids <- sort(unique(rows$accession_id))
  stopifnot(is.character(ids))
  keys <- names(panel$groups)[vapply(panel$groups, function(group) identical(sort(group), ids), logical(1))]
  stopifnot(length(keys) == 1L)
  keys
}

arab_low <- function(rows, scores) {
  temperature <- 2 * rows$high_temperature - 1
  pcs <- scores[rows$accession_id, , drop = FALSE]
  interaction <- pcs * temperature
  colnames(interaction) <- paste0(colnames(pcs), "_temperature")
  cbind("(Intercept)" = 1, temperature = temperature, pcs, interaction)
}

arab_scores <- function(ids, panel, recipe) {
  bed <- panel$beds[[recipe$panel_key]]
  result <- matrix(0, length(ids), ncol(recipe$loadings), dimnames = list(ids, colnames(recipe$loadings)))
  for (start in seq.int(1L, length(recipe$markers), by = 4096L)) {
    j <- start:min(start + 4095L, length(recipe$markers))
    G <- bed[ids, recipe$markers[j], drop = FALSE]
    stopifnot(!anyNA(G))
    G <- sweep(sweep(G, 2L, recipe$center[j], "-"), 2L, recipe$scale[j], "/")
    result <- result + G %*% recipe$loadings[j, , drop = FALSE]
  }
  result
}

arab_compute_recipe <- function(rows, panel, config) {
  ids <- sort(unique(rows$accession_id))
  key <- arab_key(rows, panel)
  bed <- panel$beds[[key]]
  markers <- colnames(bed)
  n <- length(ids)
  p <- length(markers)
  center <- scale <- numeric(p)
  K <- matrix(0, n, n)
  stopifnot(all(rows[, .(fraction = mean(high_temperature)), by = accession_id]$fraction == 0.5))
  for (start in seq.int(1L, p, by = 4096L)) {
    j <- start:min(start + 4095L, p)
    G <- bed[ids, markers[j], drop = FALSE]
    stopifnot(!anyNA(G))
    center[j] <- colMeans(G)
    G <- sweep(G, 2L, center[j], "-")
    scale[j] <- sqrt(colMeans(G^2))
    stopifnot(all(scale[j] > 1e-8))
    K <- K + tcrossprod(sweep(G, 2L, scale[j], "/")) / p
  }
  eig <- eigen(K, symmetric = TRUE)
  count <- as.integer(config$population_components)
  stopifnot(all(eig$values[seq_len(count)] > 1e-8))
  loadings <- matrix(NA_real_, p, count, dimnames = list(markers, paste0("PC", seq_len(count))))
  for (start in seq.int(1L, p, by = 4096L)) {
    j <- start:min(start + 4095L, p)
    G <- bed[ids, markers[j], drop = FALSE]
    G <- sweep(sweep(G, 2L, center[j], "-"), 2L, scale[j], "/")
    loadings[j, ] <- sweep(crossprod(G, eig$vectors[, seq_len(count), drop = FALSE]),
                           2L, eig$values[seq_len(count)] * p, "/") * sqrt(n - 1L)
  }
  scores <- eig$vectors[, seq_len(count), drop = FALSE] * sqrt(n - 1L)
  dimnames(scores) <- list(ids, colnames(loadings))
  low <- arab_low(rows, scores)
  decomposition <- qr(low)
  low_columns <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  weight <- 1 / as.numeric(table(rows$accession_id)[rows$accession_id])
  y_center <- sum(weight * rows$y) / n
  y_scale <- sqrt(sum(weight * (rows$y - y_center)^2) / n)
  stopifnot(y_scale > 0, low_columns[1L] == 1L)
  list(panel_key = key, markers = markers, center = center, scale = scale, loadings = loadings,
       low_columns = low_columns, y_center = y_center, y_scale = y_scale,
       search_dimension = 2L * p, training_accessions = ids)
}

arab_recipe <- function(rows, panel, config) {
  key <- arab_key(rows, panel)
  recipe <- panel$recipes[[key]]
  stopifnot(identical(recipe$training_accessions, sort(unique(rows$accession_id))))
  recipe
}

arab_features <- function(rows, panel, recipe) {
  stopifnot(is.character(rows$accession_id))
  ids <- sort(unique(rows$accession_id))
  scores <- arab_scores(ids, panel, recipe)
  low <- arab_low(rows, scores)[, recipe$low_columns, drop = FALSE]
  p <- length(recipe$markers)
  d <- ncol(low)
  X <- matrix(NA_real_, nrow(rows), d + 2L * p)
  X[, seq_len(d)] <- low
  temperature <- 2 * rows$high_temperature - 1
  bed <- panel$beds[[recipe$panel_key]]
  for (start in seq.int(1L, p, by = 4096L)) {
    j <- start:min(start + 4095L, p)
    G <- bed[rows$accession_id, recipe$markers[j], drop = FALSE]
    stopifnot(!anyNA(G))
    G <- sweep(sweep(G, 2L, recipe$center[j], "-"), 2L, recipe$scale[j], "/")
    X[, d + j] <- G
    X[, d + p + j] <- G * temperature
  }
  colnames(X) <- c(colnames(low), paste0("G_", recipe$markers), paste0("GT_", recipe$markers))
  stopifnot(all(is.finite(X)), all(X[, 1L] == 1), !anyDuplicated(colnames(X)))
  list(X = X, Z = matrix(1, nrow(rows), 1L), y = (rows$y - recipe$y_center) / recipe$y_scale,
       raw_y = rows$y, cluster_id = rows$accession_id, n = length(ids),
       weight = 1 / as.numeric(table(rows$accession_id)[rows$accession_id]),
       penalty_factor = c(rep(0, d), rep(1, 2L * p)), recipe = recipe)
}

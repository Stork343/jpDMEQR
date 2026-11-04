rice_panel <- function(root) {
  groups <- jsonlite::read_json(file.path(root, "training_sets.json"), simplifyVector = TRUE)
  beds <- lapply(names(groups), function(key) BEDMatrix::BEDMatrix(file.path(root, paste0(key, "_panel")), simple_names = TRUE))
  names(beds) <- names(groups)
  list(groups = groups, beds = beds, root = root)
}

rice_low <- function(rows, scores) {
  awd <- as.numeric(rows$irrigation == "AWD") - 0.5
  year <- as.numeric(rows$year == 2014) - 0.5
  pcs <- scores[rows$accession_id, , drop = FALSE]
  interaction <- pcs * awd
  colnames(interaction) <- paste0(colnames(pcs), "_AWD")
  cbind("(Intercept)" = 1, AWD = awd, year2014 = year, AWD_year = awd * year, pcs, interaction)
}

rice_scores <- function(ids, panel, recipe) {
  bed <- panel$beds[[recipe$panel_key]]
  result <- matrix(0, length(ids), ncol(recipe$loadings), dimnames = list(ids, colnames(recipe$loadings)))
  for (start in seq.int(1L, length(recipe$markers), by = 4096L)) {
    j <- start:min(start + 4095L, length(recipe$markers))
    G <- bed[ids, recipe$markers[j], drop = FALSE]
    stopifnot(!anyNA(G))
    G <- sweep(sweep(G, 2L, recipe$center[j], "-"), 2L, recipe$scale[j], "/")
    stopifnot(all(is.finite(G)))
    result <- result + G %*% recipe$loadings[j, , drop = FALSE]
  }
  result
}

rice_recipe <- function(rows, panel, config) {
  ids <- sort(unique(rows$accession_id))
  keys <- names(panel$groups)[vapply(panel$groups, function(group) identical(sort(group), ids), logical(1))]
  stopifnot(length(keys) == 1L)
  bed <- panel$beds[[keys]]
  markers <- colnames(bed)
  n <- length(ids)
  p <- length(markers)
  center <- scale <- interaction_center <- numeric(p)
  K <- matrix(0, n, n)
  awd_average <- rows[, .(a = mean(as.numeric(irrigation == "AWD") - 0.5)), by = accession_id]
  a <- awd_average$a[match(ids, awd_average$accession_id)]
  for (start in seq.int(1L, p, by = 4096L)) {
    j <- start:min(start + 4095L, p)
    G <- bed[ids, markers[j], drop = FALSE]
    stopifnot(!anyNA(G))
    center[j] <- colMeans(G)
    G <- sweep(G, 2L, center[j], "-")
    scale[j] <- sqrt(colMeans(G^2))
    stopifnot(all(scale[j] > 1e-8))
    G <- sweep(G, 2L, scale[j], "/")
    K <- K + tcrossprod(G) / p
    interaction_center[j] <- colMeans(G * a)
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
  low <- rice_low(rows, scores)
  decomposition <- qr(low)
  low_columns <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  weight <- 1 / as.numeric(table(rows$accession_id)[rows$accession_id])
  y_center <- sum(rows$y * weight) / n
  y_scale <- sqrt(sum((rows$y - y_center)^2 * weight) / n)
  interaction_scale <- sqrt(0.25 - interaction_center^2)
  stopifnot(is.finite(y_scale), y_scale > 0, all(interaction_scale > 1e-8), low_columns[1L] == 1L)
  list(panel_key = keys, markers = markers, search_dimension = 2L * p, center = center, scale = scale,
       interaction_center = interaction_center, interaction_scale = interaction_scale,
       loadings = loadings, low_columns = low_columns, y_center = y_center,
       y_scale = y_scale, training_accessions = ids)
}

rice_features <- function(rows, panel, recipe) {
  ids <- sort(unique(rows$accession_id))
  scores <- rice_scores(ids, panel, recipe)
  low <- rice_low(rows, scores)[, recipe$low_columns, drop = FALSE]
  p <- length(recipe$markers)
  d <- ncol(low)
  X <- matrix(NA_real_, nrow(rows), d + 2L * p)
  X[, seq_len(d)] <- low
  awd <- as.numeric(rows$irrigation == "AWD") - 0.5
  bed <- panel$beds[[recipe$panel_key]]
  for (start in seq.int(1L, p, by = 4096L)) {
    j <- start:min(start + 4095L, p)
    G <- bed[rows$accession_id, recipe$markers[j], drop = FALSE]
    stopifnot(!anyNA(G))
    G <- sweep(sweep(G, 2L, recipe$center[j], "-"), 2L, recipe$scale[j], "/")
    X[, d + j] <- G
    X[, d + p + j] <- sweep(sweep(G * awd, 2L, recipe$interaction_center[j], "-"),
                             2L, recipe$interaction_scale[j], "/")
  }
  colnames(X) <- c(colnames(low), paste0("G_", recipe$markers), paste0("GW_", recipe$markers))
  stopifnot(all(is.finite(low)), all(X[, 1L] == 1), !anyDuplicated(colnames(X)))
  list(X = X, Z = cbind(1, as.numeric(rows$year == 2014) - 0.5),
       y = (rows$y - recipe$y_center) / recipe$y_scale, raw_y = rows$y,
       cluster_id = rows$accession_id, n = length(ids),
       weight = 1 / as.numeric(table(rows$accession_id)[rows$accession_id]),
       penalty_factor = c(rep(0, d), rep(1, 2L * p)), recipe = recipe)
}

y1011_panel <- function(root) {
  groups <- jsonlite::read_json(file.path(root, "training_sets.json"), simplifyVector = TRUE)
  beds <- lapply(names(groups), function(key)
    BEDMatrix::BEDMatrix(file.path(root, paste0(key, "_panel")), simple_names = TRUE))
  names(beds) <- names(groups)
  copy <- data.table::fread(file.path(root, "copy_number.csv.gz"))
  ids <- copy$strain_id
  copy <- as.matrix(copy[, -1L, with = FALSE])
  rownames(copy) <- ids
  stopifnot(all(copy[!is.na(copy)] >= 0))
  list(groups = groups, beds = beds, copy = log2(1 + copy),
       map = data.table::fread(file.path(root, "genotype_identifiers.csv")))
}

y1011_genotypes <- function(ids, panel, key, markers, center = NULL, scale = NULL) {
  iid <- panel$map$IID[match(ids, panel$map$strain_id)]
  stopifnot(!anyNA(iid))
  G <- panel$beds[[key]][iid, markers, drop = FALSE]
  if (is.null(center)) return(G)
  missing <- which(is.na(G), arr.ind = TRUE)
  if (nrow(missing)) G[missing] <- center[missing[, 2L]]
  G <- sweep(sweep(G, 2L, center, "-"), 2L, scale, "/")
  stopifnot(all(is.finite(G)))
  G
}

y1011_low <- function(rows, scores) {
  concentration <- rows$high_salt - 0.5
  base <- cbind(ploidy = rows$ploidy - 2, ploidy_squared = (rows$ploidy - 2)^2,
                aneuploid = rows$aneuploid, scores[rows$strain_id, , drop = FALSE])
  interaction <- base * concentration
  colnames(interaction) <- paste0(colnames(base), "_salt")
  cbind("(Intercept)" = 1, high_salt = concentration, base, interaction)
}

y1011_scores <- function(ids, panel, recipe) {
  result <- matrix(0, length(ids), ncol(recipe$loadings), dimnames = list(ids, colnames(recipe$loadings)))
  for (start in seq.int(1L, length(recipe$markers), by = 4096L)) {
    j <- start:min(start + 4095L, length(recipe$markers))
    G <- y1011_genotypes(ids, panel, recipe$panel_key, recipe$markers[j], recipe$g_center[j], recipe$g_scale[j])
    result <- result + G %*% recipe$loadings[j, , drop = FALSE]
  }
  result
}

y1011_recipe <- function(rows, panel, config) {
  ids <- sort(unique(rows$strain_id))
  keys <- names(panel$groups)[vapply(panel$groups, function(group) identical(sort(group), ids), logical(1))]
  stopifnot(length(keys) == 1L, all(table(rows$strain_id) == 2L),
            all(rows[, .(value = sum(high_salt)), by = strain_id]$value == 1L))
  markers <- colnames(panel$beds[[keys]])
  p <- length(markers)
  n <- length(ids)
  g_center <- g_scale <- numeric(p)
  K <- matrix(0, n, n)
  for (start in seq.int(1L, p, by = 4096L)) {
    j <- start:min(start + 4095L, p)
    G <- y1011_genotypes(ids, panel, keys, markers[j])
    g_center[j] <- colMeans(G, na.rm = TRUE)
    missing <- which(is.na(G), arr.ind = TRUE)
    if (nrow(missing)) G[missing] <- g_center[j][missing[, 2L]]
    G <- sweep(G, 2L, g_center[j], "-")
    g_scale[j] <- sqrt(colMeans(G^2))
    stopifnot(all(is.finite(g_scale[j])), all(g_scale[j] > 1e-8))
    K <- K + tcrossprod(sweep(G, 2L, g_scale[j], "/")) / p
  }
  eig <- eigen(K, symmetric = TRUE)
  count <- as.integer(config$population_components)
  stopifnot(all(eig$values[seq_len(count)] > 1e-8))
  loadings <- matrix(NA_real_, p, count, dimnames = list(markers, paste0("PC", seq_len(count))))
  for (start in seq.int(1L, p, by = 4096L)) {
    j <- start:min(start + 4095L, p)
    G <- y1011_genotypes(ids, panel, keys, markers[j], g_center[j], g_scale[j])
    loadings[j, ] <- sweep(crossprod(G, eig$vectors[, seq_len(count), drop = FALSE]),
                           2L, eig$values[seq_len(count)] * p, "/") * sqrt(n - 1L)
  }
  scores <- eig$vectors[, seq_len(count), drop = FALSE] * sqrt(n - 1L)
  dimnames(scores) <- list(ids, colnames(loadings))
  C <- panel$copy[ids, , drop = FALSE]
  observed <- colMeans(!is.na(C)) >= config$copy_observed_fraction_minimum
  C <- C[, observed, drop = FALSE]
  c_median <- apply(C, 2L, median, na.rm = TRUE)
  missing <- which(is.na(C), arr.ind = TRUE)
  if (nrow(missing)) C[missing] <- c_median[missing[, 2L]]
  c_center <- colMeans(C)
  c_scale <- sqrt(colMeans(sweep(C, 2L, c_center, "-")^2))
  variable <- is.finite(c_scale) & c_scale > 1e-8
  copy_names <- colnames(C)[variable]
  low <- y1011_low(rows, scores)
  decomposition <- qr(low)
  low_columns <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  y_center <- mean(rows$y)
  y_scale <- sqrt(mean((rows$y - y_center)^2))
  stopifnot(length(copy_names) > 0L, is.finite(y_scale), y_scale > 0, low_columns[1L] == 1L)
  list(panel_key = keys, markers = markers, copy_names = copy_names,
       search_dimension = 2L * (p + length(copy_names)), g_center = g_center, g_scale = g_scale,
       c_median = c_median[variable], c_center = c_center[variable], c_scale = c_scale[variable],
       loadings = loadings, low_columns = low_columns, y_center = y_center, y_scale = y_scale,
       training_strains = ids)
}

y1011_features <- function(rows, panel, recipe) {
  ids <- sort(unique(rows$strain_id))
  scores <- y1011_scores(ids, panel, recipe)
  low <- y1011_low(rows, scores)[, recipe$low_columns, drop = FALSE]
  pg <- length(recipe$markers)
  pc <- length(recipe$copy_names)
  p <- pg + pc
  d <- ncol(low)
  X <- matrix(NA_real_, nrow(rows), d + 2L * p)
  X[, seq_len(d)] <- low
  concentration <- 2 * rows$high_salt - 1
  for (start in seq.int(1L, pg, by = 4096L)) {
    j <- start:min(start + 4095L, pg)
    G <- y1011_genotypes(rows$strain_id, panel, recipe$panel_key, recipe$markers[j], recipe$g_center[j], recipe$g_scale[j])
    X[, d + j] <- G
    X[, d + p + j] <- G * concentration
  }
  C <- panel$copy[rows$strain_id, recipe$copy_names, drop = FALSE]
  missing <- which(is.na(C), arr.ind = TRUE)
  if (nrow(missing)) C[missing] <- recipe$c_median[missing[, 2L]]
  C <- sweep(sweep(C, 2L, recipe$c_center, "-"), 2L, recipe$c_scale, "/")
  X[, d + pg + seq_len(pc)] <- C
  X[, d + p + pg + seq_len(pc)] <- C * concentration
  names <- c(paste0("G_", recipe$markers), paste0("C_", recipe$copy_names))
  colnames(X) <- c(colnames(low), names, paste0(names, "_salt"))
  stopifnot(all(is.finite(X)), all(X[, 1L] == 1), !anyDuplicated(colnames(X)))
  list(X = X, Z = matrix(1, nrow(rows), 1L), y = (rows$y - recipe$y_center) / recipe$y_scale,
       raw_y = rows$y, cluster_id = rows$strain_id, n = length(ids),
       weight = 1 / as.numeric(table(rows$strain_id)[rows$strain_id]),
       penalty_factor = c(rep(0, d), rep(1, 2L * p)), recipe = recipe)
}

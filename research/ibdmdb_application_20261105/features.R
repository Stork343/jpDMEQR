ibd_panel <- function(root) {
  abundance <- data.table::fread(file.path(root, "metabolite_abundance.csv.gz"))
  P <- as.matrix(abundance[, -1L, with = FALSE])
  rownames(P) <- abundance$sample_id
  signal <- data.table::fread(file.path(root, "method_signal.csv"))
  S <- as.matrix(signal[, -1L, with = FALSE])
  rownames(S) <- signal$sample_id
  list(abundance = P, signal = S, annotation = data.table::fread(file.path(root, "feature_annotation.csv")))
}

ibd_covariates <- function(rows, panel) {
  site <- vapply(c("Cincinnati", "Emory", "MGH", "MGH Pediatrics"), function(value) as.numeric(rows$site == value), numeric(nrow(rows)))
  colnames(site) <- paste0("C_site_", make.names(colnames(site)))
  cbind(C_age = rows$age, C_week = rows$week_num, C_male = as.numeric(rows$sex == "Male"),
    C_CD = as.numeric(rows$diagnosis == "CD"), C_UC = as.numeric(rows$diagnosis == "UC"),
    C_antibiotics = as.numeric(rows$antibiotics == "Yes"),
    C_immunosuppressants = as.numeric(rows$immunosuppressants == "Yes"), site,
    panel$signal[rows$sample_id, , drop = FALSE])
}

ibd_transform <- function(rows, panel, genes, replacement) {
  P <- panel$abundance[rows$sample_id, genes, drop = FALSE]
  for (j in seq_along(genes)) P[P[, j] == 0, j] <- replacement[j]
  methods <- panel$annotation$Method[match(genes, panel$annotation$feature_id)]
  stopifnot(!anyNA(methods), all(P > 0))
  log2(P) - panel$signal[rows$sample_id, methods, drop = FALSE]
}

ibd_recipe <- function(rows, panel, config) {
  ids <- sort(unique(rows$strain_id))
  weight <- 1 / as.numeric(table(rows$strain_id)[rows$strain_id])
  P <- panel$abundance[rows$sample_id, , drop = FALSE]
  genes <- colnames(P)[colSums((P > 0) * weight) / length(ids) >= config$positive_fraction_minimum]
  replacement <- apply(P[, genes, drop = FALSE], 2L, function(x) min(x[x > 0]) / 2)
  values <- cbind(ibd_covariates(rows, panel), ibd_transform(rows, panel, genes, replacement))
  center <- colSums(values * weight) / length(ids)
  scale <- sqrt(colSums(sweep(values, 2L, center, "-")^2 * weight) / length(ids))
  keep <- names(scale)[is.finite(scale) & scale > 1e-8]
  y_center <- sum(rows$y * weight) / length(ids)
  y_scale <- sqrt(sum((rows$y - y_center)^2 * weight) / length(ids))
  stopifnot(length(keep) > length(ids), y_scale > 0)
  list(genes = genes, replacement = replacement, columns = keep, center = center[keep], scale = scale[keep],
    y_center = y_center, y_scale = y_scale, training_strains = ids, search_dimension = length(keep),
    contrasts = matrix(numeric(), 1L, 0L, dimnames = list("fecalcal", character())))
}

ibd_features <- function(rows, panel, recipe) {
  values <- cbind(ibd_covariates(rows, panel), ibd_transform(rows, panel, recipe$genes, recipe$replacement))
  values <- sweep(sweep(values[, recipe$columns, drop = FALSE], 2L, recipe$center, "-"), 2L, recipe$scale, "/")
  X <- cbind("(Intercept)" = 1, values)
  weight <- 1 / as.numeric(table(rows$strain_id)[rows$strain_id])
  stopifnot(all(is.finite(X)), !anyDuplicated(colnames(X)), all(X[, 1L] == 1))
  list(X = X, Z = matrix(1, nrow(rows), 1L), y = (rows$y - recipe$y_center) / recipe$y_scale,
    raw_y = rows$y, cluster_id = rows$strain_id, n = data.table::uniqueN(rows$strain_id), weight = weight,
    penalty_factor = c(0, rep(1, ncol(values))), recipe = recipe)
}

drug_clinical <- function(rows) {
  stopifnot(all(rows$expression_source %in% c("Sanger RNASeq", "Broad (CCLE) RNASeq")),
            all(rows$ASSAY %in% c("a", "s")))
  Z <- cbind("(Intercept)" = 1,
             expression_sanger = as.numeric(rows$expression_source == "Sanger RNASeq"),
             assay_s = as.numeric(rows$ASSAY == "s"))
  attr(Z, "assign") <- 0:2
  Z
}

drug_recipe <- function(rows, expression, config) {
  models <- sort(unique(rows$SANGER_MODEL_ID))
  stopifnot(all(models %in% rownames(expression$fpkm)))
  raw <- expression$fpkm[models, , drop = FALSE]
  observed <- colMeans(is.finite(raw))
  expressed <- colMeans(raw > 1, na.rm = TRUE)
  selected <- which(observed >= config$gene_filter$observed_training_fraction &
                    expressed >= config$gene_filter$training_fpkm_above_one_fraction)
  transformed <- log2(raw[, selected, drop = FALSE] + 1)
  medians <- matrixStats::colMedians(transformed, na.rm = TRUE)
  missing <- which(!is.finite(transformed), arr.ind = TRUE)
  if (nrow(missing)) transformed[missing] <- medians[missing[, 2L]]
  centers <- colMeans(transformed)
  scales <- sqrt(colSums(sweep(transformed, 2L, centers, "-")^2) / (nrow(transformed) - 1L))
  keep <- which(is.finite(scales) & scales > 1e-8)
  stopifnot(length(keep) > length(models))
  Z <- drug_clinical(rows)
  decomposition <- qr(Z)
  clinical_columns <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  stopifnot(clinical_columns[1L] == 1L)
  maximum_dose <- max(rows$CONC)
  minimum_dose <- min(rows$CONC)
  dose <- (sqrt(rows$CONC) - sqrt(minimum_dose)) / (sqrt(maximum_dose) - sqrt(minimum_dose))
  basis <- splines::bs(dose, df = 5L, degree = 3L, intercept = TRUE)
  list(genes = colnames(raw)[selected[keep]], centers = centers[keep],
       scales = scales[keep], medians = medians[keep],
       tissues = sort(unique(rows$tissue)), clinical_columns = clinical_columns,
       maximum_dose = maximum_dose, minimum_dose = minimum_dose, basis = basis,
       training_models = models, training_patients = sort(unique(rows$patient_id)))
}

drug_features <- function(rows, expression, recipe, expand = TRUE) {
  raw <- expression$fpkm[match(rows$SANGER_MODEL_ID, rownames(expression$fpkm)), recipe$genes, drop = FALSE]
  genes <- log2(raw + 1)
  missing <- which(!is.finite(genes), arr.ind = TRUE)
  if (nrow(missing)) genes[missing] <- recipe$medians[missing[, 2L]]
  genes <- sweep(sweep(genes, 2L, recipe$centers, "-"), 2L, recipe$scales, "/")
  tissue <- vapply(recipe$tissues[-1L], function(value) as.numeric(rows$tissue == value), numeric(nrow(rows)))
  if (length(recipe$tissues) > 1L) {
    colnames(tissue) <- paste0("tissue_", recipe$tissues[-1L])
    genes <- cbind(genes, tissue)
  }
  Z <- drug_clinical(rows)[, recipe$clinical_columns, drop = FALSE]
  attr(Z, "assign") <- seq_len(ncol(Z)) - 1L
  dose <- (sqrt(rows$CONC) - sqrt(recipe$minimum_dose)) /
    (sqrt(recipe$maximum_dose) - sqrt(recipe$minimum_dose))
  basis <- predict(recipe$basis, newx = dose)
  result <- list(y = rows$normalized_intensity, cluster_id = rows$patient_id,
                 dose = dose, concentration = rows$CONC, clinical = Z,
                 n = length(unique(rows$patient_id)), imputed_gene_values = nrow(missing),
                 unseen_tissue_rows = sum(!rows$tissue %in% recipe$tissues))
  if (!expand) {
    result$genes <- genes
    return(result)
  }
  low <- do.call(cbind, lapply(seq_len(ncol(Z)), function(j) basis * Z[, j]))
  colnames(low) <- as.vector(t(outer(colnames(Z), seq_len(ncol(basis)), paste, sep = "_dose")))
  X <- matrix(0, nrow(rows), ncol(low) + ncol(genes) * ncol(basis))
  X[, seq_len(ncol(low))] <- low
  names <- colnames(low)
  for (j in seq_len(ncol(basis))) {
    positions <- ncol(low) + (j - 1L) * ncol(genes) + seq_len(ncol(genes))
    X[, positions] <- genes * basis[, j]
    names <- c(names, paste0(colnames(genes), "_dose", j))
  }
  colnames(X) <- names
  stopifnot(all(is.finite(X)))
  result$X <- X
  result$Z <- matrix(1, nrow(rows), 1L)
  result$penalty_factor <- c(rep(0, ncol(low)), rep(1, ncol(X) - ncol(low)))
  result$weight <- 1 / as.numeric(table(result$cluster_id)[result$cluster_id])
  result
}

drug_loss <- function(y, prediction, tau) {
  residual <- y - prediction
  residual * (tau - (residual < 0))
}

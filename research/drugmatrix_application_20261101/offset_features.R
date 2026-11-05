dmo_panel <- function(config) {
  panel <- readRDS(config$input_panel)
  annotation <- data.table::fread(config$input_annotation)
  panel$allowed_probes <- intersect(annotation[TYPE == config$probe_type, ID], colnames(panel$X))
  stopifnot(length(panel$allowed_probes) >= 1000L)
  panel
}

dmo_recipe <- function(rows, panel, config) {
  recipe <- dm_recipe(rows, panel)
  recipe$training_strains <- recipe$training_compounds
  recipe$contrasts <- matrix(numeric(), 1L, 0L, dimnames = list("ALT", character()))
  weight <- 1 / as.numeric(table(rows$compound_id)[rows$compound_id])
  weight <- weight / sum(weight)
  clinical <- panel$C[rows$sample_id, recipe$low_columns[-1L], drop = FALSE]
  recipe$clinical_center <- colSums(clinical * weight)
  recipe$clinical_scale <- sqrt(colSums(sweep(clinical, 2L, recipe$clinical_center, "-")^2 * weight))
  stopifnot(all(is.finite(recipe$clinical_scale)), all(recipe$clinical_scale > 1e-8))
  recipe$search_dimension <- recipe$search_dimension + ncol(clinical)
  recipe
}

dmo_features <- function(rows, panel, recipe) {
  features <- dm_features(rows, panel, recipe)
  indices <- seq.int(2L, length(recipe$low_columns))
  features$X[, indices] <- sweep(sweep(features$X[, indices, drop = FALSE], 2L,
    recipe$clinical_center, "-"), 2L, recipe$clinical_scale, "/")
  features$penalty_factor <- c(0, rep(1, ncol(features$X) - 1L))
  features
}

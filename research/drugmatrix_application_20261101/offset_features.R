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
  recipe
}

dmo_features <- function(rows, panel, recipe) dm_features(rows, panel, recipe)

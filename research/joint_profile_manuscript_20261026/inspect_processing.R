source("research/score_validity_20261026/load_joint.R")
source_v2_module(".", "application_v2", .GlobalEnv)
paths <- c(
  "data/derived/GSE65391/processing_20261023_v1/official_R1_R2_normexp_training_quantile/preprocessing_cache/prediction_r01_f01_inner01/prepared_base.rds",
  "data/derived/GSE65391/processing_20261023_v1/processed_median_across_mapped_probes/preprocessing_cache/module_full_sample/prepared_base.rds",
  "results/raw/application/sle_complete_20261022_v1_exploratory/preprocessing/split01_swap0.rds")
for (path in paths) {
  cat("Source", path, "\n")
  if (!file.exists(path)) {
    cat("Absent\n")
    next
  }
  object <- readRDS(path)
  cat("Fields", names(object), "\n")
  for (name in names(object)) {
    value <- object[[name]]
    cat(name, "class", class(value), "dimensions", dim(value), "length", length(value), "\n")
    if (is.list(value) && !is.data.frame(value)) cat("Nested fields", names(value), "\n")
  }
}
for (id in c("B17", "B18", "C06", "C07", "C08", "C09", "C10", "C11")) {
  path <- profile_pair_target_asset_path_v2(".", id)
  cat("Target", id, path, "exists", file.exists(path), "\n")
  if (!file.exists(path)) next
  asset <- readRDS(path)
  cat("Fields", names(asset), "\n")
  cat("Population", asset$n_population, "repeats", asset$repeats, "target length", length(asset$beta_star_mc), "\n")
}

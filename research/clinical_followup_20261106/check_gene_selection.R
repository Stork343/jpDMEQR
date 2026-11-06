library(data.table)
library(jsonlite)
library(yaml)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
study <- read_yaml(args[1L])
panel_data <- fread(file.path(study$input_directory, "features.csv.gz"))
genes <- grep("^P_gene_", names(panel_data), value = TRUE)
records <- list()
for (fold in seq_len(study$outer_folds)) {
  outer <- file.path(study$output_root, paste0("outer_", fold))
  training <- fread(file.path(outer, "training.csv"))
  recipes <- readRDS(file.path(outer, "recipes.rds"))
  stopifnot(all(training$role == "training"))
  for (key in names(recipes)) {
    samples <- if (key == "training") training else training[inner_fold != as.integer(sub("inner_", "", key))]
    panels <- as.matrix(panel_data[match(samples$sample_id, sample_id), ..genes])
    observed_mask <- is.finite(panels)
    group <- samples$strain_id
    visits <- as.numeric(table(group)[group])
    observed <- colMeans(rowsum(observed_mask / visits, group))
    panels[!observed_mask] <- 0
    mean_value <- colMeans(rowsum(panels / visits, group)) / observed
    residual <- sweep(panels, 2L, mean_value, "-")
    residual[!observed_mask] <- 0
    variance <- colMeans(rowsum(residual^2 / visits, group)) / observed
    eligible <- genes[observed >= study$minimum_observed_fraction & is.finite(variance)]
    expected <- head(eligible[order(-variance[eligible], eligible, method = "radix")], study$maximum_high_features)
    recipe <- recipes[[key]]$direction_recipe$base_recipe
    selected <- intersect(genes, recipe$features)
    stopifnot(setequal(selected, expected), setequal(recipe$training_strains, samples$strain_id))
    records[[length(records)+1L]] <- data.table(fold, panel = key, patients = uniqueN(group), genes = length(selected))
  }
}
summary <- rbindlist(records)
fwrite(summary, file.path(study$report_directory, "gene_selection_checks.csv"))
write_json(list(status = "passed", verified_panels = nrow(summary), genes_per_panel = unique(summary$genes),
  reserved_outcomes_used = FALSE), file.path(study$report_directory, "gene_selection_verification.json"),
  auto_unbox = TRUE, pretty = TRUE)

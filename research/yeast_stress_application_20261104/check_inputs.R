library(data.table)
library(jsonlite)
library(yaml)
library(BEDMatrix)
source("research/yeast_1011_application_20261104/features.R")
source("research/yeast_stress_application_20261104/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")

config <- read_yaml("config/discovery_yeast_stress_20261104.yml")
panel <- y1011_panel(config$genomic_directory)
previous <- readRDS("results/raw/application/yeast_1011_20261104/discovery/NaCl/JOINT_prediction_tau025/fit.rds")
base <- y1011_genomic_recipe(previous$recipe$training_strains, panel, config)
shared <- intersect(names(base), names(previous$recipe))
for (key in shared) stopifnot(isTRUE(all.equal(base[[key]], previous$recipe[[key]], tolerance = 1e-12)))
records <- list()
for (assay_name in c("Anisomycin", "Caffeine")) {
  rows <- fread(file.path(config$input_directory, assay_name, "discovery.csv"))[exploration_role == "training"]
  current <- config
  current$conditions <- config$assays[[assay_name]]
  recipe <- stress_recipe(rows, panel, current)
  sample <- rows[strain_id %in% head(sort(unique(rows$strain_id)), 40L)]
  features <- stress_features(sample, panel, recipe)
  d <- length(recipe$low_columns)
  p <- recipe$base_dimension
  condition <- recipe$contrasts[sample$condition, , drop = FALSE]
  error <- 0
  for (j in seq_len(ncol(condition))) {
    error <- max(error, abs(features$X[, d + j * p + seq_len(p)] - features$X[, d + seq_len(p)] * condition[, j]))
  }
  orthogonality <- max(abs(crossprod(recipe$contrasts) / nrow(recipe$contrasts) - diag(ncol(recipe$contrasts))))
  stopifnot(error < 1e-12, orthogonality < 1e-12,
            max(abs(rowsum(features$weight, features$cluster_id) - 1)) < 1e-12)
  check <- features
  keep <- c(seq_len(d), d + seq_len(8L), d + length(recipe$markers) + match("X2730.YDR038C_NumOfGenes_3", recipe$copy_names))
  check$X <- features$X[, keep, drop = FALSE]
  check$penalty_factor <- features$penalty_factor[keep]
  fit <- ep_path(check, "JOINT", 0.25, config, 2L, list(Lambda = 1, size = 2L))
  gaps <- vapply(fit$diagnostics[[1L]], `[[`, numeric(1), "gap")
  stopifnot(max(gaps) <= config$gap_tolerance)
  records[[assay_name]] <- data.table(assay = assay_name, conditions = length(current$conditions),
    variables = ncol(features$X), interaction_error = error, contrast_error = orthogonality, maximum_gap = max(gaps))
}
output <- "results/preflight/empirical_selection_20261101/yeast_stress"
fwrite(rbindlist(records), file.path(output, "input_checks.csv"))
write_json(list(status = "passed", genomic_recipe_matches_saved_fit = TRUE,
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "input_checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(records))

library(data.table)
library(jsonlite)
library(BEDMatrix)
source("research/yeast_1011_application_20261104/features.R")

root <- "data/derived/yeast_proteome_20261104"
destination <- file.path(root, "population_controls.rds")
if (file.exists(destination)) stop("遗传背景调整数据已经存在")
roles <- fread(file.path(root, "roles.csv"))
ids <- roles[role == "discovery"]$strain_id
panel <- y1011_panel("data/derived/yeast_1011_20261104")
source_directory <- "results/raw/application/yeast_1011_20261104/discovery/NaCl/JOINT_prediction_tau025"
records <- list()
checks <- list()
for (key in c("training", "inner_1", "inner_2", "inner_3")) {
  filename <- if (key == "training") "fit.rds" else paste0(key, ".rds")
  saved <- readRDS(file.path(source_directory, filename))
  recipe <- saved$recipe[c("panel_key", "markers", "g_center", "g_scale", "loadings", "training_strains")]
  stopifnot(identical(recipe$panel_key, key))
  expected <- if (key == "training") {
    roles[exploration_role == "training"]$strain_id
  } else {
    fold <- as.integer(sub("inner_", "", key))
    roles[exploration_role == "training" & inner_fold != fold]$strain_id
  }
  stopifnot(all(expected %in% recipe$training_strains),
            !length(intersect(recipe$training_strains, roles[role == "validation"]$strain_id)))
  if (key != "training") stopifnot(!length(intersect(recipe$training_strains,
    roles[exploration_role == "training" & inner_fold == fold]$strain_id)))
  scores <- y1011_scores(ids, panel, recipe)
  records[[key]] <- list(scores = scores, eligible_training_strains = sort(expected),
                          genotype_training_strains = recipe$training_strains,
                          source_recipe = recipe, source_file = file.path(source_directory, filename))
  checks[[key]] <- data.table(panel = key, genotype_training_strains = length(recipe$training_strains),
                             proteome_training_strains = length(expected), scored_discovery_strains = length(ids))
}
saveRDS(records, destination)
output <- "results/preflight/empirical_selection_20261101/yeast_proteome"
fwrite(rbindlist(checks), file.path(output, "population_control_checks.csv"))
write_json(list(status = "passed", fitted_from_genotype_training_roles = TRUE,
                formal_validation_outcomes_used = FALSE, formal_validation_genotypes_used_for_fitting = FALSE),
           file.path(output, "population_control_checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(checks))

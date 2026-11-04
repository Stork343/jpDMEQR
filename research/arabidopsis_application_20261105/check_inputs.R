library(data.table)
library(jsonlite)
library(yaml)
library(BEDMatrix)
source("research/arabidopsis_application_20261105/features.R")

root <- "data/derived/arabidopsis_20261105"
output <- "results/preflight/empirical_selection_20261101/arabidopsis"
config <- read_yaml("config/discovery_arabidopsis_20261105.yml")
panel <- arab_panel(root, load_recipes = FALSE)
destination <- file.path(root, "training_recipes.rds")
if (file.exists(destination)) stop("训练变换已经存在")
reference <- fread(file.path(root, "genotype_reader_reference.csv"), colClasses = c(accession_id = "character"))
original <- BEDMatrix::BEDMatrix(file.path(root, "full_genotypes"), simple_names = TRUE)
expected <- as.matrix(reference[, -1L, with = FALSE])
actual <- original[reference$accession_id, names(reference)[-1L], drop = FALSE]
genotype_error <- max(abs(actual - expected))
stopifnot(genotype_error == 0)
rows <- fread(file.path(root, "discovery.csv"), colClasses = c(accession_id = "character"))[exploration_role == "training"]
recipes <- list()
checks <- list()
for (key in names(panel$groups)) {
  current <- rows[accession_id %in% panel$groups[[key]]]
  recipe <- arab_compute_recipe(current, panel, config)
  scores <- arab_scores(recipe$training_accessions, panel, recipe)
  center_error <- max(abs(colMeans(scores)))
  scale_error <- max(abs(apply(scores, 2L, sd) - 1))
  sample <- current[accession_id %in% head(recipe$training_accessions, 6L)]
  features <- arab_features(sample, panel, recipe)
  d <- length(recipe$low_columns)
  p <- length(recipe$markers)
  interaction_error <- max(abs(features$X[, d + p + seq_len(p)] -
                                features$X[, d + seq_len(p)] * (2 * sample$high_temperature - 1)))
  stopifnot(center_error < 1e-8, scale_error < 1e-8, interaction_error < 1e-12,
            max(abs(rowsum(features$weight, features$cluster_id) - 1)) < 1e-12)
  recipes[[key]] <- recipe
  checks[[key]] <- data.table(panel = key, accessions = length(recipe$training_accessions),
    markers = p, variables = d + 2L * p, genotype_error = genotype_error,
    pc_center_error = center_error, pc_scale_error = scale_error, interaction_error = interaction_error)
  cat(key, "completed\n")
  flush.console()
}
saveRDS(recipes, destination)
fwrite(rbindlist(checks), file.path(output, "input_checks.csv"))
write_json(list(status = "passed", genotype_calls_checked = length(expected), formal_validation_used = FALSE,
                manuscript_adoption = FALSE), file.path(output, "input_checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(checks))

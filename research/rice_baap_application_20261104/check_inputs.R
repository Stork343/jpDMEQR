library(data.table)
library(jsonlite)
library(yaml)
library(BEDMatrix)
source("research/rice_baap_application_20261104/features.R")

root <- "data/derived/rice_baap_20261104"
panel <- rice_panel(root)
config <- read_yaml("config/discovery_rice_baap_20261104.yml")
reference <- fread("results/scratch/rice_baap_20261104/bed_reader_reference.csv")
expected <- as.matrix(reference[, -1L, with = FALSE])
actual <- panel$beds$training[reference$accession_id, names(reference)[-1L], drop = FALSE]
genotype_error <- max(abs(actual - expected))
stopifnot(genotype_error == 0)
rows <- fread(file.path(root, "discovery.csv"))[exploration_role == "training"]
recipe <- rice_recipe(rows, panel, config)
sample <- rows[accession_id %in% head(sort(unique(rows$accession_id)), 3L)]
features <- rice_features(sample, panel, recipe)
selected <- unique(as.integer(quantile(seq_along(recipe$markers), seq(0, 1, length.out = 64L))))
G <- panel$beds$training[sample$accession_id, recipe$markers[selected], drop = FALSE]
G <- sweep(sweep(G, 2L, recipe$center[selected], "-"), 2L, recipe$scale[selected], "/")
interaction <- sweep(sweep(G * (as.numeric(sample$irrigation == "AWD") - 0.5),
                            2L, recipe$interaction_center[selected], "-"),
                      2L, recipe$interaction_scale[selected], "/")
d <- length(recipe$low_columns)
p <- length(recipe$markers)
design_error <- max(abs(features$X[, d + selected] - G), abs(features$X[, d + p + selected] - interaction))
pcs <- rice_scores(recipe$training_accessions, panel, recipe)
pc_center_error <- max(abs(colMeans(pcs)))
pc_scale_error <- max(abs(apply(pcs, 2L, sd) - 1))
stopifnot(design_error < 1e-12, pc_center_error < 1e-8, pc_scale_error < 1e-8,
          all(rows$role == "discovery"), ncol(features$Z) == 2L,
          max(abs(rowsum(features$weight, features$cluster_id) - 1)) < 1e-12)
output <- "results/preflight/empirical_selection_20261101/rice_baap"
write_json(list(status = "passed", actual_genotypes_checked = length(expected),
                genotype_reader_error = genotype_error, design_error = design_error,
                population_component_center_error = pc_center_error,
                population_component_scale_error = pc_scale_error,
                training_markers = p, prediction_columns = ncol(features$X),
                nuisance_dimensions = ncol(features$Z), formal_validation_used = FALSE),
           file.path(output, "input_checks.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("输入检查完成\n")

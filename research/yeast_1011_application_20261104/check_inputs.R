library(data.table)
library(jsonlite)
library(yaml)
library(BEDMatrix)
source("research/yeast_1011_application_20261104/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")

root <- "data/derived/yeast_1011_20261104"
panel <- y1011_panel(root)
config <- read_yaml("config/discovery_yeast_1011_20261104.yml")
reference <- fread(file.path(root, "genotype_reader_reference.csv"))
original <- BEDMatrix::BEDMatrix("data/raw/yeast_1011_20261104/gwas/1011GWAS_matrix", simple_names = TRUE)
expected <- as.matrix(reference[, -1L, with = FALSE])
actual <- original[reference$IID, names(reference)[-1L], drop = FALSE]
stopifnot(identical(unname(is.na(actual)), unname(is.na(expected))))
genotype_error <- max(abs(actual - expected), na.rm = TRUE)
stopifnot(genotype_error == 0)
rows <- fread(file.path(root, "discovery.csv"))[exploration_role == "training"]
roles <- fread(file.path(root, "roles.csv"))
stopifnot(!anyDuplicated(roles$strain_id), all(rows$role == "discovery"),
          !length(intersect(rows$strain_id, roles[role == "validation"]$strain_id)))
recipe <- y1011_recipe(rows, panel, config)
sample <- rows[strain_id %in% head(sort(unique(rows$strain_id)), 40L)]
features <- y1011_features(sample, panel, recipe)
pcs <- y1011_scores(recipe$training_strains, panel, recipe)
pc_center_error <- max(abs(colMeans(pcs)))
pc_scale_error <- max(abs(apply(pcs, 2L, sd) - 1))
d <- length(recipe$low_columns)
p <- length(recipe$markers) + length(recipe$copy_names)
interaction_error <- max(abs(features$X[, d + p + seq_len(p)] -
                              features$X[, d + seq_len(p)] * (2 * sample$high_salt - 1)))
ena <- match("X2730.YDR038C_NumOfGenes_3", recipe$copy_names)
stopifnot(!is.na(ena))
raw_copy <- fread("data/raw/yeast_1011_20261104/genesMatrix_CopyNumber.tab.gz")
source_ids <- sub("^SACE_", "", raw_copy[[1L]])
source_values <- raw_copy[[recipe$copy_names[ena]]][match(sample$strain_id, source_ids)]
expected_copy <- (log2(1 + source_values) - recipe$c_center[ena]) / recipe$c_scale[ena]
copy_error <- max(abs(features$X[, d + length(recipe$markers) + ena] - expected_copy))
stopifnot(interaction_error < 1e-12, copy_error < 1e-12,
          pc_center_error < 1e-8, pc_scale_error < 1e-8,
          max(abs(rowsum(features$weight, features$cluster_id) - 1)) < 1e-12)
check_features <- features
keep <- c(seq_len(d), d + seq_len(8L), d + length(recipe$markers) + ena)
check_features$X <- features$X[, keep, drop = FALSE]
check_features$penalty_factor <- features$penalty_factor[keep]
fit <- ep_path(check_features, "JOINT", 0.25, config, 2L, list(Lambda = 1, size = 2L))
stopifnot(all(is.finite(fit$beta)))
output <- "results/preflight/empirical_selection_20261101/yeast_1011"
saveRDS(list(fit = fit, sample_ids = sample$sample_id), file.path(output, "numerical_input_fit.rds"))
write_json(list(status = "passed", genotypes_checked = length(expected), genotype_error = genotype_error,
                interaction_error = interaction_error, copy_transform_error = copy_error,
                population_component_center_error = pc_center_error, population_component_scale_error = pc_scale_error,
                training_markers = length(recipe$markers), training_copy_families = length(recipe$copy_names),
                prediction_columns = d + 2L * p, adjustment_columns = d,
                formal_validation_used = FALSE),
           file.path(output, "input_checks.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(fit$diagnostics)

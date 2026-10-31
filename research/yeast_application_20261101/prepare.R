library(data.table)
library(jsonlite)
library(yaml)

directory <- "data/derived/yeast_20261101"
audit_dir <- "results/preflight/empirical_selection_20261101/yeast_discovery"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(directory, "panel.rds"))) stop("数据对象已经存在")
config <- read_yaml("config/discovery_yeast_20261101.yml")
original <- fread("data/raw/yeast_bloom_2013/BYxRM_PhenoData.txt")
ids <- original$V1
stopifnot(length(ids) == 1008L, !anyDuplicated(ids))
objects <- new.env()
load("data/raw/yeast_bloom_2013/cross_original.RData", objects)
identity_error <- max(abs(as.matrix(original[, -1L, with = FALSE]) - as.matrix(objects$cross$pheno)), na.rm = TRUE)
stopifnot(identity_error < 1e-10)
genotypes <- 2 * do.call(cbind, lapply(objects$cross$geno, `[[`, "data")) - 3
rownames(genotypes) <- ids
colnames(genotypes) <- paste0("G_", colnames(genotypes))
load("data/raw/yeast_bloom_2013/pheno_raw.RData", objects)
pheno_raw <- objects$pheno_raw
source_dir <- "results/scratch/empirical_selection_20261101/eqtl_source/RData"
load(file.path(source_dir, "log2_t.tpm.matrix.RData"), objects)
load(file.path(source_dir, "covariates.OD.RData"), objects)
expression <- objects$t.tpm.matrix
covariates <- objects$covariates.OD
expression_ids <- sub("-.*$", "", rownames(expression))
covariate_ids <- sub("-.*$", "", rownames(covariates))
stopifnot(!anyDuplicated(expression_ids), !anyDuplicated(covariate_ids),
          identical(expression_ids, covariate_ids), all(is.finite(expression)),
          all(is.finite(covariates)))
rownames(expression) <- expression_ids
rownames(covariates) <- covariate_ids
available <- ids[ids %in% expression_ids]
X <- cbind(expression[available, , drop = FALSE], genotypes[available, , drop = FALSE])
C <- covariates[available, , drop = FALSE]
stopifnot(all(C[, 1L] == 1), all(is.finite(X)), !anyDuplicated(colnames(X)))
set.seed(config$role_seed)
roles <- data.table(strain_id = sort(ids))
n_discovery <- floor(nrow(roles) * config$discovery_fraction)
roles[, role := sample(rep(c("discovery", "validation"), c(n_discovery, .N - n_discovery)))]
discovery_n <- sum(roles$role == "discovery")
training_n <- floor(discovery_n * config$exploration_training_fraction)
roles[role == "discovery", exploration_role := sample(rep(c("training", "evaluation"), c(training_n, .N - training_n)))]
roles[exploration_role == "training", inner_fold := sample(rep(seq_len(config$inner_folds), length.out = .N))]
outcomes <- rbindlist(lapply(names(pheno_raw), function(trait) {
  measurements <- pheno_raw[[trait]]
  rbindlist(lapply(names(measurements), function(id) {
    values <- measurements[[id]]
    data.table(strain_id = id, trait = trait, replicate = seq_along(values), y = as.numeric(values))
  }))
}))
outcomes <- merge(outcomes, roles, by = "strain_id", all.x = TRUE)
source_counts <- outcomes[, .(rows = .N, strains = uniqueN(strain_id),
                              finite_rows = sum(is.finite(y)),
                              usable_rows = sum(is.finite(y) & strain_id %in% available & !is.na(role))), by = trait]
fwrite(source_counts, file.path(audit_dir, "measurement_counts.csv"))
outcomes <- outcomes[is.finite(y) & strain_id %in% available & !is.na(role)]
setorder(outcomes, trait, strain_id, replicate)
counts <- outcomes[, .(rows = .N, strains = uniqueN(strain_id)), by = .(trait, role, exploration_role)]
fwrite(counts, file.path(audit_dir, "cohort_counts.csv"))
fwrite(roles, file.path(directory, "roles.csv"))
saveRDS(list(X = X, C = C, expression_columns = colnames(expression),
             genotype_columns = colnames(genotypes), roles = roles),
        file.path(directory, "panel.rds"), compress = FALSE)
saveRDS(outcomes[role == "discovery"], file.path(directory, "growth_discovery.rds"))
saveRDS(outcomes[role == "validation"], file.path(directory, "growth_validation_reserved.rds"))
write_json(list(status = "prepared", original_strains = length(ids),
                 expression_matched_strains = length(available), expression_variables = ncol(expression),
                 genotype_variables = ncol(genotypes), source_mean_identity_error = identity_error,
                 grouping = "one_segregant_per_tetrad", phenotype_values_imputed = FALSE,
                 formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(audit_dir, "preparation.json"), auto_unbox = TRUE, pretty = TRUE)
print(counts[trait == config$initial_trait])

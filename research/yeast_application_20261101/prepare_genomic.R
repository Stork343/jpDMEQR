library(data.table)
library(jsonlite)

directory <- "data/derived/yeast_genomic_20261101"
audit_dir <- "results/preflight/empirical_selection_20261101/yeast_genomic"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(directory, "panel.rds"))) stop("遗传标记数据对象已经存在")
original <- fread("data/raw/yeast_bloom_2013/BYxRM_PhenoData.txt")
objects <- new.env()
load("data/raw/yeast_bloom_2013/cross_original.RData", objects)
stopifnot(max(abs(as.matrix(original[, -1L, with = FALSE]) - as.matrix(objects$cross$pheno)), na.rm = TRUE) < 1e-10)
X <- 2 * do.call(cbind, lapply(objects$cross$geno, `[[`, "data")) - 3
rownames(X) <- original$V1
colnames(X) <- paste0("G_", colnames(X))
C <- matrix(1, nrow(X), 1L, dimnames = list(rownames(X), "(Intercept)"))
roles <- fread("data/derived/yeast_20261101/roles.csv")
load("data/raw/yeast_bloom_2013/pheno_raw.RData", objects)
outcomes <- rbindlist(lapply(names(objects$pheno_raw), function(trait_name) {
  measurements <- objects$pheno_raw[[trait_name]]
  rbindlist(lapply(names(measurements), function(id) data.table(
    strain_id = id, trait = trait_name, replicate = seq_along(measurements[[id]]), y = as.numeric(measurements[[id]]))))
}))
outcomes <- merge(outcomes, roles, by = "strain_id", all.x = TRUE)
outcomes <- outcomes[is.finite(y) & strain_id %in% rownames(X) & !is.na(role)]
setorder(outcomes, trait, strain_id, replicate)
saveRDS(list(X = X, C = C, roles = roles), file.path(directory, "panel.rds"), compress = FALSE)
saveRDS(outcomes[role == "discovery"], file.path(directory, "growth_discovery.rds"))
saveRDS(outcomes[role == "validation"], file.path(directory, "growth_validation_reserved.rds"))
counts <- outcomes[, .(rows = .N, strains = uniqueN(strain_id)), by = .(trait, role, exploration_role)]
fwrite(counts, file.path(audit_dir, "cohort_counts.csv"))
write_json(list(status = "prepared", strains = nrow(X), genetic_markers = ncol(X),
                grouping = "one_segregant_per_tetrad", role_assignments_preserved = TRUE,
                phenotype_values_imputed = FALSE, formal_validation_used = FALSE,
                manuscript_adoption = FALSE), file.path(audit_dir, "preparation.json"), auto_unbox = TRUE, pretty = TRUE)
print(counts[trait == "Ethanol"])

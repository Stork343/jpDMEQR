root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())
source_dir <- file.path(root, "data", "raw", "public_modules", "banchereau2016_20261022")
rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
expression <- readRDS(file.path(root, "data", "interim", "GSE65391", "expression_processed_probe.rds"))
metadata <- read.csv(file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv"),
                     stringsAsFactors = FALSE)
audit <- audit_primary_outcome_v2(metadata)
metadata <- metadata[audit$eligible_visit, , drop = FALSE]
rownames(metadata) <- metadata$sample_id
expression <- expression[, rownames(metadata), drop = FALSE]
train <- metadata$subject_id %in% sort(unique(metadata$subject_id))[seq_len(60L)]
scores <- score_frozen_module_representations_v2(expression, annotation, signatures, members, train, rule)
for (id in rule$module_order) {
  rows <- members[members$module_id == id, ]
  if (rows$member_type[1] == "ILMN_probe") {
    direct <- colMeans(expression[rows$member_id, , drop = FALSE])
  } else {
    map <- scores$gene_probe_map[match(rows$member_id, scores$gene_probe_map$gene_symbol), ]
    direct <- colMeans(expression[map$probe_id, , drop = FALSE])
  }
  if (max(abs(direct - scores$raw_scores[, id])) > 1e-10) stop("The weighted module mean differs from direct measurement means.")
}
dimensions <- do.call(rbind, lapply(c(FALSE, TRUE), function(include_sledai) {
  design <- prepare_gse65391_lowdim_design_v2(
    metadata, expression, annotation, signatures, rep(TRUE, nrow(metadata)), rule,
    include_sledai = include_sledai, module_score_members = members
  )
  data.frame(model = if (include_sledai) "B" else "A", n_subjects = length(unique(metadata$subject_id)),
             n_visits = nrow(metadata), design_dimension = ncol(design$X),
             design_rank = qr(design$X)$rank,
             dimension_limit = floor(length(unique(metadata$subject_id)) /
                                       log(length(unique(metadata$subject_id)))))
}))
write.csv(scores$audit, file.path(source_dir, "score_representation_coverage.csv"), row.names = FALSE)
write.csv(dimensions, file.path(source_dir, "clinical_module_design_dimensions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(source_dir, "score_check_sessionInfo.txt"))
print(scores$audit)
print(dimensions)

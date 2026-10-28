root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())
args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args)) args[1] else paste0("workflow_preflight_", format(Sys.time(), "%Y%m%dT%H%M%S"))
out_dir <- file.path(root, "results", "raw", "application", run_id)
if (file.exists(out_dir) || dir.exists(out_dir)) stop("The research workflow output already exists.")
rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
expression <- readRDS(file.path(root, "data", "interim", "GSE65391", "expression_processed_probe.rds"))
metadata <- read.csv(file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv"),
                     stringsAsFactors = FALSE)
audit <- audit_primary_outcome_v2(metadata)
eligible <- metadata[audit$eligible_visit, , drop = FALSE]
subjects <- unlist(lapply(split(eligible$subject_id, eligible$batch), function(ids) {
  utils::head(sort(unique(ids)), 6L)
}), use.names = FALSE)
if (length(subjects) != 12L || anyDuplicated(subjects)) {
  stop("The declared research check requires six distinct subjects from each source batch.")
}
metadata <- metadata[audit$eligible_visit & metadata$subject_id %in% subjects, , drop = FALSE]
rownames(metadata) <- metadata$sample_id
expression <- expression[, rownames(metadata), drop = FALSE]
# 固定坐标用于研究计算检查，全部module preprocessing仍逐次重新计算。
coordinates <- c("(Intercept)", "module__M1.2")
result <- lowdim_profile_pair_workflow_inference_v2(
  metadata, expression, annotation, signatures, members, rule,
  tau = 0.5, h = 12^(-3 / 10), coordinate_subset = coordinates,
  checkpoint_dir = file.path(out_dir, "workflow_checkpoints")
)
checkpoint_files <- list.files(file.path(out_dir, "workflow_checkpoints"), full.names = TRUE)
if (sum(grepl("_preprocessing[.]rds$", checkpoint_files)) != 13L ||
    sum(grepl("(full_sample|delete_subject_[0-9]+)[.]rds$", checkpoint_files)) != 13L) {
  stop("Every workflow stage must have both preprocessing and complete-fit checkpoints.")
}
full_map <- result$preprocessing$gene_preprocessing$gene_probe_map
probe_changes <- vapply(result$jackknife$leave_one, function(stage) {
  map <- stage$preprocessing$gene_preprocessing$gene_probe_map
  sum(full_map$probe_id != map$probe_id[match(full_map$gene_symbol, map$gene_symbol)])
}, integer(1))
saveRDS(result, file.path(out_dir, "workflow_result.rds"))
write.csv(result$table, file.path(out_dir, "research_coordinate_computation.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
jsonlite::write_json(list(
  run_id = run_id, research_check = TRUE, coordinate_subset = coordinates,
  subject_selection = "first_six_distinct_eligible_subject_ids_per_source_batch",
  n_clusters = 12L, n_visits = nrow(metadata), jackknife_completed = nrow(result$jackknife$estimates),
  h = 12^(-3 / 10), h_half = 12^(-3 / 10) / 2,
  method_version = result$fit$method_version,
  parameterisation = result$jackknife$parameterisation,
  jackknife_scope = result$jackknife$jackknife_scope,
  preprocessing_checkpoints = 13L, complete_fit_checkpoints = 13L,
  max_changed_gene_probes = max(probe_changes),
  max_refit_gradient = max(result$jackknife$diagnostics$refit_gradient),
  max_identity_error = max(result$jackknife$diagnostics$profile_identity_error),
  min_paired_eigenvalue = min(result$jackknife$diagnostics$hessian_min_eigenvalue),
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(out_dir, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(result$jackknife$diagnostics)
message("Research workflow check: ", normalizePath(out_dir))

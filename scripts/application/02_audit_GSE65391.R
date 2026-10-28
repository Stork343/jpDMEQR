full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/application/02_audit_GSE65391.R"
root <- if (file.exists("R/application_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
metadata_only <- identical(get_arg("metadata-only", "false"), "true")
interim_dir <- file.path(root, "data", "interim", "GSE65391")
legacy_dir <- file.path(root, "data", "derived", "GSE65391")
run_id <- get_arg("run-id", paste0("audit_", format(Sys.time(), "%Y%m%dT%H%M%S")))
derived_dir <- get_arg("output-dir", file.path(legacy_dir, run_id))
if (dir.exists(derived_dir)) stop("AUDIT_OUTPUT_EXISTS: choose a new audit output directory.")
meta_path <- file.path(interim_dir, "phenotype_parsed.csv")
expr_path <- file.path(interim_dir, "expression_processed_probe.rds")
if (!file.exists(meta_path)) stop("Run 01_build_GSE65391.R before auditing.")
metadata <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
rownames(metadata) <- metadata$geo_accession
expression_sample_ids <- NULL
source_metadata <- NULL
if (!metadata_only) {
  eset_path <- file.path(interim_dir, "gse65391_eset.rds")
  if (!file.exists(eset_path)) stop("Original GEO ExpressionSet is required for the full audit.")
  eset <- readRDS(eset_path)
  source_metadata <- Biobase::pData(eset)
  expression_sample_ids <- Biobase::sampleNames(eset)
  if (!identical(as.character(source_metadata[metadata$geo_accession, "batch_replicate:ch1"]),
                 toupper(as.character(metadata$batch_replicate_flag)))) {
    stop("Original GEO replicate flags differ from the parsed metadata.")
  }
}
config_path <- file.path(root, "config", "application.yml")
config <- yaml::read_yaml(config_path)
curvature_path <- get_arg("curvature-rule")
curvature_rule <- if (is.null(curvature_path)) NULL else yaml::read_yaml(curvature_path)
audit <- do.call(audit_primary_outcome_v2, c(
  list(metadata = metadata, expression_sample_ids = expression_sample_ids,
       minimum_visits_per_subject = config$cluster$minimum_visits,
       primary_quantiles = unlist(config$quantiles$primary),
       sensitivity_quantiles = unlist(config$quantiles$sensitivity),
       bandwidth_multipliers = unlist(config$method$h_multiplier_grid),
       curvature_rule = curvature_rule),
  config$primary_outcome$audit_gates
))
replicates <- audit_gse65391_replicate_relations_v2(metadata, source_metadata)
exclusions <- gse65391_exclusion_reasons_v2(metadata, audit$eligible_visit)
exclusion_summary <- as.data.frame(table(exclusions$reason))
names(exclusion_summary) <- c("reason", "samples")
dir.create(derived_dir, recursive = TRUE)
write.csv(audit$gates, file.path(derived_dir, "outcome_audit_gates.csv"), row.names = FALSE)
write.csv(audit$summary, file.path(derived_dir, "outcome_audit_summary.csv"), row.names = FALSE)
write.csv(audit$identifiers$conflicts, file.path(derived_dir, "identifier_conflicts.csv"), row.names = FALSE)
write.csv(audit$identifiers$summary, file.path(derived_dir, "identifier_summary.csv"), row.names = FALSE)
write.csv(audit$identifiers$canonical_keys, file.path(derived_dir, "canonical_visit_keys.csv"), row.names = FALSE)
write.csv(audit$batch$counts, file.path(derived_dir, "batch_outcome_counts.csv"), row.names = FALSE)
write.csv(audit$batch$summary, file.path(derived_dir, "batch_outcome_summary.csv"), row.names = FALSE)
write.csv(audit$curvature$diagnostics, file.path(derived_dir, "local_curvature.csv"), row.names = FALSE)
if (!is.null(audit$curvature$bootstrap_means)) {
  saveRDS(list(means = audit$curvature$bootstrap_means,
               rule = audit$curvature$rule, rng_kind = audit$curvature$rng_kind),
          file.path(derived_dir, "curvature_bootstrap.rds"))
}
if (!is.null(curvature_path)) {
  file.copy(curvature_path, file.path(derived_dir, "curvature_rule.yml"), overwrite = FALSE)
}
write.csv(exclusions, file.path(derived_dir, "cohort_exclusions.csv"), row.names = FALSE)
write.csv(exclusion_summary, file.path(derived_dir, "cohort_flow.csv"), row.names = FALSE)
write.csv(replicates$records, file.path(derived_dir, "replicate_record_relations.csv"), row.names = FALSE)
write.csv(replicates$summary, file.path(derived_dir, "replicate_summary.csv"), row.names = FALSE)
writeLines(replicates$source_processing, file.path(derived_dir, "source_data_processing.txt"))
legacy_summary_path <- file.path(legacy_dir, "outcome_audit_summary.csv")
if (file.exists(legacy_summary_path)) {
  previous <- read.csv(legacy_summary_path)
  fields <- intersect(names(previous), names(audit$summary))
  comparison <- data.frame(
    quantity = fields, previous = as.character(previous[1, fields]),
    current = as.character(audit$summary[1, fields]), stringsAsFactors = FALSE
  )
  write.csv(comparison, file.path(derived_dir, "cohort_comparison.csv"), row.names = FALSE)
}

visit_counts <- as.data.frame(table(metadata$subject_id[audit$eligible_visit]),
                              stringsAsFactors = FALSE)
names(visit_counts) <- c("subject_id", "eligible_visits")
write.csv(visit_counts, file.path(derived_dir, "visit_distribution.csv"), row.names = FALSE)

missingness <- data.frame(
  variable = names(metadata),
  missing_n = vapply(metadata, function(x) sum(is.na(x) | (is.character(x) & !nzchar(x))), integer(1)),
  missing_fraction = vapply(metadata, function(x) mean(is.na(x) | (is.character(x) & !nzchar(x))), numeric(1)),
  stringsAsFactors = FALSE
)
write.csv(missingness, file.path(derived_dir, "missingness_by_variable.csv"), row.names = FALSE)

eligible_metadata <- metadata[audit$eligible_visit, , drop = FALSE]
eligible_metadata$log_c3 <- audit$transformed_outcome[audit$eligible_visit]
eligible_metadata$canonical_visit_key <- audit$identifiers$canonical_keys$canonical_visit_key[audit$eligible_visit]
eligible_metadata$key_source <- audit$identifiers$canonical_keys$key_source[audit$eligible_visit]
folds <- if (isTRUE(audit$summary$all_gates_pass)) make_repeated_subject_folds_v2(
  metadata, audit$eligible_visit, k = config$validation$outer_folds,
  repetitions = config$validation$repetitions, seed = config$validation$base_seed
) else data.frame()
write.csv(folds, file.path(derived_dir, "subject_folds.csv"), row.names = FALSE)
write.csv(eligible_metadata, file.path(derived_dir, "eligible_visits.csv"), row.names = FALSE)
saveRDS(list(metadata = metadata, audit = audit, eligible_metadata = eligible_metadata,
             expression_file = normalizePath(expr_path, mustWork = TRUE),
             subject_folds = folds, replicates = replicates, config = config),
        file.path(derived_dir, "analysis_cohort.rds"))
capture.output(sessionInfo(), file = file.path(derived_dir, "audit_sessionInfo.txt"))

print(audit$gates)
print(audit$summary)
message("Audit outputs: ", normalizePath(derived_dir))
if (!isTRUE(audit$summary$all_gates_pass)) {
  stop("EMPIRICAL_GATE_INCOMPLETE: inspect the versioned audit report before fitting.")
}

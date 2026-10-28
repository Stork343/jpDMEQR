root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())
args <- commandArgs(trailingOnly = TRUE)
cohort_dir <- if (length(args)) args[1] else file.path(root, "data", "derived", "GSE65391", "audit_20261022_bootstrap_v1")
output <- file.path(root, "data", "derived", "GSE65391", "subject_design_20261022_v1")
if (dir.exists(output)) stop("The versioned subject design already exists.")
cohort <- readRDS(file.path(cohort_dir, "analysis_cohort.rds"))
metadata <- cohort$eligible_metadata
folds <- cohort$subject_folds
if (!isTRUE(cohort$audit$summary$all_gates_pass) ||
    nrow(metadata) != 867L || length(unique(metadata$subject_id)) != 129L) {
  stop("The complete frozen 129-subject cohort and its C3 audit are required.")
}
summary <- validate_gse65391_subject_folds_v2(metadata, folds)
dir.create(output, recursive = TRUE)
write.csv(folds, file.path(output, "subject_folds.csv"), row.names = FALSE)
write.csv(summary, file.path(output, "outer_fold_summary.csv"), row.names = FALSE)
time <- do.call(rbind, lapply(split(metadata, metadata$subject_id), function(rows) {
  data.frame(subject_id = as.character(rows$subject_id[1]), visits = nrow(rows),
             distinct_time = length(unique(rows$cumulative_time)),
             time_range = diff(range(rows$cumulative_time)), stringsAsFactors = FALSE)
}))
write.csv(time, file.path(output, "subject_time_diagnostics.csv"), row.names = FALSE)
saveRDS(cohort, file.path(output, "analysis_inputs.rds"))
jsonlite::write_json(list(
  cohort_dir = normalizePath(cohort_dir), n_subjects = 129L, n_visits = 867L,
  outer_folds = 5L, repetitions = 5L, base_seed = 20260817L,
  subject_fold_md5 = unname(tools::md5sum(file.path(output, "subject_folds.csv"))),
  subjects_ge3_visits = sum(time$visits >= 3), subjects_with_time_variation = sum(time$distinct_time > 1L),
  matching_method_gate = "required_before_association_fitting",
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "subject_design_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(summary)

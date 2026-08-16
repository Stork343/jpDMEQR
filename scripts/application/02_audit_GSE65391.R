#!/usr/bin/env Rscript
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "scripts/application/02_audit_GSE65391.R"
root <- if (file.exists("R/application_v2.R")) "." else normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())

interim_dir <- file.path(root, "data", "interim", "GSE65391")
derived_dir <- file.path(root, "data", "derived", "GSE65391")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
meta_path <- file.path(interim_dir, "phenotype_parsed.csv")
expr_path <- file.path(interim_dir, "expression_processed_probe.rds")
if (!file.exists(meta_path) || !file.exists(expr_path)) stop("Run 01_build_GSE65391.R first.")
metadata <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
rownames(metadata) <- metadata$geo_accession
expression <- readRDS(expr_path)
metadata <- metadata[colnames(expression), , drop = FALSE]

audit <- audit_primary_outcome_v2(metadata)
write.csv(audit$gates, file.path(derived_dir, "outcome_audit_gates.csv"), row.names = FALSE)
write.csv(audit$summary, file.path(derived_dir, "outcome_audit_summary.csv"), row.names = FALSE)
write.csv(cohort_flow_gse65391_v2(metadata, audit$eligible_visit),
          file.path(derived_dir, "cohort_flow.csv"), row.names = FALSE)

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
eligible_metadata$log_c3 <- log(eligible_metadata$c3)
folds <- if (audit$summary$all_gates_pass) {
  make_repeated_subject_folds_v2(metadata, audit$eligible_visit,
                                 k = 5, repetitions = 5, seed = 20260817)
} else data.frame()
write.csv(folds, file.path(derived_dir, "subject_folds.csv"), row.names = FALSE)
write.csv(eligible_metadata, file.path(derived_dir, "eligible_visits.csv"), row.names = FALSE)
saveRDS(list(metadata = metadata, expression = expression, audit = audit,
             eligible_metadata = eligible_metadata, subject_folds = folds),
        file.path(derived_dir, "analysis_cohort.rds"), compress = "xz")
capture.output(sessionInfo(), file = file.path(derived_dir, "audit_sessionInfo.txt"))

print(audit$gates)
print(audit$summary)
if (!isTRUE(audit$summary$all_gates_pass)) {
  stop("OUTCOME_GATE_FAILED: C3 analysis is blocked. Inspect audit files; do not auto-switch outcomes.")
}
message("Outcome gate passed. Frozen folds and cohort files are under: ", derived_dir)

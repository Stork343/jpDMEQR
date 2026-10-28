args <- commandArgs(trailingOnly = TRUE)
geometry_dir <- if (length(args)) args[[1]] else file.path(
  "results", "scratch", "profile_pair_geometry_development_20261022", "strict_development05"
)
run_id <- if (length(args) >= 2L) args[[2]] else "integration_20261022_05"
geometry <- jsonlite::read_json(file.path(geometry_dir, "manifest.json"), simplifyVector = TRUE)
run_dir <- file.path("results", "raw", "simulation", run_id)
metrics <- read.csv(file.path(run_dir, "replication_metrics.csv"), stringsAsFactors = FALSE)
coordinates <- read.csv(file.path(run_dir, "coordinate_metrics.csv"), stringsAsFactors = FALSE)
if (!isTRUE(geometry$pass) || geometry$method_version != "profile_pair_20261022" ||
    any(!metrics$status %in% c("ok", "warning")) ||
    nrow(metrics) != 18L || nrow(coordinates) != 64L) {
  stop("Matching successful paired verification outputs are required.")
}
out_dir <- file.path("manuscript", "generated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
verification <- data.frame(
  quantity = c("Independent numerical seeds", "Differential geometry cases",
               "Complete fitting and inference cases", "Largest score finite difference error",
               "Largest Jacobian finite difference error", "Largest Schur identity error",
               "Largest paired variance identity error", "Registry method evaluations",
               "Registry coordinate records"),
  value = c(geometry$seed_count, geometry$case_count, geometry$fit_inference_count,
            geometry$maxima$gradient_error, geometry$maxima$hessian_error,
            geometry$maxima$profile_identity_error, geometry$maxima$sandwich_identity_error,
            nrow(metrics), nrow(coordinates)), stringsAsFactors = FALSE
)
values <- vapply(seq_len(nrow(verification)), function(i) {
  if (i %in% 4:7) paste0("\\texttt{", formatC(verification$value[i], format = "e", digits = 2), "}")
  else as.character(as.integer(verification$value[i]))
}, character(1))
writeLines(c(
  "\\begin{table}[t]",
  "\\caption{Numerical verification of the paired profile construction}\\label{tab:pair-verification}",
  "\\begin{tabularx}{\\textwidth}{@{}Xr@{}}", "\\toprule",
  "Calculation & Value \\\\", "\\midrule",
  paste0(verification$quantity, " & ", values, " \\\\"),
  "\\bottomrule", "\\end{tabularx}", "\\end{table}"
), file.path(out_dir, "tab_profile_pair_verification.tex"))
write.csv(verification, file.path(out_dir, "profile_pair_verification_source.csv"), row.names = FALSE)

cohort_dir <- file.path("data", "derived", "GSE65391", "audit_20261022_bootstrap_v1")
cohort <- read.csv(file.path(cohort_dir, "outcome_audit_summary.csv"))
gates <- read.csv(file.path(cohort_dir, "outcome_audit_gates.csv"))
module_dir <- file.path("data", "raw", "public_modules", "banchereau2016_20261022")
members <- read.csv(file.path(module_dir, "score_representation_coverage.csv"), stringsAsFactors = FALSE)
dimensions <- read.csv(file.path(module_dir, "clinical_module_design_dimensions.csv"))
definitions <- read.csv(file.path("config", "module_signatures.csv"), stringsAsFactors = FALSE)
if (nrow(cohort) != 1L || !cohort$all_gates_pass || !all(gates$pass) ||
    nrow(members) != 5L || !all(members$eligible) ||
    any(members$source_members != members$observed_members) ||
    any(dimensions$design_dimension != dimensions$design_rank)) {
  stop("Complete cohort and module verification is required.")
}
writeLines(c(
  sprintf("\\newcommand{\\SLEEligibleSubjects}{%d}", cohort$final_subjects),
  sprintf("\\newcommand{\\SLEEligibleVisits}{%d}", cohort$final_visits),
  sprintf("\\newcommand{\\SLEModelADimension}{%d}", dimensions$design_dimension[dimensions$model == "A"]),
  sprintf("\\newcommand{\\SLEModelBDimension}{%d}", dimensions$design_dimension[dimensions$model == "B"])
), file.path(out_dir, "profile_pair_application_values.tex"))
cohort_rows <- data.frame(
  quantity = c("Source SLE subjects", "Source SLE expression records", "SLE visits with observed C3",
               "Eligible subjects with repeated visits", "Eligible visits", "Distinct observed C3 values",
               "C3 missingness among source SLE visits", "C3 audit conditions passed"),
  value = c(as.character(cohort$source_sle_subjects), as.character(cohort$source_sle_visits),
            as.character(cohort$outcome_observed_visits), as.character(cohort$final_subjects),
            as.character(cohort$final_visits), as.character(cohort$distinct_values),
            sprintf("%.2f\\%%", 100 * cohort$missing_fraction),
            paste0(sum(gates$pass), "/", nrow(gates))), stringsAsFactors = FALSE
)
writeLines(c(
  "\\begin{table}[t]", "\\caption{Observed cohort construction for the C3 analysis}\\label{tab:data-audit}",
  "\\begin{tabularx}{\\textwidth}{@{}Xr@{}}", "\\toprule", "Quantity & Value \\\\", "\\midrule",
  paste0(cohort_rows$quantity, " & ", cohort_rows$value, " \\\\"),
  "\\bottomrule", "\\end{tabularx}", "\\end{table}"
), file.path(out_dir, "tab_profile_pair_cohort.tex"))
write.csv(cohort_rows, file.path(out_dir, "profile_pair_cohort_source.csv"), row.names = FALSE)
labels <- vapply(members$module_id, function(id) {
  value <- unique(definitions$module_label[definitions$module_id == id])
  if (length(value) != 1L) stop("Module biological label is not unique.")
  value
}, character(1))
units <- ifelse(members$member_type == "gene_symbol", "Genes", "ILMN probes")
writeLines(c(
  "\\begin{table}[t]", "\\caption{Complete measured membership of the prespecified immune modules}\\label{tab:immune-modules}",
  "\\begin{tabularx}{\\textwidth}{@{}lXlr@{}}", "\\toprule",
  "Module & Biological programme & Measurement unit & Members \\\\", "\\midrule",
  paste0(members$module_id, " & ", labels, " & ", units, " & ", members$observed_members, " \\\\"),
  "\\bottomrule", "\\end{tabularx}", "\\end{table}"
), file.path(out_dir, "tab_profile_pair_modules.tex"))
jsonlite::write_json(list(
  method_version = geometry$method_version, score_definition = geometry$score_definition,
  geometry_manifest = file.path(geometry_dir, "manifest.json"), simulation_run = run_dir,
  cohort_dir = cohort_dir, module_coverage = file.path(module_dir, "score_representation_coverage.csv"),
  verification_only = TRUE, final_monte_carlo_results = FALSE,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(out_dir, "profile_pair_source_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
message("Paired verification and observed cohort tables generated.")

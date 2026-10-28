source("research/submission_completion_20261027/load_comparators.R")
methods <- commandArgs(trailingOnly = TRUE)
if (!length(methods)) stop("At least one comparator identifier is required.")
registry <- submission_registry("main")
root <- "results/scratch/submission_completion_20261027/comparator_preflight"
dir.create(root, recursive = TRUE, showWarnings = FALSE)
for (method in methods) {
  admitted <- which(vapply(strsplit(registry$methods, "|", fixed = TRUE), function(x) method %in% x, logical(1)))
  if (!length(admitted)) stop("The comparator has no registered setting.")
  config <- config_row_to_list_v2(registry[admitted[1L], ])
  config$methods <- method
  path <- file.path(root, paste0(method, ".rds"))
  if (file.exists(path)) stop("The comparator preflight already exists.")
  result <- run_one_replication_v2(".", config, 1L, base_seed = 202610272L,
    final = TRUE, allow_missing_benchmarks = FALSE, run_id = "submission_comparator_preflight_20261028")
  saveRDS(result, path)
  print(result$metrics[c("experiment_id", "method_id", "status", "converged", "runtime_sec", "failure_stage", "failure_message")])
  validate_paired_failed_coordinate_accounting_v2(result$metrics, result$coordinates)
}

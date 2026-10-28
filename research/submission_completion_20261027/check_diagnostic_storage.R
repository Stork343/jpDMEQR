source("research/submission_completion_20261027/diagnostic_storage.R")
root <- "results/raw/simulation/submission_20261028_main_sqr/replications"
paths <- head(list.files(root, "^result[.]rds$", recursive = TRUE, full.names = TRUE), 20L)
if (length(paths) < 10L) stop("Completed SQR records are required for storage verification.")
output <- "results/scratch/submission_completion_20261027/diagnostic_storage_check"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
records <- list()
for (i in seq_along(paths)) {
  original <- readRDS(paths[i])
  stored <- submission_compact_diagnostics(original)
  for (field in setdiff(names(original), "method_diagnostics")) {
    if (!identical(original[[field]], stored[[field]])) stop("Statistical outputs changed during diagnostic storage.")
  }
  left <- original$method_diagnostics[["SQR-DEBIASED-IID"]]
  right <- stored$method_diagnostics[["SQR-DEBIASED-IID"]]
  fields <- c("beta_hat", "directions", "source_precision", "score", "requested_target_coordinates", "native_warning_records")
  for (field in fields) if (!identical(left[[field]], right[[field]]))
    stop("An inferential or source diagnostic changed during storage.")
  if (!identical(left$tuning$candidates, right$tuning$candidates) ||
      !identical(left$tuning$gamma_cv$candidates, right$tuning$gamma_cv$candidates))
    stop("Tuning candidates changed during storage.")
  path <- file.path(output, sprintf("replication_%02d.rds", i))
  saveRDS(stored, path)
  if (!identical(stored, readRDS(path))) stop("The compact result did not reproduce after serialization.")
  records[[i]] <- data.frame(source = paths[i], original_bytes = file.info(paths[i])$size,
    stored_bytes = file.info(path)$size, unchanged_statistical_outputs = TRUE)
}
records <- do.call(rbind, records)
write.csv(records, file.path(output, "checks.csv"), row.names = FALSE)
print(colSums(records[c("original_bytes", "stored_bytes")]))

library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide one completed simulation queue.")
root <- args[1L]
if (!file.exists(file.path(root, "completed.rds"))) stop("The simulation queue is incomplete.")
request <- readRDS(file.path(root, "request.rds"))
jobs <- as.data.table(request$jobs)
registry <- as.data.table(request$registry)
methods <- rbindlist(lapply(seq_len(nrow(registry)), function(i) {
  names <- strsplit(registry$methods[i], "|", fixed = TRUE)[[1L]]
  if (identical(request$family, "joint")) names <- names[startsWith(names, "JOINT-PROFILE")]
  if (!is.null(request$methods)) names <- intersect(names, request$methods)
  data.table(experiment_id = registry$experiment_id[i], method_id = names)
}))
expected <- merge(jobs[, .(experiment_id, replicate)], methods, by = "experiment_id", allow.cartesian = TRUE)
metrics <- coordinates <- theory <- vector("list", nrow(jobs))
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i]
  path <- file.path(root, "replications", job$experiment_id, sprintf("replicate_%05d", job$replicate), "result.rds")
  if (!file.exists(path)) stop("Missing registered simulation record: ", path)
  result <- readRDS(path)
  metric <- as.data.table(result$metrics)
  if (any(metric$experiment_id != job$experiment_id) || any(metric$replicate != job$replicate) ||
      !setequal(metric$method_id, methods[experiment_id == job$experiment_id, method_id]))
    stop("A simulation record differs from its registered identity.")
  metric[, component_nonconvergence_warnings := vapply(method_id, function(method) {
    warning_records <- result$method_diagnostics[[method]]$native_warning_records
    messages <- vapply(warning_records, function(item) item$message, character(1))
    sum(grepl("^(Lower|Upper) loop did not converge", messages, ignore.case = TRUE))
  }, integer(1))]
  metric[, optimizer_converged := converged & component_nonconvergence_warnings == 0L]
  metrics[[i]] <- metric
  coordinates[[i]] <- as.data.table(result$coordinates)
  theory[[i]] <- as.data.table(result$theory)
}
metrics <- rbindlist(metrics, fill = TRUE)
coordinates <- rbindlist(coordinates, fill = TRUE)
theory <- rbindlist(theory, fill = TRUE)
keys <- c("experiment_id", "replicate", "method_id")
if (anyDuplicated(metrics[, ..keys]) || nrow(metrics) != nrow(expected) ||
    nrow(fsetdiff(expected[, ..keys], metrics[, ..keys])) || nrow(fsetdiff(metrics[, ..keys], expected[, ..keys])))
  stop("Simulation method accounting is incomplete or duplicated.")
coordinate_keys <- c(keys, "coordinate")
if (anyDuplicated(coordinates[, ..coordinate_keys])) stop("Duplicated simulation coordinate records.")
setorderv(metrics, keys)
setorderv(coordinates, coordinate_keys)
setorderv(theory, keys)
output <- file.path("results/tables/submission_20261028/collections", basename(root))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(metrics = metrics, coordinates = coordinates, theory = theory,
  registered = expected, registry = registry, source = root,
  implementation_commit = unique(metrics$implementation_commit)), file.path(output, "records.rds"))
accounting <- metrics[, .(registered = .N, converged_count = sum(optimizer_converged),
  numerical_failures = sum(!optimizer_converged & failure_stage != "qif_covariance_rank"),
  component_nonconvergence = sum(component_nonconvergence_warnings > 0L),
  domain_unavailable = sum(failure_stage == "qif_covariance_rank"), warnings = sum(status == "warning"),
  unavailable = sum(status %in% c("failed", "not_implemented"))), by = .(experiment_id, method_id)]
fwrite(accounting, file.path(output, "accounting.csv"))
jsonlite::write_json(list(source = root, datasets = nrow(jobs), registered_method_records = nrow(expected),
  observed_method_records = nrow(metrics), coordinate_records = nrow(coordinates),
  implementation_commit = unique(metrics$implementation_commit), completed = as.character(Sys.time())),
  file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Collected", nrow(jobs), "datasets", nrow(metrics), "method records and", nrow(coordinates), "coordinates\n")

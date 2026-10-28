source("research/submission_completion_20261027/load_comparators.R")
source("research/submission_completion_20261027/dispatch_state.R")
source("research/submission_completion_20261027/diagnostic_storage.R")
args <- commandArgs(trailingOnly = TRUE)
destination <- args[1L]
batch <- as.integer(args[2L])
request <- readRDS(file.path(destination, "request.rds"))
jobs <- request$jobs[request$jobs$batch == batch, ]
if (!nrow(jobs)) stop("Empty comparator batch.")
current <- NULL
globalCallingHandlers(error = function(condition) {
  if (!is.null(current)) saveRDS(list(condition = condition, finished = Sys.time()), file.path(current, "runner_failure.rds"))
})
completed <- 0L
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  current <- file.path(destination, "replications", job$experiment_id, sprintf("replicate_%05d", job$replicate))
  if (file.exists(file.path(current, "result.rds"))) next
  if (file.exists(file.path(current, "runner_failure.rds"))) next
  if (dir.exists(current)) stop("An unfinished comparator replication requires explicit recovery.")
  dir.create(current, recursive = TRUE)
  config <- config_row_to_list_v2(request$registry[request$registry$experiment_id == job$experiment_id, ])
  config$methods <- intersect(config$methods, request$methods)
  if (!length(config$methods)) stop("The batch selected no registered comparator.")
  submission_atomic_rds(list(config = config, replicate = job$replicate, base_seed = request$base_seed,
    numerical_specification = request$numerical_specification, started = Sys.time()), file.path(current, "request.rds"))
  result <- run_one_replication_v2(".", config, job$replicate,
    base_seed = request$base_seed, allow_missing_benchmarks = FALSE, final = TRUE,
    run_id = basename(destination), implementation_commit = request$implementation_commit)
  validate_paired_failed_coordinate_accounting_v2(result$metrics, result$coordinates)
  if (!setequal(result$metrics$method_id, config$methods)) stop("Missing registered comparator record.")
  result <- submission_compact_diagnostics(result)
  submission_atomic_rds(result, file.path(current, "result.rds"))
  completed <- completed + 1L
  write.csv(data.frame(experiment_id = job$experiment_id, replicate = job$replicate, completed = completed,
    time = as.character(Sys.time())), file.path(destination, sprintf("batch_%03d_progress.csv", batch)), row.names = FALSE)
}
saveRDS(list(completed = completed, finished = Sys.time()), file.path(destination, sprintf("batch_%03d_completed.rds", batch)))
cat("Completed comparator batch", batch, "replications", completed, "\n")

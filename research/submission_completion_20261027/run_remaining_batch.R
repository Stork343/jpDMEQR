source("research/submission_completion_20261027/load_comparators.R")
source("research/submission_completion_20261027/dispatch_state.R")
source("research/submission_completion_20261027/diagnostic_storage.R")
args <- commandArgs(trailingOnly = TRUE)
root <- args[1L]
batch <- as.integer(args[2L])
request <- readRDS(file.path(root, "request.rds"))
jobs <- request$jobs[request$jobs$batch == batch, ]
if (!nrow(jobs)) stop("Empty remaining simulation batch.")
commit <- current_commit_v2()
paths <- file.path(root, "replications", jobs$experiment_id, sprintf("replicate_%05d", jobs$replicate))
completed <- sum(file.exists(file.path(paths, "result.rds")))
current <- NULL
globalCallingHandlers(error = function(condition) {
  if (!is.null(current)) saveRDS(list(condition = condition, implementation_commit = commit, finished = Sys.time()),
    file.path(current, paste0("resume_failure_", Sys.getpid(), ".rds")))
})
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  current <- paths[i]
  if (file.exists(file.path(current, "result.rds"))) next
  if (file.exists(file.path(current, "runner_failure.rds"))) stop("A recorded runner error requires investigation: ", current)
  config <- config_row_to_list_v2(request$registry[request$registry$experiment_id == job$experiment_id, ])
  if (identical(request$family, "joint")) config$methods <- config$methods[startsWith(config$methods, "JOINT-PROFILE")]
  if (!is.null(request$methods)) config$methods <- intersect(config$methods, request$methods)
  if (!length(config$methods)) stop("The remaining batch has no registered method.")
  if (dir.exists(current)) {
    prior <- readRDS(file.path(current, "request.rds"))
    if (!isTRUE(all.equal(prior$config, config)) || prior$replicate != job$replicate || prior$base_seed != request$base_seed)
      stop("An interrupted replication has a different statistical definition.")
  } else {
    dir.create(current, recursive = TRUE)
    submission_atomic_rds(list(config = config, replicate = job$replicate, base_seed = request$base_seed,
      implementation_commit = commit, started = Sys.time()), file.path(current, "request.rds"))
  }
  saveRDS(list(config = config, replicate = job$replicate, base_seed = request$base_seed,
    implementation_commit = commit, numerical_specification = yaml::read_yaml("config/submission_comparator_numerics_20261028.yml"),
    started = Sys.time()), file.path(current, paste0("completion_request_", Sys.getpid(), ".rds")))
  result <- run_one_replication_v2(".", config, job$replicate, base_seed = request$base_seed,
    allow_missing_benchmarks = FALSE, final = TRUE, run_id = basename(root), implementation_commit = commit)
  validate_paired_failed_coordinate_accounting_v2(result$metrics, result$coordinates)
  if (!setequal(result$metrics$method_id, config$methods)) stop("A remaining simulation method is missing.")
  result <- submission_compact_diagnostics(result)
  temporary <- file.path(current, paste0("result_", Sys.getpid(), ".writing"))
  saveRDS(result, temporary, compress = "xz")
  if (!file.rename(temporary, file.path(current, "result.rds"))) stop("Unable to publish the complete simulation record.")
  completed <- completed + 1L
  submission_atomic_csv(data.frame(experiment_id = job$experiment_id, replicate = job$replicate,
    completed = completed, time = as.character(Sys.time())), file.path(root, sprintf("batch_%03d_progress.csv", batch)))
}
if (completed != nrow(jobs)) stop("The remaining batch accounting differs.")
saveRDS(list(completed = completed, implementation_commit = commit, finished = Sys.time()),
  file.path(root, sprintf("batch_%03d_completed.rds", batch)))
cat("Completed registered simulation batch", batch, completed, "\n")

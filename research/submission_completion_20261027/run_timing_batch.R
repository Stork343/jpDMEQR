source("research/submission_completion_20261027/load_comparators.R")
args <- commandArgs(trailingOnly = TRUE)
root <- args[1L]
experiment <- args[2L]
worker <- as.integer(args[3L])
request <- readRDS(file.path(root, "request.rds"))
jobs <- request$jobs[request$jobs$experiment_id == experiment & request$jobs$worker == worker, ]
if (!nrow(jobs)) stop("Empty timing batch.")
started <- proc.time()[[3L]]
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  configuration <- config_row_to_list_v2(request$registry[request$registry$experiment_id == experiment, ])
  fitting <- configuration
  if (experiment %in% c("G06", "G07")) fitting$experiment_id <- "G_parallel"
  destination <- file.path(root, "replications", experiment, sprintf("replicate_%05d", job$replicate))
  if (dir.exists(destination)) stop("The timing replication already exists.")
  dir.create(destination, recursive = TRUE)
  saveRDS(list(config = configuration, fitting_config = fitting, replicate = job$replicate,
    base_seed = request$base_seed, started = Sys.time()), file.path(destination, "request.rds"))
  result <- run_one_replication_v2(".", fitting, job$replicate, base_seed = request$base_seed,
    allow_missing_benchmarks = FALSE, final = TRUE, run_id = basename(root),
    implementation_commit = request$implementation_commit)
  result$config <- configuration
  result$computational_seed_id <- fitting$experiment_id
  for (field in c("metrics", "coordinates", "theory")) result[[field]]$experiment_id <- experiment
  validate_paired_failed_coordinate_accounting_v2(result$metrics, result$coordinates)
  if (!setequal(result$metrics$method_id, configuration$methods)) stop("A timing method record is missing.")
  saveRDS(result, file.path(destination, "result.rds"))
  write.csv(data.frame(completed = i, registered = nrow(jobs), elapsed = proc.time()[[3L]] - started,
    time = as.character(Sys.time())), file.path(root, paste0(experiment, "_worker_", worker, "_progress.csv")), row.names = FALSE)
}
cat("Completed timing batch", experiment, worker, nrow(jobs), "\n")

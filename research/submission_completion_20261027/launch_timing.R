source("research/submission_completion_20261027/load_comparators.R")
library(processx)
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[1L] else "final"
if (!phase %in% c("probe", "final")) stop("Unknown timing phase.")
registry <- submission_registry("main")
registry <- registry[registry$module == "G_computation", ]
if (phase == "probe") {
  registry <- registry[registry$experiment_id == "G01", ]
  registry$replications <- 1L
}
jobs <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) data.frame(
  experiment_id = registry$experiment_id[i], replicate = seq_len(registry$replications[i]),
  worker = (seq_len(registry$replications[i]) - 1L) %% registry$workers[i] + 1L)))
root <- if (phase == "probe") "results/scratch/submission_completion_20261027/timing_probe" else
  "results/raw/simulation/submission_20261028_computation"
if (dir.exists(root)) stop("The timing dispatch already exists.")
dir.create(root, recursive = TRUE)
saveRDS(list(registry = registry, jobs = jobs, base_seed = 20261027L,
  implementation_commit = current_commit_v2(), numerical_specification = submission_numerical_specification,
  session = sessionInfo(), system = Sys.info(),
  shared_parallel_fitting_id = "G_parallel", started = Sys.time()), file.path(root, "request.rds"))
write.csv(jobs, file.path(root, "scheduled_replications.csv"), row.names = FALSE)
write.csv(registry, file.path(root, "config_snapshot.csv"), row.names = FALSE)
timing <- list()
for (i in seq_len(nrow(registry))) {
  id <- registry$experiment_id[i]
  clock_started <- proc.time()[[3L]]
  processes <- lapply(seq_len(registry$workers[i]), function(worker) processx::process$new(
    file.path(R.home("bin"), "Rscript"), c("research/submission_completion_20261027/run_timing_batch.R", root, id, as.character(worker)),
    stdout = file.path(root, paste0(id, "_worker_", worker, ".log")), stderr = "2>&1"))
  repeat {
    alive <- vapply(processes, function(process) process$is_alive(), logical(1))
    if (!any(alive)) break
    Sys.sleep(.5)
  }
  codes <- vapply(processes, function(process) process$get_exit_status(), integer(1))
  if (any(codes != 0L)) stop("A timing worker did not complete: ", id)
  timing[[i]] <- data.frame(experiment_id = id, workers = registry$workers[i],
    datasets = registry$replications[i], wall_seconds = proc.time()[[3L]] - clock_started,
    definition = "complete_batch_including_worker_startup_and_data_generation")
  write.csv(do.call(rbind, timing), file.path(root, "batch_wall_times.csv"), row.names = FALSE)
}
paths <- list.files(file.path(root, "replications"), "^result[.]rds$", recursive = TRUE)
if (length(paths) != nrow(jobs)) stop("The complete timing dataset accounting differs.")
saveRDS(list(completed_replications = length(paths), finished = Sys.time()), file.path(root, "completed.rds"))
cat("Completed", length(paths), "timing datasets\n")

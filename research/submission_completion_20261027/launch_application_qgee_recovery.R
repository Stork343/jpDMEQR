library(processx)
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[1L] else "final"
if (!phase %in% c("probe", "final")) stop("Unknown application recovery phase.")
source <- "results/raw/application_comparison/submission_20261028_qgee"
if (phase == "final") while (!file.exists(file.path(source, "completed.rds"))) Sys.sleep(15)
original <- readRDS(file.path(source, "request.rds"))
status <- read.csv(file.path(source, "status.csv"))
jobs <- original$jobs[, c("method", "model", "tau", "processing", "repetition", "fold", "job_id")]
if (phase == "probe") jobs <- jobs[jobs$job_id == 267L, ] else {
  unfinished_record <- vapply(jobs$job_id, function(id) {
    path <- file.path(source, sprintf("job_%04d", id))
    !file.exists(file.path(path, "result.rds")) && !file.exists(file.path(path, "failure.rds"))
  }, logical(1))
  jobs <- jobs[unfinished_record | status$exit_status[match(jobs$job_id, status$job_id)] < 0L, ]
}
root <- if (phase == "probe") "results/scratch/submission_completion_20261027/qgee_recovery_probe" else
  "results/raw/application_comparison/submission_20261028_qgee_recovery"
if (dir.exists(root)) stop("The application recovery dispatch already exists.")
dir.create(root, recursive = TRUE)
jobs$state <- rep("pending", nrow(jobs))
jobs$exit_status <- rep(NA_integer_, nrow(jobs))
saveRDS(list(jobs = jobs, workers = 2L, started = Sys.time(), source_attempt = source,
  source_root = original$source_root, numerical_backend = "cluster_batched_products_LAPACK_SVD",
  penalty_failure_policy = "source_HBIC_individual_penalties",
  implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)), file.path(root, "request.rds"))
active <- list()
repeat {
  for (key in names(active)) {
    process <- active[[key]]
    if (process$is_alive()) next
    i <- as.integer(key)
    jobs$state[i] <- "completed"
    jobs$exit_status[i] <- process$get_exit_status()
    active[[key]] <- NULL
    write.csv(jobs, file.path(root, "status.csv"), row.names = FALSE)
  }
  if (all(jobs$state == "completed")) break
  for (i in head(which(jobs$state == "pending"), max(0L, 2L - length(active)))) {
    active[[as.character(i)]] <- processx::process$new(file.path(R.home("bin"), "Rscript"),
      c("research/submission_completion_20261027/run_application_qgee.R", root, as.character(jobs$job_id[i])),
      stdout = file.path(root, sprintf("job_%04d.log", jobs$job_id[i])), stderr = "2>&1")
    jobs$state[i] <- "running"
  }
  Sys.sleep(.5)
}
write.csv(jobs, file.path(root, "status.csv"), row.names = FALSE)
for (id in jobs$job_id) {
  path <- file.path(root, sprintf("job_%04d", id))
  if (!file.exists(file.path(path, "result.rds")) && !file.exists(file.path(path, "failure.rds")))
    stop("An application recovery record is incomplete.")
}
saveRDS(list(completed_jobs = nrow(jobs), successful = sum(jobs$exit_status == 0L),
  failed = sum(jobs$exit_status != 0L), finished = Sys.time()), file.path(root, "completed.rds"))
cat("Completed application recovery records", nrow(jobs), "\n")

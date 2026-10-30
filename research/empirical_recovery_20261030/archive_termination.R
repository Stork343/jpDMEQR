library(data.table)
library(filelock)
source("research/submission_completion_20261027/dispatch_state.R")
source("research/empirical_recovery_20261030/queue.R")
record <- read.csv(file.path(empirical_queue_root, "active_dispatcher.csv"))
if (any(record$pid %in% ps::ps_pids())) stop("The empirical dispatcher is still running.")
status_path <- file.path(empirical_queue_root, "status.csv")
status <- fread(status_path)
if (any(status$pid[status$state == "running"] %in% ps::ps_pids()))
  stop("An empirical worker is still running.")
snapshot <- file.path(empirical_queue_root, "status_at_termination.csv")
if (file.exists(snapshot)) stop("The termination snapshot already exists.")
if (!file.copy(status_path, snapshot)) stop("Unable to preserve the dispatch snapshot.")
jobs <- empirical_registry()$jobs
for (i in which(status$state == "running")) {
  status$state[i] <- if (empirical_job_complete(jobs[i, ])) "completed" else "terminated_by_user"
  status$exit_status[i] <- if (status$state[i] == "completed") 0L else NA_integer_
}
status[state == "pending", state := "cancelled_by_user"]
submission_atomic_csv(status, status_path)
counts <- status[, .N, by = .(phase, state)]
record <- list(state = "terminated_by_user", recorded = as.character(Sys.time()),
  completed_outputs = "preserved", remaining_processes = 0L, dispatch_counts = counts,
  implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
submission_atomic_rds(record, file.path(empirical_queue_root, "termination.rds"))
jsonlite::write_json(record, "results/preflight/empirical_recovery_20261030/termination_manifest.json",
  auto_unbox = TRUE, pretty = TRUE)
empirical_backup(jobs)
destination <- file.path(path.expand("~"), "jpDMEQR-backups", "empirical_recovery_20261030", empirical_queue_root)
if (!all(file.copy(file.path(empirical_queue_root, c("termination.rds", "status_at_termination.csv")),
    destination, overwrite = TRUE))) stop("Unable to back up the termination records.")
print(counts)

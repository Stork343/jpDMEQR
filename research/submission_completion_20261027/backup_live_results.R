library(filelock)
source("research/submission_completion_20261027/dispatch_state.R")
base <- Sys.getenv("JPDMEQR_BACKUP_DIR", unset = file.path(path.expand("~"), "jpDMEQR-backups"))
dir.create(base, recursive = TRUE, showWarnings = FALSE)
lock <- filelock::lock(file.path(base, "live_backup.lock"), timeout = 0)
if (is.null(lock)) stop("The project output backup is already running.")
queues <- file.path("results/raw/simulation", c("submission_20261027_main_joint",
  paste0("submission_20261028_main_", c("mixed", "pooled", "qgee", "bias_adjusted", "sqr"))))
snapshots <- list()
repeat {
  active <- queues[dir.exists(queues)]
  for (root in active) {
    progress <- list.files(root, "^batch_[0-9]+_progress[.]csv$", full.names = TRUE)
    markers <- c(progress, file.path(root, c("request.rds", "completed.rds")))
    markers <- markers[file.exists(markers)]
    state <- setNames(as.numeric(file.info(markers)$mtime), markers)
    if (!identical(snapshots[[root]], state)) {
      if (isTRUE(submission_backup_run(root))) snapshots[[root]] <- state
    }
  }
  completed <- file.exists(file.path(queues, "completed.rds"))
  submission_atomic_csv(data.frame(queue = basename(queues), completed = completed,
    inspected = as.character(Sys.time())), file.path(base, "live_backup_status.csv"))
  if (all(completed)) break
  Sys.sleep(60)
}
for (root in queues) submission_backup_run(root, force = TRUE)
filelock::unlock(lock)

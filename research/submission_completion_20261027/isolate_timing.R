library(ps)
library(filelock)
library(processx)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || anyNA(as.integer(args))) stop("Provide the explicitly selected project controller process IDs.")
controllers <- lapply(as.integer(args), ps::ps_handle)
root <- normalizePath(".")
for (handle in controllers) {
  if (!identical(ps::ps_cwd(handle), root) ||
      !any(grepl("research/submission_completion_20261027/launch_", ps::ps_cmdline(handle), fixed = TRUE)))
    stop("An isolation target is not a simulation or application controller in this project.")
}
run <- function() {
  suspended <- list()
  locks <- list()
  output <- "results/scratch/submission_completion_20261027/timing_isolation"
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  resume <- function() {
    for (lock in locks) filelock::unlock(lock)
    for (handle in rev(suspended)) if (ps::ps_is_running(handle)) ps::ps_resume(handle)
    saveRDS(list(resumed = Sys.time()), file.path(output, "resumed.rds"))
  }
  on.exit(resume(), add = TRUE)
  suspend <- function(handle) {
    if (!ps::ps_is_running(handle)) return(invisible(NULL))
    ps::ps_suspend(handle)
    suspended[[length(suspended) + 1L]] <<- handle
  }
  for (handle in controllers) suspend(handle)
  lock_paths <- c("results/scratch/score_validity_20261026/compile.lock",
    "results/scratch/submission_completion_20261027/compile.lock",
    "results/scratch/submission_completion_20261027/qgee_compile.lock")
  for (path in lock_paths) {
    lock <- filelock::lock(path, timeout = 60000)
    if (is.null(lock)) stop("A compilation lock is unavailable for timing isolation.")
    locks[[length(locks) + 1L]] <- lock
  }
  suspend_children <- function(handle) {
    for (child in ps::ps_children(handle)) {
      suspend(child)
      suspend_children(child)
    }
  }
  for (handle in controllers) suspend_children(handle)
  for (lock in locks) filelock::unlock(lock)
  locks <- list()
  records <- lapply(suspended, function(handle) data.frame(pid = ps::ps_pid(handle),
    created = as.character(ps::ps_create_time(handle)), command = paste(ps::ps_cmdline(handle), collapse = " ")))
  write.csv(do.call(rbind, records), file.path(output, "suspended.csv"), row.names = FALSE)
  process <- processx::run(file.path(R.home("bin"), "Rscript"),
    c("research/submission_completion_20261027/launch_timing.R", "final"),
    error_on_status = FALSE, stdout = file.path(output, "timing.log"), stderr = "2>&1")
  if (process$status != 0L) stop("The isolated timing study did not complete.")
}
run()

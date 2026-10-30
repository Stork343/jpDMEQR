library(processx)
source("research/empirical_recovery_20261030/execution_control.R")
empirical_require_active()
source("research/submission_completion_20261027/dispatch_state.R")
Sys.setenv(TMPDIR = normalizePath("results/scratch/runtime"))
root <- "results/scratch/empirical_recovery_20261030"
dir.create(root, recursive = TRUE, showWarnings = FALSE)
records <- list()
previous_path <- file.path(root, "active_downloads.csv")
previous <- if (file.exists(previous_path)) read.csv(previous_path) else NULL
for (release in c("R1", "R2")) {
  if (file.exists(file.path("data/raw/GSE65391", paste0("download_recovery_20261030_", release, ".csv")))) {
    records[[release]] <- data.frame(release = release, pid = NA_integer_, started = as.character(Sys.time()))
    next
  }
  old <- if (is.null(previous)) NULL else previous[previous$release == release, ]
  if (!is.null(old) && nrow(old) == 1L && !is.na(old$pid) && old$pid %in% ps::ps_pids()) {
    handle <- ps::ps_handle(old$pid)
    if (!"--file=research/empirical_recovery_20261030/download_inputs.R" %in% ps::ps_cmdline(handle))
      stop("A recorded acquisition process has a different identity.")
    records[[release]] <- old
    next
  }
  process <- processx::process$new(file.path(R.home("bin"), "Rscript"),
    c("research/empirical_recovery_20261030/download_inputs.R", release), cleanup = FALSE,
    stdout = file.path(root, paste0("download_", release, ".log")), stderr = "2>&1")
  records[[release]] <- data.frame(release = release, pid = process$get_pid(), started = as.character(Sys.time()))
}
submission_atomic_csv(do.call(rbind, records), previous_path)
print(do.call(rbind, records))

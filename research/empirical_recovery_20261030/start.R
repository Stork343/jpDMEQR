library(processx)
Sys.setenv(TMPDIR = normalizePath("results/scratch/runtime"),
  VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
root <- "results/raw/application/empirical_recovery_20261030"
dir.create(root, recursive = TRUE, showWarnings = FALSE)
process <- processx::process$new(file.path(R.home("bin"), "Rscript"),
  "research/empirical_recovery_20261030/dispatch.R", cleanup = FALSE,
  stdout = file.path(root, "dispatcher.log"), stderr = "2>&1")
pid <- process$get_pid()
processx::process$new("/usr/bin/caffeinate", c("-s", "-w", as.character(pid)), cleanup = FALSE)
write.csv(data.frame(pid = pid, started = as.character(Sys.time()), workers = 2L),
  file.path(root, "active_dispatcher.csv"), row.names = FALSE)
cat("Started empirical dispatcher", pid, "with two workers\n")

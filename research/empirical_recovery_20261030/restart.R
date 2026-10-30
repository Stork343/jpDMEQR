library(ps)
root <- "results/raw/application/empirical_recovery_20261030"
record <- read.csv(file.path(root, "active_dispatcher.csv"))
if (nrow(record) != 1L || !ps::ps_pid_exists(record$pid)) stop("The recorded empirical dispatcher is unavailable.")
handle <- ps::ps_handle(record$pid)
if (!"--file=research/empirical_recovery_20261030/dispatch.R" %in% ps::ps_cmdline(handle) ||
    !identical(ps::ps_cwd(handle), normalizePath("."))) stop("The dispatcher process identity differs.")
ps::ps_terminate(handle)
deadline <- Sys.time() + 10
while (ps::ps_is_running(handle)) {
  if (Sys.time() > deadline) stop("The dispatcher did not terminate.")
  Sys.sleep(.1)
}
status <- read.csv(file.path(root, "status.csv"), stringsAsFactors = FALSE)
retry <- status$job_id[status$state == "failed" & status$phase == "processing_modules"]
code <- system2(file.path(R.home("bin"), "Rscript"),
  c("research/empirical_recovery_20261030/start.R", as.character(retry)))
if (code != 0L) stop("The empirical dispatcher did not restart.")

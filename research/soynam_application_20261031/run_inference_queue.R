library(data.table)
library(yaml)
library(processx)
library(filelock)
config <- read_yaml("config/application_soynam_20261031.yml")
root <- file.path(config$results_directory, "inference", config$execution_version)
dir.create(root, recursive = TRUE, showWarnings = FALSE)
handle <- filelock::lock(file.path(root, "queue.lock"), timeout = 0)
if (is.null(handle)) stop("The inference queue is already active.")
jobs <- fread(file.path(config$results_directory, "registry.csv"))
jobs[, `:=`(state = "pending", pid = NA_integer_, exit_status = NA_integer_)]
for (i in seq_len(nrow(jobs))) {
  directory <- file.path(root, jobs$case_id[i])
  if (file.exists(file.path(directory, "completion.json"))) {
    jobs$state[i] <- "complete"
    jobs$exit_status[i] <- 0L
  } else if (dir.exists(directory)) stop("An existing incomplete inference requires inspection: ", jobs$case_id[i])
}
active <- list()
repeat {
  for (key in names(active)) {
    process <- active[[key]]
    if (process$is_alive()) next
    i <- as.integer(key)
    jobs$exit_status[i] <- process$get_exit_status()
    complete <- file.exists(file.path(root, jobs$case_id[i], "completion.json"))
    jobs$state[i] <- if (jobs$exit_status[i] == 0L && complete) "complete" else "failed"
    active[[key]] <- NULL
    fwrite(jobs, file.path(root, "status.csv"))
    cat(format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), jobs$case_id[i], jobs$state[i], "\n")
    flush.console()
    if (jobs$state[i] == "failed") stop("An inference case failed: ", jobs$case_id[i])
  }
  for (i in which(jobs$state == "pending")) {
    if (length(active) >= 2L) break
    process <- processx::process$new(file.path(R.home("bin"), "Rscript"),
      c("research/soynam_application_20261031/run_case.R", jobs$case_id[i]), cleanup = FALSE,
      stdout = file.path(root, paste0(jobs$case_id[i], ".log")), stderr = "2>&1")
    active[[as.character(i)]] <- process
    jobs$state[i] <- "running"
    jobs$pid[i] <- process$get_pid()
    cat(format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), jobs$case_id[i], "started\n")
    flush.console()
  }
  fwrite(jobs, file.path(root, "status.csv"))
  if (all(jobs$state == "complete")) break
  Sys.sleep(2)
}
filelock::unlock(handle)

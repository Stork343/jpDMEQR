library(data.table)
library(yaml)
library(processx)
library(filelock)
config <- read_yaml("config/application_milk_20261031.yml")
run_directory <- file.path(config$results_directory, config$execution_version)
dir.create(run_directory, recursive = TRUE, showWarnings = FALSE)
handle <- filelock::lock(file.path(run_directory, "dispatch.lock"), timeout = 0)
if (is.null(handle)) stop("The empirical queue is already active.")
registry <- fread(file.path(config$results_directory, "registry.csv"))
registry[, `:=`(state = "pending", pid = NA_integer_, exit_status = NA_integer_)]
for (i in seq_len(nrow(registry))) {
  directory <- file.path(run_directory, registry$case_id[i])
  if (file.exists(file.path(directory, "completion.json"))) {
    registry$state[i] <- "complete"
    registry$exit_status[i] <- 0L
  } else if (dir.exists(directory)) {
    stop("An existing case requires inspection before queue continuation: ", registry$case_id[i])
  }
}
active <- list()
repeat {
  for (key in names(active)) {
    process <- active[[key]]
    if (process$is_alive()) next
    i <- as.integer(key)
    registry$exit_status[i] <- process$get_exit_status()
    completion <- file.path(run_directory, registry$case_id[i], "completion.json")
    registry$state[i] <- if (registry$exit_status[i] == 0L && file.exists(completion)) "complete" else "failed"
    active[[key]] <- NULL
    fwrite(registry, file.path(run_directory, "status.csv"))
    cat(format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), registry$case_id[i], registry$state[i], "\n")
    flush.console()
    if (registry$state[i] == "failed") stop("The empirical case failed: ", registry$case_id[i])
  }
  for (i in which(registry$state == "pending")) {
    if (length(active) >= 2L) break
    process <- processx::process$new(file.path(R.home("bin"), "Rscript"),
      c("research/milk_application_20261031/run_case.R", registry$case_id[i]),
      stdout = file.path(run_directory, paste0(registry$case_id[i], ".log")),
      stderr = "2>&1", cleanup = FALSE)
    active[[as.character(i)]] <- process
    registry$state[i] <- "running"
    registry$pid[i] <- process$get_pid()
    cat(format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), registry$case_id[i], "started\n")
    flush.console()
  }
  fwrite(registry, file.path(run_directory, "status.csv"))
  if (all(registry$state == "complete")) break
  Sys.sleep(2)
}
filelock::unlock(handle)
cat("All registered empirical cases completed.\n")

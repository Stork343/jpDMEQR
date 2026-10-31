library(data.table)
library(yaml)
library(processx)
library(ps)
library(filelock)
config <- read_yaml("config/application_soynam_20261031.yml")
root <- file.path(config$results_directory, "prediction", config$execution_version)
dir.create(root, recursive = TRUE, showWarnings = FALSE)
handle <- filelock::lock(file.path(root, "queue.lock"), timeout = 0)
if (is.null(handle)) stop("The prediction queue is already active.")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) %% 2L) stop("Existing jobs require case and process identifier pairs.")
adopted <- if (length(arguments)) setNames(as.integer(arguments[seq(2L, length(arguments), 2L)]),
  arguments[seq(1L, length(arguments), 2L)]) else integer()
jobs <- CJ(method = unlist(config$prediction$methods), tau = unlist(config$quantiles), sorted = FALSE)
jobs[, `:=`(case_id = paste0(method, "_q", sprintf("%02d", round(100 * tau))),
  state = "pending", pid = NA_integer_, exit_status = NA_integer_)]
jobs[, priority := match(tau, c(.5, .25, .75))]
setorder(jobs, priority, method)
active <- list()
for (i in seq_len(nrow(jobs))) {
  directory <- file.path(root, jobs$case_id[i])
  if (file.exists(file.path(directory, "completion.json"))) {
    jobs$state[i] <- "complete"
    jobs$exit_status[i] <- 0L
  } else if (dir.exists(directory)) {
    if (!jobs$case_id[i] %in% names(adopted)) stop("An incomplete prediction requires inspection: ", jobs$case_id[i])
    process <- ps::ps_handle(adopted[[jobs$case_id[i]]])
    command <- ps::ps_cmdline(process)
    if (!ps::ps_is_running(process) || !jobs$method[i] %in% command ||
        !any(grepl("run_prediction.R", command, fixed = TRUE))) stop("The supplied process does not match its prediction case.")
    active[[as.character(i)]] <- list(handle = process, process = NULL)
    jobs$state[i] <- "running"
    jobs$pid[i] <- adopted[[jobs$case_id[i]]]
  }
}
repeat {
  for (key in names(active)) {
    current <- active[[key]]
    if (ps::ps_is_running(current$handle)) next
    i <- as.integer(key)
    complete <- file.exists(file.path(root, jobs$case_id[i], "completion.json"))
    status <- if (is.null(current$process)) as.integer(!complete) else current$process$get_exit_status()
    jobs$exit_status[i] <- status
    jobs$state[i] <- if (status == 0L && complete) "complete" else "failed"
    active[[key]] <- NULL
    fwrite(jobs, file.path(root, "status.csv"))
    cat(format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), jobs$case_id[i], jobs$state[i], "\n")
    flush.console()
    if (jobs$state[i] == "failed") stop("A prediction case failed: ", jobs$case_id[i])
  }
  for (i in which(jobs$state == "pending")) {
    if (length(active) >= 3L) break
    process <- processx::process$new(file.path(R.home("bin"), "Rscript"),
      c("-e", "source('research/soynam_application_20261031/run_prediction.R')", jobs$method[i], as.character(jobs$tau[i])),
      cleanup = FALSE, stdout = file.path(root, paste0(jobs$case_id[i], ".log")), stderr = "2>&1")
    active[[as.character(i)]] <- list(handle = ps::ps_handle(process$get_pid()), process = process)
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

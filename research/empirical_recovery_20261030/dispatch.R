library(processx)
library(filelock)
source("research/submission_completion_20261027/dispatch_state.R")
source("research/empirical_recovery_20261030/queue.R")
run <- function() {
  dir.create(empirical_queue_root, recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(file.path(empirical_queue_root, "dispatch.lock"), timeout = 0)
  if (is.null(lock)) stop("The empirical queue is already running.")
  on.exit(filelock::unlock(lock), add = TRUE)
  registry <- empirical_registry()
  request_path <- file.path(empirical_queue_root, "request.rds")
  if (file.exists(request_path)) {
    previous <- readRDS(request_path)
    if (!identical(previous$jobs, registry$jobs)) stop("The empirical job definitions changed.")
  } else submission_atomic_rds(list(jobs = registry$jobs, workers = 2L,
    implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    session = sessionInfo(), started = Sys.time()), request_path)
  comparator_path <- file.path(empirical_comparator_root, "request.rds")
  if (!file.exists(comparator_path)) submission_atomic_rds(list(jobs = registry$comparator,
    started = Sys.time(), source = request_path), comparator_path)
  jobs <- registry$jobs
  jobs$state <- "pending"
  jobs$exit_status <- NA_integer_
  jobs$pid <- NA_integer_
  status_path <- file.path(empirical_queue_root, "status.csv")
  if (file.exists(status_path)) {
    prior <- read.csv(status_path, stringsAsFactors = FALSE)
    if (!identical(prior$job_id, jobs$job_id)) stop("The empirical status identities changed.")
    for (i in which(prior$state == "running")) {
      if (!is.na(prior$pid[i]) && ps::ps_pid_exists(prior$pid[i]))
        stop("A recorded empirical worker is still running: ", prior$pid[i])
    }
    jobs$state[prior$state == "failed"] <- "failed"
    jobs$exit_status[prior$state == "failed"] <- prior$exit_status[prior$state == "failed"]
  }
  for (i in seq_len(nrow(jobs))) if (empirical_job_complete(jobs[i, ])) {
    jobs$state[i] <- "completed"
    jobs$exit_status[i] <- 0L
  }
  active <- list()
  last_backup <- Sys.time() - 61
  repeat {
    for (key in names(active)) {
      process <- active[[key]]
      if (process$is_alive()) next
      i <- as.integer(key)
      code <- process$get_exit_status()
      if (code == 0L && !empirical_job_complete(jobs[i, ]))
        stop("An empirical worker returned without its complete results: ", i)
      jobs$state[i] <- if (code == 0L) "completed" else "failed"
      jobs$exit_status[i] <- code
      active[[key]] <- NULL
    }
    for (i in which(jobs$state == "pending")) {
      if (length(active) >= 2L) break
      if (!empirical_job_ready(jobs[i, ])) next
      process <- processx::process$new(file.path(R.home("bin"), "Rscript"),
        c(jobs$script[i], jobs$arguments[[i]]), cleanup = FALSE,
        stdout = file.path(empirical_queue_root, sprintf("job_%04d.log", i)), stderr = "2>&1")
      active[[as.character(i)]] <- process
      jobs$state[i] <- "running"
      jobs$pid[i] <- process$get_pid()
    }
    submission_atomic_csv(jobs[, setdiff(names(jobs), "arguments")], status_path)
    if (as.numeric(difftime(Sys.time(), last_backup, units = "secs")) >= 60) {
      empirical_backup(jobs)
      last_backup <- Sys.time()
    }
    if (all(jobs$state %in% c("completed", "failed"))) break
    Sys.sleep(5)
  }
  empirical_backup(jobs)
  submission_atomic_rds(list(successful = sum(jobs$state == "completed"),
    failed = sum(jobs$state == "failed"), finished = Sys.time()),
    file.path(empirical_queue_root, "completed.rds"))
}
run()

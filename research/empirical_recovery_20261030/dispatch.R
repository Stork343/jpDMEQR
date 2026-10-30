library(processx)
library(filelock)
source("research/empirical_recovery_20261030/execution_control.R")
empirical_require_active()
source("research/submission_completion_20261027/dispatch_state.R")
source("research/empirical_recovery_20261030/queue.R")
run <- function() {
  retry_ids <- as.integer(commandArgs(trailingOnly = TRUE))
  if (anyNA(retry_ids)) stop("Supply numeric job identifiers for an explicit continuation.")
  dir.create(empirical_queue_root, recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(file.path(empirical_queue_root, "dispatch.lock"), timeout = 0)
  if (is.null(lock)) stop("The empirical queue is already running.")
  on.exit(filelock::unlock(lock), add = TRUE)
  registry <- empirical_registry()
  request_path <- file.path(empirical_queue_root, "request.rds")
  if (file.exists(request_path)) {
    previous <- readRDS(request_path)
    identity <- setdiff(names(registry$jobs), "arguments")
    if (!identical(previous$jobs[, identity], registry$jobs[, identity]))
      stop("The empirical job definitions changed.")
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
  active <- list()
  status_path <- file.path(empirical_queue_root, "status.csv")
  if (file.exists(status_path)) {
    prior <- read.csv(status_path, stringsAsFactors = FALSE)
    if (!identical(prior$job_id, jobs$job_id)) stop("The empirical status identities changed.")
    for (i in which(prior$state == "running")) {
      if (!is.na(prior$pid[i]) && prior$pid[i] %in% ps::ps_pids()) {
        handle <- ps::ps_handle(prior$pid[i])
        arguments <- ps::ps_cmdline(handle)
        position <- match("--args", arguments)
        if (!paste0("--file=", jobs$script[i]) %in% arguments ||
            !identical(ps::ps_cwd(handle), normalizePath(".")) ||
            (!is.na(position) && !identical(arguments[-seq_len(position)], jobs$arguments[[i]])))
          stop("A recorded empirical process has a different identity.")
        active[[as.character(i)]] <- list(process = NULL, handle = handle)
        jobs$state[i] <- "running"
        jobs$pid[i] <- prior$pid[i]
      }
    }
    jobs$state[prior$state == "failed"] <- "failed"
    jobs$exit_status[prior$state == "failed"] <- prior$exit_status[prior$state == "failed"]
    for (i in retry_ids) {
      if (!i %in% jobs$job_id || jobs$state[i] != "failed") stop("A requested continuation is not a failed job.")
      log <- file.path(empirical_queue_root, sprintf("job_%04d.log", i))
      retained <- file.path(empirical_queue_root, "attempt_logs", paste0(sprintf("job_%04d", i), "_", Sys.getpid(), ".log"))
      dir.create(dirname(retained), recursive = TRUE, showWarnings = FALSE)
      if (!file.rename(log, retained)) stop("Unable to preserve the completed process log.")
      jobs$state[i] <- "pending"
      jobs$exit_status[i] <- NA_integer_
    }
  }
  for (i in seq_len(nrow(jobs))) if (empirical_job_complete(jobs[i, ])) {
    jobs$state[i] <- "completed"
    jobs$exit_status[i] <- 0L
  }
  submission_atomic_rds(list(jobs = registry$jobs, adopted = names(active), retry_ids = retry_ids,
    implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), started = Sys.time()),
    file.path(empirical_queue_root, paste0("execution_request_", Sys.getpid(), ".rds")))
  last_backup <- Sys.time() - 61
  repeat {
    for (key in names(active)) {
      worker <- active[[key]]
      if (submission_worker_alive(worker)) next
      i <- as.integer(key)
      code <- if (!is.null(worker$process)) worker$process$get_exit_status() else
        if (empirical_job_complete(jobs[i, ])) 0L else 1L
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
      active[[as.character(i)]] <- list(process = process, handle = ps::ps_handle(process$get_pid()))
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
  builders <- list(
    c("research/empirical_recovery_20261030/audit_progress.R"),
    c("research/empirical_recovery_20261030/collect_comparators.R", "final"),
    c("research/submission_completion_20261027/summarise_application.R"),
    c("research/submission_completion_20261027/build_application_tables.R"))
  for (arguments in builders) {
    status <- system2(file.path(R.home("bin"), "Rscript"), arguments)
    if (status != 0L) stop("The empirical result collection did not complete: ", arguments[1L])
  }
}
run()

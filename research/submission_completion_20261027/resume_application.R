library(ps)
library(processx)
library(filelock)
source("research/submission_completion_20261027/dispatch_state.R")
source("research/submission_completion_20261027/application_state.R")
args <- commandArgs(trailingOnly = TRUE)
family <- args[1L]
workers <- as.integer(args[2L])
if (!family %in% c("mixed", "qgee") || !workers %in% 1:8) stop("Invalid application continuation.")
root <- if (family == "mixed") "results/raw/application_comparison/submission_20261027_v2" else
  "results/raw/application_comparison/submission_20261028_qgee"
script <- if (family == "mixed") "run_application_comparator.R" else "run_application_qgee.R"
run <- function() {
  lock <- filelock::lock(file.path(root, "continuation.lock"), timeout = 0)
  if (is.null(lock)) stop("This application continuation is already running.")
  on.exit(filelock::unlock(lock), add = TRUE)
  if (file.exists(file.path(root, "completed.rds"))) return(invisible(NULL))
  request <- readRDS(file.path(root, "request.rds"))
  jobs <- request$jobs
  if (!"job_id" %in% names(jobs)) jobs$job_id <- seq_len(nrow(jobs))
  fields <- c("method", "model", "tau", "processing", "repetition", "fold")
  jobs$state <- "pending"
  jobs$exit_status <- NA_integer_
  active <- submission_processes(root, script)
  inspect <- function(i) {
    path <- file.path(root, sprintf("job_%04d", jobs$job_id[i]))
    result_file <- file.path(path, "result.rds")
    failure_file <- file.path(path, "failure.rds")
    if (file.exists(result_file) && file.exists(failure_file)) stop("An application job has conflicting terminal records.")
    if (file.exists(result_file) || file.exists(failure_file)) {
      success <- file.exists(result_file)
      result <- readRDS(if (success) result_file else failure_file)
      if (!isTRUE(all.equal(result$job[, fields], jobs[i, fields], check.attributes = FALSE)))
        stop("An application continuation identity differs: ", path)
      if (success && (!isTRUE(result$fit$converged) || !nrow(result$predictions) ||
          any(!is.finite(result$predictions$predicted_quantile)))) stop("An application result is incomplete.")
      if (!success && (is.null(result$condition) || is.null(result$stage))) stop("An application failure is incomplete.")
      if (!file.exists(file.path(path, "completion.rds"))) saveRDS(list(job = jobs[i, fields],
        exit_status = if (success) 0L else 1L, checked = Sys.time()), file.path(path, "completion.rds"))
      return(if (success) 0L else 1L)
    }
    log <- file.path(root, sprintf("job_%04d.log", jobs$job_id[i]))
    if (file.exists(log) && any(grepl("caught segfault|irrecoverable exception", readLines(log, warn = FALSE)))) return(-11L)
    NA_integer_
  }
  for (i in seq_len(nrow(jobs))) {
    key <- as.character(jobs$job_id[i])
    if (!is.null(active[[key]])) jobs$state[i] <- "running" else {
      code <- inspect(i)
      if (!is.na(code)) {
        jobs$state[i] <- "completed"
        jobs$exit_status[i] <- code
      }
    }
  }
  saveRDS(list(started = Sys.time(), workers = workers, adopted_jobs = names(active),
    implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), session = sessionInfo()),
    file.path(root, paste0("continuation_", Sys.getpid(), ".rds")))
  written_status <- NULL
  repeat {
    for (key in names(active)) {
      if (submission_worker_alive(active[[key]])) next
      i <- match(as.integer(key), jobs$job_id)
      code <- inspect(i)
      if (is.na(code)) stop("An application worker ended without a terminal record: ", key)
      jobs$state[i] <- "completed"
      jobs$exit_status[i] <- code
      active[[key]] <- NULL
    }
    if (!identical(written_status, jobs)) {
      submission_atomic_csv(jobs, file.path(root, "status.csv"))
      written_status <- jobs
    }
    if (all(jobs$state == "completed")) break
    pending <- which(jobs$state == "pending")
    if (family == "qgee") pending <- pending[vapply(jobs$job_id[pending],
      function(id) submission_application_completed(request$source_root, id), logical(1))]
    for (i in head(pending, max(0L, workers - length(active)))) {
      id <- jobs$job_id[i]
      key <- as.character(id)
      directory <- file.path(root, sprintf("job_%04d", id))
      log <- file.path(root, sprintf("job_%04d.log", id))
      if (dir.exists(directory)) {
        prior <- readRDS(file.path(directory, "request.rds"))
        if (!isTRUE(all.equal(prior[, fields], jobs[i, fields], check.attributes = FALSE)))
          stop("The interrupted application request differs.")
        retained <- file.path(root, "incomplete_jobs", paste0(sprintf("job_%04d", id), "_", Sys.getpid()))
        dir.create(dirname(retained), recursive = TRUE, showWarnings = FALSE)
        if (dir.exists(retained) || !file.rename(directory, retained)) stop("Unable to preserve an interrupted application calculation.")
        if (file.exists(log) && !file.rename(log, file.path(retained, "process.log"))) stop("Unable to preserve its process log.")
      }
      process <- processx::process$new(file.path(R.home("bin"), "Rscript"),
        c(file.path("research/submission_completion_20261027", script), root, key), cleanup = FALSE,
        stdout = log, stderr = "2>&1")
      active[[key]] <- list(process = process, handle = ps::ps_handle(process$get_pid()))
      jobs$state[i] <- "running"
    }
    Sys.sleep(1)
  }
  saveRDS(list(completed_jobs = nrow(jobs), successful = sum(jobs$exit_status == 0L),
    failed = sum(jobs$exit_status != 0L), finished = Sys.time()), file.path(root, "completed.rds"))
  cat("Completed application queue", family, nrow(jobs), "\n")
}
run()

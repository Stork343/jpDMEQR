library(ps)
library(processx)
library(filelock)
source("research/submission_completion_20261027/dispatch_state.R")
args <- commandArgs(trailingOnly = TRUE)
root <- args[1L]
workers <- as.integer(args[2L])
if (!workers %in% 1:8 || !file.exists(file.path(root, "request.rds"))) stop("Invalid simulation continuation.")
run <- function() {
  if (length(submission_launchers_alive(root))) stop("A live dispatcher already owns this simulation queue.")
  lock <- filelock::lock(file.path(root, "continuation.lock"), timeout = 0)
  if (is.null(lock)) stop("This simulation continuation is already running.")
  on.exit(filelock::unlock(lock), add = TRUE)
  if (file.exists(file.path(root, "completed.rds"))) return(invisible(NULL))
  request <- readRDS(file.path(root, "request.rds"))
  batches <- data.frame(batch = sort(unique(request$jobs$batch)), state = "pending", exit_status = NA_integer_)
  active <- submission_processes(root, c("run_simulation_batch.R", "run_comparator_batch.R", "run_remaining_batch.R"))
  complete <- function(batch) {
    jobs <- request$jobs[request$jobs$batch == batch, ]
    paths <- file.path(root, "replications", jobs$experiment_id, sprintf("replicate_%05d", jobs$replicate), "result.rds")
    all(file.exists(paths))
  }
  for (i in seq_len(nrow(batches))) {
    key <- as.character(batches$batch[i])
    if (!is.null(active[[key]])) batches$state[i] <- "running" else if (complete(batches$batch[i])) {
      batches$state[i] <- "completed"
      batches$exit_status[i] <- 0L
    }
  }
  saveRDS(list(started = Sys.time(), workers = workers, adopted_batches = names(active),
    implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), session = sessionInfo()),
    file.path(root, paste0("continuation_", Sys.getpid(), ".rds")))
  written_status <- NULL
  repeat {
    for (key in names(active)) {
      if (submission_worker_alive(active[[key]])) next
      i <- match(as.integer(key), batches$batch)
      if (!complete(batches$batch[i])) stop("A simulation worker ended before completing its registered batch: ", key)
      batches$state[i] <- "completed"
      batches$exit_status[i] <- 0L
      active[[key]] <- NULL
    }
    if (!identical(written_status, batches)) {
      submission_atomic_csv(batches, file.path(root, "batch_status.csv"))
      written_status <- batches
    }
    if (all(batches$state == "completed")) break
    capacity <- if (identical(request$family, "joint")) workers else submission_comparator_capacity(request$family, workers)
    for (i in head(which(batches$state == "pending"), max(0L, capacity - length(active)))) {
      key <- as.character(batches$batch[i])
      process <- processx::process$new(file.path(R.home("bin"), "Rscript"),
        c("research/submission_completion_20261027/run_remaining_batch.R", root, key), cleanup = FALSE,
        stdout = file.path(root, sprintf("batch_%03d_continuation_%d.log", batches$batch[i], Sys.getpid())), stderr = "2>&1")
      active[[key]] <- list(process = process, handle = ps::ps_handle(process$get_pid()))
      batches$state[i] <- "running"
    }
    Sys.sleep(2)
  }
  results <- list.files(file.path(root, "replications"), "^result[.]rds$", recursive = TRUE)
  if (length(results) != nrow(request$jobs)) stop("The complete simulation record count differs.")
  saveRDS(list(completed_replications = length(results), finished = Sys.time()), file.path(root, "completed.rds"))
  cat("Completed simulation queue", basename(root), length(results), "\n")
}
run()

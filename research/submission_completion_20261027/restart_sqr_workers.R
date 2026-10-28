library(ps)
root <- "results/raw/simulation/submission_20261028_main_sqr"
request <- readRDS(file.path(root, "request.rds"))
args <- as.integer(commandArgs(trailingOnly = TRUE))
if (!length(args) || anyNA(args)) stop("Provide the current SQR worker process IDs.")
handles <- lapply(args, ps::ps_handle)
kill_worker <- ps::ps_kill
batches <- integer(length(handles))
for (i in seq_along(handles)) {
  command <- ps::ps_cmdline(handles[[i]])
  marker <- match("--args", command)
  if (!"--file=research/submission_completion_20261027/run_comparator_batch.R" %in% command ||
      is.na(marker) || command[marker + 1L] != root || ps::ps_cwd(handles[[i]]) != normalizePath("."))
    stop("A selected process is not an original SQR worker.")
  batches[i] <- as.integer(command[marker + 2L])
}
transition <- function() {
  suspended <- list()
  on.exit(for (handle in suspended) if (isTRUE(ps::ps_is_running(handle))) ps::ps_resume(handle), add = TRUE)
  for (handle in handles) {
    ps::ps_suspend(handle)
    suspended[[length(suspended) + 1L]] <- handle
  }
  records <- list()
  for (i in seq_along(batches)) {
    jobs <- request$jobs[request$jobs$batch == batches[i], ]
    paths <- file.path(root, "replications", jobs$experiment_id, sprintf("replicate_%05d", jobs$replicate))
    existing <- which(dir.exists(paths))
    complete <- file.exists(file.path(paths, "result.rds"))
    incomplete <- existing[!complete[existing]]
    if (length(incomplete) != 1L) stop("The SQR worker is not within a single incomplete replication.")
    latest <- tail(which(complete), 1L)
    if (length(latest)) {
      value <- readRDS(file.path(paths[latest], "result.rds"))
      if (nrow(value$metrics) != 1L || value$metrics$method_id != "SQR-DEBIASED-IID" ||
          value$metrics$replicate != jobs$replicate[latest] || value$metrics$experiment_id != jobs$experiment_id[latest])
        stop("The most recent SQR record is incomplete.")
    }
    current <- readRDS(file.path(paths[incomplete], "request.rds"))
    if (current$replicate != jobs$replicate[incomplete] || current$config$experiment_id != jobs$experiment_id[incomplete])
      stop("The current SQR request does not match its registered data.")
    records[[i]] <- list(batch = batches[i], pid = args[i], incomplete_directory = paths[incomplete],
      completed_replications = sum(complete), request = current)
  }
  saveRDS(list(workers = records, storage_policy = "SQR_unrounded_end_conditions_and_candidate_coefficients", created = Sys.time()),
    file.path(root, paste0("storage_transition_", Sys.getpid(), ".rds")))
  for (handle in handles) kill_worker(handle)
  cat("Preserved completed SQR records and prepared", length(records), "registered replications for continuation\n")
}
transition()

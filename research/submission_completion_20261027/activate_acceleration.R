library(ps)
library(processx)
library(filelock)
source("research/submission_completion_20261027/dispatch_state.R")
accepted <- jsonlite::read_json("results/preflight/sqr_acceleration_20261030/acceptance.json")
if (!isTRUE(accepted$passed) || accepted$complete_recorded_datasets != 4L ||
    accepted$largest_verified_dimension != 2000L) stop("The numerical release checks have not passed.")
if (!is.character(accepted$validated_code_commit) || length(accepted$validated_code_commit) != 1L ||
    system2("git", c("diff", "--quiet", accepted$validated_code_commit, "--", unlist(accepted$validated_paths))) != 0L)
  stop("The current numerical sources differ from the verified implementation.")
families <- c("sqr", "pooled", "qgee")
roots <- setNames(file.path("results/raw/simulation", paste0("submission_20261028_main_", families)), families)
roots <- roots[!file.exists(file.path(roots, "completed.rds"))]
if (!length(roots)) stop("No unfinished comparator queue requires activation.")
processes <- ps::ps()
processes <- processes[processes$name == "R" & processes$pid != Sys.getpid(), ]
managers <- setNames(vector("list", length(roots)), names(roots))
for (i in seq_len(nrow(processes))) {
  handle <- processes$ps_handle[[i]]
  if (is.null(handle) || !ps::ps_is_running(handle)) next
  arguments <- ps::ps_cmdline(handle)
  if (!identical(ps::ps_cwd(handle), normalizePath("."))) next
  marker <- match("--args", arguments)
  if (is.na(marker)) next
  for (family in names(roots)) {
    continuation <- "--file=research/submission_completion_20261027/resume_simulation.R" %in% arguments &&
      arguments[marker + 1L] == roots[[family]]
    original <- "--file=research/submission_completion_20261027/launch_comparators.R" %in% arguments &&
      arguments[marker + 1L] == family
    if (continuation || original) {
      if (!is.null(managers[[family]])) stop("Two dispatchers own the same comparator queue.")
      managers[[family]] <- handle
    }
  }
}
workers <- list()
for (family in names(roots)) {
  group <- submission_processes(roots[[family]], c("run_comparator_batch.R", "run_remaining_batch.R"))
  for (key in names(group)) workers[[paste0(family, "_", key)]] <-
    c(group[[key]], list(family = family, batch = as.integer(key)))
}
cat("Current comparator workers:", length(workers), "\n")
print(data.frame(queue = names(roots), dispatcher = vapply(managers,
  function(x) if (is.null(x)) NA_integer_ else ps::ps_pid(x), integer(1))))
if (!"apply" %in% commandArgs(trailingOnly = TRUE)) quit(status = 0L)
changes <- system2("git", c("status", "--porcelain", "--untracked-files=normal"), stdout = TRUE)
if (length(changes)) stop("Commit the implementation before activating its calculation processes.")

activate <- function() {
  held <- list()
  on.exit(for (handle in held) if (isTRUE(ps::ps_is_running(handle))) ps::ps_resume(handle), add = TRUE)
  for (handle in managers) if (!is.null(handle)) {
    ps::ps_suspend(handle)
    held[[length(held) + 1L]] <- handle
  }
  for (worker in workers) {
    ps::ps_suspend(worker$handle)
    held[[length(held) + 1L]] <- worker$handle
  }
  snapshots <- list()
  if (length(workers)) {
    requests <- lapply(roots, function(root) readRDS(file.path(root, "request.rds")))
    methods <- c(sqr = "SQR-DEBIASED-IID", pooled = "POOLED-QR-LASSO", qgee = "QGEE-SCAD")
    for (key in names(workers)) {
      worker <- workers[[key]]
      family <- worker$family
      root <- roots[[family]]
      request <- requests[[family]]
      jobs <- request$jobs[request$jobs$batch == worker$batch, ]
      paths <- file.path(root, "replications", jobs$experiment_id, sprintf("replicate_%05d", jobs$replicate))
      completed <- file.exists(file.path(paths, "result.rds"))
      unfinished <- which(dir.exists(paths) & !completed)
      if (length(unfinished) > 1L) stop("A worker has multiple incomplete registered datasets.")
      if (length(unfinished)) {
        saved <- readRDS(file.path(paths[unfinished], "request.rds"))
        if (saved$replicate != jobs$replicate[unfinished] ||
            saved$config$experiment_id != jobs$experiment_id[unfinished] || saved$base_seed != request$base_seed)
          stop("An unfinished dataset differs from its registered identity.")
      }
      last <- tail(which(completed), 1L)
      if (length(last)) {
        value <- readRDS(file.path(paths[last], "result.rds"))
        if (nrow(value$metrics) != 1L || value$metrics$method_id != methods[[family]] ||
            value$metrics$replicate != jobs$replicate[last] || value$metrics$experiment_id != jobs$experiment_id[last])
          stop("A completed dataset differs from its registered identity.")
      }
      snapshots[[key]] <- list(family = family, batch = worker$batch, pid = ps::ps_pid(worker$handle),
        completed_datasets = sum(completed), unfinished_directories = paths[unfinished])
    }
    for (family in names(roots)) {
      root <- roots[[family]]
      path <- file.path(root, paste0("numerical_transition_", Sys.getpid(), ".rds"))
      saved_workers <- Filter(function(worker) worker$family == family, snapshots)
      submission_atomic_rds(list(created = Sys.time(), workers = saved_workers,
        numerical_specification = accepted$numerical_specification,
        implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)), path)
      submission_backup_run(root, force = TRUE)
    }
  }
  for (worker in workers) ps::ps_kill(worker$handle)
  for (handle in managers) if (!is.null(handle) && isTRUE(ps::ps_is_running(handle))) ps::ps_kill(handle)
  launched <- list()
  for (family in names(roots)) {
    log <- file.path("results/logs/recovery_20261028", paste0("accelerated_", family, "_", Sys.getpid(), ".log"))
    process <- processx::process$new(file.path(R.home("bin"), "Rscript"),
      c("research/submission_completion_20261027/resume_simulation.R", roots[[family]], "12"),
      cleanup = FALSE, stdout = log, stderr = "2>&1")
    pid <- process$get_pid()
    if (Sys.info()[["sysname"]] == "Darwin") processx::process$new("/usr/bin/caffeinate",
      c("-s", "-w", as.character(pid)), cleanup = FALSE, stdout = paste0(log, ".power"), stderr = "2>&1")
    launched[[family]] <- data.frame(queue = family, pid = pid, log = log, started = as.character(Sys.time()))
  }
  record <- do.call(rbind, launched)
  write.csv(record, "results/scratch/sqr_acceleration_20261030/active_dispatchers.csv", row.names = FALSE)
  print(record)
}
activate()

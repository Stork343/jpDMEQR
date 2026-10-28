submission_processes <- function(root, scripts) {
  listing <- ps::ps()
  listing <- listing[listing$name == "R" & listing$pid != Sys.getpid(), ]
  found <- list()
  expected <- paste0("--file=research/submission_completion_20261027/", scripts)
  for (i in seq_len(nrow(listing))) {
    handle <- listing$ps_handle[[i]]
    if (is.null(handle) || !ps::ps_is_running(handle)) next
    arguments <- ps::ps_cmdline(handle)
    if (!any(arguments %in% expected) || !identical(ps::ps_cwd(handle), normalizePath("."))) next
    marker <- match("--args", arguments)
    if (is.na(marker) || length(arguments) < marker + 2L || arguments[marker + 1L] != root) next
    identifier <- as.integer(arguments[marker + 2L])
    if (is.na(identifier) || identifier < 1L) stop("A running submission worker has an invalid identity.")
    key <- as.character(identifier)
    if (!is.null(found[[key]])) stop("Two running workers have the same registered identity.")
    found[[key]] <- list(handle = handle, process = NULL)
  }
  found
}

submission_atomic_csv <- function(value, path) {
  temporary <- paste0(path, ".", Sys.getpid(), ".writing")
  write.csv(value, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop("Unable to publish the dispatch status.")
}

submission_worker_alive <- function(worker) {
  if (!is.null(worker$process)) isTRUE(worker$process$is_alive()) else
    isTRUE(ps::ps_is_running(worker$handle))
}

submission_launchers_alive <- function(root) {
  listing <- ps::ps()
  listing <- listing[listing$name == "R" & listing$pid != Sys.getpid(), ]
  found <- integer()
  for (i in seq_len(nrow(listing))) {
    handle <- listing$ps_handle[[i]]
    if (is.null(handle) || !ps::ps_is_running(handle)) next
    arguments <- ps::ps_cmdline(handle)
    simulation <- "--file=research/submission_completion_20261027/launch_simulation.R" %in% arguments
    comparator <- "--file=research/submission_completion_20261027/launch_comparators.R" %in% arguments
    if (!simulation && !comparator) next
    if (!identical(ps::ps_cwd(handle), normalizePath("."))) next
    marker <- match("--args", arguments)
    if (is.na(marker) || length(arguments) < marker + 2L) next
    destination <- if (simulation) file.path("results/raw/simulation",
      paste("submission_20261027", arguments[marker + 1L], arguments[marker + 2L], sep = "_")) else
      file.path("results/raw/simulation", paste0("submission_20261028_main_", arguments[marker + 1L]))
    if (identical(normalizePath(destination, mustWork = FALSE), normalizePath(root, mustWork = FALSE)))
      found <- c(found, listing$pid[i])
  }
  found
}

submission_comparator_capacity <- function(family, workers) {
  families <- c("mixed", "pooled", "sqr", "qgee", "bias_adjusted")
  if (!family %in% families || !workers %in% 1:8) stop("Invalid comparator capacity request.")
  if (!file.exists("results/raw/simulation/submission_20261027_main_joint/completed.rds"))
    return(min(workers, 1L))
  if (family != "sqr") return(min(workers, 2L))
  other <- setdiff(families, "sqr")
  remaining <- sum(!file.exists(file.path("results/raw/simulation",
    paste0("submission_20261028_main_", other), "completed.rds")))
  min(workers, 12L - 2L * remaining)
}

submission_atomic_rds <- function(value, path, compress = "xz") {
  if (file.exists(path)) stop("The completed output already exists: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".", Sys.getpid(), ".writing")
  saveRDS(value, temporary, compress = compress)
  if (!file.rename(temporary, path)) stop("Unable to publish the complete output: ", path)
  invisible(path)
}

submission_backup_run <- function(root, force = FALSE) {
  base <- Sys.getenv("JPDMEQR_BACKUP_DIR", unset = file.path(path.expand("~"), "jpDMEQR-backups"))
  destination <- file.path(base, basename(dirname(root)), basename(root))
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(file.path(destination, "backup.lock"), timeout = if (force) 60000 else 0)
  if (is.null(lock)) {
    if (force) stop("Unable to acquire the output backup lock.")
    return(invisible(FALSE))
  }
  on.exit(filelock::unlock(lock), add = TRUE)
  outputs <- list.files(root, "^(result|fit_[0-9]+)[.]rds$", recursive = TRUE)
  metadata <- c("request.rds", "dispatcher.rds", "config_snapshot.csv", "scheduled_replications.csv")
  metadata <- metadata[file.exists(file.path(root, metadata))]
  files <- sort(unique(c(outputs, metadata)))
  if (!length(files)) stop("No completed files are available for backup.")
  manifest <- file.path(root, paste0("backup_files_", Sys.getpid(), ".txt"))
  writeLines(files, manifest)
  status <- system2("rsync", c("-a", shQuote(paste0("--files-from=", normalizePath(manifest))),
    shQuote(paste0(normalizePath(root), "/")), shQuote(paste0(normalizePath(destination), "/"))))
  if (status != 0L) stop("The completed output backup failed.")
  submission_atomic_csv(data.frame(source = normalizePath(root), files = length(files),
    completed_outputs = length(outputs), saved = as.character(Sys.time())), file.path(destination, "backup_status.csv"))
  invisible(TRUE)
}

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

submission_atomic_rds <- function(value, path, compress = "xz") {
  if (file.exists(path)) stop("The completed output already exists: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".", Sys.getpid(), ".writing")
  saveRDS(value, temporary, compress = compress)
  if (!file.rename(temporary, path)) stop("Unable to publish the complete output: ", path)
  invisible(path)
}

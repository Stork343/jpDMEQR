submission_application_completed <- function(root, job_id) {
  if (is.null(root)) return(FALSE)
  directory <- file.path(root, sprintf("job_%04d", job_id))
  if (file.exists(file.path(directory, "completion.rds"))) return(TRUE)
  log <- file.path(root, sprintf("job_%04d.log", job_id))
  file.exists(log) && any(grepl("^Completed comparator job|^Completed source QGEE job|^Execution halted$",
    tail(readLines(log, warn = FALSE), 12L)))
}

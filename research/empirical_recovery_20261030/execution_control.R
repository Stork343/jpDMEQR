empirical_require_active <- function() {
  path <- "research/empirical_recovery_20261030/execution_control.json"
  if (file.exists(path) && jsonlite::read_json(path)$state != "active")
    stop("The user has terminated this empirical study.")
}

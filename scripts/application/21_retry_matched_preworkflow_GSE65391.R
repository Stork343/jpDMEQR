root <- normalizePath(".", mustWork = TRUE)
library(processx)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Supply original run ID, worker count and a new source-retry run ID.")
original <- file.path(root, "results", "raw", "application_matched", args[1])
workers <- as.integer(args[2])
output <- file.path(root, "results", "raw", "application_matched", args[3])
if (!is.finite(workers) || workers < 1L || dir.exists(output)) stop("Invalid or existing retry execution settings.")
source_identity <- validate_application_numeric_source_identity_v2(root, "078559a")
scheduled <- read.csv(file.path(original, "scheduled_replicates.csv"), stringsAsFactors = FALSE)
processes <- read.csv(file.path(original, "job_results.csv"), stringsAsFactors = FALSE)
failed <- processes[processes$exit_status != 0L, ]
failed_directory <- file.path(original, failed$case_id, sprintf("replicate_%04d", failed$replicate))
before_workflow <- !file.exists(file.path(failed_directory, "run_request.rds"))
failed <- failed[before_workflow, ]
if (!nrow(failed)) stop("No actual pre-workflow source failures require a retry.")
logs <- file.path(original, "console", sprintf("%s_rep%04d.log", failed$case_id, failed$replicate))
reasons <- character(length(logs))
for (i in seq_along(logs)) {
  lines <- readLines(logs[i], warn = FALSE)
  if (any(grepl("unexpected input", lines, fixed = TRUE)) &&
      any(grepl("application_numeric_git_expressions_v2", lines, fixed = TRUE))) {
    reasons[i] <- "source_parser_input_line_endings_before_workflow"
  } else if (any(grepl("APPLICATION_NUMERIC_SOURCE_FAILED: committed numerical dependencies must match the declared ancestor build.",
                       lines, fixed = TRUE))) {
    reasons[i] <- "source_identity_file_transition_before_workflow"
  } else {
    stop("A pre-workflow failure has a different cause and requires its own diagnosis.")
  }
}
key <- function(table) paste(table$case_id, table$replicate, sep = ":")
jobs <- scheduled[match(key(failed), key(scheduled)), , drop = FALSE]
if (anyNA(jobs$case_id) || anyDuplicated(key(jobs))) stop("The retry keys do not match the original frozen schedule.")
config <- yaml::read_yaml(file.path(root, "config", "application_matched_prospective_20261022.yml"))
jobs$seed <- config$base_seed + (jobs$q + 2L * (match(jobs$error_family, config$error_families) - 1L)) * 1000000L + jobs$replicate
jobs$original_console <- logs
jobs$retry_reason <- reasons
jobs$state <- "pending"
jobs$exit_status <- NA_integer_
jobs$actual_source_commit <- NA_character_
dir.create(output, recursive = TRUE)
dir.create(file.path(output, "console"))
write.csv(jobs, file.path(output, "scheduled_replicates.csv"), row.names = FALSE)
jsonlite::write_json(list(
  run_id = args[3], original_run = args[1], scheduled = nrow(jobs), workers = workers,
  source_identity = source_identity, solver_configuration = solver_configuration_v2(),
  eligibility = "actual_source_identity_termination_without_workflow_request_or_fit",
  reason_counts = as.list(table(reasons)),
  execution_order = "after_original_coordinator_completion", retry_definition = "reuse_original_frozen_seed_and_target_asset",
  registered_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
while (!file.exists(file.path(original, "completed_manifest.json"))) Sys.sleep(15)
execution_identity <- validate_application_numeric_source_identity_v2(root, "078559a")
jsonlite::write_json(list(source_identity = execution_identity,
                          started_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                     file.path(output, "execution_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
active <- list()
while (any(jobs$state != "finished")) {
  for (name in names(active)) {
    child <- active[[name]]
    if (!child$is_alive()) {
      index <- as.integer(name)
      jobs$state[index] <- "finished"
      jobs$exit_status[index] <- child$get_exit_status()
      request_path <- file.path(output, jobs$case_id[index], sprintf("replicate_%04d", jobs$replicate[index]), "run_request.rds")
      if (!file.exists(request_path)) {
        write.csv(jobs, file.path(output, "source_failure_snapshot.csv"), row.names = FALSE)
        stop("A source retry did not reach the workflow; inspect its actual console before dispatching further jobs.")
      }
      request <- readRDS(request_path)
      if (request$seed != jobs$seed[index] || !identical(request$target_asset, jobs$target_asset[index])) {
        stop("A retry changed its original seed or target asset.")
      }
      jobs$actual_source_commit[index] <- request$source_identity$actual_source_commit
      utils::write.table(jobs[index, , drop = FALSE], file.path(output, "job_results.csv"), sep = ",",
                          row.names = FALSE, append = file.exists(file.path(output, "job_results.csv")),
                          col.names = !file.exists(file.path(output, "job_results.csv")))
      active[[name]] <- NULL
    }
  }
  available <- which(jobs$state == "pending")
  for (index in head(available, max(0L, workers - length(active)))) {
    job <- jobs[index, ]
    jobs$started_source_commit[index] <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
    arguments <- c("scripts/application/13_run_matched_GSE65391.R", paste0("--case-id=", job$case_id),
                   paste0("--replicate=", job$replicate), paste0("--target-asset=", job$target_asset),
                   paste0("--run-id=", args[3]), "--application-build-commit=078559a")
    log <- file.path(output, "console", sprintf("%s_rep%04d.log", job$case_id, job$replicate))
    active[[as.character(index)]] <- processx::process$new(
      file.path(R.home("bin"), "Rscript"), arguments, stdout = log, stderr = log,
      wd = root, cleanup = FALSE, cleanup_tree = FALSE
    )
    jobs$state[index] <- "running"
  }
  write.csv(jobs, file.path(output, "status_snapshot.csv"), row.names = FALSE)
  Sys.sleep(1)
}
write.csv(jobs, file.path(output, "all_job_results.csv"), row.names = FALSE)
jsonlite::write_json(list(
  run_id = args[3], original_run = args[1], scheduled = nrow(jobs),
  successful_processes = sum(jobs$exit_status == 0L), failed_workflow_processes = sum(jobs$exit_status != 0L),
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

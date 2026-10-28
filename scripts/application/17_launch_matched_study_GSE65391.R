root <- normalizePath(".", mustWork = TRUE)
library(processx)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
args <- commandArgs(trailingOnly = TRUE)
target_run <- if (length(args)) args[1] else "matched_targets_20261022_v1"
workers <- if (length(args) > 1L) as.integer(args[2]) else 2L
run_id <- if (length(args) > 2L) args[3] else "matched_study_20261022_v1"
build_commit <- if (length(args) > 3L) args[4] else NULL
source_identity <- validate_application_numeric_source_identity_v2(root, build_commit)
Sys.setenv(JPDMEQR_APPLICATION_BUILD_COMMIT = source_identity$application_build_commit)
if (!is.finite(workers) || workers < 1L ||
    !grepl("^[A-Za-z0-9_.-]+$", run_id)) stop("Invalid matched-study launcher settings.")
target_directory <- file.path(root, "results", "raw", "application_matched", target_run)
targets <- read.csv(file.path(target_directory, "scheduled_target_jobs.csv"), stringsAsFactors = FALSE)
registry <- read.csv(file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1", "matched_registry.csv"),
                      stringsAsFactors = FALSE)
jobs <- registry[rep(seq_len(nrow(registry)), times = 200L), , drop = FALSE]
jobs$replicate <- rep(seq_len(200L), each = nrow(registry))
target_case <- sprintf("AM_%s_q%d_%s", jobs$model, jobs$q, jobs$error_family)
target_run_id <- targets$run_id[match(target_case, targets$case_id)]
if (anyNA(target_run_id) || nrow(registry) != 24L || any(registry$replications != 200L)) {
  stop("All eight targets and the complete 24 by 200 study matrix are required.")
}
jobs$target_asset <- file.path(root, "results", "raw", "application_matched", target_run_id,
                                "matched_target_asset.rds")
jobs$state <- "pending"
jobs$exit_status <- NA_integer_
jobs$started_source_commit <- NA_character_
jobs$actual_source_commit <- NA_character_
jobs$application_build_commit <- source_identity$application_build_commit
output <- file.path(root, "results", "raw", "application_matched", run_id)
if (dir.exists(output)) stop("The versioned matched-study directory already exists.")
dir.create(output, recursive = TRUE)
dir.create(file.path(output, "console"))
write.csv(jobs, file.path(output, "scheduled_replicates.csv"), row.names = FALSE)
jsonlite::write_json(
  list(run_id = run_id, workers = workers, scheduled = nrow(jobs),
       implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       source_identity = source_identity,
       solver_configuration = solver_configuration_v2(),
       started_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
  file.path(output, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE
)
active <- list()
last_snapshot <- Sys.time()
while (any(jobs$state != "finished")) {
  for (name in names(active)) {
    child <- active[[name]]
    if (!child$is_alive()) {
      index <- as.integer(name)
      jobs$state[index] <- "finished"
      jobs$exit_status[index] <- child$get_exit_status()
      request_path <- file.path(output, jobs$case_id[index], sprintf("replicate_%04d", jobs$replicate[index]), "run_request.rds")
      if (file.exists(request_path)) {
        jobs$actual_source_commit[index] <- readRDS(request_path)$source_identity$actual_source_commit
      }
      row <- jobs[index, , drop = FALSE]
      utils::write.table(row, file.path(output, "job_results.csv"), sep = ",", row.names = FALSE,
                          append = file.exists(file.path(output, "job_results.csv")),
                          col.names = !file.exists(file.path(output, "job_results.csv")))
      active[[name]] <- NULL
    }
  }
  available <- which(jobs$state == "pending" & file.exists(jobs$target_asset) &
                       file.exists(file.path(dirname(jobs$target_asset), "completed_manifest.json")))
  free <- max(0L, workers - length(active))
  for (index in utils::head(available, free)) {
    job <- jobs[index, ]
    jobs$started_source_commit[index] <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
    log <- file.path(output, "console", sprintf("%s_rep%04d.log", job$case_id, job$replicate))
    arguments <- c("scripts/application/13_run_matched_GSE65391.R", paste0("--case-id=", job$case_id),
                   paste0("--replicate=", job$replicate), paste0("--target-asset=", job$target_asset),
                   paste0("--run-id=", run_id))
    child <- processx::process$new(file.path(R.home("bin"), "Rscript"), arguments,
                                    stdout = log, stderr = log, wd = root,
                                    cleanup = FALSE, cleanup_tree = FALSE)
    jobs$state[index] <- "running"
    active[[as.character(index)]] <- child
  }
  if (as.numeric(difftime(Sys.time(), last_snapshot, units = "secs")) >= 30) {
    write.csv(jobs, file.path(output, "status_snapshot.csv"), row.names = FALSE)
    last_snapshot <- Sys.time()
  }
  target_status_path <- file.path(target_directory, "target_job_results.csv")
  if (!length(active) && !length(available) && any(jobs$state == "pending") && file.exists(target_status_path)) {
    target_status <- read.csv(target_status_path, stringsAsFactors = FALSE)
    if (any(target_status$exit_status != 0L)) stop("A required population target failed; pending study jobs remain recorded.")
  }
  Sys.sleep(1)
}
write.csv(jobs, file.path(output, "all_job_results.csv"), row.names = FALSE)
jsonlite::write_json(list(run_id = run_id, scheduled = nrow(jobs),
                          completed_processes = sum(jobs$state == "finished"),
                          successful_processes = sum(jobs$exit_status == 0L),
                          failed_processes = sum(jobs$exit_status != 0L),
                          target_run = target_run,
                          completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                     file.path(output, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

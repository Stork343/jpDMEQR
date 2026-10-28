root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
library(processx)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
source(file.path(root, "scripts", "application", "_paired_method_gate_20261022.R"))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
workload_directory <- get_arg("workload-directory")
checks_path <- get_arg("checks-manifest")
phase <- get_arg("phase")
workers <- as.integer(get_arg("workers", "2"))
run_id <- get_arg("run-id")
if (is.null(workload_directory) || is.null(checks_path) || is.null(run_id) || is.null(phase) ||
    !phase %in% c("module", "validation", "exploratory") || !is.finite(workers) || workers < 1L) {
  stop("Supply a frozen workload directory, actual checks, phase, positive workers and a new dispatcher run ID.")
}
implementation <- validate_application_paired_gate_20261022(root, checks_path, inference = TRUE)
commands <- readRDS(file.path(workload_directory, "commands.rds"))[[phase]]
expected <- c(module = 40L, validation = 1050L, exploratory = 120L)
if (length(commands) != expected[[phase]]) stop("The complete prespecified phase workload is required.")
if (phase == "validation") {
  ordering <- do.call(rbind, lapply(commands, function(command) command$job))
  index <- order(ordering$repetition, ordering$fold, ordering$model, ordering$tau, ordering$specification_id)
  commands <- commands[index]
}
output <- file.path(root, "results", "raw", "application_dispatch", run_id)
if (dir.exists(output)) stop("The application dispatcher run already exists.")
dir.create(output, recursive = TRUE)
dir.create(file.path(output, "console"))
jobs <- do.call(rbind, lapply(commands, function(command) command$job))
rownames(jobs) <- NULL
jobs$state <- "pending"
jobs$exit_status <- NA_integer_
jobs$actual_source_commit <- NA_character_
jobs$failure_stage <- NA_character_
write.csv(jobs, file.path(output, "scheduled_jobs.csv"), row.names = FALSE)
saveRDS(list(implementation = implementation, commands = commands,
             workload_directory = normalizePath(workload_directory)), file.path(output, "run_request.rds"))
if (phase == "exploratory") {
  arguments <- c("scripts/application/09_fit_exploratory_pair_GSE65391.R", "--prepare-only=true",
                 paste0("--run-id=", jobs$run_id[1]), paste0("--checks-manifest=", normalizePath(checks_path)))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(arguments),
                     stdout = file.path(output, "preprocessing_console.log"),
                     stderr = file.path(output, "preprocessing_console.log"))
  if (status != 0L) stop("The complete frozen exploratory preprocessing did not finish.")
}
active <- list()
last_snapshot <- Sys.time()
while (any(jobs$state != "finished")) {
  for (name in names(active)) {
    child <- active[[name]]
    if (!child$is_alive()) {
      index <- as.integer(name)
      jobs$state[index] <- "finished"
      jobs$exit_status[index] <- child$get_exit_status()
      completed <- file.path(root, jobs$expected_manifest[index])
      failed <- file.path(dirname(completed), "failed_manifest.json")
      if (jobs$exit_status[index] == 0L) {
        if (!file.exists(completed)) stop("A successful application process lacks its actual completion record.")
        manifest <- jsonlite::read_json(completed, simplifyVector = TRUE)
        actual <- if (phase == "exploratory") manifest$implementation$actual_source_commit else manifest$actual_source_commit
        jobs$actual_source_commit[index] <- actual
      } else {
        if (!file.exists(failed)) {
          write.csv(jobs, file.path(output, "source_failure_snapshot.csv"), row.names = FALSE)
          stop("An application process failed before its workflow failure record; inspect its actual console before further dispatch.")
        }
        manifest <- jsonlite::read_json(failed, simplifyVector = TRUE)
        jobs$failure_stage[index] <- manifest$failure_stage
        jobs$actual_source_commit[index] <- manifest$implementation$actual_source_commit
      }
      utils::write.table(jobs[index, , drop = FALSE], file.path(output, "job_results.csv"), sep = ",",
                          row.names = FALSE, append = file.exists(file.path(output, "job_results.csv")),
                          col.names = !file.exists(file.path(output, "job_results.csv")))
      active[[name]] <- NULL
    }
  }
  available <- which(jobs$state == "pending")
  for (index in head(available, max(0L, workers - length(active)))) {
    command <- commands[[index]]
    arguments <- c(command$script, command$arguments, paste0("--checks-manifest=", normalizePath(checks_path)))
    log <- file.path(output, "console", sprintf("job%04d.log", jobs$job_id[index]))
    active[[as.character(index)]] <- processx::process$new(
      file.path(R.home("bin"), "Rscript"), arguments, wd = root, stdout = log, stderr = log,
      cleanup = FALSE, cleanup_tree = FALSE
    )
    jobs$state[index] <- "running"
  }
  if (as.numeric(difftime(Sys.time(), last_snapshot, units = "secs")) >= 30) {
    write.csv(jobs, file.path(output, "status_snapshot.csv"), row.names = FALSE)
    last_snapshot <- Sys.time()
  }
  Sys.sleep(1)
}
if (phase == "exploratory") {
  arguments <- c("scripts/application/09_fit_exploratory_pair_GSE65391.R", "--summarise-only=true",
                 paste0("--run-id=", jobs$run_id[1]), paste0("--checks-manifest=", normalizePath(checks_path)))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(arguments),
                     stdout = file.path(output, "summary_console.log"), stderr = file.path(output, "summary_console.log"))
  if (status != 0L) stop("The complete exploratory result summary did not finish.")
}
write.csv(jobs, file.path(output, "all_job_results.csv"), row.names = FALSE)
jsonlite::write_json(list(
  run_id = run_id, phase = phase, scheduled = nrow(jobs), successful = sum(jobs$exit_status == 0L),
  failed_workflows = sum(jobs$exit_status != 0L), numeric_certificate_commit = implementation$numeric_certificate_commit,
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

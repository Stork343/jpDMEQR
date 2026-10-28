library(filelock)
library(ps)
source("research/submission_completion_20261027/dispatch_state.R")
directory <- "results/scratch/submission_completion_20261027"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
lock <- filelock::lock(file.path(directory, "sqr_queue.lock"), timeout = 0)
if (is.null(lock)) stop("The SQR continuation is already running.")
output <- file.path(directory, "sqr_procedure_check")
checks <- file.path(output, "checks.csv")
status_path <- file.path(directory, "sqr_continuation_status.csv")
submission_atomic_csv(data.frame(state = "waiting_numerical_validation", changed = as.character(Sys.time())), status_path)
while (!file.exists(checks)) Sys.sleep(15)
validation <- read.csv(checks, stringsAsFactors = FALSE)
if (nrow(validation) != 1L || !all(c("maximum_difference", "tuning_difference") %in% names(validation)) ||
    !is.finite(validation$maximum_difference) || validation$maximum_difference > 2e-5 ||
    !is.finite(validation$tuning_difference) || validation$tuning_difference > 1e-12)
  stop("The complete SQR numerical comparison did not pass.")
candidate_path <- file.path(output, "candidate.rds")
candidate <- readRDS(candidate_path)
if (nrow(candidate$metrics) != 1L || candidate$metrics$method_id != "SQR-DEBIASED-IID" ||
    candidate$metrics$experiment_id != "A01" || !isTRUE(candidate$metrics$converged))
  stop("The completed SQR validation candidate is unavailable.")
preflight_path <- file.path(directory, "comparator_preflight", "SQR-DEBIASED-IID.rds")
if (!file.exists(preflight_path)) submission_atomic_rds(candidate, preflight_path)
jsonlite::write_json(list(source = candidate_path, numerical_checks = checks,
  preflight_reuses_completed_candidate = TRUE), paste0(preflight_path, ".source.json"),
  auto_unbox = TRUE, pretty = TRUE)
root <- "results/raw/simulation/submission_20261028_main_sqr"
if (file.exists(file.path(root, "completed.rds"))) quit(status = 0L)
if (length(submission_launchers_alive(root))) quit(status = 0L)
arguments <- if (file.exists(file.path(root, "request.rds")))
  c("research/submission_completion_20261027/resume_simulation.R", root, "8") else
  c("research/submission_completion_20261027/launch_comparators.R", "sqr", "8")
submission_atomic_csv(data.frame(state = "running", changed = as.character(Sys.time())), status_path)
result <- system2(file.path(R.home("bin"), "Rscript"), shQuote(arguments))
if (result != 0L || !file.exists(file.path(root, "completed.rds")))
  stop("The SQR queue ended before its registered calculations completed.")
submission_atomic_csv(data.frame(state = "completed", changed = as.character(Sys.time())), status_path)
filelock::unlock(lock)

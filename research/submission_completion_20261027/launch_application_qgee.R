library(processx)
source("research/submission_completion_20261027/application_state.R")
args <- commandArgs(trailingOnly = TRUE)
phase <- args[1L]
workers <- as.integer(args[2L])
reuse <- if (length(args) >= 3L) args[3L] else "reuse"
if (!phase %in% c("probe", "final") || !workers %in% 1:8) stop("Invalid source QGEE dispatch.")
if (!reuse %in% c("reuse", "fresh")) stop("Invalid QGEE cache policy.")
source_root <- if (reuse == "reuse") "results/raw/application_comparison/submission_20261027_v2" else NULL
if (!is.null(source_root) && !dir.exists(source_root)) stop("The requested QGEE cache is unavailable.")
checks <- read.csv("results/scratch/submission_completion_20261027/qgee_penalty_policy_check/checks.csv")
if (checks$maximum_coefficient_difference > 1e-5) stop("The QGEE penalty-path reference check did not pass.")
all <- expand.grid(method = c("LQMM", "POOLED-QR-LASSO", "QGEE-SCAD"), model = c("A", "B"),
  tau = c(.25, .5, .75), processing = c("reference", "official_R1_R2_normexp_training_quantile", "processed_median_across_mapped_probes"),
  repetition = 1:5, fold = 1:5, stringsAsFactors = FALSE)
all$job_id <- seq_len(nrow(all))
jobs <- all[all$method == "QGEE-SCAD", ]
if (phase == "probe") jobs <- jobs[1L, ]
root <- if (phase == "probe") "results/scratch/submission_completion_20261027/qgee_application_job_check" else
  "results/raw/application_comparison/submission_20261028_qgee"
if (dir.exists(root)) stop("The source QGEE dispatch already exists.")
dir.create(root, recursive = TRUE)
jobs$state <- "pending"
jobs$exit_status <- NA_integer_
saveRDS(list(jobs = jobs, workers = workers, started = Sys.time(),
  source_root = source_root,
  penalty_failure_policy = "source_HBIC_individual_penalties", reference_check = checks,
  implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), session = sessionInfo()), file.path(root, "request.rds"))
active <- list()
source_ready <- function(job_id) {
  if (is.null(source_root)) return(TRUE)
  submission_application_completed(source_root, job_id)
}
repeat {
  for (key in names(active)) {
    process <- active[[key]]
    if (process$is_alive()) next
    index <- as.integer(key)
    jobs$state[index] <- "completed"
    jobs$exit_status[index] <- process$get_exit_status()
    active[[key]] <- NULL
    write.csv(jobs, file.path(root, "status.csv"), row.names = FALSE)
  }
  if (all(jobs$state == "completed")) break
  pending <- head(which(jobs$state == "pending"), 40L)
  ready <- pending[vapply(jobs$job_id[pending], source_ready, logical(1))]
  for (index in head(ready, max(0L, workers - length(active)))) {
    active[[as.character(index)]] <- processx::process$new(file.path(R.home("bin"), "Rscript"),
      c("research/submission_completion_20261027/run_application_qgee.R", root, as.character(jobs$job_id[index])),
      stdout = file.path(root, sprintf("job_%04d.log", jobs$job_id[index])), stderr = "2>&1")
    jobs$state[index] <- "running"
  }
  Sys.sleep(.5)
}
write.csv(jobs, file.path(root, "status.csv"), row.names = FALSE)
saveRDS(list(successful = sum(jobs$exit_status == 0L), failed = sum(jobs$exit_status != 0L),
  completed_jobs = nrow(jobs), finished = Sys.time()), file.path(root, "completed.rds"))
print(table(jobs$exit_status))

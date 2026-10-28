library(processx)
args <- commandArgs(trailingOnly = TRUE)
phase <- args[1L]
workers <- as.integer(args[2L])
if (!phase %in% c("prepare", "comparator_probe", "comparators") || !workers %in% 1:8) stop("Invalid application dispatch.")
if (phase == "prepare") {
  jobs <- expand.grid(repetition = 1:5, fold = 1:5)
  jobs <- jobs[!(jobs$repetition == 1L & jobs$fold == 1L), ]
  root <- "results/scratch/submission_application_20261027/prepare_dispatch"
  script <- "research/submission_completion_20261027/prepare_application.R"
  arguments <- lapply(seq_len(nrow(jobs)), function(i) as.character(unlist(jobs[i, ])))
} else {
  checks <- read.csv("results/scratch/submission_application_20261027/checks/clinical_units.csv")
  if (checks$value[1L] > checks$bound[1L]) stop("Clinical coordinate verification failed.")
  jobs <- expand.grid(method = c("LQMM", "POOLED-QR-LASSO", "QGEE-SCAD"), model = c("A", "B"),
    tau = c(.25, .5, .75), processing = c("reference", "official_R1_R2_normexp_training_quantile", "processed_median_across_mapped_probes"),
    repetition = 1:5, fold = 1:5, stringsAsFactors = FALSE)
  if (phase == "comparator_probe") jobs <- jobs[1:3, ]
  root <- if (phase == "comparator_probe") "results/scratch/submission_application_20261027/comparator_probe" else
    "results/raw/application_comparison/submission_20261027_v2"
  script <- "research/submission_completion_20261027/run_application_comparator.R"
  arguments <- lapply(seq_len(nrow(jobs)), function(i) c(root, as.character(i)))
}
if (dir.exists(root)) stop("The application dispatch already exists.")
dir.create(root, recursive = TRUE)
jobs$state <- "pending"
jobs$exit_status <- NA_integer_
saveRDS(list(jobs = jobs, script = script, arguments = arguments, workers = workers,
             started = Sys.time(), session = sessionInfo()), file.path(root, "request.rds"))
active <- list()
repeat {
  for (id in names(active)) {
    process <- active[[id]]
    if (process$is_alive()) next
    i <- as.integer(id)
    jobs$state[i] <- "completed"
    jobs$exit_status[i] <- process$get_exit_status()
    active[[id]] <- NULL
    write.csv(jobs, file.path(root, "status.csv"), row.names = FALSE)
  }
  if (all(jobs$state == "completed")) break
  for (i in head(which(jobs$state == "pending"), max(0L, workers - length(active)))) {
    active[[as.character(i)]] <- processx::process$new(file.path(R.home("bin"), "Rscript"), c(script, arguments[[i]]),
      stdout = file.path(root, sprintf("job_%04d.log", i)), stderr = "2>&1")
    jobs$state[i] <- "running"
  }
  Sys.sleep(.5)
}
write.csv(jobs, file.path(root, "status.csv"), row.names = FALSE)
saveRDS(list(successful = sum(jobs$exit_status == 0L), failed = sum(jobs$exit_status != 0L), finished = Sys.time()), file.path(root, "completed.rds"))
print(table(jobs$exit_status))

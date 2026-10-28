library(processx)
args <- commandArgs(trailingOnly = TRUE)
phase <- args[1L]
if (!phase %in% c("processing_modules", "prediction", "slope_prediction", "exploratory")) stop("Unknown manuscript phase.")
processings <- c("reference", "official_R1_R2_normexp_training_quantile", "processed_median_across_mapped_probes")
if (phase == "processing_modules") {
  jobs <- expand.grid(processing = processings[-1L], model = c("A", "B"), tau = c(.25, .5, .75), stringsAsFactors = FALSE)
  script <- "research/joint_profile_manuscript_20261026/run_processing_module.R"
  arguments <- lapply(seq_len(nrow(jobs)), function(i) c(jobs$processing[i], jobs$model[i], as.character(jobs$tau[i])))
} else if (phase == "slope_prediction") {
  jobs <- expand.grid(processing = "reference", repetition = 1:5, fold = 1:5, stringsAsFactors = FALSE)
  script <- "research/joint_profile_manuscript_20261026/run_prediction_fold.R"
  arguments <- lapply(seq_len(nrow(jobs)), function(i) c(jobs$processing[i], as.character(jobs$repetition[i]), as.character(jobs$fold[i]), "slope_only"))
} else if (phase == "prediction") {
  jobs <- expand.grid(processing = processings, repetition = 1:5, fold = 1:5, stringsAsFactors = FALSE)
  jobs <- jobs[!(jobs$processing == "reference" & jobs$repetition == 1L & jobs$fold == 1L), ]
  script <- "research/joint_profile_manuscript_20261026/run_prediction_fold.R"
  arguments <- lapply(seq_len(nrow(jobs)), function(i) c(jobs$processing[i], as.character(jobs$repetition[i]), as.character(jobs$fold[i])))
} else {
  roles <- unlist(lapply(1:10, function(split) sprintf("split%02d_swap%d", split, 0:1)))
  jobs <- expand.grid(processing = processings, role = roles, stringsAsFactors = FALSE)
  script <- "research/joint_profile_manuscript_20261026/run_exploratory_role.R"
  arguments <- lapply(seq_len(nrow(jobs)), function(i) c(jobs$processing[i], jobs$role[i]))
}
output <- file.path("results/scratch/joint_profile_manuscript_20261026", paste0(phase, "_dispatch"))
if (dir.exists(output)) stop("This manuscript dispatcher already exists.")
dir.create(output, recursive = TRUE)
jobs$state <- "pending"
jobs$exit_status <- NA_integer_
saveRDS(list(jobs = jobs, script = script, arguments = arguments, started = Sys.time()), file.path(output, "request.rds"))
active <- list()
repeat {
  for (id in names(active)) {
    process <- active[[id]]
    if (process$is_alive()) next
    i <- as.integer(id)
    jobs$state[i] <- "completed"
    jobs$exit_status[i] <- process$get_exit_status()
    active[[id]] <- NULL
    write.csv(jobs, file.path(output, "status.csv"), row.names = FALSE)
  }
  if (all(jobs$state == "completed")) break
  for (i in head(which(jobs$state == "pending"), max(0L, 4L - length(active)))) {
    active[[as.character(i)]] <- processx::process$new(file.path(R.home("bin"), "Rscript"), c(script, arguments[[i]]),
                                                       stdout = file.path(output, sprintf("job_%03d.log", i)), stderr = "2>&1")
    jobs$state[i] <- "running"
  }
  Sys.sleep(.5)
}
write.csv(jobs, file.path(output, "status.csv"), row.names = FALSE)
saveRDS(list(successful = sum(jobs$exit_status == 0L), failed = sum(jobs$exit_status != 0L),
             finished = Sys.time()), file.path(output, "completed.rds"))
print(table(jobs$exit_status))

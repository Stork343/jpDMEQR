library(processx)
input <- "results/scratch/joint_profile_manuscript_20261026/module_inputs"
jobs <- read.csv(file.path(input, "jobs.csv"), stringsAsFactors = FALSE)
output <- "results/raw/application/joint_profile_20261026"
if (file.exists(file.path(output, "dispatcher.rds")) && !"resume" %in% commandArgs(trailingOnly = TRUE)) stop("The module dispatcher already exists.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(output, "dispatcher.rds"))) {
  saveRDS(list(jobs = jobs, workers = 4L, started = Sys.time()), file.path(output, "dispatcher.rds"))
} else if (!identical(readRDS(file.path(output, "dispatcher.rds"))$jobs, jobs)) stop("The resumed module request differs.")
jobs$state <- "pending"
jobs$exit_status <- NA_integer_
active <- list()
repeat {
  for (id in names(active)) {
    process <- active[[id]]
    if (process$is_alive()) next
    i <- as.integer(id)
    jobs$state[i] <- "completed"
    jobs$exit_status[i] <- process$get_exit_status()
    active[[id]] <- NULL
    write.csv(jobs, file.path(output, "job_status.csv"), row.names = FALSE)
  }
  if (all(jobs$state == "completed")) break
  available <- which(jobs$state == "pending")
  for (i in head(available, max(0L, 4L - length(active)))) {
    active[[as.character(i)]] <- processx::process$new(
      file.path(R.home("bin"), "Rscript"),
      c("research/joint_profile_manuscript_20261026/run_module_job.R", as.character(jobs$job_id[i])),
      cleanup = FALSE,
      stdout = file.path(output, sprintf("job_%03d.log", jobs$job_id[i])), stderr = "2>&1")
    jobs$state[i] <- "running"
  }
  Sys.sleep(0.5)
}
write.csv(jobs, file.path(output, "job_status.csv"), row.names = FALSE)
saveRDS(list(completed = nrow(jobs), successful = sum(jobs$exit_status == 0L),
             failed = sum(jobs$exit_status != 0L), finished = Sys.time()), file.path(output, "completed.rds"))
print(table(jobs$exit_status))

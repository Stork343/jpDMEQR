source("research/submission_completion_20261027/load.R")
source("research/submission_completion_20261027/gaussian_curvature.R")
args <- commandArgs(trailingOnly = TRUE)
worker <- as.integer(args[1L])
workers <- as.integer(args[2L])
if (!workers %in% 1:8 || !worker %in% seq_len(workers)) stop("Invalid Gaussian diagnostic worker.")
root <- "results/raw/simulation/submission_20261027_main_joint"
request <- readRDS(file.path(root, "request.rds"))
registry <- request$registry[request$registry$module == "A_theory_scaling", ]
jobs <- request$jobs[request$jobs$experiment_id %in% registry$experiment_id, ]
jobs <- jobs[(seq_len(nrow(jobs)) - 1L) %% workers + 1L == worker, ]
output <- "results/raw/simulation/submission_20261028_gaussian_diagnostics"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
curvature <- submission_gaussian_curvature(.5, 3:8, 1, 1)
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  input <- file.path(root, "replications", job$experiment_id, sprintf("replicate_%05d", job$replicate), "result.rds")
  destination <- file.path(output, job$experiment_id, sprintf("replicate_%05d.rds", job$replicate))
  if (!file.exists(input)) stop("The registered inference calculation is incomplete: ", input)
  saved <- if (file.exists(destination)) readRDS(destination) else NULL
  compatible <- !is.null(saved) && identical(saved$curvature$definition, curvature$definition)
  if (!compatible) {
    result <- readRDS(input)
    diagnostic <- submission_gaussian_diagnostics(result, curvature)
    diagnostic$source_file <- input
    diagnostic$experiment_id <- job$experiment_id
    diagnostic$replicate <- job$replicate
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    saveRDS(diagnostic, destination)
  }
  write.csv(data.frame(completed = i, registered = nrow(jobs), time = as.character(Sys.time())),
    file.path(output, sprintf("worker_%02d_progress.csv", worker)), row.names = FALSE)
}
saveRDS(list(completed = nrow(jobs), workers = workers, finished = Sys.time()),
  file.path(output, sprintf("worker_%02d_completed.rds", worker)))
cat("Completed", nrow(jobs), "registered Gaussian diagnostic datasets\n")

source("research/submission_completion_20261027/load_comparators.R")
library(processx)
args <- commandArgs(trailingOnly = TRUE)
family <- args[1L]
workers <- as.integer(args[2L])
if (!workers %in% 1:8) stop("Invalid comparator worker count.")
methods <- switch(family, mixed = c("DOUBLE-PEN-QLMM", "LQMM", "QIF-SEE"),
  pooled = "POOLED-QR-LASSO", sqr = "SQR-DEBIASED-IID", qgee = "QGEE-SCAD", bias_adjusted = "BIAS-ADJ-LQMM",
  stop("Unknown comparator queue."))
for (method in methods) {
  path <- file.path("results/scratch/submission_completion_20261027/comparator_preflight", paste0(method, ".rds"))
  if (!file.exists(path)) stop("The comparator preflight is missing: ", method)
  result <- readRDS(path)
  if (!all(result$metrics$method_id == method) || any(result$metrics$status == "not_implemented"))
    stop("Invalid comparator preflight identity.")
}
reference <- jsonlite::read_json("results/preflight/submission_reference_20261027/reference_fidelity.json")
if (!length(reference$checks) || !all(vapply(reference$checks, function(x) isTRUE(x$passed), logical(1))))
  stop("The published-method reference validation did not pass.")
checks <- read.csv("results/scratch/submission_completion_20261027/comparator_controls_checks.csv")
if (any(checks$difference > 1e-6)) stop("Comparator numerical control verification failed.")
if (family == "pooled" && !file.exists("results/scratch/submission_completion_20261027/pooled_procedure_check/checks.csv"))
  stop("The complete pooled numerical comparison is missing.")
if (family == "sqr" && !file.exists("results/scratch/submission_comparator_probes_20261027/SQR_A01_parametric_warm/comparison.csv"))
  stop("The complete SQR numerical comparison is missing.")
registry <- submission_registry("main")
registry <- registry[registry$module != "G_computation", ]
registry <- registry[vapply(strsplit(registry$methods, "|", fixed = TRUE), function(x) any(x %in% methods), logical(1)), ]
jobs <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) data.frame(
  experiment_id = registry$experiment_id[i], replicate = seq_len(registry$replications[i]))))
jobs$batch <- (seq_len(nrow(jobs)) - 1L) %% (workers * 4L) + 1L
destination <- file.path("results/raw/simulation", paste0("submission_20261028_main_", family))
if (dir.exists(destination)) stop("The comparator queue already exists.")
dir.create(destination, recursive = TRUE)
saveRDS(list(registry = registry, jobs = jobs, methods = methods, workers = workers, base_seed = 20261027L,
  implementation_commit = current_commit_v2(), session = sessionInfo(), system = Sys.info(),
  numerical_specification = yaml::read_yaml("config/submission_comparator_numerics_20261028.yml"),
  reference_fidelity = reference, numerical_checks = checks, started = Sys.time()), file.path(destination, "request.rds"))
write.csv(registry, file.path(destination, "config_snapshot.csv"), row.names = FALSE)
write.csv(jobs, file.path(destination, "scheduled_replications.csv"), row.names = FALSE)
batches <- data.frame(batch = sort(unique(jobs$batch)), state = "pending", exit_status = NA_integer_)
active <- list()
repeat {
  for (key in names(active)) {
    process <- active[[key]]
    if (process$is_alive()) next
    index <- as.integer(key)
    batches$state[index] <- "completed"
    batches$exit_status[index] <- process$get_exit_status()
    active[[key]] <- NULL
    write.csv(batches, file.path(destination, "batch_status.csv"), row.names = FALSE)
  }
  if (all(batches$state == "completed")) break
  for (index in head(which(batches$state == "pending"), max(0L, workers - length(active)))) {
    active[[as.character(index)]] <- processx::process$new(file.path(R.home("bin"), "Rscript"),
      c("research/submission_completion_20261027/run_comparator_batch.R", destination, as.character(batches$batch[index])),
      stdout = file.path(destination, sprintf("batch_%03d.log", index)), stderr = "2>&1")
    batches$state[index] <- "running"
  }
  Sys.sleep(1)
}
write.csv(batches, file.path(destination, "batch_status.csv"), row.names = FALSE)
if (any(batches$exit_status != 0L)) stop("A comparator batch requires investigation.")
results <- list.files(file.path(destination, "replications"), "^result[.]rds$", recursive = TRUE)
if (length(results) != nrow(jobs)) stop("The comparator replication accounting is incomplete.")
saveRDS(list(completed_replications = length(results), finished = Sys.time()), file.path(destination, "completed.rds"))
cat("Completed", length(results), "registered comparator datasets\n")

source("research/submission_completion_20261027/load.R")
library(processx)
args <- commandArgs(trailingOnly = TRUE)
phase <- args[1L]
family <- args[2L]
workers <- as.integer(args[3L])
if (!phase %in% c("pilot", "main", "preflight") || !workers %in% 1:8) stop("Invalid simulation dispatch.")
registry <- submission_registry(phase)
if (phase == "main") registry <- registry[registry$module != "G_computation", ]
selected <- vapply(strsplit(registry$methods, "|", fixed = TRUE), function(methods) {
  if (family == "joint") any(startsWith(methods, "JOINT-PROFILE")) else
    if (family == "comparators") any(!startsWith(methods, "JOINT-PROFILE")) else family %in% methods
}, logical(1))
registry <- registry[selected, ]
if (!nrow(registry)) stop("No registered simulation rows.")
jobs <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) data.frame(
  experiment_id = registry$experiment_id[i], replicate = seq_len(if (phase == "preflight") 1L else registry$replications[i]))))
jobs$batch <- (seq_len(nrow(jobs)) - 1L) %% (workers * 4L) + 1L
destination <- file.path("results/raw/simulation", paste("submission_20261027", phase, family, sep = "_"))
if (dir.exists(destination)) stop("The simulation dispatch already exists.")
dir.create(destination, recursive = TRUE)
dispatch_lock <- filelock::lock(file.path(destination, "continuation.lock"), timeout = 0)
if (is.null(dispatch_lock)) stop("This simulation queue already has a dispatcher.")
saveRDS(list(registry = registry, jobs = jobs, phase = phase, family = family, workers = workers,
  base_seed = if (phase == "pilot") 202610271L else if (phase == "preflight") 202610272L else 20261027L,
  implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  session = sessionInfo(), system = Sys.info(), started = Sys.time()), file.path(destination, "request.rds"))
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
      c("research/submission_completion_20261027/run_simulation_batch.R", destination, as.character(batches$batch[index])),
      cleanup = FALSE,
      stdout = file.path(destination, sprintf("batch_%03d.log", index)), stderr = "2>&1")
    batches$state[index] <- "running"
  }
  Sys.sleep(1)
}
write.csv(batches, file.path(destination, "batch_status.csv"), row.names = FALSE)
if (any(batches$exit_status != 0L)) stop("A simulation batch requires investigation.")
results <- list.files(file.path(destination, "replications"), "^result[.]rds$", recursive = TRUE)
if (length(results) != nrow(jobs)) stop("The simulation replication accounting is incomplete.")
saveRDS(list(completed_replications = length(results), finished = Sys.time()), file.path(destination, "completed.rds"))
cat("Completed", length(results), "registered replications\n")

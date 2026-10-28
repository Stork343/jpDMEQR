options(warn = 1, width = 160)
source("research/score_validity_20261026/load_joint.R")
args <- commandArgs(trailingOnly = TRUE)
worker <- as.integer(args[1L])
workers <- as.integer(args[2L])
run_id <- if (length(args) >= 3L) args[3L] else "design_independent"
seed_base <- if (length(args) >= 4L) as.integer(args[4L]) else 1800000000L
resume <- length(args) >= 5L && identical(args[5L], "resume")
if (!is.finite(worker) || worker < 1L || worker > workers || !is.finite(seed_base)) stop("Invalid evaluation request.")
registry_path <- if (length(args) >= 6L) args[6L] else "research/score_validity_20261026/evaluation.csv"
registry <- read.csv(registry_path, stringsAsFactors = FALSE)
output <- file.path("results/scratch/score_validity_20261026", run_id)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
request_file <- file.path(output, paste0("worker_request_", worker, ".rds"))
if (file.exists(request_file)) {
  if (!resume) stop("This worker request already exists.")
  previous <- readRDS(request_file)
  if (!identical(previous$registry, registry) || previous$seed_base != seed_base || previous$workers != workers) stop("Different resumed design.")
  request_file <- file.path(output, paste0("worker_resume_", worker, "_", format(Sys.time(), "%Y%m%dT%H%M%S"), ".rds"))
}
saveRDS(list(registry = registry, worker = worker, workers = workers, seed_base = seed_base,
             method_version = "joint_profile_design_20261026", started = Sys.time(), session = sessionInfo()), request_file)
current_directory <- NULL
current_request <- NULL
globalCallingHandlers(error = function(condition) {
  if (!is.null(current_directory)) saveRDS(list(request = current_request, condition = condition),
                                          file.path(current_directory, "failure.rds"))
})
completed <- 0L
for (case_index in seq_len(nrow(registry))) {
  cfg <- registry[case_index, ]
  for (replicate in seq_len(cfg$replications)) {
    if ((replicate - 1L) %% workers != worker - 1L) next
    seed <- seed_base + 10000L * case_index + replicate
    current_request <- list(config = cfg, replicate = replicate, seed = seed)
    current_directory <- file.path(output, cfg$case_id, sprintf("replicate_%03d", replicate))
    if (dir.exists(current_directory)) {
      if (!resume) stop("An evaluation unit already exists.")
      completed_files <- list.files(current_directory, pattern = "^result[.]rds$", recursive = TRUE)
      if (length(completed_files)) {
        if (length(completed_files) != 1L) stop("Duplicate completed evaluation unit.")
        completed <- completed + 1L
        next
      }
      current_directory <- file.path(current_directory, sprintf("attempt_%02d", 2L + length(list.dirs(current_directory, recursive = FALSE))))
    }
    dir.create(current_directory, recursive = TRUE)
    saveRDS(current_request, file.path(current_directory, "request.rds"))
    dat <- generate_profile_qr_data_v2(n_clusters = cfg$n_clusters, p = cfg$p, s = cfg$s,
                                       tau = cfg$tau, q = cfg$q, error_dist = cfg$error_dist,
                                       signal = cfg$signal, rho_x = cfg$rho_x, seed = seed)
    coordinates <- c(1L, 2L, 6L, 7L)
    fit <- ds_fit_profile(dat$y, dat$X, dat$Z, dat$cluster_id, cfg$tau, coordinates)
    table <- fit$table
    table$case_id <- cfg$case_id
    table$replicate <- replicate
    table$truth <- dat$beta0[coordinates]
    table$covered <- table$lower <= table$truth & table$upper >= table$truth
    table$reject_zero <- table$lower > 0 | table$upper < 0
    table$z_truth <- (table$beta_tilde - table$truth) / table$se
    fit$request <- current_request
    fit$table <- table
    saveRDS(fit, file.path(current_directory, "result.rds"))
    write.csv(table, file.path(current_directory, "coordinates.csv"), row.names = FALSE)
    completed <- completed + 1L
  }
  cat("Worker", worker, "completed", cfg$case_id, "datasets", completed, "\n")
}
saveRDS(list(worker = worker, completed = completed, finished = Sys.time()),
         file.path(output, paste0("worker_complete_", worker, ".rds")))

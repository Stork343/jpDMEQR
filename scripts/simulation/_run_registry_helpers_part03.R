run_registry_v2 <- function(root, config_path, run_id, jobs = 1L,
                            max_reps = Inf, allow_missing_benchmarks = TRUE,
                            base_seed = 20260817L) {
  cfg_df <- read_simulation_registry_v2(config_path)
  tasks <- list()
  tcount <- 0L
  for (i in seq_len(nrow(cfg_df))) {
    cfg <- config_row_to_list_v2(cfg_df[i, , drop = FALSE])
    B <- if (is.finite(max_reps)) min(as.integer(cfg$replications), as.integer(max_reps)) else as.integer(cfg$replications)
    for (r in seq_len(B)) {
      tcount <- tcount + 1L
      tasks[[tcount]] <- list(config = cfg, replicate = r)
    }
  }

  runner <- function(task) {
    tryCatch(
      run_one_replication_v2(root, task$config, task$replicate, base_seed,
                             allow_missing_benchmarks),
      error = function(e) list(error = conditionMessage(e), config = task$config,
                               replicate = task$replicate)
    )
  }
  if (jobs > 1L && .Platform$OS.type != "windows") {
    output <- parallel::mclapply(tasks, runner, mc.cores = jobs, mc.preschedule = FALSE)
  } else {
    output <- lapply(tasks, runner)
  }

  run_dir <- file.path(root, "results", "raw", "simulation", run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(config_path, file.path(run_dir, "config_snapshot.csv"), overwrite = TRUE)
  writeLines(current_commit_v2(root), file.path(run_dir, "implementation_commit.txt"))
  capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))
  saveRDS(output, file.path(run_dir, "raw_objects.rds"))

  good <- output[!vapply(output, function(x) !is.null(x$error), logical(1))]
  errors <- output[vapply(output, function(x) !is.null(x$error), logical(1))]
  metrics <- if (length(good)) do.call(rbind, lapply(good, `[[`, "metrics")) else data.frame()
  coords <- if (length(good)) do.call(rbind, lapply(good, function(x) x$coordinates)) else data.frame()
  if (nrow(metrics)) write.csv(metrics, file.path(run_dir, "replication_metrics.csv"), row.names = FALSE)
  if (nrow(coords)) write.csv(coords, file.path(run_dir, "coordinate_metrics.csv"), row.names = FALSE)
  if (length(errors)) {
    err_df <- do.call(rbind, lapply(errors, function(x) data.frame(
      experiment_id = x$config$experiment_id,
      replicate = x$replicate,
      error = x$error,
      stringsAsFactors = FALSE
    )))
    write.csv(err_df, file.path(run_dir, "runner_failures.csv"), row.names = FALSE)
  }
  invisible(list(run_dir = run_dir, metrics = metrics, coordinates = coords, errors = errors))
}

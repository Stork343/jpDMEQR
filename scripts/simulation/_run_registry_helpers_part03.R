# Registry execution, deterministic parallelisation, and immutable raw outputs.

screening_record_v2 <- function(x) {
  s <- x$screening
  data.frame(
    experiment_id = x$config$experiment_id,
    replicate = x$config$replicate %||% NA_integer_,
    seed = x$seed,
    screen_method = s$method,
    screen_engine = s$engine,
    screen_cluster_count = length(s$screen_clusters),
    fit_cluster_count = length(s$fit_clusters),
    overlap_count = s$overlap_count,
    requested_dimension = s$screen_dimension_requested,
    realised_dimension = s$screen_dimension_realised,
    sure_inclusion = s$sure_inclusion,
    forced_target_count = length(s$forced_target_names),
    selected_names = paste(s$selected_names, collapse = ";"),
    forced_target_names = paste(s$forced_target_names, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

sort_result_table_v2 <- function(x, columns) {
  if (!nrow(x)) return(x)
  columns <- intersect(columns, names(x))
  if (!length(columns)) return(x)
  ord <- do.call(order, c(x[columns], list(na.last = TRUE)))
  x[ord, , drop = FALSE]
}

run_registry_v2 <- function(root,
                            config_path,
                            run_id,
                            jobs = 1L,
                            max_reps = Inf,
                            allow_missing_benchmarks = TRUE,
                            base_seed = 20260817L,
                            experiment_ids = character(),
                            final = FALSE) {
  cfg_df <- read_simulation_registry_v2(config_path)
  if (length(experiment_ids)) {
    cfg_df <- cfg_df[cfg_df$experiment_id %in% experiment_ids, , drop = FALSE]
    if (!nrow(cfg_df)) stop("No requested experiment_id was found in the registry.")
  }
  jobs <- max(1L, as.integer(jobs))
  tasks <- list()
  tcount <- 0L
  for (i in seq_len(nrow(cfg_df))) {
    cfg <- config_row_to_list_v2(cfg_df[i, , drop = FALSE])
    B <- if (is.finite(max_reps)) {
      min(as.integer(cfg$replications), as.integer(max_reps))
    } else as.integer(cfg$replications)
    for (r in seq_len(B)) {
      tcount <- tcount + 1L
      cfg$replicate <- r
      tasks[[tcount]] <- list(config = cfg, replicate = r)
    }
  }
  if (!length(tasks)) stop("The selected registry contains no tasks.")

  run_dir <- file.path(root, "results", "raw", "simulation", run_id)
  if (dir.exists(run_dir)) stop("Run directory already exists and is immutable: ", run_dir)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  lock_path <- file.path(run_dir, "RUNNING.lock")
  writeLines(c(
    paste0("pid=", Sys.getpid()),
    paste0("started_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("commit=", current_commit_v2(root))
  ), lock_path)
  completed <- FALSE
  on.exit({
    if (file.exists(lock_path)) unlink(lock_path)
    if (!completed) {
      writeLines(
        paste0("Run did not complete: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
        file.path(run_dir, "INCOMPLETE.txt")
      )
    }
  }, add = TRUE)

  file.copy(config_path, file.path(run_dir, "config_snapshot.csv"), overwrite = FALSE)
  writeLines(current_commit_v2(root), file.path(run_dir, "implementation_commit.txt"))
  writeLines(current_branch_v2(root), file.path(run_dir, "implementation_branch.txt"))
  writeLines(file_checksum_v2(config_path), file.path(run_dir, "config_sha256.txt"))
  write_json_v2(hardware_manifest_v2(), file.path(run_dir, "hardware.json"))
  capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))
  write_json_v2(list(
    run_id = run_id,
    task_count = length(tasks),
    jobs = jobs,
    max_reps = if (is.finite(max_reps)) max_reps else "Inf",
    allow_missing_benchmarks = allow_missing_benchmarks,
    final = final,
    base_seed = base_seed,
    experiment_ids = cfg_df$experiment_id,
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ), file.path(run_dir, "run_request.json"))

  runner <- function(task) {
    tryCatch(
      run_one_replication_v2(
        root, task$config, task$replicate,
        base_seed = base_seed,
        allow_missing_benchmarks = allow_missing_benchmarks,
        final = final
      ),
      error = function(e) list(
        error = conditionMessage(e),
        config = task$config,
        replicate = task$replicate,
        seed = seed_from_id_v2(base_seed, task$config$experiment_id,
                               task$replicate, "training")
      )
    )
  }

  if (jobs > 1L && .Platform$OS.type != "windows") {
    output <- parallel::mclapply(
      tasks, runner, mc.cores = jobs,
      mc.preschedule = FALSE, mc.set.seed = FALSE
    )
  } else {
    output <- lapply(tasks, runner)
  }

  saveRDS(output, file.path(run_dir, "raw_objects.rds"), compress = "xz")
  is_error <- vapply(output, function(x) !is.null(x$error), logical(1))
  good <- output[!is_error]
  errors <- output[is_error]
  metrics <- if (length(good)) do.call(rbind, lapply(good, `[[`, "metrics")) else data.frame()
  coords <- if (length(good)) do.call(rbind, lapply(good, function(x) x$coordinates)) else data.frame()
  theory <- if (length(good)) do.call(rbind, lapply(good, function(x) x$theory)) else data.frame()
  screening <- if (length(good)) {
    do.call(rbind, lapply(good, screening_record_v2))
  } else data.frame()

  metrics <- sort_result_table_v2(metrics, c("experiment_id", "replicate", "method_id"))
  coords <- sort_result_table_v2(coords, c("experiment_id", "replicate", "method_id", "coordinate"))
  theory <- sort_result_table_v2(theory, c("experiment_id", "replicate", "method_id"))
  screening <- sort_result_table_v2(screening, c("experiment_id", "replicate"))

  if (nrow(metrics)) write_atomic_csv_v2(metrics, file.path(run_dir, "replication_metrics.csv"))
  if (nrow(coords)) write_atomic_csv_v2(coords, file.path(run_dir, "coordinate_metrics.csv"))
  if (nrow(theory)) write_atomic_csv_v2(theory, file.path(run_dir, "theory_diagnostics.csv"))
  if (nrow(screening)) write_atomic_csv_v2(screening, file.path(run_dir, "screening_records.csv"))

  if (length(errors)) {
    err_df <- do.call(rbind, lapply(errors, function(x) data.frame(
      experiment_id = x$config$experiment_id,
      replicate = x$replicate,
      seed = x$seed,
      error = x$error,
      stringsAsFactors = FALSE
    )))
    err_df <- sort_result_table_v2(err_df, c("experiment_id", "replicate"))
    write_atomic_csv_v2(err_df, file.path(run_dir, "runner_failures.csv"))
  }

  manifest <- list(
    run_id = run_id,
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    implementation_commit = current_commit_v2(root),
    implementation_branch = current_branch_v2(root),
    config_sha256 = file_checksum_v2(config_path),
    task_count = length(tasks),
    successful_task_objects = length(good),
    runner_failure_count = length(errors),
    method_rows = nrow(metrics),
    coordinate_rows = nrow(coords),
    theory_rows = nrow(theory),
    jobs = jobs,
    final = final,
    raw_objects_sha256 = file_checksum_v2(file.path(run_dir, "raw_objects.rds"))
  )
  write_json_v2(manifest, file.path(run_dir, "completed_manifest.json"))
  completed <- TRUE
  if (file.exists(file.path(run_dir, "INCOMPLETE.txt"))) {
    unlink(file.path(run_dir, "INCOMPLETE.txt"))
  }
  invisible(list(
    run_dir = run_dir,
    metrics = metrics,
    coordinates = coords,
    theory = theory,
    screening = screening,
    errors = errors,
    manifest = manifest
  ))
}

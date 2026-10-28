# Registry execution, deterministic parallelisation, and immutable raw outputs.

screening_record_v2 <- function(x) {
  if (!is.null(x$screening_table)) return(x$screening_table)
  s <- x$screening
  data.frame(
    run_id = x$run_id %||% NA_character_,
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
                            final = FALSE,
                            replicate_subset = NULL,
                            seed_experiment = NULL,
                            resume = FALSE,
                            storage_mode = NULL) {
  storage <- registry_storage_configuration_v2(storage_mode)
  cfg_df <- read_simulation_registry_v2(config_path)
  phase <- attr(cfg_df, "registry_phase") %||% registry_type_for_path_v2(config_path)
  if (length(experiment_ids)) {
    cfg_df <- cfg_df[cfg_df$experiment_id %in% experiment_ids, , drop = FALSE]
    if (!nrow(cfg_df)) stop("No requested experiment_id was found in the registry.")
  }
  if (isTRUE(final) && "method_version" %in% names(cfg_df) &&
      any(cfg_df$method_version == "profile_pair_20261022") && identical(phase, "main")) {
    validate_profile_pair_final_gate_v2(root, config_path)
  }
  jobs <- max(1L, as.integer(jobs))
  tasks <- list()
  tcount <- 0L
  for (i in seq_len(nrow(cfg_df))) {
    cfg <- config_row_to_list_v2(cfg_df[i, , drop = FALSE])
    B <- if (is.finite(max_reps)) {
      min(as.integer(cfg$replications), as.integer(max_reps))
    } else as.integer(cfg$replications)
    reps <- seq_len(B)
    if (!is.null(replicate_subset)) {
      reps <- reps[reps %in% as.integer(replicate_subset)]
    }
    for (r in reps) {
      tcount <- tcount + 1L
      cfg$replicate <- r
      tasks[[tcount]] <- list(config = cfg, replicate = r)
    }
  }
  if (!length(tasks)) stop("The selected registry contains no tasks.")

  run_dir <- file.path(root, "results", "raw", "simulation", run_id)
  objects_dir <- file.path(run_dir, "replication_objects")
  lock_path <- file.path(run_dir, "RUNNING.lock")
  task_key <- function(cfg, replicate) {
    sprintf("%s_%05d", cfg$experiment_id, as.integer(replicate))
  }
  completed_results <- list()
  if (dir.exists(run_dir)) {
    if (!isTRUE(resume)) stop("Run directory already exists and is immutable: ", run_dir)
    if (file.exists(file.path(run_dir, "completed_manifest.json"))) {
      stop("Run directory is already complete; start a new run_id: ", run_dir)
    }
    request_path <- file.path(run_dir, "run_request.json")
    if (!file.exists(request_path)) {
      stop("Cannot resume without run_request.json in ", run_dir)
    }
    previous <- read_json_v2(request_path)
    previous_storage <- previous$storage$mode %||% "memory"
    if (!identical(previous_storage, storage$mode)) stop("Cannot resume with a different result-storage mode.")
    if (storage$mode == "streaming" &&
        previous$storage$matrix_inline_entries != storage$matrix_inline_entries) {
      stop("Cannot resume with a different diagnostic-storage configuration.")
    }
    if (!identical(isTRUE(previous$final), isTRUE(final))) {
      stop("Cannot resume: the recorded final flag differs from this request.")
    }
    if (!identical(isTRUE(previous$allow_missing_benchmarks),
                   isTRUE(allow_missing_benchmarks))) {
      stop("Cannot resume: the recorded benchmark policy differs from this request.")
    }
    if (!identical(as.numeric(previous$base_seed %||% NA_real_), as.numeric(base_seed))) {
      stop("Cannot resume: the recorded base seed differs from this request.")
    }
    recorded_max <- previous$max_reps %||% NA
    if (!identical(as.character(recorded_max), as.character(max_reps))) {
      stop("Cannot resume: the recorded replication cap differs from this request.")
    }
    if (!setequal(as.character(previous$experiment_ids %||% character()),
                  cfg_df$experiment_id)) {
      stop("Cannot resume: the recorded experiment set differs from this request.")
    }
    recorded_sha <- trimws(readLines(file.path(run_dir, "config_sha256.txt"),
                                     warn = FALSE)[1] %||% "")
    if (!identical(recorded_sha, file_checksum_v2(config_path))) {
      stop("Cannot resume: the registry checksum differs from the recorded run.")
    }
    recorded_commit <- trimws(readLines(file.path(run_dir, "implementation_commit.txt"),
                                        warn = FALSE)[1] %||% "")
    if (!identical(recorded_commit, current_commit_v2(root))) {
      stop("Cannot resume: the implementation commit differs from the recorded run.")
    }
    if (file.exists(lock_path)) {
      lock_lines <- readLines(lock_path, warn = FALSE)
      pid_line <- grep("^pid=", lock_lines, value = TRUE)
      lock_pid <- if (length(pid_line)) {
        suppressWarnings(as.integer(sub("^pid=", "", pid_line[1])))
      } else NA_integer_
      if (isTRUE(pid_running_v2(lock_pid))) {
        stop("Cannot resume: shard lock is held by a running process (pid ",
             lock_pid, "). Remove ", lock_path, " only if that process is gone.")
      }
      message("Resuming shard with a stale lock from pid ", lock_pid, ".")
    }
    dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
    object_files <- sort(list.files(objects_dir, pattern = "\\.rds$", full.names = TRUE))
    completed_results <- lapply(object_files, function(path) {
      value <- readRDS(path)
      if (storage$mode == "streaming") registry_table_record_v2(value) else value
    })
    if (length(completed_results)) {
      schema <- list(
        metrics = contract_schema_columns_v2(root, 2L),
        coordinates = contract_schema_columns_v2(root, 3L),
        theory = contract_schema_columns_v2(root, 4L),
        screening = contract_schema_columns_v2(root, 5L)
      )
      observed <- lapply(completed_results[vapply(completed_results, function(x) is.null(x$error), logical(1))], function(x) {
        list(metrics = names(x$metrics),
             coordinates = names(x$coordinates),
             theory = names(x$theory),
             screening = names(screening_record_v2(x)))
      })
      for (tb in if (length(observed)) names(schema) else character()) {
        columns <- lapply(observed, function(o) o[[tb]])
        reference <- columns[[1]]
        for (k in seq_along(columns)) {
          if (!identical(columns[[k]], reference)) {
            stop("Cannot resume: completed objects use different column sets in ", tb, ".")
          }
        }
        missing <- setdiff(schema[[tb]], reference)
        if (length(missing)) {
          stop("Cannot resume: completed objects predate the current record layer in ", tb,
               ", missing columns: ", paste(missing, collapse = ", "))
        }
      }
    }
  } else {
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
  }
  done_keys <- sub("\\.rds$", "", basename(list.files(objects_dir, pattern = "\\.rds$")))
  tasks <- tasks[!vapply(tasks, function(task) {
    task_key(task$config, task$replicate) %in% done_keys
  }, logical(1))]
  writeLines(c(
    paste0("pid=", Sys.getpid()),
    paste0("started_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("commit=", current_commit_v2(root)),
    paste0("resumed_objects=", length(completed_results))
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

  utils::write.csv(cfg_df, file.path(run_dir, "config_snapshot.csv"), row.names = FALSE)
  if (tolower(tools::file_ext(config_path)) %in% c("yml", "yaml")) {
    file.copy(config_path, file.path(run_dir, "config_specification.yml"), overwrite = FALSE)
    source_registry <- attr(cfg_df, "registry_source_file")
    if (!is.null(source_registry)) {
      file.copy(source_registry, file.path(run_dir, "config_source_registry.csv"), overwrite = FALSE)
    }
  }
  writeLines(current_commit_v2(root), file.path(run_dir, "implementation_commit.txt"))
  writeLines(current_branch_v2(root), file.path(run_dir, "implementation_branch.txt"))
  writeLines(file_checksum_v2(config_path), file.path(run_dir, "config_sha256.txt"))
  write_json_v2(hardware_manifest_v2(), file.path(run_dir, "hardware.json"))
  write_json_v2(solver_configuration_v2(), file.path(run_dir, "solver_configuration.json"))
  capture.output(sessionInfo(), file = file.path(run_dir, "sessionInfo.txt"))
  write_json_v2(list(
    run_id = run_id,
    task_count = length(tasks),
    jobs = jobs,
    max_reps = if (is.finite(max_reps)) max_reps else "Inf",
    allow_missing_benchmarks = allow_missing_benchmarks,
    final = final,
    phase = phase,
    base_seed = base_seed,
    seed_experiment = seed_experiment %||% "default",
    storage = storage,
    experiment_ids = cfg_df$experiment_id,
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ), file.path(run_dir, "run_request.json"))

  runner <- function(task) {
    result <- tryCatch(
      run_one_replication_v2(
        root, task$config, task$replicate,
        base_seed = base_seed,
        allow_missing_benchmarks = allow_missing_benchmarks,
        final = final,
        seed_experiment = seed_experiment,
        run_id = run_id
      ),
      error = function(e) list(
        error = conditionMessage(e),
        error_class = class(e), failure_stage = e$failure_stage %||% "unknown",
        failure_diagnostics = e$diagnostics, partial_fit = e$fit %||% e$partial_fit,
        condition = e,
        config = task$config,
        replicate = task$replicate,
        seed = seed_from_id_v2(base_seed, task$config$experiment_id,
                               task$replicate, "training")
      )
    )
    if (storage$mode == "streaming") {
      result <- compact_registry_replication_v2(result, storage, root)
      key <- task_key(task$config, task$replicate)
      saveRDS(result, file.path(objects_dir, paste0(key, ".rds")), compress = "xz")
      return(registry_table_record_v2(result))
    }
    result
  }

  # 逐重复写出对象文件：每个重复完成即写入 replication_objects/，
  # 长跑中断后可以按已完成的对象恢复，最终表格仍按确定性顺序重建。
  chunk_size <- max(1L, jobs)
  new_results <- list()
  if (length(tasks)) {
    for (start in seq(1L, length(tasks), by = chunk_size)) {
      idx <- start:min(start + chunk_size - 1L, length(tasks))
      chunk <- tasks[idx]
      chunk_out <- if (jobs > 1L && .Platform$OS.type != "windows") {
        parallel::mclapply(
          chunk, runner, mc.cores = jobs,
          mc.preschedule = FALSE, mc.set.seed = FALSE
        )
      } else {
        lapply(chunk, runner)
      }
      for (k in seq_along(chunk)) {
        task <- chunk[[k]]
        key <- task_key(task$config, task$replicate)
        if (storage$mode == "memory") {
          saveRDS(chunk_out[[k]], file.path(objects_dir, paste0(key, ".rds")), compress = "xz")
        }
        new_results[[length(new_results) + 1L]] <- chunk_out[[k]]
        message(sprintf("recorded %s (%d of %d new tasks)", key,
                        length(new_results), length(tasks)))
      }
    }
  }
  output <- c(completed_results, new_results)

  raw_objects <- if (storage$mode == "streaming") registry_replication_index_v2(run_dir, run_id) else output
  saveRDS(raw_objects, file.path(run_dir, "raw_objects.rds"), compress = "xz")
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
      failure_stage = x$failure_stage %||% "unknown",
      error_class = paste(x$error_class %||% character(), collapse = ";"),
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
    task_count = length(tasks) + length(completed_results),
    resumed_task_objects = length(completed_results),
    new_task_objects = length(new_results),
    successful_task_objects = length(good),
    runner_failure_count = length(errors),
    method_rows = nrow(metrics),
    coordinate_rows = nrow(coords),
    theory_rows = nrow(theory),
    jobs = jobs,
    final = final,
    storage = storage,
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

registry_storage_configuration_v2 <- function(storage_mode = NULL) {
  mode <- storage_mode %||% Sys.getenv("JPDMEQR_RESULT_STORAGE", unset = "memory")
  if (!mode %in% c("memory", "streaming")) stop("Unknown registry result-storage mode.")
  entries <- suppressWarnings(as.integer(Sys.getenv("JPDMEQR_DIAGNOSTIC_MATRIX_ENTRIES", unset = "4096")))
  if (length(entries) != 1L || !is.finite(entries) || entries < 1L) stop("Invalid diagnostic matrix storage limit.")
  list(mode = mode, diagnostics = if (mode == "streaming") "compact_diagnostics_v1" else "full",
       matrix_inline_entries = entries,
       raw_objects_format = if (mode == "streaming") "replication_file_index_v1" else "replication_object_list")
}

registry_numeric_summary_v2 <- function(x) {
  values <- as.numeric(x)
  finite <- values[is.finite(values)]
  scale <- if (length(finite)) max(abs(finite)) else NA_real_
  list(length = length(values), finite = length(finite), missing = sum(is.na(values)),
       positive_infinite = sum(values == Inf, na.rm = TRUE),
       negative_infinite = sum(values == -Inf, na.rm = TRUE),
       finite_minimum = if (length(finite)) min(finite) else NA_real_,
       finite_maximum = if (length(finite)) max(finite) else NA_real_,
       finite_mean = if (length(finite)) mean(finite) else NA_real_,
       finite_sd = if (length(finite) > 1L) stats::sd(finite) else NA_real_,
       finite_l2 = if (!length(finite)) NA_real_ else if (scale == 0) 0 else
         scale * sqrt(sum((finite / scale)^2)))
}

compact_registry_diagnostics_v2 <- function(x, matrix_inline_entries = 4096L, path = character()) {
  field <- if (length(path)) tail(path, 1L) else ""
  if (is.environment(x)) return(list(storage = "environment_reference", class = class(x), names = ls(x)))
  if (is.function(x)) return(list(storage = "function_reference", formals = formals(x), body = deparse(body(x))))
  if (inherits(x, "formula")) return(list(storage = "formula_reference", formula = deparse(x)))
  common <- field %in% c("X", "Z", "y", "cluster_id", "X_i", "Z_i", "X_tilde", "Xt")
  matrix <- is.matrix(x) || is.array(x)
  direction <- any(grepl("omega|direction", path, ignore.case = TRUE)) && matrix &&
    length(dim(x)) == 2L && min(dim(x)) <= 64L
  if ((common && (is.atomic(x) || matrix)) ||
      (matrix && length(x) > matrix_inline_entries && !direction)) {
    summary <- list(storage = "reconstructible_numeric_summary_v1", class = class(x),
                    dimensions = dim(x), names = if (is.null(dim(x))) names(x) else dimnames(x),
                    source_field = paste(path, collapse = "$"))
    if (is.numeric(x) || is.logical(x)) {
      summary$values <- registry_numeric_summary_v2(x)
      if (is.matrix(x)) {
        summary$column_means <- colMeans(x)
        if (nrow(x) == ncol(x)) summary$diagonal <- diag(x)
      }
    } else {
      summary$length <- length(x)
      summary$unique_values <- unique(as.character(x))
    }
    return(summary)
  }
  if (is.data.frame(x)) return(x)
  if (inherits(x, "condition")) {
    fields <- as.list(x)
    result <- lapply(seq_along(fields), function(i) {
      compact_registry_diagnostics_v2(fields[[i]], matrix_inline_entries, c(path, names(fields)[i]))
    })
    names(result) <- names(fields)
    return(list(storage = "condition_record_v1", class = class(x), message = conditionMessage(x), fields = result))
  }
  if (is.list(x)) {
    result <- lapply(seq_along(x), function(i) {
      name <- if (!is.null(names(x)) && nzchar(names(x)[i])) names(x)[i] else as.character(i)
      compact_registry_diagnostics_v2(x[[i]], matrix_inline_entries, c(path, name))
    })
    attributes(result) <- attributes(x)
    return(result)
  }
  x
}

compact_registry_replication_v2 <- function(x, configuration, root) {
  result <- x
  commits <- unique(x$metrics$implementation_commit)
  source_commit <- if (length(commits) == 1L && !is.na(commits)) commits else current_commit_v2(root)
  for (field in c("method_diagnostics", "failure_diagnostics", "partial_fit", "condition")) {
    if (!is.null(x[[field]])) result[[field]] <- compact_registry_diagnostics_v2(
      x[[field]], configuration$matrix_inline_entries, field
    )
  }
  result$storage <- list(format = "compact_diagnostics_v1", source_seed = x$seed,
                         source_config_field = "config", implementation_commit = source_commit,
                         numerical_configuration = solver_configuration_v2(),
                         matrix_inline_entries = configuration$matrix_inline_entries,
                         reconstruction = "run_one_replication_v2_with_recorded_config_seed_and_source")
  result
}

registry_table_record_v2 <- function(x) {
  if (!is.null(x$error)) return(list(error = x$error, error_class = x$error_class,
                                    failure_stage = x$failure_stage, config = x$config,
                                    replicate = x$replicate, seed = x$seed))
  list(metrics = x$metrics, coordinates = x$coordinates, theory = x$theory,
       screening_table = screening_record_v2(x),
       config = list(experiment_id = x$config$experiment_id, replicate = x$config$replicate),
       seed = x$seed, run_id = x$run_id)
}

registry_replication_index_v2 <- function(run_dir, run_id) {
  files <- sort(list.files(file.path(run_dir, "replication_objects"), pattern = "\\.rds$", full.names = FALSE))
  keys <- sub("\\.rds$", "", files)
  list(format = "replication_file_index_v1", run_id = run_id,
       base_directory = normalizePath(run_dir),
       files = data.frame(key = keys, relative_path = file.path("replication_objects", files),
                          stringsAsFactors = FALSE))
}

read_registry_replication_v2 <- function(raw_objects, key, base_directory = NULL) {
  index_path <- NULL
  if (is.character(raw_objects) && length(raw_objects) == 1L) {
    index_path <- normalizePath(raw_objects, mustWork = TRUE)
    raw_objects <- readRDS(index_path)
  }
  if (identical(raw_objects$format, "replication_file_index_v1")) {
    matches <- which(raw_objects$files$key == key)
    if (length(matches) != 1L) stop("The requested replication key is absent or duplicated.")
    directory <- base_directory %||% if (!is.null(index_path)) dirname(index_path) else raw_objects$base_directory
    return(readRDS(file.path(directory, raw_objects$files$relative_path[matches])))
  }
  if (!is.list(raw_objects)) stop("Unknown raw replication storage format.")
  keys <- vapply(raw_objects, function(x) sprintf("%s_%05d", x$config$experiment_id,
                                                 x$config$replicate %||% x$replicate), character(1))
  matches <- which(keys == key)
  if (length(matches) != 1L) stop("The requested replication key is absent or duplicated.")
  raw_objects[[matches]]
}

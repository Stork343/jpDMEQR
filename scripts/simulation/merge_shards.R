#!/usr/bin/env Rscript

# 合并分片运行目录。每个分片目录必须带有 completed_manifest.json，
# 合并结果按实验、重复、方法与坐标排序，并核查覆盖范围与重复键。

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/merge_shards.R"
root <- if (file.exists("R/metrics_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "simulation_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
run_id <- cli$`run-id` %||% stop("Supply --run-id=<shard prefix>.")
n_shards <- as.integer(cli$`n-shards` %||% stop("Supply --n-shards=<count>."))
require_complete <- as_bool_cli_v2(cli$`require-complete`, TRUE)
merged_id <- cli$`merged-id` %||% run_id
raw_dir <- file.path(root, "results", "raw", "simulation")
shard_dirs <- file.path(raw_dir, paste0(run_id, "_s", seq_len(n_shards)))
missing_dirs <- shard_dirs[!dir.exists(shard_dirs)]
if (length(missing_dirs)) stop("Missing shard directories: ", paste(missing_dirs, collapse = ", "))
if (isTRUE(require_complete)) {
  incomplete <- shard_dirs[!file.exists(file.path(shard_dirs, "completed_manifest.json"))]
  if (length(incomplete)) {
    stop("Shards without a completed manifest: ", paste(basename(incomplete), collapse = ", "))
  }
}
merged <- file.path(raw_dir, merged_id)
if (dir.exists(merged)) stop("Merged run directory already exists: ", merged)

csvs <- c("replication_metrics.csv", "coordinate_metrics.csv", "theory_diagnostics.csv",
          "screening_records.csv", "runner_failures.csv")
key_of <- function(tab) {
  columns <- intersect(c("experiment_id", "replicate", "method_id", "coordinate"), names(tab))
  if (!length(columns)) stop("Merged table has no identifying columns.")
  do.call(paste, c(tab[columns], list(sep = "|")))
}
merged_tables <- list()
for (cc in csvs) {
  parts <- lapply(shard_dirs, function(d) {
    path <- file.path(d, cc)
    if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) next
  tab <- do.call(rbind, parts)
  if (nrow(tab)) {
    keys <- key_of(tab)
    if (anyDuplicated(keys)) {
      stop("Duplicated rows across shards in ", cc, ": ",
           sum(duplicated(keys)), " duplicates.")
    }
    ord <- do.call(order, c(tab[intersect(c("experiment_id", "replicate", "method_id",
                                            "coordinate"), names(tab))],
                            list(na.last = TRUE)))
    tab <- tab[ord, , drop = FALSE]
  }
  merged_tables[[cc]] <- tab
}
if (!length(merged_tables)) stop("No shard table was found to merge.")

dir.create(merged, recursive = TRUE, showWarnings = FALSE)
for (cc in names(merged_tables)) {
  write.csv(merged_tables[[cc]], file.path(merged, cc), row.names = FALSE)
  cat(cc, ": rows =", nrow(merged_tables[[cc]]), "\n")
}

# 共享的初始化文件从第一个分片复制；运行请求重新写出，因为各分片的
# 实验集合不同，直接复制某一个分片的请求会给出错误的来源信息。
shared <- c("config_snapshot.csv", "config_sha256.txt", "implementation_commit.txt",
            "implementation_branch.txt", "sessionInfo.txt",
            "hardware.json", "solver_configuration.json")
for (name in shared) {
  source_path <- file.path(shard_dirs[1], name)
  if (file.exists(source_path)) file.copy(source_path, file.path(merged, name),
                                          overwrite = TRUE)
}
first_request_path <- file.path(shard_dirs[1], "run_request.json")
first_request <- if (file.exists(first_request_path)) {
  read_json_v2(first_request_path)
} else {
  stop("Missing run_request.json in ", shard_dirs[1])
}
registry_snapshot <- read.csv(file.path(merged, "config_snapshot.csv"),
                              stringsAsFactors = FALSE)
write_json_v2(list(
  run_id = merged_id,
  merged = TRUE,
  shard_count = n_shards,
  shard_ids = paste0(run_id, "_s", seq_len(n_shards)),
  task_count = sum(vapply(shard_dirs, function(d) {
    manifest_path <- file.path(d, "completed_manifest.json")
    if (!file.exists(manifest_path)) return(0L)
    as.integer(read_json_v2(manifest_path)$task_count %||% 0L)
  }, integer(1))),
  jobs = as.integer(first_request$jobs %||% 1L),
  max_reps = first_request$max_reps %||% "Inf",
  allow_missing_benchmarks = isTRUE(first_request$allow_missing_benchmarks),
  final = isTRUE(first_request$final),
  base_seed = as.numeric(first_request$base_seed %||% NA_real_),
  experiment_ids = registry_snapshot$experiment_id,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(merged, "run_request.json"))

registry_path <- file.path(merged, "config_snapshot.csv")
if (!file.exists(registry_path)) stop("Missing ", registry_path)
registry <- read.csv(registry_path, stringsAsFactors = FALSE)
expected <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) {
  data.frame(experiment_id = registry$experiment_id[i],
             method_id = parse_pipe_character_v2(registry$methods[i]),
             scheduled_replications = as.integer(registry$replications[i]),
             stringsAsFactors = FALSE)
}))
metrics <- merged_tables[["replication_metrics.csv"]]
observed_key <- unique(paste(metrics$experiment_id, metrics$method_id, sep = "|"))
expected_key <- paste(expected$experiment_id, expected$method_id, sep = "|")
uncovered <- setdiff(expected_key, observed_key)
if (length(uncovered) && isTRUE(require_complete)) {
  stop("Merged run does not cover ", length(uncovered), " registered cells, for example: ",
       paste(head(uncovered, 5L), collapse = "; "))
}

meta <- lapply(seq_len(n_shards), function(s) {
  d <- shard_dirs[s]
  manifest_path <- file.path(d, "completed_manifest.json")
  manifest <- if (file.exists(manifest_path)) read_json_v2(manifest_path) else list()
  list(shard = s,
       run_dir = d,
       completed = file.exists(manifest_path),
       replication_rows = if (file.exists(file.path(d, "replication_metrics.csv"))) {
         length(readLines(file.path(d, "replication_metrics.csv"), warn = FALSE)) - 1L
       } else 0L,
       commit = trimws(readLines(file.path(d, "implementation_commit.txt"),
                                 warn = FALSE)[1] %||% NA_character_),
       config_sha = trimws(readLines(file.path(d, "config_sha256.txt"),
                                     warn = FALSE)[1] %||% NA_character_),
       method_rows = manifest$method_rows %||% NA_integer_,
       coordinate_rows = manifest$coordinate_rows %||% NA_integer_,
       runner_failure_count = manifest$runner_failure_count %||% NA_integer_)
})
merge_manifest <- list(
  merged_run_id = merged_id,
  shard_prefix = run_id,
  n_shards = n_shards,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  implementation_commit = trimws(readLines(file.path(merged, "implementation_commit.txt"),
                                           warn = FALSE)[1] %||% NA_character_),
  config_sha256 = trimws(readLines(file.path(merged, "config_sha256.txt"),
                                   warn = FALSE)[1] %||% NA_character_),
  registered_cells = length(expected_key),
  observed_cells = length(observed_key),
  uncovered_cells = length(uncovered),
  replication_rows = nrow(metrics),
  coordinate_rows = nrow(merged_tables[["coordinate_metrics.csv"]]),
  shards = meta
)
writeLines(jsonlite::toJSON(merge_manifest, auto_unbox = TRUE, pretty = TRUE, null = "null"),
           file.path(merged, "merge_manifest.json"))
writeLines(jsonlite::toJSON(list(
  run_id = merged_id,
  merged = TRUE,
  n_shards = n_shards,
  method_rows = nrow(metrics),
  coordinate_rows = nrow(merged_tables[["coordinate_metrics.csv"]]),
  registered_cells = length(expected_key),
  uncovered_cells = length(uncovered),
  implementation_commit = merge_manifest$implementation_commit,
  final = isTRUE(first_request$final)
), auto_unbox = TRUE, pretty = TRUE), file.path(merged, "completed_manifest.json"))
cat("MERGED:", merged, "\n")
if (length(uncovered)) {
  cat("UNCOVERED CELLS:", length(uncovered), "\n")
}
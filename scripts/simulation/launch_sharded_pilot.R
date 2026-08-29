#!/usr/bin/env Rscript
# Round-3 pilot launcher: cluster self-normalised first-stage lambda regime.
args <- commandArgs(trailingOnly = TRUE)
root <- "D:/OneDrive/paper/Jointly Penalised and Debiased Inference for High-Dimensional Mixed-Effects Quantile Regression/jpDMEQR"
setwd(root)
source("scripts/00_source_v2.R", local = FALSE)
for (mod in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, mod, envir = .GlobalEnv)
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

config_path <- if (length(args) > 7 && nzchar(args[8])) {
  file.path(root, "config", args[8])
} else file.path(root, "config", "simulation_pilot_v5.csv")
run_id <- args[1] %||% stop("run_id required")
experiments <- if (length(args) > 1 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else character()
max_reps <- if (length(args) > 2 && nzchar(args[3])) as.numeric(args[3]) else Inf
jobs <- if (length(args) > 3 && nzchar(args[4])) as.integer(args[4]) else 1L
n_shards <- if (length(args) > 4 && nzchar(args[5])) as.integer(args[5]) else 1L
shard_index <- if (length(args) > 5 && nzchar(args[6])) as.integer(args[6]) else 1L
seed_exp <- if (length(args) > 6 && nzchar(args[7])) args[7] else NULL
options(jpDMEQR.cluster_cores = jobs)

rep_subset <- NULL
if (n_shards > 1L) {
  rep_subset <- seq(from = shard_index, to = max_reps, by = n_shards)
}
run_id <- if (n_shards > 1L) paste0(run_id, "_s", shard_index) else run_id

cat("run_id:", run_id, "| experiments:", paste(experiments, collapse = ","),
    "| max_reps:", max_reps, "| jobs:", jobs, "| shard:", shard_index, "/", n_shards,
    "| seed_exp:", seed_exp %||% "default", "\n")
ans <- run_registry_v2(
  root, config_path, run_id,
  jobs = jobs,
  max_reps = max_reps,
  allow_missing_benchmarks = TRUE,
  final = FALSE,
  experiment_ids = experiments,
  replicate_subset = rep_subset,
  seed_experiment = seed_exp
)
cat("DONE run_dir:", ans$run_dir, "\n")
cat("method_rows:", nrow(ans$metrics), "| coordinate_rows:", nrow(ans$coordinates),
    "| errors:", length(ans$errors), "\n")

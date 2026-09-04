#!/usr/bin/env Rscript
# Sharded pilot launcher (portable): root is auto-detected from the script
# location or overridden by the JPDMEQR_ROOT environment variable.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/launch_sharded_pilot.R"
root <- Sys.getenv("JPDMEQR_ROOT", unset = NA)
if (is.na(root) || !nzchar(root)) {
  root <- if (file.exists("R/profile_v2.R")) "." else
    normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
}
setwd(root)
args <- commandArgs(trailingOnly = TRUE)
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
# Cluster-level parallelism inside each task: default follows jobs; may be
# set independently via JPDMEQR_CLUSTER_CORES (Linux mclapply; results are
# deterministic - cluster contributions are order-preserved and additive).
options(jpDMEQR.cluster_cores = as.integer(Sys.getenv("JPDMEQR_CLUSTER_CORES", unset = as.character(jobs))))
# mu-CV Dantzig parallelism inside each task (the ~75-80% cost term): default
# 1 (serial, bit-identical to the frozen reference). May be raised via
# JPDMEQR_CV_CORES. Each (fold, mu) solve is independent; results are
# aggregated in fold order so wall-clock parallelism does not change the
# selected mu. On Linux this uses mclapply inside the task process, so keep
# `jobs` at 1 (use shards for task-level parallelism) to avoid nested-fork
# memory blowup; on Windows it uses a deterministic PSOCK cluster.
options(jpDMEQR.cv_cores = as.integer(Sys.getenv("JPDMEQR_CV_CORES", unset = "1")))

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

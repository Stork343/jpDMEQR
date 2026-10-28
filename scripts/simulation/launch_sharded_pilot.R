#!/usr/bin/env Rscript

# 分片运行试点的启动器：按实验集合切分，与主运行使用相同的门禁与最终模式标记，
# 每个分片写一个独立的运行目录，随后用 merge_shards.R 合并并运行试点门禁。

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/launch_sharded_pilot.R"
root <- Sys.getenv("JPDMEQR_ROOT", unset = NA)
if (is.na(root) || !nzchar(root)) {
  root <- if (file.exists("R/profile_v2.R")) "." else
    normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
}
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source_v2_module(root, "benchmark_adapters_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
config <- cli$config %||% file.path(root, "config", "simulation_pilot.csv")
run_id <- cli$`run-id` %||% stop("Supply --run-id=<shard prefix>.")
experiment_ids <- cli$`experiment-id` %||% stop("Supply --experiment-id=A,B,...")
experiment_ids <- trimws(strsplit(experiment_ids, ",", fixed = TRUE)[[1]])
jobs <- as.integer(cli$jobs %||% 1L)
max_reps <- as.numeric(cli$`max-reps` %||% Inf)
development <- as_bool_cli_v2(cli$development, FALSE)
allow_missing <- as_bool_cli_v2(cli$`allow-missing-benchmarks`, development)
resume <- as_bool_cli_v2(cli$resume, FALSE)
shard <- suppressWarnings(as.integer(cli$shard %||% NA_integer_))
if (is.finite(shard)) run_id <- paste0(run_id, "_s", shard)
# mu-CV Dantzig parallelism inside each task (default 1 = serial, bit-identical
# to the frozen reference; shards supply the process-level parallelism).
options(jpDMEQR.cv_cores = as.integer(Sys.getenv("JPDMEQR_CV_CORES", unset = "1")))

if (!file.exists(config)) stop("Configuration file is missing: ", config)
registry <- validate_registry_contract_v2(root, config, "pilot")
if (!registry$pass) {
  stop("Pilot registry validation failed: ", paste(registry$problems, collapse = "; "))
}
geometry <- validate_geometry_gate_v2(root, final = !development)
if (!geometry$pass) stop("Geometry gate must pass before the pilot: ", geometry$reason)
if (!development && allow_missing) stop("Non-development pilot cannot allow missing benchmarks.")

message("Pilot shard ", run_id, " covers ", paste(experiment_ids, collapse = ","))
ans <- run_registry_v2(
  root, config, run_id,
  jobs = jobs,
  max_reps = max_reps,
  allow_missing_benchmarks = allow_missing,
  final = !development,
  experiment_ids = experiment_ids,
  resume = resume
)
message("Pilot shard written to: ", ans$run_dir)
message("Method rows: ", nrow(ans$metrics), "; coordinate rows: ", nrow(ans$coordinates),
        "; runner failures: ", length(ans$errors))
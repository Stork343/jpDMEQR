#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/02_run_main.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source_v2_module(root, "benchmark_adapters_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
config <- cli$config %||% file.path(root, "config", "simulation_main.csv")
jobs <- as.integer(cli$jobs %||% 1L)
max_reps <- as.numeric(cli$`max-reps` %||% Inf)
run_id <- cli$`run-id` %||% paste0("main_", format(Sys.time(), "%Y%m%dT%H%M%S"))
development <- as_bool_cli_v2(cli$development, FALSE)
allow_missing <- as_bool_cli_v2(cli$`allow-missing-benchmarks`, development)
experiment_ids <- cli$`experiment-id` %||% ""
experiment_ids <- if (nzchar(experiment_ids)) {
  trimws(strsplit(experiment_ids, ",", fixed = TRUE)[[1]])
} else character()
# mu-CV Dantzig parallelism inside each task (default 1 = serial, bit-identical
# to the frozen reference; see launch_sharded_pilot.R for the nested-fork note).
options(jpDMEQR.cv_cores = as.integer(Sys.getenv("JPDMEQR_CV_CORES", unset = "1")))

if (!file.exists(config)) stop("Configuration file is missing: ", config)
registry_type <- if (identical(basename(config), "simulation_main.csv")) "main" else "preflight"
registry_check <- validate_registry_contract_v2(root, config, registry_type)
if (!registry_check$pass) {
  stop("Registry validation failed: ", paste(registry_check$problems, collapse = "; "))
}
if (length(registry_check$warnings)) {
  warning(paste(registry_check$warnings, collapse = "; "))
}

if (!development) {
  if (allow_missing) stop("Final execution cannot allow missing benchmarks.")
  if (is.finite(max_reps)) stop("Final execution cannot use --max-reps.")
  verify_freeze_manifest_v2(root, config)
} else {
  message("DEVELOPMENT MODE: outputs cannot be used for manuscript claims.")
}

ans <- run_registry_v2(
  root = root,
  config_path = config,
  run_id = run_id,
  jobs = jobs,
  max_reps = max_reps,
  allow_missing_benchmarks = allow_missing,
  experiment_ids = experiment_ids,
  final = !development
)
if (identical(basename(config), "simulation_preflight.csv")) {
  micro_pass <- length(ans$errors) == 0L && nrow(ans$metrics) > 0L &&
    all(ans$metrics$status %in% c("ok", "warning"))
  write_json_v2(list(
    pass = micro_pass,
    commit_sha = current_commit_v2(root),
    config_sha256 = file_checksum_v2(config),
    run_dir = normalizePath(ans$run_dir, mustWork = TRUE),
    method_rows = nrow(ans$metrics),
    coordinate_rows = nrow(ans$coordinates),
    runner_failure_count = length(ans$errors),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ), micro_preflight_manifest_path_v2(root))
  if (!micro_pass) stop("Micro-preflight run failed; inspect ", ans$run_dir)
}
message("Simulation run written to: ", ans$run_dir)

#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/02_run_main.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source_v2_module(root, "benchmark_adapters_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
config <- cli$config %||% file.path(root, "config", "simulation_main_b200.csv")
jobs <- as.integer(cli$jobs %||% 1L)
max_reps <- as.numeric(cli$`max-reps` %||% Inf)
run_id <- cli$`run-id` %||% paste0("main_", format(Sys.time(), "%Y%m%dT%H%M%S"))
development <- as_bool_cli_v2(cli$development, FALSE)
allow_missing <- as_bool_cli_v2(cli$`allow-missing-benchmarks`, development)
experiment_ids <- cli$`experiment-id` %||% ""
resume <- as_bool_cli_v2(cli$resume, FALSE)
experiment_ids <- if (nzchar(experiment_ids)) {
  trimws(strsplit(experiment_ids, ",", fixed = TRUE)[[1]])
} else character()
# mu-CV Dantzig parallelism inside each task (default 1 = serial, bit-identical
# to the frozen reference; see launch_sharded_pilot.R for the nested-fork note).
options(jpDMEQR.cv_cores = as.integer(Sys.getenv("JPDMEQR_CV_CORES", unset = "1")))

if (!file.exists(config)) stop("Configuration file is missing: ", config)
registry_type <- registry_type_for_path_v2(config)
method_version <- registry_method_version_v2(read_simulation_registry_v2(config))
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
  if (requires_freeze_manifest_v2(config)) verify_freeze_manifest_v2(root, config)
  if (identical(method_version, "profile_pair_20261022") &&
      identical(registry_type, "preflight")) {
    preflight_registry <- read_simulation_registry_v2(config)
    geometry <- validate_geometry_gate_v2(root, final = TRUE, method_version = method_version)
    benchmarks <- validate_benchmark_gate_v2(root, preflight_registry, final = TRUE)
    targets <- validate_target_gate_v2(root, preflight_registry, final = TRUE)
    if (!geometry$pass || !benchmarks$pass || !targets$pass) {
      stop("The preflight numerical prerequisites have not passed.")
    }
  }
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
  final = !development,
  resume = resume
)
if (identical(registry_type, "preflight")) {
  if (development) {
    message("Development-mode preflight run: no micro-preflight gate manifest is written.")
  } else {
    preflight_registry <- read_simulation_registry_v2(config)
    absent <- setdiff(preflight_registry$experiment_id, unique(ans$metrics$experiment_id))
    if (length(absent)) {
      stop("The micro-preflight gate needs every registered row; absent: ",
           paste(absent, collapse = ", "))
    }
    execution <- validate_micro_registry_execution_v2(root, ans, preflight_registry)
    micro_pass <- execution$pass
    manifest <- list(
      pass = micro_pass,
      final = TRUE,
      method_version = method_version,
      score_definition = if (identical(method_version, "profile_pair_20261022"))
        "profile_bandwidth_pair_2_half_minus_full" else NULL,
      scope = "registered_dgp_dispatch_and_record_accounting",
      commit_sha = current_commit_v2(root),
      config_sha256 = file_checksum_v2(config),
      run_dir = normalizePath(ans$run_dir, mustWork = TRUE),
      registry_rows = nrow(preflight_registry),
      method_rows = nrow(ans$metrics),
      coordinate_rows = nrow(ans$coordinates),
      runner_failure_count = length(ans$errors),
      execution_validation = execution,
      solver_configuration = solver_configuration_v2(),
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    path <- micro_preflight_manifest_path_v2(root, method_version)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write_json_v2(manifest, path)
    if (!micro_pass) stop("Micro-preflight run failed; inspect ", ans$run_dir)
  }
}
message("Simulation run written to: ", ans$run_dir)

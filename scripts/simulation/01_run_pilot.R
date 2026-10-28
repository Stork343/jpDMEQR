#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/01_run_pilot.R"
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
config <- cli$config %||% file.path(root, "config", "simulation_pilot.csv")
jobs <- as.integer(cli$jobs %||% 1L)
max_reps <- as.numeric(cli$`max-reps` %||% Inf)
run_id <- cli$`run-id` %||% paste0("pilot_", format(Sys.time(), "%Y%m%dT%H%M%S"))
development <- as_bool_cli_v2(cli$development, FALSE)
allow_missing <- as_bool_cli_v2(cli$`allow-missing-benchmarks`, development)

registry <- validate_registry_contract_v2(root, config, "pilot")
if (!registry$pass) stop("Pilot registry validation failed: ", paste(registry$problems, collapse = "; "))
method_version <- registry_method_version_v2(read_simulation_registry_v2(config))
geometry <- validate_geometry_gate_v2(root, final = !development, method_version = method_version)
if (!geometry$pass) stop("Geometry gate must pass before the pilot: ", geometry$reason)
if (!development) {
  pilot_registry <- read_simulation_registry_v2(config)
  benchmarks <- validate_benchmark_gate_v2(root, pilot_registry, final = TRUE)
  targets <- validate_target_gate_v2(root, pilot_registry, final = TRUE)
  if (!benchmarks$pass || !targets$pass) stop("The pilot benchmark or target prerequisites have not passed.")
}
if (!development && allow_missing) stop("Non-development pilot cannot allow missing benchmarks.")

ans <- run_registry_v2(
  root, config, run_id,
  jobs = jobs,
  max_reps = max_reps,
  allow_missing_benchmarks = allow_missing,
  final = !development
)
message("Pilot run written to: ", ans$run_dir)
message("Run scripts/simulation/03_summarise.R --run-id=", run_id,
        " --pilot-gate=true to evaluate the frozen pilot gate.")

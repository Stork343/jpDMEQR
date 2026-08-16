#!/usr/bin/env Rscript
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "scripts/simulation/01_run_pilot.R"
root <- if (file.exists("R/profile_v2.R")) "." else normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
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

message("Run geometry validation before interpreting this pilot.")
ans <- run_registry_v2(root, config, run_id, jobs = jobs, max_reps = max_reps,
                       allow_missing_benchmarks = TRUE)
message("Pilot run written to: ", ans$run_dir)

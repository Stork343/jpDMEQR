#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/05_preflight_freeze.R"
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
development <- as_bool_cli_v2(cli$development, FALSE)
manifest <- freeze_preflight_v2(
  root = root,
  config_path = config,
  final = !development,
  write_manifest = TRUE
)
cat("\nFreeze preflight summary\n")
cat("commit:", manifest$commit_sha, "\n")
cat("config sha256:", manifest$config_sha256, "\n")
cat("governance:", manifest$governance$pass, "\n")
cat("registry:", manifest$registry$pass, "\n")
cat("benchmark:", manifest$benchmark$pass, "\n")
cat("geometry:", manifest$geometry$pass, "\n")
cat("targets:", manifest$targets$pass, "\n")
cat("micro preflight:", manifest$micro_preflight$pass, "\n")
cat("pilot:", manifest$pilot$pass, "\n")
cat("overall:", manifest$pass, "\n")
if (!manifest$benchmark$pass) print(manifest$benchmark$table)
if (!manifest$targets$pass) print(manifest$targets$table)
if (!manifest$pass) {
  stop("Freeze preflight failed. Final-scale simulation remains blocked.")
}
message("Freeze preflight passed. The manifest now authorises the matching commit/config only.")

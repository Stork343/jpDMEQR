#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/00_validate_registry.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
config <- cli$config %||% file.path(root, "config", "simulation_main.csv")
type <- cli$type %||% if (identical(basename(config), "simulation_main.csv")) {
  "main"
} else if (identical(basename(config), "simulation_pilot.csv")) {
  "pilot"
} else "preflight"
ans <- validate_registry_contract_v2(root, config, type)
print(ans)
if (!ans$pass) stop("Registry validation failed: ", paste(ans$problems, collapse = "; "))
message("Registry validation passed: ", config)

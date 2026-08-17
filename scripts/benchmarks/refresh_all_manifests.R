#!/usr/bin/env Rscript
# Regenerate ALL acceptance manifests on the current commit (except the
# SQR-DEBIASED-IID reproduction cell, which is generated separately by
# accept_sqr_debiased_iid.R because it runs a B=200 diagnostic).
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/benchmarks/refresh_all_manifests.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)

source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source_v2_module(root, "benchmark_adapters_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

message("Refreshing manifests for the remaining adapters on commit ",
        current_commit_v2(root), " ...")
sys.source(file.path(root, "scripts", "benchmarks",
                     "accept_patch_adapters.R"), envir = environment())
sys.source(file.path(root, "scripts", "benchmarks",
                     "accept_remaining_benchmarks.R"), envir = environment())
sys.source(file.path(root, "scripts", "benchmarks",
                     "accept_oracle_adapters.R"), envir = environment())
sys.source(file.path(root, "scripts", "benchmarks",
                     "accept_sqr_debiased_iid.R"), envir = environment())
message("All acceptance manifests refreshed on commit ",
        current_commit_v2(root), ".")

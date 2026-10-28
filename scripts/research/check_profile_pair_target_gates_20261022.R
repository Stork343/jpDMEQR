root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2")) source_v2_module(root, module)
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args)) args[[1]] else format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
out <- file.path(root, "results", "scratch", "profile_pair_target_gates_20261022", run_id)
if (dir.exists(out)) stop("The paired target-gate verification directory already exists.")
dir.create(out, recursive = TRUE)
paths <- c(main = "config/simulation_pair_main_20261022.yml",
           pilot = "config/simulation_pair_pilot_20261022.yml",
           preflight = "config/simulation_pair_preflight_20261022.yml")
expected_counts <- c(main = 33L, pilot = 3L, preflight = 4L)
mode <- if (length(args) > 1L) args[[2]] else "all"
if (!mode %in% c("all", "root_context")) stop("Unknown paired gate verification mode.")
if (mode == "root_context") paths <- paths["preflight"]
checks <- list()
for (phase in names(paths)) {
  cfg <- read_simulation_registry_v2(file.path(root, paths[[phase]]))
  result <- validate_target_gate_v2(root, cfg, final = TRUE)
  frozen <- profile_pair_freeze_asset_paths_v2(root, cfg)
  stopifnot(isTRUE(result$pass), all(result$table$pass),
            identical(names(result$table), c("experiment_id", "asset", "pass", "reason")),
            nrow(result$table) == expected_counts[[phase]],
            nrow(frozen) == expected_counts[[phase]], !anyDuplicated(frozen$label),
            all(file.exists(frozen$path)),
            all(grepl("profile_pair_20261022/construction_v3/", frozen$path, fixed = TRUE)))
  saveRDS(result, file.path(out, paste0(phase, "_gate.rds")))
  utils::write.csv(result$table, file.path(out, paste0(phase, "_gate.csv")), row.names = FALSE)
  utils::write.csv(frozen, file.path(out, paste0(phase, "_freeze_paths.csv")), row.names = FALSE)
  checks[[phase]] <- list(pass = result$pass, rows = nrow(result$table),
                          profile_targets = sum(result$table$asset == "profile_target"),
                          population_directions = sum(result$table$asset == "population_direction"),
                          structural_audits = sum(result$table$asset == "structural_target_audit"))
  message("Actual paired target gate accepted: ", phase)
}
cfg <- read_simulation_registry_v2(file.path(root, "config", "simulation_pair_main_20261022.yml"))
empty <- cfg[cfg$experiment_id == "B01", , drop = FALSE]
stopifnot(nrow(empty) == 1L, nrow(profile_pair_freeze_asset_paths_v2(root, empty)) == 0L,
          isTRUE(validate_profile_pair_target_gate_v2(root, empty, final = TRUE)$pass))
legacy <- read_simulation_registry_v2(file.path(root, "config", "simulation_main.csv"))
testthat::expect_error(profile_pair_freeze_asset_paths_v2(root, legacy), "paired registry")
wrong <- read_simulation_registry_v2(file.path(root, paths[["preflight"]]))
wrong$lambda_gamma[wrong$experiment_id == "PF04"] <- 2
testthat::expect_error(validate_profile_pair_target_gate_v2(root, wrong, final = TRUE), "identity")
preflight <- read_simulation_registry_v2(file.path(root, paths[["preflight"]]))
previous_root <- getOption("jpDMEQR.source_root")
outside <- withr::with_dir(dirname(root), validate_profile_pair_target_gate_v2(root, preflight, final = TRUE))
stopifnot(isTRUE(outside$pass), nrow(outside$table) == 4L,
          identical(getOption("jpDMEQR.source_root"), previous_root))
saveRDS(outside, file.path(out, "outside_repository_cwd_gate.rds"))
manifest <- list(pass = all(vapply(checks, `[[`, logical(1), "pass")),
                 validation_commit = current_commit_v2(root),
                 actual_population_objects_checked = sum(vapply(checks, `[[`, integer(1), "rows")),
                 checks = checks, empty_asset_group = TRUE,
                 nonpaired_registry_rejected = TRUE, mismatched_working_penalty_rejected = TRUE,
                 outside_repository_cwd_accepted = TRUE, source_root_option_restored = TRUE)
jsonlite::write_json(manifest, file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines("COMPLETE", file.path(out, "COMPLETE.txt"))
print(manifest)

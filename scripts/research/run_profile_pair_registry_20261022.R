args <- commandArgs(trailingOnly = TRUE)
options(warn = 1)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
setwd(root)
source(file.path("scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(".", module, environment())
}
source(file.path("scripts", "simulation", "_run_registry_helpers.R"))
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "newton",
           JPDMEQR_DANTZIG_SOLVER = "clarabel_native",
           JPDMEQR_ASSEMBLY = "sequential", JPDMEQR_CV_CORES = "1",
           JPDMEQR_HESSIAN_REFRESH = "1", JPDMEQR_SHARED_FIRST_STAGE = "FALSE")
config_path <- file.path("config", "simulation_preflight_profile_pair_20261022.csv")
check <- validate_registry_contract_v2(".", config_path, "preflight")
if (!check$pass) stop(paste(check$problems, collapse = "; "))
registry <- read_simulation_registry_v2(config_path)
cfg <- config_row_to_list_v2(registry[registry$experiment_id == "PAIR01", , drop = FALSE])
target_path <- profile_pair_target_asset_path_v2(".", cfg$experiment_id)
direction_path <- profile_pair_direction_asset_path_v2(".", cfg$experiment_id)
if (file.exists(target_path) || file.exists(direction_path)) {
  if (!file.exists(target_path) || !file.exists(direction_path)) {
    stop("Both matching population assets are required.")
  }
  target <- readRDS(target_path)
  direction <- readRDS(direction_path)
  validate_profile_pair_target_asset_v2(target, cfg, final = FALSE)
  validate_profile_pair_direction_asset_v2(direction, cfg, final = FALSE)
} else {
  target <- build_profile_pair_target_asset_v2(
    cfg, n_population = 160L, repeats = 2L, development = TRUE
  )
  direction <- build_profile_pair_direction_asset_v2(
    cfg, target, n_population = 160L, repeats = 2L, development = TRUE
  )
  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(direction_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(target, target_path)
  saveRDS(direction, direction_path)
}
run_id <- if (length(args) >= 2L) args[[2]] else
  paste0("profile_pair_preflight_", format(Sys.time(), "%Y%m%dT%H%M%S"))
result <- run_registry_v2(
  root = ".", config_path = config_path, run_id = run_id, jobs = 1L,
  max_reps = 2L, allow_missing_benchmarks = FALSE, final = FALSE
)
if (length(result$errors) || !nrow(result$metrics) || !nrow(result$coordinates) ||
    any(!result$metrics$status %in% c("ok", "warning"))) {
  stop("Paired profile preflight contains failed methods; inspect ", result$run_dir)
}
expected_methods <- sum(vapply(seq_len(nrow(registry)), function(i) {
  cfg_i <- config_row_to_list_v2(registry[i, , drop = FALSE])
  length(cfg_i$methods) * min(cfg_i$replications, 2L)
}, integer(1)))
expected_coordinates <- sum(vapply(seq_len(nrow(registry)), function(i) {
  cfg_i <- config_row_to_list_v2(registry[i, , drop = FALSE])
  sum(cfg_i$methods != "POOLED-QR-LASSO") * cfg_i$target_coordinate_count *
    min(cfg_i$replications, 2L)
}, integer(1)))
if (nrow(result$metrics) != expected_methods ||
    nrow(result$coordinates) != expected_coordinates) {
  stop("Paired profile preflight is missing requested result rows.")
}
for (entry in list(list("metrics", 2L), list("coordinates", 3L),
                   list("theory", 4L), list("screening", 5L))) {
  if (!identical(names(result[[entry[[1]]]]), contract_schema_columns_v2(".", entry[[2]]))) {
    stop("Result schema differs from RESULTS_CONTRACT.md: ", entry[[1]])
  }
}
if (anyDuplicated(result$metrics[c("experiment_id", "replicate", "method_id")]) ||
    anyDuplicated(result$coordinates[c("experiment_id", "replicate", "method_id", "coordinate")])) {
  stop("Paired profile preflight contains duplicate result rows.")
}
message("Paired profile registry preflight completed: ", result$run_dir)

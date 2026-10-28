root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(normalizePath(Sys.getenv("JPDMEQR_NATIVE_LIBRARY", unset = file.path(
  root, "results", "scratch", "profile_pair_package_20261022", "library")), mustWork = TRUE), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2")) source_v2_module(root, module)
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
config_paths <- parse_pipe_character_v2(cli$configs %||% stop("Supply --configs=<versioned registry paths separated by |> ."))
implementation_commit <- cli$commit %||% stop("Supply the final implementation commit.")
implementation_commit <- profile_pair_asset_commit_v2(implementation_commit, FALSE)
if (!all(vapply(c("streamed", "conditional_gaussian"), function(engine) {
  profile_pair_population_commit_matches_v2(implementation_commit, current_commit_v2(root), engine, root)
}, logical(1)))) stop("Population construction dependencies differ from the declared implementation commit.")
kind <- cli$kind %||% "all"
if (!kind %in% c("all", "targets", "directions")) stop("Unknown population build kind.")
workers <- as.integer(cli$workers %||% "2")
if (!is.finite(workers) || workers < 1L || workers > 3L) stop("Population workers must be between 1 and 3.")
run_id <- cli$`run-id` %||% stop("Supply a new population build run-id.")
out <- file.path(root, "results", "scratch", "profile_pair_population_final_20261022", run_id)
if (dir.exists(out)) stop("The population run directory already exists.")
dir.create(out, recursive = TRUE)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "newton", JPDMEQR_DANTZIG_SOLVER = "clarabel_native",
           JPDMEQR_ASSEMBLY = "sequential", JPDMEQR_CV_CORES = "1",
           JPDMEQR_HESSIAN_REFRESH = "1", JPDMEQR_SHARED_FIRST_STAGE = "FALSE",
           OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
profile_pair_population_eigen_init_v2(root)
options(jpDMEQR.population_nuisance_backend = "compiled")
registries <- lapply(config_paths, read_simulation_registry_v2)
configurations <- unlist(lapply(registries, function(registry) {
  lapply(seq_len(nrow(registry)), function(i) config_row_to_list_v2(registry[i, , drop = FALSE]))
}), recursive = FALSE)
ids <- vapply(configurations, `[[`, character(1), "experiment_id")
if (anyDuplicated(ids)) stop("Population registry paths contain duplicate experiment identities.")
needed <- vapply(configurations, function(config) {
  identical(config$target_mode, "profile_mc") || "PROFILE-DQR-POP-H" %in% config$methods
}, logical(1))
configurations <- configurations[needed]
if (kind == "targets") {
  configurations <- configurations[vapply(configurations, function(config) {
    identical(config$target_mode, "profile_mc")
  }, logical(1))]
}
requested <- parse_pipe_character_v2(cli$ids %||% "")
if (length(requested)) {
  available <- vapply(configurations, `[[`, character(1), "experiment_id")
  if (!all(requested %in% available)) stop("A requested population experiment is unavailable.")
  configurations <- configurations[match(requested, available)]
}
if (!length(configurations)) stop("No population assets were selected.")
if (any(!vapply(configurations, function(config) identical(config$method_version, "profile_pair_20261022"), logical(1)))) {
  stop("Final population construction requires versioned paired registries.")
}
saveRDS(configurations, file.path(out, "configurations.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
build_case <- function(config) {
  id <- config$experiment_id
  case_dir <- file.path(out, id)
  dir.create(case_dir, recursive = TRUE)
  stage <- "input_validation"
  withCallingHandlers({
  started <- proc.time()
  checkpoint <- function(label, object) {
    stage <<- label
    saveRDS(object, file.path(case_dir, paste0(label, ".rds")), compress = FALSE)
  }
  target <- NULL
  target_path <- profile_pair_target_asset_path_v2(root, id)
  direction_path <- profile_pair_direction_asset_path_v2(root, id)
  if (identical(config$target_mode, "profile_mc")) {
    stage <- "target"
    if (file.exists(target_path)) {
      target <- readRDS(target_path)
      validate_profile_pair_target_asset_v2(target, config, final = TRUE, expected_commit = implementation_commit)
    } else {
      target <- build_profile_pair_target_asset_v2(
        config, n_population = 100000L, repeats = 4L, development = FALSE,
        implementation_commit = implementation_commit,
        control = list(population_engine = "streamed", batch_clusters = 1000L), checkpoint = checkpoint)
      dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(target, target_path, compress = FALSE)
    }
  }
  direction <- NULL
  if (kind != "targets" && "PROFILE-DQR-POP-H" %in% config$methods) {
    stage <- "direction"
    if (file.exists(direction_path)) {
      direction <- readRDS(direction_path)
      validate_profile_pair_direction_asset_v2(direction, config, final = TRUE, expected_commit = implementation_commit)
    } else {
      engine <- if (identical(config$target_mode, "structural")) "conditional_gaussian" else "streamed"
      if (engine == "conditional_gaussian") profile_pair_gaussian_row_conditions_v2(config)
      direction <- build_profile_pair_direction_asset_v2(
        config, target, n_population = 100000L, repeats = 4L, development = FALSE,
        implementation_commit = implementation_commit,
        control = list(population_engine = engine, batch_clusters = 1000L), checkpoint = checkpoint)
      dir.create(dirname(direction_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(direction, direction_path, compress = FALSE)
    }
  }
  elapsed <- proc.time() - started
  manifest <- list(
    experiment_id = id, implementation_commit = implementation_commit,
    construction_version = profile_pair_asset_version_v2()$construction_version,
    p = config$p, n_population = 100000L, repeats = 4L, development = FALSE,
    target_path = if (is.null(target)) NULL else target_path,
    direction_path = if (is.null(direction)) NULL else direction_path,
    target_source = direction$target_source %||% target$target_source,
    evaluation_point = direction$evaluation_point %||% target$evaluation_point,
    population_engine = direction$acceptance$population_engine %||% target$acceptance$population_engine,
    population_matrix_backend = getOption("jpDMEQR.population_matrix_backend", "R_BLAS"),
    population_nuisance_backend = getOption("jpDMEQR.population_nuisance_backend", "vectorised_R"),
    target_max_kkt = target$max_kkt_residual,
    max_nuisance_gradient = max(c(target$max_nuisance_gradient, direction$max_nuisance_gradient)),
    identity_error = max(c(target$profile_identity_error, direction$profile_identity_error)),
    runtime_sec = unname(elapsed[["elapsed"]]),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  jsonlite::write_json(manifest, file.path(case_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
  writeLines("COMPLETE", file.path(case_dir, "COMPLETE.txt"))
  message("Population assets completed: ", id)
  manifest
  }, error = function(error) {
    saveRDS(list(condition = error, stage = stage, config = config), file.path(case_dir, "failure.rds"))
    writeLines(c(stage, conditionMessage(error)), file.path(case_dir, "FAILED.txt"))
  })
}
result <- if (workers == 1L) lapply(configurations, build_case) else {
  parallel::mclapply(configurations, build_case, mc.cores = workers, mc.preschedule = FALSE)
}
saveRDS(result, file.path(out, "build_results.rds"))
if (any(vapply(result, inherits, logical(1), "try-error"))) stop("A population case failed; inspect its saved condition.")
jsonlite::write_json(list(pass = TRUE, implementation_commit = implementation_commit,
                          method_version = "profile_pair_20261022", cases = result,
                          workers = workers, n_population = 100000L, repeats = 4L),
                     file.path(out, "build_manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines("COMPLETE", file.path(out, "COMPLETE.txt"))

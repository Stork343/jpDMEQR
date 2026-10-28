root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(normalizePath(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library")), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2")
source_v2_module(root, "simulation_v2")
args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args)) args[[1]] else format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
out <- file.path(root, "results", "scratch", "profile_pair_population_performance_20261022", "native_stress", run_id)
if (dir.exists(out)) stop("The native population benchmark directory already exists.")
dir.create(out, recursive = TRUE)
stage <- "source"
globalCallingHandlers(error = function(error) {
  saveRDS(list(condition = error, stage = stage), file.path(out, "failure.rds"))
  writeLines(c(stage, conditionMessage(error)), file.path(out, "FAILED.txt"))
})
options(jpDMEQR.population_nuisance_backend = "compiled")
registry <- read_simulation_registry_v2(file.path(root, "config", "simulation_pilot.csv"))
config <- config_row_to_list_v2(registry[registry$experiment_id == "P01", , drop = FALSE])
seed <- seed_from_id_v2(97001L, "P01", 1L, "paired_population_direction")
started <- proc.time()[[3]]
stage <- "data_generation"
dat <- profile_pair_conditional_gaussian_data_v2(config, 100000L, seed)
generation_sec <- proc.time()[[3]] - started
started <- proc.time()[[3]]
stage <- "paired_profiles"
components <- profile_pair_conditional_gaussian_components_v2(
  dat, dat$beta0, config$tau, config$n_clusters^(-3 / 10), 1)
component_sec <- proc.time()[[3]] - started
saveRDS(components, file.path(out, "native_P01_full_scale_components.rds"))
jsonlite::write_json(list(
  development = TRUE, n_population = 100000L, p = 200L, seed = seed,
  generation_sec = generation_sec, component_sec = component_sec,
  nuisance_backend = components$population_nuisance_backend,
  maximum_nuisance_gradient = components$max_nuisance_gradient,
  profile_identity_error = components$profile_identity_error,
  covariance_identity_error = components$population_moments$covariance_identity_error),
  file.path(out, "native_P01_full_scale_benchmark.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
print(c(generation_sec = generation_sec, component_sec = component_sec,
        maximum_nuisance_gradient = components$max_nuisance_gradient))

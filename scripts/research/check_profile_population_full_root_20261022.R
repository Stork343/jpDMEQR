root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(normalizePath(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library")), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2")
source_v2_module(root, "simulation_v2")
args <- commandArgs(trailingOnly = TRUE)
id <- if (length(args)) args[[1]] else "C14"
run_id <- if (length(args) > 1L) args[[2]] else format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
out <- file.path(root, "results", "scratch", "profile_pair_population_performance_20261022", "full_root", run_id)
if (dir.exists(out)) stop("The full-root development directory already exists.")
dir.create(out, recursive = TRUE)
stage <- "source"
globalCallingHandlers(error = function(error) {
  saveRDS(list(condition = error, stage = stage), file.path(out, "failure.rds"))
  writeLines(c(stage, conditionMessage(error)), file.path(out, "FAILED.txt"))
})
profile_pair_population_eigen_init_v2(root)
options(jpDMEQR.population_nuisance_backend = "compiled")
registry <- read_simulation_registry_v2(file.path(root, "config", "simulation_pair_main_20261022.yml"))
config <- config_row_to_list_v2(registry[registry$experiment_id == id, , drop = FALSE])
seed <- seed_from_id_v2(87001L, id, 1L, "paired_population_target")
stage <- "data_generation"
started <- proc.time()[[3]]
dat <- profile_pair_population_data_v2(config, 100000L, seed)
dat$population_layout <- profile_pair_population_layout_v2(dat)
generation_sec <- proc.time()[[3]] - started
h <- min(config$h_multiplier * config$n_clusters^(-3 / 10) / 4, 100000^(-1 / 3))
stage <- "full_design_population_root"
started <- proc.time()[[3]]
fit <- solve_profile_pair_target_v2(
  dat, config$tau, h, config$lambda_gamma,
  control = list(population_engine = "streamed", batch_clusters = 1000L))
root_sec <- proc.time()[[3]] - started
saveRDS(fit, file.path(out, "fit.rds"))
jsonlite::write_json(list(
  development = TRUE, final_eligible = FALSE, experiment_id = id,
  n_population = 100000L, repeats = 1L, p = config$p, seed = seed,
  generation_sec = generation_sec, root_sec = root_sec,
  h_target = h, kkt_residual = fit$kkt_residual,
  maximum_nuisance_gradient = fit$max_nuisance_gradient,
  profile_identity_error = fit$profile_identity_error,
  evaluations = fit$component_evaluations, curvature = fit$curvature,
  population_matrix_backend = fit$population_matrix_backend,
  population_nuisance_backend = fit$population_nuisance_backend),
  file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines("COMPLETE", file.path(out, "COMPLETE.txt"))
print(c(generation_sec = generation_sec, root_sec = root_sec,
        kkt_residual = fit$kkt_residual, nuisance_gradient = fit$max_nuisance_gradient,
        evaluations = fit$component_evaluations))

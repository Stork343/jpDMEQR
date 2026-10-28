source("research/submission_completion_20261027/load.R")
args <- commandArgs(trailingOnly = TRUE)
method <- args[1L]
id <- args[2L]
registry <- read_simulation_registry_v2("config/simulation_main.csv")
config <- config_row_to_list_v2(registry[registry$experiment_id == id, ])
config$method_version <- "profile_pair_20261022"
optimizer <- if (length(args) >= 3L) args[3L] else "reference"
label <- paste(method, id, optimizer, if (length(args) >= 4L) args[4L] else "default", sep = "_")
destination <- file.path("results/scratch/submission_comparator_probes_20261027", label)
if (dir.exists(destination)) stop("The comparator probe already exists.")
dir.create(destination, recursive = TRUE)
started <- proc.time()[[3L]]
globalCallingHandlers(error = function(condition) saveRDS(list(condition = condition,
  elapsed = proc.time()[[3L]] - started), file.path(destination, "failure.rds")))
seed <- seed_from_id_v2(20261027L, id, 1L, "numerical_probe")
dat <- prepare_analysis_data_v2(simulate_from_config_v2(config, seed), config, seed)$data
coordinates <- resolve_target_coordinates_v2(config, dat)
if (method %in% c("BIAS-ADJ-LQMM", "LQMM")) {
  dat <- common_lqmm_subset_v2(dat, coordinates, config$tau, max_features = 20L,
                               seed = seed + 1L, screen_fraction = .5)$data
}
tuning <- build_tuning_for_analysis_v2(config, dat)
saveRDS(list(method = method, config = config, seed = seed, coordinates = coordinates,
             purpose = "numerical_validation", started = Sys.time()), file.path(destination, "request.rds"))
fit <- fit_benchmark_frozen_tuning_v2(method, dat, config$tau, coordinates, tuning, seed,
  context = list(config = config), control = list(ci_level = .95, B = 100L, nK = 15L,
    lambda_beta_grid = config$lambda_beta_multipliers, lambda_alpha_grid = c(.5, 1),
    solver_opts = if (length(args) >= 4L) list(basis_policy = args[4L]) else list(),
    maxit = 1000L, fit_control = list(max_iter = 10000L, max_backtrack = 50L, beta_tol = 1e-7,
                                    optimizer = optimizer)))
elapsed <- proc.time()[[3L]] - started
saveRDS(list(fit = fit, elapsed = elapsed), file.path(destination, "result.rds"))
cat(method, id, "status", fit$status, "elapsed", elapsed, "\n")

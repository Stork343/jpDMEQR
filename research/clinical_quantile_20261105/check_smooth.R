library(yaml)
library(reticulate)
source("research/clinical_quantile_20261105/engine.R")
use_python(file.path(getwd(), "results/scratch/cq_solver_env/bin/python"), required = TRUE)
Sys.setenv(NUMBA_CACHE_DIR = normalizePath("results/scratch/runtime", mustWork = TRUE))
source_python("research/clinical_quantile_20261105/smooth_path.py")
study <- read_yaml("config/application_gnhs_baseline_bmi_20261105.yml")
config <- read_yaml("config/clinical_quantile_20261105.yml")
panel <- cp_panel(study$input_directory)
outer <- file.path(study$output_root, "outer_1")
rows <- fread(file.path(outer, "training.csv"))
recipe <- readRDS(file.path(outer, "recipes.rds"))$inner_1$direction_recipe$base_recipe
features <- cq_features(rows[inner_fold != 1L], panel, recipe)
tau <- 0.25
alpha <- 0.1
lambda <- cq_lambda(features, tau, alpha, config)
clinical <- which(features$penalty_factor == 0)
initial <- rq.fit.br(features$X[, clinical, drop = FALSE] * features$weight,
  features$y * features$weight, tau = tau)
beta <- numeric(ncol(features$X))
beta[clinical] <- initial$coefficients
started <- proc.time()[["elapsed"]]
fitted <- fit_smooth_path(features$X, features$y, features$weight, features$penalty_factor,
  lambda, alpha, tau, features$n, beta)
elapsed <- proc.time()[["elapsed"]] - started
record <- data.table(rows = nrow(features$X), variables = ncol(features$X),
  models = ncol(fitted$beta), runtime_sec = elapsed, bandwidth = fitted$bandwidth,
  maximum_kkt = max(fitted$diagnostics[, 4L]))
saveRDS(list(fit = fitted, features = features, tau = tau, alpha = alpha, lambda = lambda),
  "results/scratch/runtime/cq_smooth_check.rds")
fwrite(record, "results/preflight/clinical_quantile_20261105/smooth_path_checks.csv")
print(record)

library(yaml)
source("research/clinical_quantile_20261105/engine.R")
source("research/clinical_quantile_20261105/qp_path.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
study <- read_yaml(args[1L])
config <- read_yaml("config/clinical_quantile_20261105.yml")
panel <- cp_panel(study$input_directory)
outer <- file.path(study$output_root, "outer_1")
rows <- fread(file.path(outer, "training.csv"))
recipe <- readRDS(file.path(outer, "recipes.rds"))$inner_1$direction_recipe$base_recipe
features <- cq_features(rows[inner_fold != 1L], panel, recipe)
started <- proc.time()[["elapsed"]]
path <- cq_qp_path(features, 0.25, 0, config)
elapsed <- proc.time()[["elapsed"]] - started
exact <- cq_exact(features, 0.25, path$lambda[20L], 0, config)
objective_difference <- cq_objective(features, path$beta[, 20L], 0.25, path$lambda[20L], 0) - exact$objective
stopifnot(abs(objective_difference) < 1e-5)
record <- data.table(study_id = study$study_id, rows = nrow(features$X),
  variables = ncol(features$X), models = ncol(path$beta), runtime_sec = elapsed,
  total_iterations = sum(path$iterations), maximum_iterations = max(path$iterations), objective_difference)
fwrite(record, file.path("results/preflight/clinical_quantile_20261105", paste0(study$study_id, "_qp_path.csv")))
print(record)

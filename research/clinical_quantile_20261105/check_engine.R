library(yaml)
library(jsonlite)
source("research/clinical_quantile_20261105/engine.R")
config <- read_yaml("config/clinical_quantile_20261105.yml")
study <- read_yaml("config/application_ra_westchina_20261105.yml")
panel <- cp_panel(study$input_directory)
outer <- file.path(study$output_root, "outer_1")
rows <- fread(file.path(outer, "training.csv"))
recipes <- readRDS(file.path(outer, "recipes.rds"))
recipe <- recipes$inner_1$direction_recipe$base_recipe
features <- cq_features(rows[inner_fold != 1L], panel, recipe)
clinical <- which(features$penalty_factor == 0)
protein <- which(features$penalty_factor == 1)[seq_len(12L)]
columns <- c(clinical, protein)
features$X <- features$X[, columns, drop = FALSE]
features$penalty_factor <- features$penalty_factor[columns]
checks <- list()
for (tau in unlist(config$quantiles)) for (alpha in c(0, 0.5, 1)) {
  path <- cq_path(features, tau, alpha, config)
  index <- 20L
  exact <- cq_exact(features, tau, path$lambda[index], alpha, config)
  approximate <- cq_objective(features, path$beta[, index], tau, exact$lambda, alpha)
  weight <- nrow(features$X) * features$weight / features$n
  transformed <- mean(ep_loss(features$y * weight, drop(features$X %*% exact$beta) * weight, tau))
  direct <- sum(features$weight * ep_loss(features$y, drop(features$X %*% exact$beta), tau)) / features$n
  stopifnot(abs(transformed - direct) < 1e-12, approximate - exact$objective < 0.005)
  checks[[length(checks) + 1L]] <- data.table(tau, alpha, objective = exact$objective,
    path_objective_difference = approximate - exact$objective,
    homogeneity_error = abs(transformed - direct), gap = exact$diagnostics$gap,
    stationarity = exact$diagnostics$stationarity)
}
output <- "results/preflight/clinical_quantile_20261105"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(checks), file.path(output, "engine_checks.csv"))
write_json(list(status = "passed", source = "actual_ra_discovery_training_patients",
  objective_checks = length(checks), generated_observations = 0),
  file.path(output, "engine_acceptance.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(checks))

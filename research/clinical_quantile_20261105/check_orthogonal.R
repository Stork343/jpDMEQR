library(yaml)
library(jsonlite)
source("research/clinical_quantile_20261105/engine.R")
source("research/clinical_quantile_20261105/orthogonal.R")
study <- read_yaml("config/application_gnhs_baseline_bmi_20261105.yml")
panel <- cp_panel(study$input_directory)
outer <- file.path(study$output_root, "outer_1")
rows <- fread(file.path(outer, "training.csv"))
recipe <- readRDS(file.path(outer, "recipes.rds"))$training$direction_recipe$base_recipe
features <- cq_features(rows, panel, recipe)
k <- match("C_prior_bmi", colnames(features$X))
stopifnot(!is.na(k))
columns <- c(which(features$penalty_factor == 0), which(features$penalty_factor == 1)[seq_len(12L)])
X <- features$X[, columns, drop = FALSE]
k <- match("C_prior_bmi", colnames(X))
checks <- list()
for (rho in c(0, 0.4, 0.8)) {
  tau <- 0.75
  transformed <- cq_transform_blocks(X, features$cluster_id, tau, rho)
  inverse <- cq_transform_blocks(X, features$cluster_id, tau, rho, inverse = TRUE)
  composition <- cq_transform_blocks(transformed, features$cluster_id, tau, rho)
  stopifnot(max(abs(composition - inverse)) < 1e-9)
  direct_error <- 0
  for (indices in split(seq_len(nrow(X)), features$cluster_id)) {
    m <- length(indices)
    C <- tau * (1 - tau) * ((1 - rho) * diag(m) + matrix(rho, m, m))
    direct_error <- max(direct_error, max(abs(solve(C, X[indices, , drop = FALSE]) - inverse[indices, , drop = FALSE])))
  }
  d <- transformed[, k]
  controls <- transformed[, -k, drop = FALSE]
  decomposition <- qr(controls)
  stopifnot(decomposition$rank == ncol(controls))
  residual <- drop(qr.resid(decomposition, d))
  derivative <- drop(crossprod(controls, residual)) / features$n
  J <- sum(residual * d) / features$n
  variance <- sum(residual^2) / features$n
  stopifnot(direct_error < 1e-9, max(abs(derivative)) < 1e-9, abs(J - variance) < 1e-9, J > 0)
  checks[[length(checks) + 1L]] <- data.table(rho, direct_error,
    maximum_nuisance_derivative = max(abs(derivative)), sensitivity = J,
    variance_identity_error = abs(J - variance))
}
output <- "results/preflight/clinical_quantile_20261105"
fwrite(rbindlist(checks), file.path(output, "orthogonal_identity_checks.csv"))
write_json(list(status = "passed", source = "actual_gnhs_discovery_training_design",
  checks = length(checks), clinical_outcomes_used_in_identity_test = FALSE),
  file.path(output, "orthogonal_identity_acceptance.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(checks))

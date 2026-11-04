library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
library(glmnet)
source("research/drugmatrix_application_20261101/prediction.R")
source("research/empirical_selection_20261101/prediction_core.R")

config <- read_yaml("config/discovery_rice_baap_20261104.yml")
config$elastic_net_alpha <- 0
config$elastic_net_path_length <- 3L
config$ridge_variance_ratio_max <- 100
config$ridge_variance_ratio_min <- 0.01
panel <- readRDS("data/derived/drugmatrix_20261101/panel.rds")
panel$allowed_probes <- colnames(panel$X)
rows <- readRDS("data/derived/drugmatrix_20261101/discovery.rds")[exploration_role == "training"]
rows <- rows[compound_id %in% head(sort(unique(rows$compound_id)), 50L)]
features <- dm_features(rows, panel, dm_recipe(rows, panel))
columns <- c(which(features$penalty_factor == 0), head(which(features$penalty_factor > 0), 30L))
features$X <- features$X[, columns, drop = FALSE]
features$penalty_factor <- features$penalty_factor[columns]
fit <- ep_path(features, "ENET", NA_real_, config, 0L)
ratios <- exp(seq(log(100), log(0.01), length.out = 3L))
gram <- crossprod(features$X, features$X * features$weight) / features$n
rhs <- crossprod(features$X, features$y * features$weight) / features$n
checks <- rbindlist(lapply(seq_along(ratios), function(j) {
  penalty <- sum(features$penalty_factor > 0) / features$n * ratios[j]
  expected <- drop(solve(gram + diag(penalty * features$penalty_factor), rhs))
  beta <- fit$beta[, j]
  objective <- function(value) {
    sum(features$weight * (features$y - drop(features$X %*% value))^2) / (2 * features$n) +
      penalty * sum(features$penalty_factor * value^2) / 2
  }
  data.table(variance_ratio = ratios[j], objective_error = abs(objective(beta) - objective(expected)),
             maximum_prediction_error = max(abs(features$X %*% (beta - expected))),
             intercept_native = beta[1L], intercept_reference = expected[1L])
}))
print(checks)
stopifnot(max(checks$objective_error) < 1e-8, max(checks$maximum_prediction_error) < 1e-5)
output <- "results/preflight/empirical_selection_20261101/prediction_core"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(checks, file.path(output, "ridge_checks.csv"))

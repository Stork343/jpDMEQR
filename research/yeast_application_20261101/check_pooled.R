library(data.table)
library(jsonlite)
source("research/score_validity_20261026/load_joint.R")
source("research/soynam_application_20261031/load_pooled_lp.R")
source("research/yeast_application_20261101/prediction.R")
source("research/yeast_application_20261101/comparators.R")

data <- readRDS("data/derived/yeast_genomic_20261101/growth_discovery.rds")
rows <- data[trait == "Magnesium_Sulfate" & exploration_role == "training"]
rows <- rows[strain_id %in% head(unique(rows$strain_id), 60L)]
panel <- readRDS("data/derived/yeast_genomic_20261101/panel.rds")
recipe <- yeast_recipe(rows, panel)
features <- yeast_features(rows, panel, recipe)
chromosome <- sub(".*_(chr[0-9]+)_.*", "\\1", colnames(features$X))
selected <- c(which(features$penalty_factor == 0), which(!duplicated(chromosome) & grepl("^chr", chromosome)))
features$X <- features$X[, selected, drop = FALSE]
features$penalty_factor <- features$penalty_factor[selected]
row_weight <- nrow(features$X) * features$weight / features$n
lambda <- c(0.1, 0.05, 0)
checks <- rbindlist(lapply(c(0.25, 0.5, 0.75), function(tau_value) {
  fit <- yeast_pooled_path(features, tau_value, lambda = lambda)
  rbindlist(lapply(seq_along(lambda), function(j) {
    reference <- if (lambda[j] == 0) {
      quantreg::rq.fit.fnb(features$X * row_weight, features$y * row_weight, tau = tau_value, eps = 1e-9)
    } else {
      quantreg::rq.fit.lasso(features$X * row_weight, features$y * row_weight,
        tau = tau_value, lambda = 2 * nrow(features$X) * lambda[j] * features$penalty_factor, eps = 1e-9)
    }
    objective <- function(beta) {
      residual <- features$y - drop(features$X %*% beta)
      sum(features$weight * residual * (tau_value - (residual < 0))) / features$n +
        lambda[j] * sum(features$penalty_factor * abs(beta))
    }
    data.table(tau = tau_value, lambda = lambda[j],
      objective_error = abs(objective(fit$beta[, j]) - objective(reference$coefficients)),
      duality_gap = fit$primal_dual_gap[j])
  }))
}))
stopifnot(max(checks$objective_error) < 1e-7, max(abs(checks$duality_gap)) < 1e-7)
directory <- "results/preflight/empirical_selection_20261101/yeast_pooled"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
fwrite(checks, file.path(directory, "checks.csv"))
print(checks)

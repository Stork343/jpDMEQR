.libPaths(c(normalizePath("results/scratch/dependencies/library"), .libPaths()))
library(data.table)
library(yaml)
library(highs)
library(quantreg)
source("research/score_validity_20261026/load_joint.R")
source("research/soynam_application_20261031/data.R")
source("research/soynam_application_20261031/load_pooled_lp.R")
source("research/soynam_application_20261031/pooled_lp.R")
config <- read_yaml("config/application_soynam_20261031.yml")
input <- readRDS(file.path(config$results_directory, "input.rds"))
eligible <- input$line_folds[prediction_fold != config$prediction$heldout_fold]
ids <- eligible[, head(ID, 10L), by = Family]$V1
full <- soy_design(input, prediction = TRUE)
data <- soy_subset(full, which(full$observations$ID %in% ids))
small <- data
selected <- c(match(c("ACC10_2013", "ACC10_2014", paste0("Block", 1:8), "check_surface"), colnames(data$X)),
  which(data$penalty_factor > 0)[c(1L, 2L)])
small$X <- data$X[, selected]
small$penalty_factor <- data$penalty_factor[selected]
checks <- list()
for (tau in c(.25, .5, .75)) {
  candidate <- soy_pooled_compiled_path(small, tau, lambda = c(.1, .01, 0))
  repeated <- soy_pooled_compiled_path(small, tau, lambda = c(.1, .01, 0))
  if (max(abs(candidate$beta - repeated$beta)) > 1e-12) stop("The pooled path is not deterministic.")
  for (j in 1:3) {
    scale <- apply(small$X, 2L, function(x) sqrt(mean((x - mean(x))^2)))
    penalty <- length(small$y) * candidate$lambda[j] * small$penalty_factor * scale
    positive <- which(penalty > 0)
    P <- diag(penalty)[positive, , drop = FALSE]
    reference <- quantreg::rq.fit.fnb(rbind(small$X, P, -P), c(small$y, rep(0, 2L * length(positive))), tau = tau)
    objective <- function(beta) {
      residual <- small$y - drop(small$X %*% beta)
      (sum(residual * (tau - (residual < 0))) + sum(penalty * abs(beta))) / length(small$y)
    }
    difference <- abs(objective(candidate$beta[, j]) - objective(reference$coefficients))
    if (difference > 2e-6) stop("The compiled pooled fit differs from independent quantreg.")
    checks[[length(checks) + 1L]] <- data.table(tau, lambda = candidate$lambda[j], objective_difference = difference,
      primal_dual_gap = candidate$diagnostics$primal_dual_gap[j],
      dual_violation = candidate$diagnostics$dual_violation[j], mapping_error = candidate$diagnostics$mapping_error[j])
  }
}
fwrite(rbindlist(checks), "results/preflight/soynam_20261031/pooled_lp_checks.csv")
print(rbindlist(checks))
started <- proc.time()[3L]
fitted <- soy_pooled_compiled_path(data, .5, progress = TRUE)
seconds <- proc.time()[3L] - started
reference <- soy_pooled_lp_path(data, .5, lambda = fitted$lambda, progress = TRUE)
penalty_scale <- apply(data$X, 2L, function(x) sqrt(mean((x - mean(x))^2))) * data$penalty_factor
objective <- function(beta, lambda) {
  residual <- data$y - drop(data$X %*% beta)
  mean(abs(residual)) / 2 + lambda * sum(penalty_scale * abs(beta))
}
discrepancy <- vapply(seq_along(fitted$lambda), function(j)
  abs(objective(fitted$beta[, j], fitted$lambda[j]) - objective(reference$beta[, j], reference$lambda[j])), numeric(1))
if (max(abs(fitted$lambda - reference$lambda)) > 1e-9 || max(discrepancy) > 2e-6)
  stop("The persistent pooled LP path differs from the independent rebuilt path.")
saveRDS(fitted, "results/scratch/soynam_20261031/pooled_lp_compiled_path.rds")
benchmark <- data.table(lines = data$n, p = ncol(data$X), path_points = length(fitted$lambda), seconds,
  maximum_objective_difference = max(discrepancy), maximum_gap = max(abs(fitted$diagnostics$primal_dual_gap)),
  maximum_dual_violation = max(fitted$diagnostics$dual_violation), mapping_error = max(fitted$diagnostics$mapping_error))
fwrite(benchmark, "results/preflight/soynam_20261031/pooled_lp_path_checks.csv")
print(benchmark)

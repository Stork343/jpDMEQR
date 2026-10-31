library(data.table)
library(yaml)
library(hqreg)
library(quantreg)
library(CVXR)
source("research/score_validity_20261026/load_joint.R")
source("research/soynam_application_20261031/prediction.R")
config <- read_yaml("config/application_soynam_20261031.yml")
input <- readRDS(file.path(config$results_directory, "input.rds"))
ids <- input$line_folds[, head(ID, 5L), by = Family]$V1
rows <- which(input$observations$ID %in% ids)
observed <- input$observations[rows]
fixed_columns <- c("ACC10_2013", "ACC10_2014", paste0("Block", 1:8), "check_surface")
fixed <- input$fixed_design[rows, fixed_columns, drop = FALSE]
G <- input$genotypes[match(observed$ID, rownames(input$genotypes)), 1:2, drop = FALSE]
X <- cbind(fixed, G)
data <- sb_data(observed$YLD / 10, X, matrix(1, nrow(X), 1L), observed$cluster_id)
data$penalty_factor <- c(rep(0, ncol(fixed)), rep(1, ncol(G)))
options <- config$prediction
options$hqreg_message <- TRUE
checks <- list()
for (tau in c(.25, .5, .75)) {
  printed <- capture.output(fitted <- soy_pooled_path(data, tau, options, lambda = c(.1, .01, 0)))
  invisible(capture.output(repeated <- soy_pooled_path(data, tau, options, lambda = c(.1, .01, 0))))
  gamma_lines <- grep("Gamma =", printed, value = TRUE)
  gammas <- as.numeric(sub(".*Gamma = ([0-9.]+).*", "\\1", gamma_lines))
  if (length(gammas) != 3L || any(!is.finite(gammas))) stop("The native smoothing diagnostics are incomplete.")
  reference <- quantreg::rq.fit.fnb(X, data$y, tau = tau)
  residual <- data$y - X %*% fitted$beta[, 3L]
  original <- data$y - X %*% reference$coefficients
  loss_difference <- abs(mean(residual * (tau - (residual < 0))) - mean(original * (tau - (original < 0))))
  deterministic_error <- max(abs(fitted$beta - repeated$beta))
  print(data.table(tau = tau, loss_difference = loss_difference,
    deterministic_error = deterministic_error, mapping_error = fitted$mapping_error,
    maximum_iterations = max(fitted$iterations)))
  if (loss_difference > gammas[3L] / 4 + 2e-6 || deterministic_error > 1e-12 || fitted$mapping_error > 1e-9)
    stop("The actual-data pooled benchmark check failed.")
  scales <- sqrt(colMeans(sweep(X, 2L, colMeans(X), "-")^2))
  for (j in 1:3) {
    gamma <- gammas[j]
    beta <- CVXR::Variable(ncol(X))
    residual <- data$y - X %*% beta
    objective <- sum(CVXR::huber(residual, gamma) / (4 * gamma) + (tau - .5) * residual) / length(data$y) +
      fitted$lambda[j] * sum(CVXR::multiply(data$penalty_factor * scales, abs(beta)))
    invisible(capture.output(solved <- CVXR::solve(CVXR::Problem(CVXR::Minimize(objective)),
      solver = "CLARABEL", tol_gap_abs = 1e-10, tol_gap_rel = 1e-10, tol_feas = 1e-10)))
    if (solved$status != "optimal") stop("The independent reference solver failed.")
    r <- data$y - drop(X %*% fitted$beta[, j])
    native_objective <- mean(ifelse(abs(r) <= gamma, r^2 / (4 * gamma), abs(r) / 2 - gamma / 4) +
      (tau - .5) * r) + fitted$lambda[j] * sum(data$penalty_factor * scales * abs(fitted$beta[, j]))
    objective_difference <- abs(native_objective - solved$value)
    if (objective_difference > 2e-6) stop("The native quantile continuation differs from the independent convex reference.")
    checks[[length(checks) + 1L]] <- data.table(tau = tau, lambda = fitted$lambda[j],
      actual_lines = data$n, actual_observations = length(data$y), native_gamma = gamma,
      convex_objective_difference = objective_difference, zero_penalty_loss_difference = loss_difference,
      deterministic_error = deterministic_error, mapping_error = fitted$mapping_error,
      maximum_iterations = max(fitted$iterations))
  }
}
fwrite(rbindlist(checks), "results/preflight/soynam_20261031/pooled_benchmark_checks.csv")
print(rbindlist(checks))

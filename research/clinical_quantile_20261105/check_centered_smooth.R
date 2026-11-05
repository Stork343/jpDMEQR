library(yaml)
library(jsonlite)
library(reticulate)
source("research/clinical_quantile_20261105/engine.R")
use_python(file.path(getwd(), "results/scratch/cq_solver_env/bin/python"), required = TRUE)
Sys.setenv(NUMBA_CACHE_DIR = normalizePath("results/scratch/runtime", mustWork = TRUE))
source_python("research/clinical_quantile_20261105/smooth_path.py")
config <- read_yaml("config/clinical_quantile_20261105.yml")
study <- read_yaml("config/application_ra_westchina_20261105.yml")
panel <- cp_panel(study$input_directory)
outer <- file.path(study$output_root, "outer_1")
rows <- fread(file.path(outer, "training.csv"))
recipe <- readRDS(file.path(outer, "recipes.rds"))$inner_1$direction_recipe$base_recipe
features <- cq_features(rows[inner_fold != 1L], panel, recipe)
columns <- c(which(features$penalty_factor == 0), which(features$penalty_factor == 1)[seq_len(12L)])
features$X <- features$X[, columns, drop = FALSE]
features$penalty_factor <- features$penalty_factor[columns]
N <- nrow(features$X)
p <- ncol(features$X)
penalized <- which(features$penalty_factor == 1)
clinical <- which(features$penalty_factor == 0)
d <- length(penalized)
a <- N * features$weight / features$n
X <- features$X * a
Q <- qr.Q(qr(X[, clinical, drop = FALSE]))
numerical_gram <- crossprod(crossprod(Q, X)) / N + diag(features$penalty_factor)
delta <- features$n^(-0.2)
select <- sparseMatrix(i = seq_len(d), j = penalized, x = 1, dims = c(d, p))
equality <- cbind(X, Diagonal(N), -Diagonal(N), Diagonal(N), -Diagonal(N), Matrix(0, N, d, sparse = TRUE))
bounds <- cbind(Matrix(0, 4L * N, p, sparse = TRUE), -Diagonal(4L * N), Matrix(0, 4L * N, d, sparse = TRUE))
positive <- cbind(select, Matrix(0, d, 4L * N, sparse = TRUE), -Diagonal(d))
negative <- cbind(-select, Matrix(0, d, 4L * N, sparse = TRUE), -Diagonal(d))
A <- cq_sparse(rbind(equality, bounds, positive, negative))
checks <- list()
for (tau in unlist(config$quantiles)) for (alpha in unlist(config$alpha)) {
  lambda <- cq_lambda(features, tau, alpha, config)[20L]
  initial <- rq.fit.br(features$X[, clinical, drop = FALSE] * features$weight,
    features$y * features$weight, tau = tau)
  beta <- numeric(p)
  beta[clinical] <- initial$coefficients
  fit <- fit_smooth_path(features$X, features$y, features$weight, features$penalty_factor,
    c(lambda), alpha, tau, features$n, beta)
  diagonal <- c(lambda * (1 - alpha) * features$penalty_factor,
    rep(tau / (delta * N), N), rep((1 - tau) / (delta * N), N), numeric(2L * N + d))
  P <- Matrix::Matrix(Diagonal(length(diagonal), diagonal), sparse = TRUE)
  P[seq_len(p), seq_len(p)] <- P[seq_len(p), seq_len(p)] + numerical_gram / features$n
  P <- cq_sparse(P)
  linear <- c(numeric(p + 2L * N), rep(tau / N, N), rep((1 - tau) / N, N), rep(lambda * alpha, d))
  shifted_y <- features$y * a + (2 * tau - 1) * delta / 2
  reference <- clarabel(A = A, b = c(shifted_y, numeric(4L * N + 2L * d)), P = P, q = linear,
    cones = list(z = as.integer(N), l = as.integer(4L * N + 2L * d)),
    control = clarabel_control(max_iter = 300L, verbose = FALSE,
      tol_gap_abs = 1e-10, tol_gap_rel = 1e-10, tol_feas = 1e-10))
  stopifnot(names(solver_status_descriptions())[reference$status] == "Solved")
  objective_error <- abs(fit$diagnostics[1L, 5L] - reference$obj_val)
  prediction_rmse <- sqrt(mean((X %*% fit$beta[, 1L] - X %*% reference$x[seq_len(p)])^2))
  if (objective_error >= 1e-6 || prediction_rmse >= 1e-3)
    stop(sprintf("tau=%g alpha=%g objective=%g prediction=%g", tau, alpha, objective_error, prediction_rmse))
  checks[[length(checks) + 1L]] <- data.table(tau, alpha, objective_error,
    prediction_rmse, maximum_kkt = fit$diagnostics[1L, 4L])
}
output <- "results/preflight/clinical_quantile_20261105"
fwrite(rbindlist(checks), file.path(output, "centered_smooth_reference.csv"))
write_json(list(status = "passed", source = "actual_ra_training_data", checks = length(checks),
  maximum_objective_error = max(rbindlist(checks)$objective_error),
  maximum_prediction_rmse = max(rbindlist(checks)$prediction_rmse)),
  file.path(output, "centered_smooth_acceptance.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(checks))

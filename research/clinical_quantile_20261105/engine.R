library(data.table)
library(Matrix)
library(hqreg)
library(clarabel)
library(quantreg)
library(matrixStats)
source("research/clinical_search_20261105/panel.R")
source("research/empirical_selection_20261101/prediction_core.R")

cq_features <- function(rows, panel, recipe, columns = NULL) {
  features <- cp_features(rows, panel, recipe)
  clinical <- which(colnames(features$X) == "(Intercept)" |
    startsWith(colnames(features$X), "C_") | startsWith(colnames(features$X), "M_C_"))
  if (is.null(columns)) {
    original <- features$X[, clinical, drop = FALSE] * sqrt(features$weight)
    decomposition <- qr(original, tol = 1e-10)
    chosen <- sort(decomposition$pivot[seq_len(decomposition$rank)])
    preserved <- qr.fitted(qr(original[, chosen, drop = FALSE], tol = 1e-10), original)
    stopifnot(max(abs(preserved - original)) < 1e-8)
    columns <- colnames(features$X)[c(clinical[chosen], setdiff(seq_len(ncol(features$X)), clinical))]
  }
  stopifnot(all(columns %in% colnames(features$X)))
  features$X <- features$X[, columns, drop = FALSE]
  features$penalty_factor <- as.numeric(!(columns == "(Intercept)" |
    startsWith(columns, "C_") | startsWith(columns, "M_C_")))
  features
}

cq_objective <- function(features, beta, tau, lambda, alpha) {
  sum(features$weight * ep_loss(features$y, drop(features$X %*% beta), tau)) / features$n +
    lambda * sum(features$penalty_factor * (alpha * abs(beta) + (1 - alpha) * beta^2 / 2))
}

cq_lambda <- function(features, tau, alpha, config) {
  clinical <- which(features$penalty_factor == 0)
  weight <- features$weight / features$n
  fit <- quantreg::rq.fit.br(features$X[, clinical, drop = FALSE] * weight, features$y * weight, tau = tau)
  dual <- weight * (fit$dual - (1 - tau))
  stopifnot(max(abs(crossprod(features$X[, clinical, drop = FALSE], dual))) < 1e-7)
  maximum <- max(abs(crossprod(features$X[, features$penalty_factor == 1, drop = FALSE], dual)))
  stopifnot(is.finite(maximum), maximum > 0)
  maximum / max(alpha, config$ridge_lambda_reference_alpha) *
    exp(seq(0, log(config$lambda_min_ratio), length.out = config$path_length))
}

cq_path <- function(features, tau, alpha, config) {
  lambda <- cq_lambda(features, tau, alpha, config)
  weight <- nrow(features$X) * features$weight / features$n
  fitted <- hqreg::hqreg_raw(features$X * weight, features$y * weight,
    method = "quantile", tau = tau, alpha = alpha, lambda = lambda, intercept = FALSE,
    penalty.factor = features$penalty_factor, screen = "ASR", dfmax = ncol(features$X) + 1L,
    max.iter = config$path_max_iterations, eps = config$path_tolerance)
  stopifnot(ncol(fitted$beta) == length(lambda), all(is.finite(fitted$beta)),
    all(fitted$iter < config$path_max_iterations), !as.logical(fitted$saturated),
    max(abs(fitted$lambda - lambda)) < 1e-12)
  list(beta = fitted$beta, lambda = lambda, alpha = alpha, iterations = fitted$iter)
}

cq_sparse <- function(value) as(as(Matrix::Matrix(value, sparse = TRUE), "generalMatrix"), "CsparseMatrix")

cq_exact <- function(features, tau, lambda, alpha, config) {
  X <- features$X
  N <- nrow(X)
  p <- ncol(X)
  penalized <- which(features$penalty_factor == 1)
  d <- length(penalized)
  weight <- features$weight / features$n
  select <- sparseMatrix(i = seq_len(d), j = penalized, x = 1, dims = c(d, p))
  equality <- cbind(cq_sparse(X), Diagonal(N), -Diagonal(N), Matrix(0, N, d, sparse = TRUE))
  nonnegative <- cbind(Matrix(0, 2L * N, p, sparse = TRUE), -Diagonal(2L * N),
    Matrix(0, 2L * N, d, sparse = TRUE))
  positive <- cbind(select, Matrix(0, d, 2L * N, sparse = TRUE), -Diagonal(d))
  negative <- cbind(-select, Matrix(0, d, 2L * N, sparse = TRUE), -Diagonal(d))
  A <- cq_sparse(rbind(equality, nonnegative, positive, negative))
  P <- cq_sparse(Diagonal(p + 2L * N + d,
    c(lambda * (1 - alpha) * features$penalty_factor, numeric(2L * N + d))))
  linear <- c(numeric(p), tau * weight, (1 - tau) * weight, rep(lambda * alpha, d))
  settings <- clarabel_control(max_iter = 300L, verbose = FALSE,
    tol_gap_abs = 1e-10, tol_gap_rel = 1e-10, tol_feas = 1e-10)
  messages <- capture.output(fitted <- clarabel(A = A, b = c(features$y, numeric(2L * N + 2L * d)),
    q = linear, P = P, cones = list(z = as.integer(N), l = as.integer(2L * N + 2L * d)), control = settings))
  status <- names(solver_status_descriptions())[fitted$status]
  stopifnot(status == "Solved")
  beta <- setNames(fitted$x[seq_len(p)], colnames(X))
  dual <- -fitted$z[seq_len(N)]
  gradient <- drop(crossprod(X, dual))
  stationarity <- drop(P %*% fitted$x + linear + crossprod(A, fitted$z))
  box_error <- max(0, dual - tau * weight, (tau - 1) * weight - dual)
  residual_error <- max(abs(drop(equality %*% fitted$x) - features$y))
  objective <- cq_objective(features, beta, tau, lambda, alpha)
  primal_error <- abs(objective - fitted$obj_val)
  gap <- fitted$obj_val - fitted$obj_val_dual
  stopifnot(max(abs(stationarity), box_error, residual_error, primal_error, abs(gap)) <
    config$exact_fit_tolerance * max(1, abs(objective)))
  list(beta = beta, lambda = lambda, alpha = alpha, tau = tau, objective = objective,
    dual = dual, diagnostics = list(status = status, iterations = fitted$iterations,
      gap = gap, stationarity = max(abs(stationarity)), box_error = box_error,
      residual_error = residual_error, primal_error = primal_error, messages = messages))
}

cq_fit <- function(rows, evaluation, panel, recipes, study, config, tau, directory) {
  stopifnot(!length(intersect(rows$strain_id, evaluation$strain_id)))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  grid <- data.table::CJ(alpha = unlist(config$alpha), index = seq_len(config$path_length), sorted = FALSE)
  total_loss <- numeric(nrow(grid))
  oof <- matrix(NA_real_, nrow(rows), nrow(grid))
  for (fold in sort(unique(rows$inner_fold))) {
    recipe <- recipes[[paste0("inner_", fold)]]$direction_recipe$base_recipe
    train_rows <- rows[inner_fold != fold]
    test_rows <- rows[inner_fold == fold]
    stopifnot(identical(sort(unique(train_rows$strain_id)), recipe$training_strains))
    training <- cq_features(train_rows, panel, recipe)
    testing <- cq_features(test_rows, panel, recipe, colnames(training$X))
    for (alpha in unlist(config$alpha)) {
      cat(study$study_id, tau, "inner", fold, "alpha", alpha, "\n")
      flush.console()
      fit <- cq_path(training, tau, alpha, config)
      predictions <- ep_predict(testing, fit$beta)
      losses <- ep_loss(test_rows$y, predictions, tau)
      indices <- which(grid$alpha == alpha)
      total_loss[indices] <- total_loss[indices] + colSums(losses * testing$weight)
      oof[rows$inner_fold == fold, indices] <- predictions
      saveRDS(list(path = fit, columns = colnames(training$X), recipe = recipe),
        file.path(directory, sprintf("inner_%d_alpha_%s.rds", fold, alpha)))
    }
  }
  stopifnot(all(is.finite(oof)))
  choice <- which.min(total_loss)
  setting <- grid[choice]
  recipe <- recipes$training$direction_recipe$base_recipe
  stopifnot(identical(sort(unique(rows$strain_id)), recipe$training_strains))
  training <- cq_features(rows, panel, recipe)
  testing <- cq_features(evaluation, panel, recipe, colnames(training$X))
  path <- cq_path(training, tau, setting$alpha, config)
  final <- cq_exact(training, tau, path$lambda[setting$index], setting$alpha, config)
  path_objective <- cq_objective(training, path$beta[, setting$index], tau, final$lambda, final$alpha)
  objective_difference <- path_objective - final$objective
  stopifnot(objective_difference >= -1e-8,
    objective_difference <= config$maximum_path_objective_difference * max(1, abs(final$objective)))
  prediction <- drop(ep_predict(testing, final$beta))
  fit <- c(final, list(recipe = recipe, columns = colnames(training$X), grid = grid,
    inner_loss = total_loss / uniqueN(rows$strain_id), choice = choice, oof = oof[, choice],
    prediction = prediction, path_objective_difference = objective_difference,
    training_patients = uniqueN(rows$strain_id), test_patients = uniqueN(evaluation$strain_id),
    clinical_columns = names(final$beta)[training$penalty_factor == 0]))
  saveRDS(fit, file.path(directory, "fit.rds"))
  fit
}

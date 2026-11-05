library(glmnet)

cq_density <- function(features, lower_beta, upper_beta, delta, config) {
  difference <- drop(features$X %*% (upper_beta - lower_beta))
  limits <- unlist(config$density_range)
  floor <- 2 * delta / limits[2L]
  ceiling <- 2 * delta / limits[1L]
  density <- 2 * delta / pmax(floor, pmin(ceiling, difference))
  list(value = density, lower_bound_count = sum(difference <= floor),
    upper_bound_count = sum(difference >= ceiling), crossing_count = sum(difference <= 0))
}

cq_transform_blocks <- function(values, cluster, tau, rho, inverse = FALSE) {
  stopifnot(rho >= 0, rho < 1, tau > 0, tau < 1)
  values <- as.matrix(values)
  result <- matrix(NA_real_, nrow(values), ncol(values), dimnames = dimnames(values))
  blocks <- split(seq_along(cluster), cluster)
  variance <- tau * (1 - tau)
  power <- if (inverse) 1 else 0.5
  for (indices in blocks) {
    current <- values[indices, , drop = FALSE]
    average <- matrix(colMeans(current), length(indices), ncol(current), byrow = TRUE)
    result[indices, ] <- (current - average) / (variance * (1 - rho))^power +
      average / (variance * (1 - rho + length(indices) * rho))^power
  }
  stopifnot(all(is.finite(result)))
  result
}

cq_score_correlation <- function(rows, prediction, tau, config) {
  score <- tau - as.numeric(rows$y <= prediction)
  blocks <- split(seq_len(nrow(rows)), rows$strain_id)
  repeated <- blocks[lengths(blocks) > 1L]
  stopifnot(length(repeated) >= 20L)
  products <- vapply(repeated, function(indices) {
    values <- score[indices]
    (sum(values)^2 - sum(values^2)) / (length(values) * (length(values) - 1L))
  }, numeric(1))
  raw <- mean(products) / (tau * (1 - tau))
  limits <- unlist(config$score_correlation_range)
  list(raw = raw, working = min(limits[2L], max(limits[1L], raw)),
    repeated_patients = length(repeated))
}

cq_projection <- function(features, density, target, rho, tau, fold_id, config) {
  k <- match(target, colnames(features$X))
  stopifnot(!is.na(k), all(is.finite(density)), all(density > 0))
  transformed <- cq_transform_blocks(features$X * density, features$cluster_id, tau, rho)
  d <- transformed[, k]
  X <- transformed[, -k, drop = FALSE]
  penalty <- features$penalty_factor[-k]
  clinical <- which(penalty == 0)
  penalized <- which(penalty == 1)
  decomposition <- qr(X[, clinical, drop = FALSE], tol = 1e-10)
  stopifnot(decomposition$rank == length(clinical))
  target_residual <- drop(qr.resid(decomposition, d))
  predictor_residual <- qr.resid(decomposition, X[, penalized, drop = FALSE])
  active <- which(colSums(predictor_residual^2) / nrow(X) > 1e-12)
  stopifnot(length(active) >= 2L)
  glmnet.control(fdev = 0, devmax = 1)
  fit <- cv.glmnet(predictor_residual[, active, drop = FALSE], target_residual,
    family = "gaussian", alpha = 1, intercept = FALSE, standardize = FALSE,
    nlambda = config$projection_path_length, lambda.min.ratio = config$projection_lambda_min_ratio,
    foldid = fold_id, type.measure = "mse", thresh = 1e-12, maxit = 1000000L,
    grouped = FALSE, parallel = FALSE)
  stopifnot(fit$glmnet.fit$jerr == 0L, length(fit$lambda) == config$projection_path_length)
  selected <- which.min(fit$cvm)
  gamma <- setNames(numeric(ncol(X)), colnames(X))
  gamma[penalized[active]] <- as.numeric(fit$glmnet.fit$beta[, selected])
  gamma[clinical] <- drop(qr.coef(decomposition, d - drop(X[, penalized, drop = FALSE] %*% gamma[penalized])))
  residual <- d - drop(X %*% gamma)
  sensitivity <- sum(residual * d) / features$n
  quadratic <- sum(residual^2) / features$n
  balance <- drop(crossprod(X, residual)) / features$n
  stopifnot(is.finite(sensitivity), sensitivity > 1e-8,
    max(abs(balance[clinical])) < 1e-7)
  list(gamma = gamma, sensitivity = sensitivity, residual_quadratic = quadratic,
    maximum_balance = max(abs(balance)), clinical_balance = max(abs(balance[clinical])),
    lambda = fit$lambda[selected], selected_index = selected,
    lambda_grid = fit$lambda, cv_loss = fit$cvm, rho = rho, target = target, tau = tau)
}

cq_influence <- function(features, beta, density, projection) {
  k <- match(projection$target, colnames(features$X))
  nuisance <- features$X[, names(projection$gamma), drop = FALSE]
  residual <- features$X[, k] - drop(nuisance %*% projection$gamma)
  score <- projection$tau - as.numeric(features$y <= drop(features$X %*% beta))
  adjusted <- drop(cq_transform_blocks(score, features$cluster_id, projection$tau,
    projection$rho, inverse = TRUE))
  contribution <- residual * density * adjusted
  influence <- rowsum(contribution, features$cluster_id, reorder = TRUE)[, 1L] / projection$sensitivity
  multiplier <- features$recipe$y_scale / features$recipe$scale[projection$target]
  list(initial = unname(beta[k] * multiplier), influence = influence * multiplier,
    patient = names(influence), multiplier = unname(multiplier))
}

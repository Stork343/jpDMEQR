ep_ridge_press <- function(features, config) {
  low_columns <- which(features$penalty_factor == 0)
  high_columns <- which(features$penalty_factor == 1)
  stopifnot(identical(low_columns, seq_along(low_columns)), all(features$X[, 1L] == 1))
  N <- nrow(features$X)
  p <- length(high_columns)
  root_weight <- sqrt(N * features$weight / sum(features$weight))
  low <- features$X[, low_columns, drop = FALSE] * root_weight
  decomposition <- qr(low)
  stopifnot(decomposition$rank == length(low_columns))
  Q <- qr.Q(decomposition)
  weighted_X <- features$X[, high_columns, drop = FALSE] * root_weight
  low_X <- qr.coef(decomposition, weighted_X)
  weighted_y <- features$y * root_weight
  low_y <- drop(qr.coef(decomposition, weighted_y))
  residual_X <- qr.resid(decomposition, weighted_X)
  residual_y <- drop(qr.resid(decomposition, weighted_y))
  kernel <- tcrossprod(residual_X) / N
  eig <- eigen(kernel, symmetric = TRUE)
  tolerance <- 10 * max(N, p) * .Machine$double.eps * max(1, eig$values[1L])
  stopifnot(min(eig$values) >= -tolerance)
  values <- pmax(eig$values, 0)
  lambda <- p / features$n * exp(seq(log(config$ridge_variance_ratio_max),
    log(config$ridge_variance_ratio_min), length.out = config$elastic_net_path_length))
  filters <- sweep(1 / outer(values, lambda, "+"), 1L, values, "*")
  projected_y <- drop(crossprod(eig$vectors, residual_y))
  weighted_prediction <- drop(qr.fitted(decomposition, weighted_y)) +
    eig$vectors %*% sweep(filters, 1L, projected_y, "*")
  residual <- weighted_y - weighted_prediction
  groups <- split(seq_len(N), features$cluster_id)
  pairs <- do.call(rbind, lapply(groups, function(indices) expand.grid(row = indices, column = indices)))
  pair_low <- rowSums(Q[pairs$row, , drop = FALSE] * Q[pairs$column, , drop = FALSE])
  pair_eigen <- eig$vectors[pairs$row, , drop = FALSE] * eig$vectors[pairs$column, , drop = FALSE]
  hat_pairs <- pair_low + pair_eigen %*% filters
  leave_out <- matrix(NA_real_, N, length(lambda))
  cursor <- 0L
  minimum_margin <- Inf
  for (indices in groups) {
    m <- length(indices)
    stopifnot(max(root_weight[indices]) - min(root_weight[indices]) < 1e-12)
    positions <- cursor + seq_len(m * m)
    for (j in seq_along(lambda)) {
      complement <- diag(m) - matrix(hat_pairs[positions, j], m, m)
      minimum_margin <- min(minimum_margin, min(eigen(complement, symmetric = TRUE, only.values = TRUE)$values))
      leave_out[indices, j] <- features$y[indices] -
        solve(complement, residual[indices, j]) / root_weight[indices]
    }
    cursor <- cursor + m * m
  }
  stopifnot(all(is.finite(leave_out)), minimum_margin > 0)
  loss <- colSums((features$y - leave_out)^2 * features$weight) / features$n
  choice <- which.min(loss)
  dual <- drop(eig$vectors %*% (projected_y / (values + lambda[choice])))
  high_beta <- drop(crossprod(residual_X, dual)) / N
  beta <- c(low_y - drop(low_X %*% high_beta), high_beta)
  names(beta) <- colnames(features$X)
  error <- max(abs(kernel %*% dual + lambda[choice] * dual - residual_y))
  stopifnot(error < 1e-8 * max(1, max(abs(residual_y))))
  list(beta = beta, leave_out_prediction = leave_out[, choice], choice = choice, lambda = lambda,
       loss = loss, leave_out_path = leave_out, minimum_margin = minimum_margin, system_error = error,
       penalty_scale = "p_over_cluster_count_times_variance_ratio")
}

ep_direction_recipe <- function(rows, panel, config, base_recipe_function, base_feature_function) {
  base_recipe <- base_recipe_function(rows, panel, config)
  base <- base_feature_function(rows, panel, base_recipe)
  background <- ep_ridge_press(base, config)
  condition <- cbind(mean = 1, base_recipe$contrasts[rows$condition, , drop = FALSE])
  directions <- condition * background$leave_out_prediction
  low_columns <- which(base$penalty_factor == 0)
  root_weight <- sqrt(base$weight)
  decomposition <- qr(base$X[, low_columns, drop = FALSE] * root_weight)
  projection <- qr.coef(decomposition, directions * root_weight)
  residual <- directions - base$X[, low_columns, drop = FALSE] %*% projection
  scale <- sqrt(colSums(residual^2 * base$weight) / base$n)
  stopifnot(all(is.finite(scale)), all(scale > 1e-8))
  names(background$leave_out_prediction) <- rows$sample_id
  background$leave_out_path <- NULL
  list(base_recipe = base_recipe, background = background, projection = projection, direction_scale = scale,
       low_columns = low_columns, y_center = base_recipe$y_center, y_scale = base_recipe$y_scale,
       training_strains = base_recipe$training_strains,
       search_dimension = base_recipe$search_dimension + ncol(directions))
}

ep_direction_features <- function(rows, panel, recipe, base_feature_function) {
  features <- base_feature_function(rows, panel, recipe$base_recipe)
  in_training <- rows$sample_id %in% names(recipe$background$leave_out_prediction)
  stopifnot(all(in_training) || !any(in_training))
  score <- if (all(in_training)) recipe$background$leave_out_prediction[rows$sample_id] else
    drop(features$X %*% recipe$background$beta)
  condition <- cbind(mean = 1, recipe$base_recipe$contrasts[rows$condition, , drop = FALSE])
  directions <- condition * score - features$X[, recipe$low_columns, drop = FALSE] %*% recipe$projection
  directions <- sweep(directions, 2L, recipe$direction_scale, "/")
  colnames(directions) <- paste0("RIDGE_", colnames(condition))
  features$X <- cbind(features$X, directions)
  features$penalty_factor <- c(features$penalty_factor, rep(1, ncol(directions)))
  features$recipe <- recipe
  stopifnot(all(is.finite(features$X)))
  features
}

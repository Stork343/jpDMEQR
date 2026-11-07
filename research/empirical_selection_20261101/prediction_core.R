ep_loss <- function(y, prediction, tau) {
  residual <- y - prediction
  residual * (tau - (residual < 0))
}

ep_predict <- function(features, beta) {
  features$X %*% beta * features$recipe$y_scale + features$recipe$y_center
}

ep_mean_path <- function(features, alpha, config, lambda = NULL) {
  unpenalized <- which(features$penalty_factor == 0)
  penalized <- which(features$penalty_factor == 1)
  stopifnot(identical(unpenalized, seq_along(unpenalized)), length(penalized) > 0L,
            length(unpenalized) + length(penalized) == ncol(features$X),
            all(features$X[, 1L] == 1), all(features$weight > 0))
  N <- nrow(features$X)
  p <- length(penalized)
  model_variables <- ncol(features$X) - 1L
  if (alpha > 0) {
    native_lambda <- if (is.null(lambda)) NULL else lambda * p / model_variables
    fit <- glmnet::glmnet(features$X[, -1L, drop = FALSE], features$y,
      weights = features$weight, family = "gaussian", alpha = alpha,
      lambda = native_lambda, nlambda = as.integer(config$elastic_net_path_length),
      lambda.min.ratio = config$elastic_net_lambda_min_ratio,
      penalty.factor = features$penalty_factor[-1L], intercept = TRUE, standardize = FALSE,
      thresh = 1e-13, maxit = 1000000L)
    stopifnot(fit$jerr == 0L, ncol(fit$beta) == config$elastic_net_path_length)
    return(list(beta = stats::coef(fit), lambda = fit$lambda * model_variables / p,
      diagnostics = list(solver = "glmnet", alpha = alpha, jerr = fit$jerr,
                         native_lambda = fit$lambda, passes = fit$npasses)))
  }
  stopifnot(alpha == 0)
  weight <- features$weight / sum(features$weight)
  root_weight <- sqrt(N * weight)
  low <- features$X[, unpenalized, drop = FALSE] * root_weight
  decomposition <- qr(low)
  stopifnot(decomposition$rank == length(unpenalized))
  weighted_X <- features$X[, penalized, drop = FALSE] * root_weight
  low_X <- qr.coef(decomposition, weighted_X)
  low_y <- drop(qr.coef(decomposition, features$y * root_weight))
  weighted_X <- qr.resid(decomposition, weighted_X)
  weighted_y <- drop(qr.resid(decomposition, features$y * root_weight))
  if (is.null(lambda)) {
    if (!is.null(config$ridge_variance_ratio_max)) {
      lambda <- p / features$n * exp(seq(log(config$ridge_variance_ratio_max),
        log(config$ridge_variance_ratio_min), length.out = config$elastic_net_path_length))
    } else {
      maximum <- max(abs(crossprod(weighted_X, weighted_y))) / (N * 0.001)
      stopifnot(is.finite(maximum), maximum > 0)
      lambda <- maximum * exp(seq(0, log(config$elastic_net_lambda_min_ratio),
                                   length.out = config$elastic_net_path_length))
    }
  }
  stopifnot(length(lambda) == config$elastic_net_path_length, all(is.finite(lambda)), all(lambda > 0))
  primal <- p < N
  gram <- if (primal) crossprod(weighted_X) / N else tcrossprod(weighted_X) / N
  rhs <- if (primal) drop(crossprod(weighted_X, weighted_y)) / N else weighted_y
  decomposition <- eigen(gram, symmetric = TRUE)
  eigen_tolerance <- 10 * max(N, p) * .Machine$double.eps * max(1, decomposition$values[1L])
  stopifnot(min(decomposition$values) >= -eigen_tolerance)
  values <- pmax(decomposition$values, 0)
  projected_y <- drop(crossprod(decomposition$vectors, rhs))
  inverse_scores <- outer(projected_y, rep(1, length(lambda))) / outer(values, lambda, "+")
  solution <- decomposition$vectors %*% inverse_scores
  system_error <- apply(abs(gram %*% solution + sweep(solution, 2L, lambda, "*") - rhs), 2L, max)
  stopifnot(max(system_error) < 1e-8 * max(1, max(abs(rhs))))
  high_beta <- if (primal) solution else crossprod(weighted_X, solution) / N
  low_beta <- matrix(low_y, length(low_y), length(lambda)) - low_X %*% high_beta
  beta <- rbind(Matrix::Matrix(low_beta, sparse = TRUE), Matrix::Matrix(high_beta, sparse = TRUE))
  rownames(beta) <- colnames(features$X)
  colnames(beta) <- paste0("s", seq_along(lambda) - 1L)
  list(beta = beta, lambda = lambda,
       diagnostics = list(solver = if (primal) "weighted_qr_primal_ridge" else "weighted_qr_kernel_ridge",
         alpha = alpha, lambda = lambda, linear_system_error = system_error,
         minimum_kernel_eigenvalue = if (primal) 0 else min(decomposition$values),
         system_dimension = nrow(gram), minimum_system_eigenvalue = min(decomposition$values)))
}

ep_path <- function(features, method, tau, config, maximum, setting = NULL) {
  if (method == "JOINT") {
    lambdas <- if (is.null(setting)) unlist(config$working_penalties) else setting$Lambda
    size <- if (is.null(setting)) maximum else setting$size
    fits <- lapply(lambdas, function(lambda) {
      prepared <- pc_prepare_initial(features$y, features$X, features$Z, features$cluster_id,
        tau, diag(lambda, ncol(features$Z)), features$penalty_factor, fit_always_included = FALSE)
      pc_orthogonal_path(prepared, tau, 0:size, solver = jp_joint_qp, gap_tolerance = config$gap_tolerance)
    })
    grid <- rbindlist(lapply(lambdas, function(lambda) data.table(Lambda = lambda, size = 0:size)))
    beta <- do.call(cbind, lapply(fits, `[[`, "coefficients"))
    diagnostics <- lapply(fits, `[[`, "diagnostics")
    chosen <- if (is.null(setting)) seq_len(ncol(beta)) else ncol(beta)
  } else if (method == "ENET") {
    glmnet::glmnet.control(fdev = 0, devmax = 1)
    alphas <- if (is.null(setting)) unlist(config$elastic_net_alpha) else setting$alpha
    fits <- lapply(alphas, function(alpha) ep_mean_path(features, alpha, config))
    grid <- rbindlist(lapply(alphas, function(alpha)
      data.table(alpha = alpha, penalty_index = seq_len(config$elastic_net_path_length))))
    beta <- do.call(cbind, lapply(fits, `[[`, "beta"))
    diagnostics <- lapply(fits, `[[`, "diagnostics")
    chosen <- if (is.null(setting)) seq_len(ncol(beta)) else setting$penalty_index
  } else {
    stopifnot(method == "POOLED")
    row_weight <- nrow(features$X) * features$weight / features$n
    fit <- soy_pooled_lp_cpp(features$X * row_weight, features$y * row_weight,
      features$penalty_factor, tau, 40L, 0.01, numeric(), TRUE,
      maximum_solutions = if (is.null(setting)) 0L else as.integer(setting$penalty_index))
    grid <- data.table(penalty_index = seq_len(ncol(fit$beta)))
    beta <- fit$beta
    diagnostics <- fit[c("primal_dual_gap", "dual_violation", "iterations")]
    chosen <- if (is.null(setting)) seq_len(ncol(beta)) else setting$penalty_index
  }
  list(beta = beta[, chosen, drop = FALSE], grid = grid[chosen], diagnostics = diagnostics)
}

ep_calibrate <- function(rows, prediction, tau, cluster_column, calibration_column) {
  groups <- if (is.null(calibration_column)) rep("all", nrow(rows)) else as.character(rows[[calibration_column]])
  levels <- sort(unique(groups))
  X <- vapply(levels, function(value) as.numeric(groups == value), numeric(nrow(rows)))
  ids <- rows[[cluster_column]]
  weight <- 1 / as.numeric(table(ids)[ids])
  fit <- quantreg::rq.fit(X * weight, (rows$y - prediction) * weight, tau = tau, method = "br")
  stopifnot(all(is.finite(fit$coefficients)))
  setNames(as.numeric(fit$coefficients), levels)
}

ep_fit_prediction <- function(rows, evaluation, panel, method, tau, config, progress, directory,
                              recipe_function, feature_function, cluster_column, calibration_column = NULL,
                              prediction_function = ep_predict, calibration_mode = "oof") {
  stopifnot(calibration_mode %in% c("oof", "training"))
  if (calibration_mode == "training") stopifnot(method != "ENET")
  folds <- sort(unique(rows$inner_fold))
  recipes <- lapply(folds, function(fold) recipe_function(rows[inner_fold != fold], panel, config))
  counts <- vapply(folds, function(fold) uniqueN(rows[inner_fold != fold][[cluster_column]]), integer(1))
  maximum <- min(vapply(seq_along(folds), function(j)
    min(recipes[[j]]$search_dimension, floor(counts[j] / log(counts[j]))), numeric(1)))
  if (method == "JOINT") {
    ranks <- vapply(seq_along(folds), function(j) {
      features <- feature_function(rows[inner_fold != folds[j]], panel, recipes[[j]])
      weighted <- features$X * sqrt(features$weight / features$n)
      qr(weighted, tol = 1e-5)$rank - sum(features$penalty_factor == 0)
    }, integer(1))
    stopifnot(all(ranks >= 0L))
    final_recipe <- recipe_function(rows, panel, config)
    final_design <- feature_function(rows, panel, final_recipe)
    final_rank <- qr(final_design$X * sqrt(final_design$weight/final_design$n), tol=1e-5)$rank
    maximum <- min(maximum, ranks, final_rank-sum(final_design$penalty_factor==0))
  }
  loss <- NULL
  oof <- NULL
  for (f in seq_along(folds)) {
    training <- feature_function(rows[inner_fold != folds[f]], panel, recipes[[f]])
    validation <- feature_function(rows[inner_fold == folds[f]], panel, recipes[[f]])
    progress("inner_fit", list(fold = folds[f], clusters = training$n,
      variables = ncol(training$X), additions = maximum))
    gc()
    fit <- ep_path(training, method, tau, config, maximum)
    predicted <- as.matrix(prediction_function(validation, fit$beta))
    if (calibration_mode == "training") {
      fitted_prediction <- as.matrix(prediction_function(training, fit$beta))
      validation_groups <- if (is.null(calibration_column)) rep("all", nrow(predicted)) else
        as.character(rows[inner_fold == folds[f]][[calibration_column]])
      for (j in seq_len(ncol(predicted))) {
        offset <- ep_calibrate(rows[inner_fold != folds[f]], fitted_prediction[, j], tau,
                               cluster_column, calibration_column)
        predicted[, j] <- predicted[, j] + offset[validation_groups]
      }
    }
    if (f == 1L) {
      grid <- fit$grid
      oof <- matrix(NA_real_, nrow(rows), ncol(predicted))
      loss <- numeric(ncol(predicted))
    }
    stopifnot(identical(as.list(grid), as.list(fit$grid)), ncol(predicted) == ncol(oof))
    values <- if (method == "ENET") (validation$raw_y - predicted)^2 else ep_loss(validation$raw_y, predicted, tau)
    loss <- loss + colSums(values * validation$weight)
    oof[which(rows$inner_fold == folds[f]), ] <- predicted
    saveRDS(list(path = fit, recipe = recipes[[f]], prediction = predicted),
            file.path(directory, paste0("inner_", folds[f], ".rds")))
    rm(training, validation, fit, predicted, values)
    gc()
  }
  stopifnot(all(is.finite(oof)))
  choice <- which.min(loss)
  recipe <- recipe_function(rows, panel, config)
  training <- feature_function(rows, panel, recipe)
  progress("training_fit", list(clusters = training$n, selected_index = choice))
  final <- ep_path(training, method, tau, config, maximum, as.list(grid[choice]))
  stopifnot(ncol(final$beta) == 1L)
  beta <- as.numeric(final$beta[, 1L])
  final_diagnostics <- final$diagnostics
  names(beta) <- colnames(training$X)
  quantiles <- if (method == "ENET") unlist(config$quantiles) else tau
  calibration_prediction <- if (calibration_mode == "oof") oof[, choice] else
    drop(prediction_function(training, beta))
  calibration <- lapply(quantiles, function(value)
    ep_calibrate(rows, calibration_prediction, value, cluster_column, calibration_column))
  rm(training, final)
  gc()
  test <- feature_function(evaluation, panel, recipe)
  groups <- if (is.null(calibration_column)) rep("all", nrow(evaluation)) else as.character(evaluation[[calibration_column]])
  base <- drop(prediction_function(test, beta))
  prediction <- vapply(calibration, function(offset) base + offset[groups], numeric(nrow(evaluation)))
  stopifnot(all(is.finite(prediction)))
  list(beta = beta, recipe = recipe, prediction = prediction, quantiles = quantiles,
       calibration = calibration, choice = choice, grid = grid, loss = loss,
       oof = oof[, choice], final_diagnostics = final_diagnostics, calibration_mode = calibration_mode,
       maximum_search_size = maximum)
}

ep_loss <- function(y, prediction, tau) {
  residual <- y - prediction
  residual * (tau - (residual < 0))
}

ep_predict <- function(features, beta) {
  features$X %*% beta * features$recipe$y_scale + features$recipe$y_center
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
    fits <- lapply(alphas, function(alpha) {
      lambda <- NULL
      if (alpha == 0 && !is.null(config$ridge_variance_ratio_max)) {
        penalized <- sum(features$penalty_factor > 0)
        lambda <- penalized^2 / ((ncol(features$X) - 1L) * features$n) * exp(seq(
          log(config$ridge_variance_ratio_max), log(config$ridge_variance_ratio_min),
          length.out = as.integer(config$elastic_net_path_length)))
      }
      stopifnot(all(features$X[, 1L] == 1), features$penalty_factor[1L] == 0)
      fit <- glmnet::glmnet(features$X[, -1L, drop = FALSE], features$y, weights = features$weight,
        family = "gaussian", alpha = alpha, nlambda = as.integer(config$elastic_net_path_length),
        lambda = lambda,
        lambda.min.ratio = config$elastic_net_lambda_min_ratio, penalty.factor = features$penalty_factor[-1L],
        intercept = TRUE, standardize = FALSE, thresh = 1e-13, maxit = 1000000L)
      stopifnot(fit$jerr == 0L, ncol(fit$beta) == config$elastic_net_path_length)
      fit
    })
    grid <- rbindlist(lapply(alphas, function(alpha)
      data.table(alpha = alpha, penalty_index = seq_len(config$elastic_net_path_length))))
    beta <- do.call(cbind, lapply(fits, stats::coef))
    diagnostics <- lapply(fits, function(fit) list(jerr = fit$jerr, lambda = fit$lambda))
    chosen <- if (is.null(setting)) seq_len(ncol(beta)) else setting$penalty_index
  } else {
    stopifnot(method == "POOLED")
    row_weight <- nrow(features$X) * features$weight / features$n
    fit <- soy_pooled_lp_cpp(features$X * row_weight, features$y * row_weight,
      features$penalty_factor, tau, 40L, 0.01, numeric(), TRUE)
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
                              recipe_function, feature_function, cluster_column, calibration_column = NULL) {
  folds <- sort(unique(rows$inner_fold))
  recipes <- lapply(folds, function(fold) recipe_function(rows[inner_fold != fold], panel, config))
  counts <- vapply(folds, function(fold) uniqueN(rows[inner_fold != fold][[cluster_column]]), integer(1))
  maximum <- min(vapply(seq_along(folds), function(j)
    max(0, floor(counts[j] / log(counts[j])) - length(recipes[[j]]$low_columns)), numeric(1)))
  loss <- NULL
  oof <- NULL
  for (f in seq_along(folds)) {
    training <- feature_function(rows[inner_fold != folds[f]], panel, recipes[[f]])
    validation <- feature_function(rows[inner_fold == folds[f]], panel, recipes[[f]])
    progress("inner_fit", list(fold = folds[f], clusters = training$n,
      variables = ncol(training$X), additions = maximum))
    gc()
    fit <- ep_path(training, method, tau, config, maximum)
    predicted <- as.matrix(ep_predict(validation, fit$beta))
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
  calibration <- lapply(quantiles, function(value)
    ep_calibrate(rows, oof[, choice], value, cluster_column, calibration_column))
  rm(training, final)
  gc()
  test <- feature_function(evaluation, panel, recipe)
  groups <- if (is.null(calibration_column)) rep("all", nrow(evaluation)) else as.character(evaluation[[calibration_column]])
  base <- drop(ep_predict(test, beta))
  prediction <- vapply(calibration, function(offset) base + offset[groups], numeric(nrow(evaluation)))
  stopifnot(all(is.finite(prediction)))
  list(beta = beta, recipe = recipe, prediction = prediction, quantiles = quantiles,
       calibration = calibration, choice = choice, grid = grid, loss = loss,
       oof = oof[, choice], final_diagnostics = final_diagnostics)
}

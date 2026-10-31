yeast_recipe <- function(rows, panel) {
  ids <- sort(unique(rows$strain_id))
  X <- panel$X[ids, , drop = FALSE]
  center <- colMeans(X)
  scale <- sqrt(colSums(sweep(X, 2L, center, "-")^2) / (length(ids) - 1L))
  keep <- which(is.finite(scale) & scale > 1e-8)
  levels <- sort(unique(rows$replicate))
  replicate_design <- vapply(levels[-1L], function(value) as.numeric(rows$replicate == value), numeric(nrow(rows)))
  if (length(levels) > 1L) colnames(replicate_design) <- paste0("replicate", levels[-1L])
  low <- cbind(panel$C[rows$strain_id, , drop = FALSE], replicate_design)
  decomposition <- qr(low)
  low_columns <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  weights <- 1 / as.numeric(table(rows$strain_id)[rows$strain_id])
  y_center <- sum(rows$y * weights) / length(ids)
  y_scale <- sqrt(sum((rows$y - y_center)^2 * weights) / length(ids))
  stopifnot(is.finite(y_scale), y_scale > 0, low_columns[1L] == 1L)
  list(variables = colnames(X)[keep], center = center[keep], scale = scale[keep],
       replicate_levels = levels, low_columns = low_columns, y_center = y_center,
       y_scale = y_scale, training_strains = ids)
}

yeast_features <- function(rows, panel, recipe) {
  high <- sweep(sweep(panel$X[rows$strain_id, recipe$variables, drop = FALSE],
                     2L, recipe$center, "-"), 2L, recipe$scale, "/")
  replicate_design <- vapply(recipe$replicate_levels[-1L],
    function(value) as.numeric(rows$replicate == value), numeric(nrow(rows)))
  if (length(recipe$replicate_levels) > 1L) colnames(replicate_design) <- paste0("replicate", recipe$replicate_levels[-1L])
  low <- cbind(panel$C[rows$strain_id, , drop = FALSE], replicate_design)[, recipe$low_columns, drop = FALSE]
  X <- cbind(low, high)
  stopifnot(all(is.finite(X)), all(X[, 1L] == 1), !anyDuplicated(colnames(X)))
  list(X = X, Z = matrix(1, nrow(rows), 1L), y = (rows$y - recipe$y_center) / recipe$y_scale,
       raw_y = rows$y, cluster_id = rows$strain_id, n = length(unique(rows$strain_id)),
       weight = 1 / as.numeric(table(rows$strain_id)[rows$strain_id]),
       penalty_factor = c(rep(0, ncol(low)), rep(1, ncol(high))), recipe = recipe)
}

yeast_raw_prediction <- function(features, beta) {
  features$X %*% beta * features$recipe$y_scale + features$recipe$y_center
}

yeast_loss <- function(y, prediction, tau) {
  residual <- y - prediction
  residual * (tau - (residual < 0))
}

yeast_calibrate <- function(rows, oof, tau, replicate_levels) {
  weights <- 1 / as.numeric(table(rows$strain_id)[rows$strain_id])
  X <- vapply(replicate_levels, function(value) as.numeric(rows$replicate == value), numeric(nrow(rows)))
  fit <- quantreg::rq.fit(X * weights, (rows$y - oof) * weights, tau = tau, method = "br")
  stopifnot(all(is.finite(fit$coefficients)))
  setNames(as.numeric(fit$coefficients), replicate_levels)
}

yeast_profile_prediction <- function(rows, evaluation, panel, config, tau, progress, directory) {
  folds <- sort(unique(rows$inner_fold))
  inner <- lapply(folds, function(fold) {
    recipe <- yeast_recipe(rows[inner_fold != fold], panel)
    list(train = yeast_features(rows[inner_fold != fold], panel, recipe),
         validation = yeast_features(rows[inner_fold == fold], panel, recipe))
  })
  maximum <- min(vapply(inner, function(fold) max(0, floor(fold$train$n / log(fold$train$n)) -
                                                  sum(fold$train$penalty_factor == 0)), numeric(1)))
  sizes <- 0:as.integer(maximum)
  lambdas <- unlist(config$working_penalties)
  loss <- matrix(0, length(lambdas), length(sizes))
  oof <- lapply(lambdas, function(value) matrix(NA_real_, nrow(rows), length(sizes)))
  certificates <- list()
  for (f in seq_along(folds)) {
    training <- inner[[f]]$train
    validation <- inner[[f]]$validation
    for (l in seq_along(lambdas)) {
      progress("profile_inner", list(fold = folds[f], Lambda = lambdas[l],
               strains = training$n, variables = ncol(training$X), additions = maximum))
      prepared <- pc_prepare_initial(training$y, training$X, training$Z, training$cluster_id,
        tau, matrix(lambdas[l], 1L, 1L), training$penalty_factor, fit_always_included = FALSE)
      path <- pc_orthogonal_path(prepared, tau, sizes, solver = jp_joint_qp, gap_tolerance = config$gap_tolerance)
      predicted <- yeast_raw_prediction(validation, path$coefficients)
      loss[l, ] <- loss[l, ] + colSums(yeast_loss(validation$raw_y, predicted, tau) * validation$weight)
      oof[[l]][which(rows$inner_fold == folds[f]), ] <- predicted
      certificates[[length(certificates) + 1L]] <- path$diagnostics
      saveRDS(list(path = path, recipe = training$recipe, prediction = predicted),
              file.path(directory, sprintf("inner_%d_lambda_%g.rds", folds[f], lambdas[l])))
    }
  }
  choice <- arrayInd(which.min(loss), dim(loss))
  selected_size <- sizes[choice[2L]]
  selected_lambda <- lambdas[choice[1L]]
  selected_oof <- oof[[choice[1L]]][, choice[2L]]
  stopifnot(all(is.finite(selected_oof)))
  recipe <- yeast_recipe(rows, panel)
  calibration <- yeast_calibrate(rows, selected_oof, tau, recipe$replicate_levels)
  training <- yeast_features(rows, panel, recipe)
  progress("profile_training", list(Lambda = selected_lambda, additions = selected_size))
  prepared <- pc_prepare_initial(training$y, training$X, training$Z, training$cluster_id,
    tau, matrix(selected_lambda, 1L, 1L), training$penalty_factor, fit_always_included = FALSE)
  fitted <- pc_orthogonal_path(prepared, tau, 0:selected_size, solver = jp_joint_qp,
                                gap_tolerance = config$gap_tolerance)
  beta <- fitted$coefficients[, selected_size + 1L]
  names(beta) <- colnames(training$X)
  test <- yeast_features(evaluation, panel, recipe)
  prediction <- drop(yeast_raw_prediction(test, beta)) + calibration[as.character(evaluation$replicate)]
  list(beta = beta, recipe = recipe, prediction = prediction, calibration = calibration,
       oof = selected_oof, losses = loss, selected_lambda = selected_lambda, selected_size = selected_size,
       certificates = certificates, final_gap = fitted$fitted$duality_gap)
}

yeast_enet_path <- function(features, alpha, config) {
  fit <- glmnet::glmnet(features$X, features$y, weights = features$weight, family = "gaussian",
    alpha = alpha, nlambda = as.integer(config$elastic_net_path_length),
    lambda.min.ratio = config$elastic_net_lambda_min_ratio,
    penalty.factor = features$penalty_factor, intercept = FALSE, standardize = FALSE,
    thresh = 1e-10, maxit = 1000000L)
  stopifnot(fit$jerr == 0L, ncol(fit$beta) == config$elastic_net_path_length)
  fit
}

yeast_enet_prediction <- function(rows, evaluation, panel, config, progress, directory) {
  glmnet::glmnet.control(fdev = 0, devmax = 1)
  alphas <- unlist(config$elastic_net_alpha)
  count <- config$elastic_net_path_length
  loss <- matrix(0, length(alphas), count)
  oof <- lapply(alphas, function(alpha) matrix(NA_real_, nrow(rows), count))
  for (fold in sort(unique(rows$inner_fold))) {
    recipe <- yeast_recipe(rows[inner_fold != fold], panel)
    training <- yeast_features(rows[inner_fold != fold], panel, recipe)
    validation <- yeast_features(rows[inner_fold == fold], panel, recipe)
    for (a in seq_along(alphas)) {
      progress("enet_inner", list(fold = fold, alpha = alphas[a], variables = ncol(training$X)))
      fit <- yeast_enet_path(training, alphas[a], config)
      predicted <- as.matrix(yeast_raw_prediction(validation, fit$beta))
      loss[a, ] <- loss[a, ] + colSums((validation$raw_y - predicted)^2 * validation$weight)
      oof[[a]][which(rows$inner_fold == fold), ] <- predicted
      saveRDS(list(fit = fit, recipe = recipe, prediction = predicted),
              file.path(directory, sprintf("inner_%d_alpha_%g.rds", fold, alphas[a])))
    }
  }
  choice <- arrayInd(which.min(loss), dim(loss))
  selected_oof <- oof[[choice[1L]]][, choice[2L]]
  stopifnot(all(is.finite(selected_oof)))
  recipe <- yeast_recipe(rows, panel)
  training <- yeast_features(rows, panel, recipe)
  progress("enet_training", list(alpha = alphas[choice[1L]], penalty_index = choice[2L]))
  fit <- yeast_enet_path(training, alphas[choice[1L]], config)
  beta <- as.numeric(fit$beta[, choice[2L]])
  names(beta) <- colnames(training$X)
  test <- yeast_features(evaluation, panel, recipe)
  mean_prediction <- drop(yeast_raw_prediction(test, beta))
  calibrations <- lapply(config$quantiles, function(tau) yeast_calibrate(rows, selected_oof, tau, recipe$replicate_levels))
  prediction <- vapply(calibrations, function(calibration)
    mean_prediction + calibration[as.character(evaluation$replicate)], numeric(nrow(evaluation)))
  list(beta = beta, recipe = recipe, prediction = prediction, calibration = calibrations,
       mean_oof = selected_oof, loss = loss, selected_alpha = alphas[choice[1L]],
       selected_penalty_index = choice[2L], lambda = fit$lambda[choice[2L]])
}

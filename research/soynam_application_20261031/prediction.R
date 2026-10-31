soy_joint_prediction <- function(outer, inner, tau, lambda = 1) {
  maximum <- min(vapply(inner, function(fold) min(sum(fold$train$penalty_factor > 0),
    floor(fold$train$n / log(fold$train$n))), numeric(1)))
  sizes <- 0:as.integer(maximum)
  losses <- numeric(length(sizes))
  predictions <- vector("list", length(inner))
  diagnostics <- list()
  for (f in seq_along(inner)) {
    train <- inner[[f]]$train
    validation <- inner[[f]]$validation
    soy_progress("prediction_search", list(fold = f, clusters = train$n, largest_additions = maximum))
    prepared <- pc_prepare_initial(train$y, train$X, train$Z, train$cluster_id,
      tau, matrix(lambda, 1L, 1L), train$penalty_factor, fit_always_included = FALSE)
    path <- pc_orthogonal_path(prepared, tau, sizes, gap_tolerance = 1e-8,
      max_iterations = 100000L, solver = soy_pair_qp_cpp)
    predictions[[f]] <- validation$X %*% path$coefficients
    residual <- validation$y - predictions[[f]]
    losses <- losses + colSums(residual * (tau - (residual < 0)) * validation$weight)
    diagnostics[[f]] <- data.table(fold = f,
      maximum_gap = max(vapply(path$diagnostics, `[[`, numeric(1), "gap")),
      maximum_iterations = max(vapply(path$diagnostics, `[[`, numeric(1), "iterations")))
  }
  selected <- order(losses, sizes)[1L]
  oof <- rep(NA_real_, length(outer$y))
  for (f in seq_along(inner)) oof[inner[[f]]$validation$original_rows] <- predictions[[f]][, selected]
  if (any(!is.finite(oof))) stop("The training validation predictions are incomplete.")
  soy_progress("prediction_final_fit", list(additions = sizes[selected]))
  prepared <- pc_prepare_initial(outer$y, outer$X, outer$Z, outer$cluster_id,
    tau, matrix(lambda, 1L, 1L), outer$penalty_factor, fit_always_included = FALSE)
  fitted <- pc_orthogonal_path(prepared, tau, 0:sizes[selected], gap_tolerance = 1e-8,
    max_iterations = 100000L, solver = soy_pair_qp_cpp)
  beta <- fitted$coefficients[, sizes[selected] + 1L]
  names(beta) <- colnames(outer$X)
  list(beta = beta, oof = oof, tuning = data.table(additions = sizes, loss = losses / outer$n),
    selected = sizes[selected], diagnostics = rbindlist(diagnostics),
    final_gap = fitted$fitted$duality_gap, method_version = "joint_profile_prediction_paired_qp")
}

soy_pooled_path <- function(data, tau, options, lambda = NULL) {
  omitted <- match("Block1", colnames(data$X))
  if (is.na(omitted)) stop("The pooled intercept mapping lacks its reference block.")
  keep <- setdiff(seq_len(ncol(data$X)), omitted)
  args <- list(X = data$X[, keep, drop = FALSE], y = data$y, method = "quantile", tau = tau,
    alpha = 1, nlambda = options$hqreg_nlambda, lambda.min = options$hqreg_lambda_min,
    eps = options$hqreg_eps, max.iter = options$hqreg_max_iterations,
    preprocess = "standardize", penalty.factor = data$penalty_factor[keep], screen = "ASR",
    message = isTRUE(options$hqreg_message))
  if (!is.null(lambda)) args$lambda <- lambda
  fitted <- do.call(hqreg::hqreg, args)
  expected <- if (is.null(lambda)) options$hqreg_nlambda else length(lambda)
  if (!is.matrix(fitted$beta) || ncol(fitted$beta) != expected || fitted$saturated ||
      any(fitted$iter >= options$hqreg_max_iterations) || any(!is.finite(fitted$beta)))
    stop("The native pooled quantile path did not complete.")
  beta <- matrix(0, ncol(data$X), expected, dimnames = list(colnames(data$X), NULL))
  beta[keep, ] <- fitted$beta[-1L, , drop = FALSE]
  blocks <- grep("^Block[0-9]+$", colnames(data$X))
  beta[blocks, ] <- sweep(beta[blocks, , drop = FALSE], 2L, fitted$beta[1L, ], "+")
  mapped <- data$X %*% beta
  native <- sweep(args$X %*% fitted$beta[-1L, , drop = FALSE], 2L, fitted$beta[1L, ], "+")
  if (max(abs(mapped - native)) > 1e-9) stop("The pooled coefficient mapping changed its predictions.")
  list(beta = beta, lambda = fitted$lambda, iterations = fitted$iter,
    native = fitted, mapping_error = max(abs(mapped - native)))
}

soy_pooled_prediction <- function(outer, inner, tau, options) {
  losses <- numeric(options$hqreg_nlambda)
  predictions <- vector("list", length(inner))
  diagnostics <- list()
  ratios <- exp(seq(0, log(options$hqreg_lambda_min), length.out = options$hqreg_nlambda))
  for (f in seq_along(inner)) {
    soy_progress("pooled_prediction_search", list(fold = f))
    path <- soy_pooled_path(inner[[f]]$train, tau, options)
    if (max(abs(path$lambda / path$lambda[1L] - ratios)) > 1e-8)
      stop("The native penalty ratios differ from the registered grid.")
    validation <- inner[[f]]$validation
    predictions[[f]] <- validation$X %*% path$beta
    residual <- validation$y - predictions[[f]]
    losses <- losses + colSums(residual * (tau - (residual < 0)) * validation$weight)
    diagnostics[[f]] <- data.table(fold = f, maximum_iterations = max(path$iterations), mapping_error = path$mapping_error)
  }
  selected <- which.min(losses)
  oof <- rep(NA_real_, length(outer$y))
  for (f in seq_along(inner)) oof[inner[[f]]$validation$original_rows] <- predictions[[f]][, selected]
  if (any(!is.finite(oof))) stop("The pooled validation predictions are incomplete.")
  soy_progress("pooled_prediction_final_fit", list(selected_ratio = ratios[selected]))
  fitted <- soy_pooled_path(outer, tau, options)
  list(beta = fitted$beta[, selected], oof = oof,
    tuning = data.table(ratio = ratios, loss = losses / outer$n), selected = ratios[selected],
    diagnostics = rbindlist(diagnostics), final_iterations = fitted$iterations[selected],
    method_version = "pooled_quantile_lasso_hqreg_1.4-1")
}

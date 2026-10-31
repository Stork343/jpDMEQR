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

soy_pooled_prediction <- function(outer, inner, tau, options) {
  losses <- numeric(options$pooled_nlambda)
  predictions <- vector("list", length(inner))
  diagnostics <- list()
  ratios <- exp(seq(0, log(options$pooled_lambda_min), length.out = options$pooled_nlambda))
  for (f in seq_along(inner)) {
    soy_progress("pooled_prediction_search", list(fold = f))
    path <- soy_pooled_compiled_path(inner[[f]]$train, tau, options$pooled_nlambda,
      options$pooled_lambda_min, progress = TRUE)
    saveRDS(path, file.path(soy_context$directory, paste0("training_path_", f, ".rds")), compress = FALSE)
    if (max(abs(path$lambda / path$lambda[1L] - ratios)) > 1e-8)
      stop("The pooled penalty ratios differ from the registered grid.")
    validation <- inner[[f]]$validation
    predictions[[f]] <- validation$X %*% path$beta
    residual <- validation$y - predictions[[f]]
    losses <- losses + colSums(residual * (tau - (residual < 0)) * validation$weight)
    diagnostics[[f]] <- copy(path$diagnostics)[, fold := f]
  }
  selected <- which.min(losses)
  oof <- rep(NA_real_, length(outer$y))
  for (f in seq_along(inner)) oof[inner[[f]]$validation$original_rows] <- predictions[[f]][, selected]
  if (any(!is.finite(oof))) stop("The pooled validation predictions are incomplete.")
  soy_progress("pooled_prediction_final_fit", list(selected_ratio = ratios[selected]))
  fitted <- soy_pooled_compiled_path(outer, tau, options$pooled_nlambda, options$pooled_lambda_min, progress = TRUE)
  saveRDS(fitted, file.path(soy_context$directory, "training_path_final.rds"), compress = FALSE)
  list(beta = fitted$beta[, selected], oof = oof,
    tuning = data.table(ratio = ratios, loss = losses / outer$n), selected = ratios[selected],
    diagnostics = rbindlist(diagnostics), final_diagnostics = fitted$diagnostics[selected],
    method_version = "pooled_quantile_lasso_highs_1.15.1")
}

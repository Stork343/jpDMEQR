jp_prepare_prediction_search <- function(data, tau, lambda) {
  scale <- sqrt(colSums(data$Z^2 * data$weight) / data$n)
  if (any(!is.finite(scale)) || any(scale <= 0)) stop("The nuisance design requires positive scales.")
  Z <- sweep(data$Z, 2L, scale, "/")
  Lambda <- diag(lambda, ncol(Z)) / outer(scale, scale)
  prepared <- pc_prepare_initial(data$y, data$X, Z, data$cluster_id, tau,
                                  Lambda, data$penalty_factor, max_iterations = 500000L,
                                  fit_always_included = ncol(Z) == 1L)
  prepared$nuisance_scale <- scale
  prepared
}

jp_prediction_fit <- function(inner, outer, tau, lambda_grid = c(.5, 1, 2), checkpoint = NULL) {
  maximum <- min(vapply(inner, function(fold) {
    min(sum(fold$train$penalty_factor > 0), floor(fold$train$n / log(fold$train$n)))
  }, numeric(1)))
  sizes <- 0:as.integer(maximum)
  loss <- matrix(0, length(lambda_grid), length(sizes))
  diagnostics <- list()
  for (f in seq_along(inner)) {
    train <- inner[[f]]$train
    validation <- inner[[f]]$validation
    for (l in seq_along(lambda_grid)) {
      if (!is.null(checkpoint)) checkpoint(list(stage = "preparation", fold = f, lambda = lambda_grid[l]))
      prepared <- jp_prepare_prediction_search(train, tau, lambda_grid[l])
      if (!is.null(checkpoint)) checkpoint(list(stage = "search", fold = f, lambda = lambda_grid[l]))
      engine <- if (ncol(train$Z) == 1L) joint_profile_pdal_cpp else jp_joint_qp
      path <- pc_orthogonal_path(prepared, tau, sizes, gap_tolerance = 1e-8, max_iterations = 500000L,
                                  solver = engine)
      residual <- validation$y - validation$X %*% path$coefficients
      check <- residual * (tau - (residual < 0))
      loss[l, ] <- loss[l, ] + colSums(check * validation$weight)
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(
        fold = f, lambda = lambda_grid[l],
        max_gap = if (length(path$diagnostics)) max(vapply(path$diagnostics, `[[`, numeric(1), "gap")) else prepared$null_gap,
        max_iterations = if (length(path$diagnostics)) max(vapply(path$diagnostics, `[[`, numeric(1), "iterations")) else prepared$null_iterations)
    }
  }
  candidates <- expand.grid(lambda = lambda_grid, size = sizes)
  candidates$loss <- as.vector(loss) / outer$n
  selected <- order(candidates$loss, candidates$size, candidates$lambda)[1L]
  setting <- candidates[selected, ]
  prepared <- jp_prepare_prediction_search(outer, tau, setting$lambda)
  fitted <- pc_orthogonal_path(prepared, tau, 0:setting$size,
                                gap_tolerance = 1e-8, max_iterations = 500000L,
                                solver = if (ncol(outer$Z) == 1L) joint_profile_pdal_cpp else jp_joint_qp)
  beta <- fitted$coefficients[, setting$size + 1L]
  names(beta) <- colnames(outer$X)
  list(beta = beta, selected = setting, candidates = candidates,
       diagnostics = do.call(rbind, diagnostics), selected_columns = fitted$selected_path[[setting$size + 1L]],
       final_gap = if (is.null(fitted$fitted)) prepared$null_gap else fitted$fitted$duality_gap,
       nuisance_scale = prepared$nuisance_scale,
       method_version = "joint_profile_prediction_20261026")
}

jp_prediction_data <- function(metadata, base, rows, model, filter_size, random_slope = FALSE) {
  clinical <- base$clinical[[model]]$X
  modules <- base$module_scores
  colnames(modules) <- paste0("module__", sub("^module__", "", colnames(modules)))
  genes <- base$genes[, seq_len(filter_size), drop = FALSE]
  colnames(genes) <- paste0("gene__", sub("^gene__", "", colnames(genes)))
  X <- cbind(clinical, modules, genes)
  Z <- if (random_slope) cbind(intercept = 1, centered_time = metadata$cumulative_time - base$clinical[[model]]$time_center) else matrix(1, nrow(metadata), 1L)
  data <- sb_data(metadata$log_c3[rows], X[rows, , drop = FALSE], Z[rows, , drop = FALSE], metadata$subject_id[rows])
  data$penalty_factor <- c(rep(0, ncol(clinical) + ncol(modules)), rep(1, ncol(genes)))
  data
}

drug_pooled_path <- function(data, tau, nlambda = 40L, lambda_min = 0.01, lambda = NULL,
                             progress = FALSE) {
  p <- ncol(data$X)
  penalised <- which(data$penalty_factor > 0)
  means <- numeric(p)
  means[penalised] <- colMeans(data$X[, penalised, drop = FALSE])
  X <- sweep(data$X, 2L, means, "-")
  scales <- sqrt(colMeans(X^2))
  stopifnot(all(is.finite(scales)), all(scales > 1e-10))
  X <- sweep(X, 2L, scales, "/")
  fitted <- soy_pooled_lp_cpp(X, data$y, data$penalty_factor, tau, as.integer(nlambda), lambda_min,
                              if (is.null(lambda)) numeric() else lambda, progress)
  beta <- fitted$beta / scales
  offset <- -drop(means %*% beta)
  intercept_basis <- grep("^\\(Intercept\\)_dose[1-5]$", colnames(data$X))
  stopifnot(length(intercept_basis) == 5L,
            max(abs(rowSums(data$X[, intercept_basis, drop = FALSE]) - 1)) < 1e-12)
  beta[intercept_basis, ] <- sweep(beta[intercept_basis, , drop = FALSE], 2L, offset, "+")
  mapping_error <- max(abs(data$X %*% beta - X %*% fitted$beta))
  stopifnot(mapping_error < 1e-9)
  rownames(beta) <- colnames(data$X)
  list(beta = beta, lambda = fitted$lambda, scales = scales, means = means,
       diagnostics = data.table(lambda = fitted$lambda, gap = fitted$primal_dual_gap,
         dual_violation = fitted$dual_violation, mapping_error = mapping_error,
         selected_variables = fitted$selected_variables, iterations = fitted$iterations),
       solver = paste0("HiGHS_", fitted$highs_version))
}

fit_application_comparator_v2 <- function(method, metadata, design, training_rows, tau,
                                           lambda_multiplier = 1, seed = 1L) {
  set.seed(seed)
  indices <- which(training_rows)
  indices <- indices[order(as.character(metadata$subject_id[indices]), metadata$cumulative_time[indices])]
  X <- design$X[indices, , drop = FALSE]
  y <- metadata$log_c3[indices]
  ids <- as.character(metadata$subject_id[indices])
  started <- proc.time()[[3]]
  if (method == "POOLED-QR-LASSO") {
    plugin <- quantreg::LassoLambdaHat(X, tau = tau)
    penalty <- lambda_multiplier * plugin
    penalty[design$penalty_factor == 0] <- 0
    fit <- quantreg::rq.fit.lasso(X, y, tau = tau, lambda = penalty)
    beta <- fit$coefficients
    record <- list(beta_prediction = beta, method_version = "quantreg_rq_fit_lasso",
                    target_type = "pooled_visit_weighted_quantile", package_version = as.character(utils::packageVersion("quantreg")),
                    lambda = penalty, lambda_plugin = plugin, lambda_multiplier = lambda_multiplier,
                    converged = TRUE, iterations = fit$it, status = "ok")
  } else if (method == "QGEE-SCAD") {
    if (!all(X[, 1] == 1) || colnames(X)[1] != "(Intercept)") {
      stop("The official qpgee application requires its declared intercept schema.")
    }
    attr(X, "assign") <- c(0L, seq_len(ncol(X) - 1L))
    nobs <- as.integer(table(factor(ids, levels = unique(ids))))
    kernel <- get(".fit_single_lambda", envir = asNamespace("geeVerse"))
    lambda_min_ratio <- if (length(nobs) > ncol(X)) .001 else .01
    lambda_path <- exp(seq(log(10), log(10 * lambda_min_ratio), length.out = 30L))
    path_fits <- lapply(lambda_path, function(lambda) {
      kernel(x = X, y = y, nobs = nobs, tau = tau, corstr = "exchangeable",
               lambda = lambda, penalty_factor = design$penalty_factor[-1],
               intercept = TRUE, betaint = NULL, control = geeVerse::qpgeeControl())
    })
    hbic <- vapply(path_fits, function(value) value$hbic, numeric(1))
    if (any(!is.finite(hbic))) stop("An official QGEE kernel returned nonfinite HBIC.")
    selected <- which.min(hbic)
    fit <- path_fits[[selected]]
    fit$best_lambda <- lambda_path[selected]
    fit$lambda_path <- data.frame(lambda = lambda_path, hbic = hbic,
                                   converged = vapply(path_fits, function(value) value$converged, logical(1)),
                                   iterations = vapply(path_fits, function(value) as.integer(value$iterations), integer(1)))
    beta <- fit$coefficients
    record <- list(beta_prediction = beta, method_version = "geeVerse_official_single_lambda_kernel_HBIC",
                    target_type = "marginal_longitudinal_quantile", package_version = as.character(utils::packageVersion("geeVerse")),
                    converged = fit$converged, iterations = fit$iterations,
                    lambda = fit$best_lambda, lambda_path = fit$lambda_path,
                    hbic = fit$hbic, working_correlation = "exchangeable",
                    penalty_factor_excluding_intercept = design$penalty_factor[-1],
                    implementation_source = "geeVerse_namespace_fit_single_lambda_with_original_30_point_HBIC_path",
                    status = if (isTRUE(fit$converged)) "ok" else "warning")
  } else if (method == "LQMM") {
    fixed_names <- sprintf("covariate_%03d", seq_len(ncol(X)))
    data <- data.frame(response = y, subject = factor(ids), X, check.names = FALSE)
    names(data)[seq_len(ncol(X)) + 2L] <- fixed_names
    formula <- stats::reformulate(fixed_names, response = "response", intercept = FALSE)
    fit <- lqmm::lqmm(fixed = formula, random = ~1, group = "subject", tau = tau,
                       data = data, covariance = "pdSymm", nK = 7L,
                       control = lqmm::lqmmControl(LP_tol_ll = 1e-6))
    beta <- setNames(as.numeric(stats::coef(fit)), colnames(X))
    record <- list(beta_prediction = beta, method_version = "lqmm_official_lowdimensional",
                    target_type = "conditional_linear_quantile_mixed_model",
                    package_version = as.character(utils::packageVersion("lqmm")),
                    convergence_information = fit$opt, fit_object = fit,
                    formula_coordinate_map = data.frame(fixed_name = fixed_names, coordinate = colnames(X)),
                    converged = NA, status = "returned_package_fit")
  } else stop("Unknown declared application comparator.")
  if (length(beta) != ncol(X) || any(!is.finite(beta))) stop("The comparator returned invalid coefficients.")
  names(record$beta_prediction) <- colnames(X)
  record$method <- method
  record$runtime_sec <- proc.time()[[3]] - started
  record
}

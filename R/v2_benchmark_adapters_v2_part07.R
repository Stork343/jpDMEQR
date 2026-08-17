# Fidelity implementation of QGEE-SCAD.
#
# Source method: Zu, T., Lian, H., Green, B., and Yu, Y. (2023). Ultra-High
# Dimensional Quantile Regression for Longitudinal Data: An Application to
# Blood Pressure Analysis. JASA 118(541):97-108,
# DOI 10.1080/01621459.2022.2128806.
#
# The method is the quantile penalised GEE with SCAD penalty, published with
# the maintained reference package `geeVerse` (CRAN 0.3.1, authors Zu/Green/
# Yu; GPL-3). The adapter calls the official implementation on the supplied
# training design and records the method's own marginal estimating-equation
# target. It does not relabel the proposed estimator.
#
# Source details (from the package and the paper's slides/supplement):
#   - induced-smoothing score Phi(sqrt(n) u / r_i) - (1-tau), r_i = ||X_i||;
#   - SCAD penalty with a = 3.7 via LLA/MM (c = 1e-6 perturbation);
#   - working correlation: independence / exchangeable / AR(1) / tri /
#     unstructured, moment-estimated from standardised quantile residuals;
#   - Newton updates with active-set shrinkage; HBIC or 5-fold subject-level
#     CV for lambda; SIS screening for ultra-high dimension.

fit_benchmark_qgee_scad_v2 <- function(train,
                                       tau,
                                       target_coords,
                                       tuning,
                                       seed,
                                       context = list(),
                                       control = list()) {
  if (!requireNamespace("geeVerse", quietly = TRUE)) {
    return(benchmark_failure_v2("QGEE-SCAD", "dependency",
                                "geeVerse package is required."))
  }
  set.seed(seed)
  start <- proc.time()[[3]]
  X <- as.matrix(train$X)
  y <- as.numeric(train$y)
  cluster_id <- as.character(train$cluster_id)
  if (length(y) != nrow(X) || length(cluster_id) != length(y)) {
    return(benchmark_failure_v2("QGEE-SCAD", "schema",
                                "Dimension mismatch."))
  }
  # Observation counts per cluster, aligned with the row order.
  cluster_factor <- factor(cluster_id, levels = unique(cluster_id))
  nobs <- as.integer(table(cluster_factor))
  if (sum(nobs) != length(y)) {
    return(benchmark_failure_v2("QGEE-SCAD", "schema",
                                "Cluster sizes do not reconcile."))
  }

  corstr <- control$corstr %||% "exchangeable"
  method <- control$method %||% "HBIC"
  lambda <- control$lambda %||% NULL
  intercept <- isTRUE(control$intercept %||% FALSE)
  nfold <- as.integer(control$nfold %||% 5L)
  maxit <- as.integer(control$maxit %||% 100L)
  epsilon <- control$epsilon %||% 1e-4

  # geeVerse derives the intercept flag from attr(x, "assign"): a leading 0
  # marks the intercept column. Our designs carry no intercept column, so
  # declare every column as a penalised term (assign >= 1) to avoid the
  # empty-logical bug when the attribute is absent.
  attr(X, "assign") <- rep(seq_len(ncol(X)), 1)

  fit <- tryCatch(
    geeVerse::qpgee(
      x = X, y = y, nobs = nobs, tau = tau,
      corstr = corstr, lambda = lambda, method = method,
      intercept = intercept,
      nfold = nfold,
      control = geeVerse::qpgeeControl(
        epsilon = epsilon, maxit = maxit,
        shrinkCutoff = control$shrinkCutoff %||% 1e-4
      )
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[[3]] - start
  if (inherits(fit, "error")) {
    return(benchmark_failure_v2("QGEE-SCAD", "penalised_fit",
                                conditionMessage(fit), elapsed))
  }

  beta <- tryCatch(as.numeric(fit$coefficients %||% fit$beta),
                   error = function(e) NULL)
  if (is.null(beta) || length(beta) != ncol(X)) {
    return(benchmark_failure_v2("QGEE-SCAD", "coefficient_extraction",
                                "Could not extract coefficients.", elapsed))
  }
  names(beta) <- colnames(X)
  selected <- which(abs(beta) > (control$selection_tol %||% 1e-8))

  requested <- if (is.character(target_coords)) {
    intersect(as.character(target_coords), names(beta))
  } else {
    idx <- as.integer(target_coords)
    idx <- idx[idx >= 1L & idx <= length(beta)]
    names(beta)[idx]
  }
  inf_tab <- data.frame(
    coordinate = requested %||% character(),
    index = match(requested %||% character(), names(beta)),
    feasible = rep(FALSE, length(requested %||% character())),
    beta_hat = as.numeric(beta[requested %||% character()]),
    beta_tilde = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    mu = NA_real_,
    dantzig_residual = NA_real_,
    omega_l1 = NA_real_,
    omega_l2 = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!nrow(inf_tab)) inf_tab <- data.frame()

  ans <- list(
    method_id = "QGEE-SCAD",
    status = if (isTRUE(fit$converged %||% TRUE)) "ok" else "warning",
    failure_stage = "",
    failure_message = "",
    beta_hat = beta,
    beta_tilde = NULL,
    se = NULL,
    selected = selected,
    fit_object = fit,
    inference_object = list(
      table = inf_tab,
      ci_level = NA_real_,
      target_scope = "marginal_longitudinal_quantile",
      working_correlation = corstr,
      tuning_method = method,
      lambda = if (!is.null(lambda)) lambda else
        if (!is.null(fit$best_lambda)) fit$best_lambda else NA_real_
    ),
    runtime_sec = elapsed,
    converged = isTRUE(fit$converged %||% TRUE),
    kkt_residual = NA_real_,
    warning_messages = character(),
    implementation_version = "geeVerse-0.3.1-official",
    target_scope = "marginal_longitudinal_quantile"
  )
  benchmark_add_metadata_v2(
    ans,
    reference_identifier = "Zu-Lian-Green-Yu-2023-10.1080/01621459.2022.2128806",
    fidelity_status = "official_reference_implementation"
  )
}

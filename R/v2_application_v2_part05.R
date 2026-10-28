validate_lowdim_profile_fit_v2 <- function(fit) {
  if (!isTRUE(fit$converged) || any(!is.finite(fit$beta))) {
    stop("LOWDIM_FIT_FAILED: the unpenalised profile fit did not converge.")
  }
  comp <- fit$components
  required <- c(comp$score, comp$max_nuisance_gradient, comp$profile_identity_error)
  if (length(comp$score) != length(fit$beta) || any(!is.finite(required)) ||
      max(abs(comp$score)) > 1e-7 || comp$max_nuisance_gradient > 1e-7 ||
      comp$profile_identity_error >= 1e-8) {
    stop("LOWDIM_FIT_FAILED: profile gradient, nuisance gradient or identity check failed.")
  }
  H <- (comp$hessian + t(comp$hessian)) / 2
  eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eigenvalues)) || min(eigenvalues) <= 1e-8 ||
      max(eigenvalues) / min(eigenvalues) >= 1e10) {
    stop("LOWDIM_FIT_FAILED: the reduced effective Hessian is not identified.")
  }
  list(hessian = H, minimum_eigenvalue = min(eigenvalues),
       condition_number = max(eigenvalues) / min(eigenvalues))
}

fit_unpenalised_profile_v2 <- function(y, X, Z, cluster_id, tau, h,
                                      lambda_gamma = 1, beta_init = NULL,
                                      control = list()) {
  defaults <- list(
    max_iter = 2000L, max_backtrack = 50L, beta_tol = 1e-9,
    kkt_normalized_tol = 1e-8,
    nuisance_control = list(reltol = 1e-11, maxit = 500L, grad_tol = 1e-8)
  )
  control <- utils::modifyList(defaults, control)
  fit <- fit_profile_lasso_v2(
    y, X, Z, cluster_id, tau, h, lambda_beta = 0,
    lambda_gamma = lambda_gamma, penalty_factor = rep(0, ncol(X)),
    beta_init = beta_init, control = control
  )
  fit$lowdim_diagnostics <- validate_lowdim_profile_fit_v2(fit)
  fit
}

full_subject_profile_jackknife_v2 <- function(y, X, Z, cluster_id, tau, h,
                                             lambda_gamma = 1, fit = NULL,
                                             control = list()) {
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  y <- as.numeric(y)
  if (nrow(X) != length(y) || nrow(Z) != length(y) ||
      length(cluster_id) != length(y) || anyNA(cluster_id) || any(!is.finite(y))) {
    stop("Invalid low-dimensional profile data.")
  }
  if (is.null(colnames(X)) || anyNA(colnames(X)) || any(!nzchar(colnames(X))) ||
      anyDuplicated(colnames(X))) {
    stop("Low-dimensional fixed-effect coordinates require unique names.")
  }
  clusters <- cluster_index_v2(cluster_id)
  n <- clusters$n
  if (n < 3L || ncol(X) > floor(n / log(max(n, 3)))) {
    stop("LOWDIM_DIMENSION_FAILED: d must be at most floor(n / log(max(n, 3))).")
  }
  if (is.null(fit)) {
    fit <- fit_unpenalised_profile_v2(y, X, Z, cluster_id, tau, h, lambda_gamma,
                                      control = control)
  }
  validate_lowdim_profile_fit_v2(fit)
  rows <- split(seq_along(y), as.character(cluster_id))
  estimates <- matrix(NA_real_, nrow = n, ncol = ncol(X),
                      dimnames = list(names(rows), colnames(X)))
  diagnostics <- vector("list", n)
  for (ii in seq_along(rows)) {
    keep <- setdiff(seq_along(y), rows[[ii]])
    leave_fit <- fit_unpenalised_profile_v2(
      y[keep], X[keep, , drop = FALSE], Z[keep, , drop = FALSE],
      cluster_id[keep], tau, h, lambda_gamma, beta_init = fit$beta,
      control = control
    )
    estimates[ii, ] <- leave_fit$beta
    diagnostics[[ii]] <- data.frame(
      omitted_subject = names(rows)[ii], n_clusters = n - 1L,
      converged = leave_fit$converged,
      kkt_residual = leave_fit$kkt_residual,
      nuisance_gradient = leave_fit$components$max_nuisance_gradient,
      profile_identity_error = leave_fit$components$profile_identity_error,
      hessian_min_eigenvalue = leave_fit$lowdim_diagnostics$minimum_eigenvalue,
      hessian_condition_number = leave_fit$lowdim_diagnostics$condition_number,
      stringsAsFactors = FALSE
    )
  }
  if (nrow(estimates) != n || any(!is.finite(estimates))) {
    stop("JACKKNIFE_INCOMPLETE: every subject must have a successful leave-one fit.")
  }
  centered <- sweep(estimates, 2, colMeans(estimates))
  list(
    fit = fit, estimates = estimates,
    se = sqrt((n - 1) / n * colSums(centered^2)),
    diagnostics = do.call(rbind, diagnostics), n_clusters = n
  )
}

lowdim_profile_inference_v2 <- function(y, X, Z, cluster_id, tau, h,
                                       score_source, lambda_gamma = 1,
                                       ci_level = 0.95, control = list()) {
  if (missing(score_source) || !identical(score_source, "smoothed")) {
    stop("SANDWICH_RULE_REQUIRED: the smooth criterion requires score_source='smoothed'.")
  }
  if (length(ci_level) != 1L || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("Invalid confidence level.")
  }
  jack <- full_subject_profile_jackknife_v2(
    y, X, Z, cluster_id, tau, h, lambda_gamma, control = control
  )
  field <- "cluster_scores"
  G <- jack$fit$components[[field]]
  if (is.null(G)) stop("SANDWICH_RULE_UNAVAILABLE: approved scores are absent from the fit.")
  G <- assert_numeric_matrix_v2(G, field)
  H <- jack$fit$lowdim_diagnostics$hessian
  if (nrow(G) != jack$n_clusters || ncol(G) != ncol(X)) {
    stop("Cluster-score dimensions differ from the low-dimensional fit.")
  }
  influence <- sweep(G, 2, colMeans(G)) %*% solve(H)
  sandwich_se <- sqrt(colMeans(influence^2) / jack$n_clusters)
  critical <- stats::qnorm(1 - (1 - ci_level) / 2)
  beta <- as.numeric(jack$fit$beta)
  table <- data.frame(
    coordinate = colnames(X), beta_hat = beta, score_source = score_source,
    score_definition = "gradient_of_complete_smoothed_ridge_profile_criterion",
    target_type = "finite_h_regularised_profile", h = h, tau = tau,
    lambda_gamma = lambda_gamma,
    sandwich_se = sandwich_se, sandwich_lower = beta - critical * sandwich_se,
    sandwich_upper = beta + critical * sandwich_se,
    jackknife_se = jack$se, jackknife_lower = beta - critical * jack$se,
    jackknife_upper = beta + critical * jack$se,
    ci_level = ci_level, n_clusters = jack$n_clusters,
    jackknife_completed = nrow(jack$estimates), stringsAsFactors = FALSE
  )
  list(table = table, fit = jack$fit, jackknife = jack)
}

prepare_gse65391_lowdim_design_v2 <- function(metadata, expression_probe, annotation,
                                              module_signatures, training_rows,
                                              module_rule, include_sledai = FALSE,
                                              module_score_members = NULL) {
  if (!all(c("source", "minimum_fraction") %in% names(module_rule)) ||
      !nzchar(module_rule$source) || length(module_rule$minimum_fraction) != 1L ||
      !is.finite(module_rule$minimum_fraction) || module_rule$minimum_fraction <= 0 ||
      module_rule$minimum_fraction > 1) {
    stop("MODULE_RULE_REQUIRED: supply the frozen member-proportion rule and its source.")
  }
  if (is.null(module_score_members)) {
    stop("MODULE_RULE_REQUIRED: supply the frozen measurement units, full members and weights.")
  }
  if (!identical(colnames(expression_probe), rownames(metadata))) {
    stop("MODULE_GATE_FAILED: expression columns must follow metadata sample identifiers.")
  }
  clinical <- prepare_clinical_design_split_v2(metadata, training_rows, include_sledai)
  modules <- score_frozen_module_representations_v2(
    expression_probe, annotation, module_signatures, module_score_members,
    training_rows, module_rule
  )
  if (!all(modules$audit$eligible) || any(!is.finite(modules$scores))) {
    stop("MODULE_GATE_FAILED: frozen modules have insufficient observed member genes.")
  }
  colnames(modules$scores) <- paste0("module__", colnames(modules$scores))
  X <- cbind(clinical$X, modules$scores)
  if (qr(X[training_rows, , drop = FALSE])$rank != ncol(X)) {
    stop("LOWDIM_DIMENSION_FAILED: the frozen clinical/module design is rank deficient.")
  }
  list(X = X, Z = matrix(1, nrow(X), 1L, dimnames = list(NULL, "intercept")),
       clinical = clinical, module_audit = modules$audit,
       module_rule = module_rule,
       gene_preprocessing = list(center = modules$center, scale = modules$scale,
                                 gene_probe_map = modules$gene_probe_map,
                                 score_members = modules$score_members,
                                 training_sample_ids = modules$training_sample_ids))
}

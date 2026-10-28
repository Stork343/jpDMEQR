fit_application_profile_initial_v2 <- function(y, X, Z, cluster_id, tau, h_est = NULL,
                                               lambda_gamma = 1, mean_score_multiplier = 1,
                                               lambda_multiplier = 1, penalty_factor = NULL,
                                               control = list()) {
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  n <- cluster_index_v2(cluster_id)$n
  p <- ncol(X)
  base_pf <- as.numeric(penalty_factor %||% rep(1, p))
  P <- which(base_pf > 0)
  if (!length(P) || length(base_pf) != p || any(!is.finite(base_pf)) || any(base_pf < 0) ||
      is.null(colnames(X)) || anyDuplicated(colnames(X))) {
    stop("The predictive initial estimator requires named penalised coordinates.")
  }
  if (is.null(h_est)) h_est <- (log(max(length(P), 2L)) / n)^(1 / 4)
  validate_tau_h_v2(tau, h_est)
  if (any(!is.finite(c(mean_score_multiplier, lambda_multiplier))) ||
      length(mean_score_multiplier) != 1L || length(lambda_multiplier) != 1L ||
      mean_score_multiplier <= 0 || lambda_multiplier <= 0) stop("Invalid prospective penalty multipliers.")
  allowed <- c("fit_control", "nuisance_control", "assembly", "clinical_acceleration", "clinical_cache")
  if (length(setdiff(names(control), allowed))) stop("Unknown initial-estimator control field.")
  Lambda <- make_lambda_matrix_v2(lambda_gamma, ncol(Z))
  nuisance_control <- modifyList(list(reltol = 1e-11, maxit = 400L, grad_tol = 1e-8),
                                 control$nuisance_control %||% list())
  fit_control <- modifyList(
    list(max_iter = 2000L, max_backtrack = 50L, beta_tol = 1e-7,
         kkt_normalized_tol = 1e-3, hessian_refresh = hessian_refresh_v2()),
    control$fit_control %||% list()
  )
  fit_control$nuisance_control <- nuisance_control
  acceleration <- control$clinical_acceleration %||% Sys.getenv(
    "JPDMEQR_INITIAL_ACCELERATION", unset = "none"
  )
  if (length(acceleration) != 1L || is.na(acceleration) ||
      !acceleration %in% c("none", "clinical_kkt")) stop("Unknown initial acceleration setting.")
  cache <- control$clinical_cache %||% new.env(parent = emptyenv())
  if (!is.environment(cache)) stop("clinical_cache must be a training-fold environment.")
  U <- which(base_pf == 0)
  pass_records <- list()
  solve_pass <- function(loadings, beta_init, pass_name) {
    pf <- base_pf * loadings
    penalty <- effective_lambda * pf
    record <- list(acceleration = acceleration, pass = pass_name,
                   clinical_coordinates = colnames(X)[U], cache_hit = FALSE,
                   global_certificate_accepted = FALSE,
                   coordinate_penalty = setNames(penalty, colnames(X)))
    if (acceleration == "clinical_kkt" && length(U)) {
      threshold <- min(fit_control$kkt_normalized_tol, 1e-7)
      clinical_control <- fit_control
      clinical_control$kkt_normalized_tol <- min(threshold, 1e-9)
      clinical_control$beta_tol <- min(fit_control$beta_tol, 1e-9)
      request <- list(y = as.numeric(y), X = X[, U, drop = FALSE], Z = Z,
                      cluster_id = as.character(cluster_id), tau = tau, h = h_est,
                      Lambda = Lambda, control = clinical_control,
                      assembly = control$assembly, solver_configuration = solver_configuration_v2())
      entries <- if (exists("entries", envir = cache, inherits = FALSE)) cache$entries else list()
      index <- which(vapply(entries, function(entry) identical(entry$request, request), logical(1)))
      if (length(index) > 1L) stop("The clinical request has multiple cache entries.")
      if (length(index)) {
        clinical <- entries[[index]]$fit
        record$cache_hit <- TRUE
      } else {
        clinical <- fit_profile_unpenalised_newton_v2(
          y, X[, U, drop = FALSE], Z, cluster_id, tau, h_est,
          lambda_beta = 0, lambda_gamma = Lambda, penalty_factor = rep(0, length(U)),
          beta_init = rep(0, length(U)), control = clinical_control
        )
        profile_pair_check_fit_v2(clinical, "initial_clinical_newton",
                                  clinical_control$kkt_normalized_tol, clinical_control$beta_tol)
        if (!is.finite(clinical$components$profile_identity_error) ||
            clinical$components$profile_identity_error > 1e-8) {
          profile_pair_stop_v2("initial_clinical_identity", "Clinical profile identity acceptance failed.", clinical)
        }
        entries[[length(entries) + 1L]] <- list(request = request, fit = clinical)
        cache$entries <- entries
      }
      beta_clinical <- setNames(numeric(p), colnames(X))
      beta_clinical[U] <- clinical$beta
      full <- profile_components_v2(
        y, X, Z, cluster_id, beta_clinical, tau, h_est, lambda_gamma = Lambda,
        gamma_init = clinical$gamma, need_hessian = TRUE, identity_check = TRUE,
        need_fold_sums = FALSE,
        nuisance_control = modifyList(nuisance_control, list(grad_tol = clinical$numerical_nuisance_grad_tol)),
        assembly = control$assembly
      )
      profile_pair_check_components_v2(full, "initial_full_kkt", clinical)
      if (!is.finite(full$profile_identity_error) || full$profile_identity_error > 1e-8) {
        profile_pair_stop_v2("initial_full_identity", "Full initial profile identity acceptance failed.", clinical)
      }
      residual <- setNames(pmax(0, abs(full$score) - penalty), colnames(X))
      residual[U] <- abs(full$score[U])
      global_kkt <- kkt_residual_normalized_v2(beta_clinical, full$score, penalty)
      record$global_residual <- residual
      record$global_score <- setNames(full$score, colnames(X))
      record$global_kkt_absolute <- max(residual)
      record$global_kkt_normalized <- global_kkt
      record$global_certificate_threshold <- threshold
      record$full_profile_identity_error <- full$profile_identity_error
      record$full_max_nuisance_gradient <- full$max_nuisance_gradient
      record$clinical_history <- clinical$history
      record$clinical_beta <- clinical$beta
      record$clinical_kkt <- clinical$kkt_residual
      record$clinical_beta_change <- clinical$first_stage_beta_change
      record$global_certificate_accepted <- global_kkt <= threshold && max(residual) <= threshold
      if (record$global_certificate_accepted) {
        fit <- clinical
        fit$beta <- beta_clinical
        fit$gamma <- full$gamma
        fit$components <- full
        fit$lambda_beta <- effective_lambda
        fit$penalty_factor <- pf
        fit$profile_objective <- full$objective
        fit$composite_objective <- full$objective + sum(penalty * abs(beta_clinical))
        fit$kkt_residual <- max(residual)
        fit$kkt_normalized <- global_kkt
        fit$kkt_normalized_tol <- threshold
        fit$first_stage_nonzero_count <- sum(abs(beta_clinical) > 1e-10)
        fit$optimizer <- "clinical_newton_full_kkt_certificate"
        fit$initial_hessian_fold_sums <- "not_requested"
        record$optimizer <- fit$optimizer
        pass_records[[pass_name]] <<- record
        return(fit)
      }
      beta_init <- beta_clinical
    }
    record$proximal_initial_beta <- setNames(as.numeric(beta_init), colnames(X))
    fit <- fit_profile_lasso_v2(
      y, X, Z, cluster_id, tau, h_est, effective_lambda, lambda_gamma = Lambda,
      penalty_factor = pf, beta_init = beta_init, control = fit_control
    )
    record$optimizer <- "proximal_gradient"
    record$proximal_history <- fit$history
    final_residual <- pmax(0, abs(fit$components$score) - penalty)
    active <- pf > 0 & abs(fit$beta) > 1e-10
    final_residual[active] <- abs(fit$components$score[active] + penalty[active] * sign(fit$beta[active]))
    final_residual[U] <- abs(fit$components$score[U])
    record$final_global_residual <- setNames(final_residual, colnames(X))
    record$final_global_score <- setNames(fit$components$score, colnames(X))
    record$final_beta <- fit$beta
    record$coordinate_penalty <- setNames(penalty, colnames(X))
    pass_records[[pass_name]] <<- record
    fit
  }
  stochastic <- lambda_0_n_v2(n, length(P))$lambda_0
  mean_component <- mean_score_multiplier * h_est
  effective_lambda <- lambda_multiplier * (stochastic + mean_component)
  if (!is.finite(effective_lambda) || effective_lambda <= 0) stop("Invalid prospective initial penalty.")
  component0 <- profile_components_v2(
    y, X, Z, cluster_id, rep(0, p), tau, h_est, lambda_gamma = Lambda,
    need_hessian = FALSE, nuisance_control = nuisance_control, assembly = control$assembly
  )
  profile_pair_check_components_v2(component0, "lambda_loading_pass0")
  loadings0 <- setNames(score_loadings_v2(component0$cluster_scores), colnames(X))
  check_loadings <- function(values) {
    positive <- values[P][is.finite(values[P]) & values[P] > 0]
    if (!length(positive) || any(!is.finite(values[P])) ||
        any(values[P] < 1e-8 * stats::median(positive))) {
      profile_pair_stop_v2("lambda_loading_degenerate", "Cluster score loadings are degenerate.")
    }
  }
  check_loadings(loadings0)
  preliminary <- solve_pass(loadings0, rep(0, p), "preliminary")
  profile_pair_check_fit_v2(preliminary, "penalised_fit", 1e-3, 1e-7)
  component1 <- profile_components_v2(
    y, X, Z, cluster_id, preliminary$beta, tau, h_est, lambda_gamma = Lambda,
    gamma_init = preliminary$gamma, need_hessian = FALSE,
    nuisance_control = nuisance_control, assembly = control$assembly
  )
  profile_pair_check_components_v2(component1, "lambda_loading_pass1", preliminary)
  loadings1 <- setNames(score_loadings_v2(component1$cluster_scores), colnames(X))
  final_loadings <- pmax(loadings0, loadings1)
  check_loadings(final_loadings)
  fit <- solve_pass(final_loadings, preliminary$beta, "final")
  profile_pair_check_fit_v2(fit, "penalised_fit", 1e-3, 1e-7)
  fit$beta_l1 <- fit$beta
  fit$h_est <- h_est
  fit$method_version <- "profile_pair_20261022_initial_l1"
  fit$score_definition <- "gradient_of_complete_smoothed_ridge_profile_criterion_at_h_est"
  fit$loadings <- list(pass0 = loadings0, pass1 = loadings1, final = final_loadings)
  fit$lambda_stochastic <- stochastic
  fit$lambda_mean <- mean_component
  fit$lambda_base <- stochastic + mean_component
  fit$lambda_multiplier <- lambda_multiplier
  fit$coordinate_penalty <- setNames(effective_lambda * base_pf * final_loadings, colnames(X))
  fit$initial_acceleration <- acceleration
  fit$initial_optimizer_history <- pass_records
  fit
}

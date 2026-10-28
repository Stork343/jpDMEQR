profile_pair_stop_v2 <- function(stage, message, fit = NULL, diagnostics = NULL) {
  condition <- structure(
    list(message = message, call = NULL, failure_stage = stage,
         method_version = "profile_pair_20261022",
         score_definition = "profile_bandwidth_pair_2_half_minus_full",
         fit = fit, diagnostics = diagnostics),
    class = c("profile_pair_error", "error", "condition")
  )
  stop(condition)
}

profile_pair_run_stage_v2 <- function(stage, fit = NULL, expr) {
  # 添加计算阶段与原始 condition，异常继续终止当前调用。
  withCallingHandlers(expr, error = function(condition) {
    if (!inherits(condition, "profile_pair_error")) {
      profile_pair_stop_v2(stage, conditionMessage(condition), fit,
                          list(cause = condition))
    }
  })
}

profile_pair_coordinate_index_v2 <- function(coordinates, column_names,
                                            allow_empty = FALSE) {
  if (is.null(coordinates) || !length(coordinates)) {
    if (isTRUE(allow_empty)) return(integer())
    stop("coordinates must contain at least one prespecified coordinate.")
  }
  if (is.character(coordinates)) {
    index <- match(coordinates, column_names)
    if (anyNA(index)) stop("Unknown coordinates: ",
                          paste(coordinates[is.na(index)], collapse = ", "))
  } else {
    if (!is.numeric(coordinates) || anyNA(coordinates) ||
        any(!is.finite(coordinates)) ||
        any(coordinates != as.integer(coordinates)) ||
        any(coordinates < 1L | coordinates > length(column_names))) {
      stop("coordinates must contain valid integer indices or column names.")
    }
    index <- as.integer(coordinates)
  }
  if (anyDuplicated(index)) stop("coordinates must be unique.")
  index
}

profile_pair_check_components_v2 <- function(components, stage, fit = NULL) {
  if (!isTRUE(components$nuisance_converged) ||
      !is.finite(components$max_nuisance_gradient) ||
      components$max_nuisance_gradient > 1e-7 ||
      any(!is.finite(components$score)) || !is.finite(components$objective)) {
    profile_pair_stop_v2(stage, "Nuisance profile acceptance failed.", fit,
                        list(max_nuisance_gradient =
                               components$max_nuisance_gradient))
  }
  invisible(TRUE)
}

profile_pair_check_fit_v2 <- function(fit, stage, kkt_tolerance, beta_tolerance,
                                      fit_context = NULL) {
  change <- fit$first_stage_beta_change
  accepted <- isTRUE(fit$converged) &&
    is.finite(fit$kkt_normalized) && fit$kkt_normalized <= kkt_tolerance &&
    length(change) == 1L && is.finite(change) &&
    change <= beta_tolerance * max(1, max(abs(fit$beta))) &&
    all(is.finite(fit$beta)) && is.finite(fit$composite_objective)
  if (!accepted) {
    profile_pair_stop_v2(stage, "Convex profile fit acceptance failed.",
                        fit_context %||% fit,
                        list(kkt_normalized = fit$kkt_normalized,
                             beta_change = change,
                             optimizer_failure_stage = fit$failure_stage,
                             attempted_fit = fit))
  }
  profile_pair_check_components_v2(fit$components, stage, fit_context %||% fit)
  invisible(TRUE)
}

fit_profile_pair_v2 <- function(y, X, Z, cluster_id, tau, h_inf,
                                h_est = NULL, lambda_gamma = 1,
                                mean_score_multiplier = 1,
                                lambda_multiplier = 1,
                                penalty_factor = NULL,
                                target_coordinates = NULL,
                                control = list()) {
  started <- proc.time()[[3]]
  validate_tau_h_v2(tau, h_inf)
  y <- as.numeric(y)
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  p <- ncol(X)
  if (!length(y) || any(!is.finite(y)) || p < 1L || ncol(Z) < 1L ||
      nrow(X) != length(y) || nrow(Z) != length(y) ||
      length(cluster_id) != length(y)) stop("Incompatible profile dimensions.")
  if (is.null(colnames(X)) || anyNA(colnames(X)) ||
      any(!nzchar(colnames(X))) || anyDuplicated(colnames(X))) {
    stop("X must have finite uniquely named columns.")
  }
  cluster <- cluster_index_v2(cluster_id)
  n <- cluster$n
  if (n < 2L) stop("At least two independent clusters are required.")
  if (is.matrix(lambda_gamma) &&
      max(abs(lambda_gamma - t(lambda_gamma))) >
      1e-12 * max(1, max(abs(lambda_gamma)))) {
    stop("lambda_gamma must be symmetric.")
  }
  Lambda <- make_lambda_matrix_v2(lambda_gamma, ncol(Z))
  base_pf <- as.numeric(penalty_factor %||% rep(1, p))
  if (length(base_pf) != p || any(!is.finite(base_pf)) || any(base_pf < 0)) {
    stop("penalty_factor must be a finite nonnegative vector of length p.")
  }
  for (value in list(mean_score_multiplier, lambda_multiplier)) {
    if (length(value) != 1L || !is.finite(value) || value <= 0) {
      stop("The prospective penalty multipliers must be positive scalars.")
    }
  }
  allowed <- c("fit_control", "refit_control", "nuisance_control", "assembly")
  if (!is.list(control) ||
      (length(control) && (is.null(names(control)) || any(!nzchar(names(control))))) ||
      anyDuplicated(names(control)) || length(setdiff(names(control), allowed))) {
    stop("Unknown paired fit control field.")
  }
  P <- which(base_pf > 0)
  U <- which(base_pf == 0)
  targets <- profile_pair_coordinate_index_v2(
    target_coordinates, colnames(X), allow_empty = TRUE
  )
  if (is.null(h_est)) h_est <- (log(max(length(P), 2L)) / n)^(1 / 4)
  validate_tau_h_v2(tau, h_est)
  nuisance_control <- modifyList(
    list(reltol = 1e-11, maxit = 400L, grad_tol = 1e-8),
    control$nuisance_control %||% list()
  )
  fit_control <- modifyList(
    list(max_iter = 2000L, max_backtrack = 50L, beta_tol = 1e-7,
         kkt_normalized_tol = 1e-3, hessian_refresh = hessian_refresh_v2()),
    control$fit_control %||% list()
  )
  fit_control$nuisance_control <- nuisance_control
  refit_control <- modifyList(
    list(max_iter = 1000L, max_backtrack = 50L, beta_tol = 1e-8,
         kkt_normalized_tol = 1e-7), control$refit_control %||% list()
  )
  refit_control$nuisance_control <- nuisance_control
  fit <- NULL
  preliminary <- NULL
  components_first_stage <- NULL
  loadings <- list(pass0 = NULL, pass1 = NULL, final = NULL)
  lambda_meta <- NULL
  stochastic <- mean_component <- 0
  beta_l1 <- NULL

  if (length(P)) {
    lambda_meta <- lambda_0_n_v2(n, length(P))
    stochastic <- lambda_meta$lambda_0
    mean_component <- mean_score_multiplier * h_est
    lambda_base <- stochastic + mean_component
    effective_lambda <- lambda_multiplier * lambda_base
    selection_audit <- list(
      lambda_stochastic = stochastic, lambda_mean = mean_component,
      lambda_base = lambda_base, lambda_multiplier = lambda_multiplier,
      h_est = h_est, h_inf = h_inf, Lambda = Lambda,
      method_version = "profile_pair_20261022"
    )
    component0 <- profile_pair_run_stage_v2(
      "lambda_loading_pass0", selection_audit,
      profile_components_v2(
        y, X, Z, cluster_id, rep(0, p), tau, h_est,
        lambda_gamma = Lambda, need_hessian = FALSE,
        nuisance_control = nuisance_control, assembly = control$assembly
      )
    )
    profile_pair_check_components_v2(component0, "lambda_loading_pass0")
    check_loadings <- function(values, pass) {
      positive <- values[P][is.finite(values[P]) & values[P] > 0]
      if (!length(positive) || any(!is.finite(values[P])) ||
          any(values[P] < 1e-8 * stats::median(positive))) {
        profile_pair_stop_v2("lambda_loading_degenerate",
                            "Cluster score loadings are degenerate.",
                            diagnostics = list(pass = pass, loadings = values))
      }
      invisible(TRUE)
    }
    loadings$pass0 <- setNames(score_loadings_v2(component0$cluster_scores),
                               colnames(X))
    check_loadings(loadings$pass0, "pass0")
    selection_audit$loadings <- loadings
    preliminary <- profile_pair_run_stage_v2(
      "penalised_fit", selection_audit,
      fit_profile_lasso_v2(
        y, X, Z, cluster_id, tau, h_est, effective_lambda,
        lambda_gamma = Lambda, penalty_factor = base_pf * loadings$pass0,
        beta_init = rep(0, p), control = fit_control
      )
    )
    profile_pair_check_fit_v2(preliminary, "penalised_fit", 1e-3, 1e-7,
                             selection_audit)
    selection_audit$beta_preliminary <- preliminary$beta
    component1 <- profile_pair_run_stage_v2(
      "lambda_loading_pass1", selection_audit,
      profile_components_v2(
        y, X, Z, cluster_id, preliminary$beta, tau, h_est,
        lambda_gamma = Lambda, gamma_init = preliminary$gamma,
        need_hessian = FALSE, nuisance_control = nuisance_control,
        assembly = control$assembly
      )
    )
    profile_pair_check_components_v2(component1, "lambda_loading_pass1", preliminary)
    loadings$pass1 <- setNames(score_loadings_v2(component1$cluster_scores),
                               colnames(X))
    loadings$final <- pmax(loadings$pass0, loadings$pass1)
    check_loadings(loadings$final, "final")
    selection_audit$loadings <- loadings
    fit <- profile_pair_run_stage_v2(
      "penalised_fit", selection_audit,
      fit_profile_lasso_v2(
        y, X, Z, cluster_id, tau, h_est, effective_lambda,
        lambda_gamma = Lambda, penalty_factor = base_pf * loadings$final,
        beta_init = preliminary$beta, control = fit_control
      )
    )
    profile_pair_check_fit_v2(fit, "penalised_fit", 1e-3, 1e-7, selection_audit)
    beta_l1 <- fit$beta
    components_first_stage <- fit$components
  } else {
    lambda_base <- effective_lambda <- 0
  }
  first_stage_sec <- proc.time()[[3]] - started
  selected <- if (length(P)) P[abs(beta_l1[P]) > 1e-10] else integer()
  refit_set <- sort(unique(c(selected, U, targets)))
  refit_gate <- floor(n / log(max(n, 3)))
  partial <- list(
    beta_l1 = beta_l1, selected_support = colnames(X)[selected],
    refit_set = refit_set, refit_set_names = colnames(X)[refit_set],
    refit_set_size = length(refit_set), refit_dimension_limit = refit_gate,
    loadings = loadings, lambda_stochastic = stochastic,
    lambda_mean = mean_component, lambda_base = lambda_base,
    h_est = h_est, h_inf = h_inf, Lambda = Lambda,
    score_definition = "profile_bandwidth_pair_2_half_minus_full",
    method_version = "profile_pair_20261022"
  )
  if (!length(refit_set) || length(refit_set) > refit_gate) {
    profile_pair_stop_v2("post_refit_dimension",
                        "The refit set must be nonempty and meet n/log(n).", partial)
  }
  refit_started <- proc.time()[[3]]
  refit <- profile_pair_run_stage_v2(
    "post_refit", partial,
    fit_profile_lasso_v2(
      y, X[, refit_set, drop = FALSE], Z, cluster_id, tau, h_inf,
      lambda_beta = 0, lambda_gamma = Lambda,
      penalty_factor = rep(0, length(refit_set)),
      beta_init = if (length(P)) beta_l1[refit_set] else rep(0, length(refit_set)),
      control = refit_control
    )
  )
  profile_pair_check_fit_v2(refit, "post_refit", 1e-7, 1e-8, partial)
  reduced_hessian <- refit$components$hessian
  eigenvalues <- profile_pair_run_stage_v2(
    "post_refit", partial,
    eigen((reduced_hessian + t(reduced_hessian)) / 2,
          symmetric = TRUE, only.values = TRUE)$values
  )
  condition_number <- max(eigenvalues) / min(eigenvalues)
  if (min(eigenvalues) <= 1e-8 || !is.finite(condition_number) ||
      condition_number >= 1e10) {
    profile_pair_stop_v2("post_refit", "Reduced profile curvature failed.", partial,
                        list(eigenvalues = eigenvalues,
                             condition_number = condition_number, refit = refit))
  }
  beta_refit <- setNames(numeric(p), colnames(X))
  beta_refit[refit_set] <- refit$beta
  refit_sec <- proc.time()[[3]] - refit_started
  paired_started <- proc.time()[[3]]
  partial$beta_refit <- beta_refit
  paired <- profile_pair_run_stage_v2(
    "profile_geometry", partial,
    profile_bandwidth_pair_v2(
      y, X, Z, cluster_id, beta_refit, tau, h_inf, lambda_gamma = Lambda,
      nuisance_control = nuisance_control, identity_check = TRUE,
      assembly = control$assembly, need_fold_sums = TRUE
    )
  )
  profile_pair_check_components_v2(paired, "profile_geometry", partial)
  paired$Lambda <- Lambda
  paired$beta <- beta_refit
  paired$tau <- tau
  paired$total_rows <- length(y)
  paired$p <- p
  paired$q <- ncol(Z)
  paired$identity_checked <- TRUE
  paired$nuisance_iterations <- paired$components_h$nuisance_iterations +
    paired$components_half_h$nuisance_iterations
  paired_reprofile_sec <- proc.time()[[3]] - paired_started
  if (is.null(fit)) fit <- refit
  fit$selection_profile_objective <- if (length(P)) fit$profile_objective else NA_real_
  fit$selection_composite_objective <- if (length(P)) fit$composite_objective else NA_real_
  fit$beta_l1 <- beta_l1
  fit$beta_refit <- beta_refit
  fit$beta <- beta_refit
  fit$gamma <- paired$components_h$gamma
  fit$gamma_half_h <- paired$components_half_h$gamma
  fit$gamma_definition <- "constituent_h_inf"
  fit$components_first_stage <- components_first_stage
  fit$components <- paired
  fit$profile_objective <- paired$objective
  fit$composite_objective <- paired$objective
  fit$objective_definition <- "paired_profile_at_beta_refit"
  fit$first_stage_calibrated <- length(P) > 0L
  fit$lambda_rule <- "cluster-score-loading-with-profile-mean-20261022"
  fit$lambda_stochastic <- stochastic
  fit$lambda_mean <- mean_component
  fit$lambda_base <- lambda_base
  fit$lambda_multiplier <- lambda_multiplier
  fit$mean_score_multiplier <- mean_score_multiplier
  fit$lambda_beta <- effective_lambda
  fit$lambda_gamma <- Lambda
  fit$Lambda <- Lambda
  fit$penalty_factor <- if (length(P)) base_pf * loadings$final else base_pf
  fit$loadings <- loadings
  fit$lambda_components <- list(
    stochastic = stochastic, mean = mean_component, base = lambda_base,
    multiplier = lambda_multiplier, mean_score_multiplier = mean_score_multiplier,
    coordinate = setNames(effective_lambda * fit$penalty_factor, colnames(X)),
    mean_score_constant_status = "prospective_model_class_condition"
  )
  fit$lambda_alpha <- if (length(P)) lambda_meta$alpha else NA_real_
  fit$lambda_normal_quantile <- if (length(P)) lambda_meta$normal_quantile else NA_real_
  fit$lambda_safety_constant <- if (length(P)) lambda_meta$safety_constant else NA_real_
  if (length(P)) {
    for (pass in c("pass0", "pass1")) {
      values <- loadings[[pass]][P]
      fit[[paste0("lambda_loading_", pass, "_min")]] <- min(values)
      fit[[paste0("lambda_loading_", pass, "_median")]] <- stats::median(values)
      fit[[paste0("lambda_loading_", pass, "_max")]] <- max(values)
    }
    coordinate <- fit$lambda_components$coordinate[P]
    fit$lambda_coordinate_min <- min(coordinate)
    fit$lambda_coordinate_median <- stats::median(coordinate)
    fit$lambda_coordinate_max <- max(coordinate)
    fit$zero_profile_score_max <- max(abs(component0$score))
    fit$zero_kkt_ratio_max <- max(abs(component0$score[P]) /
                                   (effective_lambda * base_pf[P] * loadings$pass0[P]))
    fit$preliminary_kkt_normalized <- preliminary$kkt_normalized
  }
  fit$final_kkt_normalized <- if (length(P)) fit$kkt_normalized else NA_real_
  fit$final_kkt_absolute <- if (length(P)) fit$kkt_residual else NA_real_
  fit$first_stage_iterations <- if (length(P)) fit$iterations else 0L
  fit$selected_support <- colnames(X)[selected]
  fit$selected_support_indices <- selected
  fit$selected_support_size <- length(selected)
  fit$refit_set <- refit_set
  fit$refit_set_names <- colnames(X)[refit_set]
  fit$refit_set_size <- length(refit_set)
  fit$refit_contains_targets <- all(targets %in% refit_set)
  fit$target_coordinates <- colnames(X)[targets]
  fit$post_refit_status <- "ok"
  fit$post_refit_iterations <- refit$iterations
  fit$post_refit_gradient_max <- max(abs(refit$components$score))
  fit$post_refit_nuisance_gradient_max <- refit$components$max_nuisance_gradient
  fit$post_refit_beta_change <- refit$first_stage_beta_change
  fit$post_refit_hessian_min_eigenvalue <- min(eigenvalues)
  fit$post_refit_hessian_condition_number <- condition_number
  fit$first_stage_sec <- first_stage_sec
  fit$refit_sec <- refit_sec
  fit$paired_reprofile_sec <- paired_reprofile_sec
  fit$nuisance_iterations <- paired$nuisance_iterations +
    refit$components$nuisance_iterations +
    if (length(P)) components_first_stage$nuisance_iterations else 0L
  fit$h_est <- h_est
  fit$h_inf <- h_inf
  fit$tau <- tau
  fit$score_definition <- paired$score_definition
  fit$method_version <- "profile_pair_20261022"
  fit$converged <- TRUE
  fit$failure_stage <- ""
  fit$diagnostics <- list(
    selection = if (length(P)) list(
      preliminary_kkt = preliminary$kkt_normalized,
      final_kkt = fit$final_kkt_normalized,
      selected_support = fit$selected_support
    ) else list(status = "unpenalised_design"),
    refit = list(gradient = fit$post_refit_gradient_max,
                 dimension = length(refit_set), dimension_limit = refit_gate,
                 min_eigenvalue = min(eigenvalues), condition_number = condition_number),
    paired = list(max_nuisance_gradient = paired$max_nuisance_gradient,
                  identity_error = paired$profile_identity_error,
                  hessian_min_eigenvalue = paired$hessian_min_eigenvalue)
  )
  fit$call <- match.call()
  class(fit) <- c("profile_pair_fit_v2", class(fit))
  fit
}

profile_pair_mu_domain_v2 <- function(mu_grid, fit = NULL) {
  if (!is.numeric(mu_grid) || !length(mu_grid) || any(!is.finite(mu_grid)) ||
      any(mu_grid <= 0)) stop("mu_grid must contain positive finite values.")
  requested <- as.numeric(mu_grid)
  unique_grid <- sort(unique(requested))
  admissible <- unique_grid[unique_grid < 1]
  excluded <- unique_grid[unique_grid >= 1]
  excluded <- data.frame(
    mu = excluded, status = rep("excluded_zero_direction", length(excluded)),
    stringsAsFactors = FALSE
  )
  domain <- list(mu_grid_requested = requested,
                 mu_grid_admissible = admissible,
                 excluded_zero_direction = excluded)
  if (!length(admissible)) {
    profile_pair_stop_v2("precision_solver",
                        "No Dantzig tolerance belongs to the parameter domain (0,1).",
                        fit, domain)
  }
  domain
}

profile_pair_solver_options_v2 <- function(solver_opts = list()) {
  backend <- profile_pair_precision_backend_v2(solver_opts)
  if (identical(backend, "dual_warm")) {
    settings <- modifyList(list(
      threads = 1L, solver = "simplex", presolve = "off", log_to_console = FALSE,
      primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8
    ), solver_opts$highs_options %||% list())
    return(list(backend = backend, control = list(static_regularization_constant = NA_real_),
                 highs_control = do.call(highs::highs_control, settings)))
  }
  modifyList(list(backend = backend, control = list(verbose = FALSE,
                                 static_regularization_constant = 1e-6)), solver_opts)
}

profile_pair_solve_row_v2 <- function(H, coordinate, mu, solver_opts = list(),
                                      fit = NULL) {
  if (length(mu) != 1L || !is.finite(mu) || mu <= 0 || mu >= 1) {
    profile_pair_stop_v2("precision_solver",
                        "A coordinate-inference tolerance must belong to (0,1).",
                        fit, list(mu = mu))
  }
  if (identical(profile_pair_precision_backend_v2(solver_opts), "dual_warm")) {
    return(profile_pair_solve_grid_v2(H, coordinate, mu, solver_opts, fit)[[1L]])
  }
  solver_opts <- profile_pair_solver_options_v2(solver_opts)
  solution <- profile_pair_run_stage_v2(
    "precision_solver", fit,
    solve_dantzig_row_native_v2(H, coordinate, mu, list(control = solver_opts$control))
  )
  solution$solver_options <- solver_opts
  statuses <- gsub("_", "", tolower(solution$attempts$status))
  valid <- c("solved", "almostsolved", "primalinfeasible", "almostprimalinfeasible")
  if (any(!statuses %in% valid)) {
    profile_pair_stop_v2("precision_solver", "The Dantzig solver failed.", fit,
                        diagnostics = solution)
  }
  if (isTRUE(solution$feasible) &&
      (!gsub("_", "", tolower(solution$solver_status)) %in%
       c("solved", "almostsolved") || solution$residual > mu + 1e-6)) {
    profile_pair_stop_v2("precision_solver", "Dantzig acceptance failed.", fit,
                        diagnostics = solution)
  }
  solution
}

profile_pair_mu_cv_v2 <- function(components, coordinate, mu_grid,
                                  solver_opts = list(), fit = NULL) {
  domain <- profile_pair_mu_domain_v2(mu_grid, fit)
  mu_grid <- domain$mu_grid_admissible
  if (!is.null(fit)) fit$inference_mu_domain <- domain
  H <- components$hessian
  p <- nrow(H)
  counts <- components$hessian_fold_counts
  sums <- components$hessian_fold_sums
  n <- components$n_clusters
  if (length(sums) != 4L || length(counts) != 4L || sum(counts) != n) {
    stop("Four deterministic paired Hessian fold buckets are required.")
  }
  n_folds <- if (n < 200L) 2L else 4L
  allocation <- rep(seq_len(n_folds), length.out = 4L)
  fold_n <- vapply(seq_len(n_folds), function(f) sum(counts[allocation == f]), numeric(1))
  if (any(fold_n <= 0 | fold_n >= n)) stop("Every CV fold must be nonempty.")
  fold_H <- lapply(seq_len(n_folds), function(f) {
    total <- Reduce(`+`, sums[allocation == f])
    (total + t(total)) / (2 * fold_n[f])
  })
  e <- numeric(p)
  e[coordinate] <- 1
  cells <- vector("list", length(mu_grid) * n_folds)
  candidates <- vector("list", length(mu_grid))
  completed_table <- function(entries) {
    completed <- entries[!vapply(entries, is.null, logical(1))]
    if (length(completed)) do.call(rbind, completed) else data.frame()
  }
  dual_grids <- if (identical(profile_pair_precision_backend_v2(solver_opts), "dual_warm")) {
    lapply(seq_len(n_folds), function(f) {
      training <- (n * H - fold_n[f] * fold_H[[f]]) / (n - fold_n[f])
      profile_pair_solve_grid_v2(training, coordinate, mu_grid, solver_opts, fit,
                                 context = list(fold = f, n_training = n - fold_n[f]))
    })
  } else NULL
  for (j in seq_along(mu_grid)) {
    defects <- quads <- rep(NA_real_, n_folds)
    for (f in seq_len(n_folds)) {
      training <- (n * H - fold_n[f] * fold_H[[f]]) / (n - fold_n[f])
      solution <- withCallingHandlers(
        if (is.null(dual_grids)) {
          profile_pair_solve_row_v2(training, coordinate, mu_grid[j], solver_opts, fit)
        } else dual_grids[[f]][[j]],
        error = function(condition) {
          context <- c(domain, list(
            coordinate = coordinate,
            coordinate_name = colnames(components$cluster_scores)[coordinate],
            mu = mu_grid[j], fold = f,
            n_training = n - fold_n[f], n_validation = fold_n[f],
            candidates = completed_table(candidates), folds = completed_table(cells),
            candidate_in_progress = list(mu = mu_grid[j], defects = defects,
                                         quads = quads),
            solver_diagnostics = condition$diagnostics,
            solver_options = profile_pair_solver_options_v2(solver_opts),
            cause = condition
          ))
          profile_pair_stop_v2("precision_solver", conditionMessage(condition), fit,
                              diagnostics = context)
        }
      )
      if (isTRUE(solution$feasible)) {
        defects[f] <- max(abs(as.numeric(fold_H[[f]] %*% solution$omega) - e))
        quads[f] <- 0.5 * drop(crossprod(solution$omega,
                                        fold_H[[f]] %*% solution$omega)) -
          solution$omega[coordinate]
      }
      cells[[(j - 1L) * n_folds + f]] <- data.frame(
        mu = mu_grid[j], fold = f, n_training = n - fold_n[f],
        n_validation = fold_n[f], feasible = solution$feasible,
        solver = solution$solver, solver_status = solution$solver_status,
        residual = solution$residual, defect = defects[f], quad = quads[f],
        static_regularization_constant =
          solution$solver_options$control$static_regularization_constant,
        stringsAsFactors = FALSE
      )
    }
    used <- sum(is.finite(defects))
    eligible <- used >= ceiling(n_folds / 2)
    candidates[[j]] <- data.frame(
      mu = mu_grid[j], eligible = eligible,
      defect_mean = if (eligible) mean(defects[is.finite(defects)]) else NA_real_,
      defect_se = if (used > 1L) stats::sd(defects[is.finite(defects)]) / sqrt(used) else NA_real_,
      quad_mean = if (eligible) mean(quads[is.finite(quads)]) else NA_real_,
      folds_used = used, infeasible_folds = n_folds - used,
      skipped_folds = 0L, stringsAsFactors = FALSE
    )
  }
  candidates <- do.call(rbind, candidates)
  cells <- do.call(rbind, cells)
  eligible <- which(candidates$eligible)
  if (!length(eligible)) {
    profile_pair_stop_v2("precision_solver", "No Dantzig CV candidate is eligible.", fit,
                        diagnostics = c(domain, list(candidates = candidates, folds = cells)))
  }
  minimum <- min(candidates$defect_mean[eligible])
  tied <- eligible[abs(candidates$defect_mean[eligible] - minimum) <=
                     1e-12 * max(1, abs(minimum))]
  chosen <- tied[which.min(candidates$mu[tied])]
  list(mu = candidates$mu[chosen], selected = TRUE,
       min_defect = minimum, defect_mean = candidates$defect_mean[chosen],
       quad_mean = candidates$quad_mean[chosen], folds_used = candidates$folds_used[chosen],
       candidates = candidates, folds = cells, n_folds = n_folds,
       mu_grid_requested = domain$mu_grid_requested,
       mu_grid_admissible = domain$mu_grid_admissible,
       excluded_zero_direction = domain$excluded_zero_direction,
       solver_options = profile_pair_solver_options_v2(solver_opts))
}

infer_profile_pair_v2 <- function(fit, coordinates, mu_grid, ci_level = 0.95,
                                  solver_preference = c("CLARABEL"),
                                  solver_opts = list(), df_correction = FALSE) {
  if (!inherits(fit, "profile_pair_fit_v2") || !isTRUE(fit$converged) ||
      !identical(fit$method_version, "profile_pair_20261022") ||
      !identical(fit$score_definition, "profile_bandwidth_pair_2_half_minus_full") ||
      !identical(fit$components$score_definition, fit$score_definition)) {
    stop("A converged profile_pair_20261022 fit is required.")
  }
  backend <- profile_pair_precision_backend_v2(solver_opts)
  if (identical(backend, "clarabel_native") &&
      (!"CLARABEL" %in% toupper(solver_preference) ||
       !requireNamespace("clarabel", quietly = TRUE) ||
       !requireNamespace("Matrix", quietly = TRUE))) {
    profile_pair_stop_v2("dependency", "Paired inference requires CLARABEL and Matrix.", fit)
  }
  if (identical(backend, "dual_warm") && !requireNamespace("highs", quietly = TRUE)) {
    profile_pair_stop_v2("dependency", "Paired dual inference requires highs.", fit)
  }
  if (length(ci_level) != 1L || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be a scalar in (0, 1).")
  }
  domain <- profile_pair_mu_domain_v2(mu_grid, fit)
  fit$inference_mu_domain <- domain
  index <- profile_pair_coordinate_index_v2(coordinates, names(fit$beta))
  if (!all(names(fit$beta)[index] %in% fit$target_coordinates)) {
    stop("Inference coordinates must be prespecified in the paired fit.")
  }
  components <- fit$components
  if (!identical(fit$beta, fit$beta_refit) ||
      !identical(components$beta, fit$beta_refit)) {
    stop("Paired components must be evaluated at the stored beta_refit.")
  }
  if (any(!is.finite(components$hessian)) ||
      max(abs(components$hessian - t(components$hessian))) > 1e-10) {
    stop("The paired Jacobian must be finite and symmetric.")
  }
  rows <- directions <- cv_results <- vector("list", length(index))
  for (j in seq_along(index)) {
    coordinate <- index[j]
    cv <- profile_pair_mu_cv_v2(components, coordinate,
                                domain$mu_grid_requested, solver_opts, fit)
    direction <- profile_pair_solve_row_v2(components$hessian, coordinate,
                                          cv$mu, solver_opts, fit)
    if (!isTRUE(direction$feasible)) {
      profile_pair_stop_v2("precision_solver", "The selected full-data row is infeasible.",
                          fit, list(coordinate = coordinate, cv = cv, direction = direction))
    }
    projection <- as.numeric(components$cluster_scores %*% direction$omega)
    variance <- mean((projection - mean(projection))^2)
    if (!is.finite(variance) || variance <= 0) {
      profile_pair_stop_v2("variance_estimation",
                          "The paired cluster-score variance must be positive.",
                          fit, list(coordinate = coordinate, cv = cv, direction = direction))
    }
    interval <- profile_bandwidth_pair_wald_v2(
      components, fit$beta_refit,
      matrix(direction$omega, ncol = 1L), coordinate, alpha = 1 - ci_level
    )
    if (isTRUE(df_correction)) {
      interval$se <- interval$se * sqrt(components$n_clusters /
                                        (components$n_clusters - 1))
      critical <- stats::qnorm((1 + ci_level) / 2)
      interval$ci_lower <- interval$beta_tilde - critical * interval$se
      interval$ci_upper <- interval$beta_tilde + critical * interval$se
    }
    constituent <- cluster_sandwich_coordinate_v2(
      components$components_h$cluster_scores, direction$omega,
      df_correction = df_correction
    )
    rows[[j]] <- data.frame(
      coordinate = interval$coordinate, index = coordinate, feasible = TRUE,
      beta_hat = fit$beta_refit[coordinate], beta_tilde = interval$beta_tilde,
      se = interval$se, se_smoothed = constituent$se,
      ci_lower = interval$ci_lower, ci_upper = interval$ci_upper,
      mu = direction$mu, mu_source = "fold_cv",
      mu_cv_min_defect = cv$min_defect, mu_cv_defect = cv$defect_mean,
      mu_cv_quad = cv$quad_mean, mu_cv_folds_used = cv$folds_used,
      dantzig_residual = direction$residual, omega_l1 = direction$l1_norm,
      omega_l2 = direction$l2_norm, adjacent_stability = NA_real_,
      solver = direction$solver, solver_status = direction$solver_status,
      solver_static_regularization_constant =
        direction$solver_options$control$static_regularization_constant,
      n_clusters = components$n_clusters, h_inf = fit$h_inf,
      h_half = fit$h_inf / 2, score_definition = fit$score_definition,
      method_version = fit$method_version, stringsAsFactors = FALSE
    )
    directions[[j]] <- direction
    cv_results[[j]] <- cv
  }
  names(directions) <- names(cv_results) <- names(fit$beta)[index]
  list(table = do.call(rbind, rows), directions = directions, cv = cv_results,
       fit = fit, ci_level = ci_level,
       mu_grid = sort(unique(domain$mu_grid_requested)),
       mu_grid_requested = domain$mu_grid_requested,
       mu_grid_admissible = domain$mu_grid_admissible,
       excluded_zero_direction = domain$excluded_zero_direction,
       solver = if (identical(backend, "dual_warm")) "HiGHS-dual-simplex" else "CLARABEL-native",
       solver_preference = solver_preference,
       solver_options = profile_pair_solver_options_v2(solver_opts),
       score_definition = fit$score_definition, method_version = fit$method_version)
}

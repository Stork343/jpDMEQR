sqr_source_tuning_stop_v2 <- function(message, diagnostics) {
  stop(structure(list(message = message, call = NULL, failure_stage = "tuning",
                       method_version = "source_20261023_CV",
                       score_definition = "iid_Gaussian_smoothed_quantile_gradient",
                       diagnostics = list(tuning = diagnostics)),
                 class = c("sqr_source_tuning_error", "error", "condition")))
}

sqr_cv_precision_record_v2 <- function(precision, coordinates) {
  matrix_record <- precision$matrix_grid_record
  certificate <- matrix_record$feasibility_certificate
  certificate$direction_matrix <- certificate$residual_matrix <- NULL
  record <- list(feasible = precision$feasible, gamma = precision$gamma,
                   requested_gamma = precision$gamma_grid_requested,
                   attempts = precision$attempts, backend = precision$backend,
                   actual_LP_attempts = matrix_record$attempts,
                   solver_control = matrix_record$solver_control,
                   basis_policy = matrix_record$basis_policy,
                   feasibility_certificate = certificate,
                   precision_definition = precision$precision_definition)
  if (!isTRUE(precision$feasible)) return(record)
  record$matrix_dimension <- nrow(precision$omega_raw)
  record$raw_row_l1 <- rowSums(abs(precision$omega_raw))
  record$symmetric_row_l1 <- rowSums(abs(precision$omega_symmetric))
  record$omega_raw_targets <- precision$omega_raw[coordinates, , drop = FALSE]
  record$omega_symmetric_targets <- precision$omega_symmetric[coordinates, , drop = FALSE]
  record$symmetric_frobenius <- sqrt(sum(precision$omega_symmetric^2))
  record$symmetry_error <- max(abs(precision$omega_symmetric - t(precision$omega_symmetric)))
  record
}

select_sqr_gamma_cv_v2 <- function(train, tau, target_coords, config, tuning, seed,
                                    lambda_selector, fold_tuning_builder, control = list()) {
  if (isTRUE(control$exact_cv_pruning)) return(submission_sqr_gamma_cv(train, tau, target_coords,
    config, tuning, seed, lambda_selector, fold_tuning_builder, control))
  if (!is.function(lambda_selector) || !is.function(fold_tuning_builder)) {
    stop("SQR gamma CV requires the declared inner lambda selector and fold-tuning builder.")
  }
  multipliers <- sort(unique(as.numeric(config$dantzig_multipliers)))
  configured <- as.numeric(config$dantzig_multipliers)
  if (!length(multipliers) || any(!is.finite(multipliers)) || any(multipliers <= 0) ||
      length(tuning$mu_grid) != length(configured)) stop("The complete frozen gamma multiplier grid is required.")
  coordinates <- if (is.character(target_coords)) match(target_coords, colnames(train$X)) else
    as.integer(target_coords)
  if (!length(coordinates) || anyNA(coordinates) || any(coordinates < 1L | coordinates > ncol(train$X))) {
    stop("Unknown SQR target coordinate.")
  }
  role_seed <- seed_from_id_v2(seed, config$experiment_id, 1L, "sqr-gamma-subject-folds")
  folds <- make_subject_folds_v2(train$cluster_id, k = 5L, seed = role_seed)
  diagnostic <- list(method_version = "source_20261023_CV",
                     selection_rule = "five_subject_folds_mean_full_matrix_inverse_defect",
                     loss_definition = "squared_Frobenius_H_validation_W_sym_minus_I_divided_by_p",
                     configured_multipliers = configured, candidate_multipliers = multipliers,
                     folds = folds, fold_seed = role_seed,
                     roles = vector("list", 5L), inner_lambda_cv = vector("list", 5L),
                     fold_initial_estimators = vector("list", 5L),
                     cells = vector("list", 5L * length(multipliers)),
                     precision_records = vector("list", 5L * length(multipliers)),
                     selected_multiplier = NA_real_, selected_gamma = NA_real_)
  row <- 0L
  for (fold in seq_len(5L)) {
    validation_subjects <- as.character(folds$cluster_id[folds$fold == fold])
    training_subjects <- as.character(folds$cluster_id[folds$fold != fold])
    diagnostic$roles[[fold]] <- list(fold = fold, training_subjects = training_subjects,
                                     validation_subjects = validation_subjects)
    training <- subset_clusters_v2(train, training_subjects)
    validation <- subset_clusters_v2(train, validation_subjects)
    fold_seed <- seed_from_id_v2(seed, config$experiment_id, fold, "sqr-gamma-inner-lambda")
    fold_tuning <- fold_tuning_builder(training)
    if (length(fold_tuning$mu_grid) != length(configured)) stop("The fold gamma grid changed its declared dimension.")
    initial_warnings <- character()
    prepared <- tryCatch(withCallingHandlers({
      selected <- lambda_selector(training, fold_seed)
      N <- length(training$y)
      p <- ncol(training$X)
      h <- max(sqrt(tau * (1 - tau)) * sqrt(log(max(p, 2L))) * N^(-3 / 10), .05)
      multiplier <- selected$tuning$lambda_beta_multiplier
      lambda <- multiplier * sqrt(tau * (1 - tau) * log(max(p, 2L)) / N)
      fit <- fit_sqr_iid_lasso_v2(training$y, training$X, tau, h, lambda,
                                  penalty_factor = control$penalty_factor %||% rep(1, p),
                                  control = control$fit_control %||% list())
      if (!isTRUE(fit$converged)) {
        condition <- simpleError("The SQR outer-training initial estimator did not converge.")
        condition$failure_stage <- "penalised_fit"
        condition$fit <- fit
        condition$diagnostics <- list(inner_lambda_cv = selected$diagnostics)
        stop(condition)
      }
      list(fit = fit, selected = selected,
           validation_hessian = sqr_iid_hessian_fn_v2(validation$y, validation$X, fit$beta, tau, h))
    }, warning = function(condition) {
      initial_warnings <<- c(initial_warnings, conditionMessage(condition))
    }), error = function(condition) list(error = condition))
    if (!is.null(prepared$error)) {
      diagnostic$inner_lambda_cv[[fold]] <- list(status = "failed",
        failure_stage = prepared$error$failure_stage %||% "inner_lambda_tuning",
        message = conditionMessage(prepared$error), diagnostics = prepared$error$diagnostics,
        condition = prepared$error, warning_messages = initial_warnings)
    } else {
      diagnostic$inner_lambda_cv[[fold]] <- prepared$selected$diagnostics
      diagnostic$fold_initial_estimators[[fold]] <- list(
        beta = prepared$fit$beta, score = prepared$fit$score, h = prepared$fit$h,
        lambda = prepared$fit$lambda, lambda_multiplier = prepared$selected$tuning$lambda_beta_multiplier,
        converged = prepared$fit$converged, iterations = prepared$fit$iterations,
        kkt_residual = prepared$fit$kkt_residual, history = prepared$fit$history,
        n_training = length(training_subjects), N_training = length(training$y),
        n_validation = length(validation_subjects), N_validation = length(validation$y),
        warning_messages = initial_warnings)
    }
    for (candidate in seq_along(multipliers)) {
      row <- row + 1L
      multiplier <- multipliers[candidate]
      gamma <- fold_tuning$mu_grid[match(multiplier, configured)]
      status <- "failed"
      failure_stage <- failure_message <- solver_status <- ""
      loss <- NA_real_
      native_warnings <- initial_warnings
      precision_record <- NULL
      if (!is.finite(gamma) || gamma <= 0 || gamma >= 1) {
        status <- "excluded_zero_direction"
        failure_stage <- "gamma_domain"
        failure_message <- "The fold gamma is outside (0,1)."
      } else if (!is.null(prepared$error)) {
        failure_stage <- prepared$error$failure_stage %||% "inner_lambda_tuning"
        failure_message <- conditionMessage(prepared$error)
      } else {
        computed <- tryCatch(withCallingHandlers({
          precision <- solve_clime_matrix_iid_v2(
            prepared$fit$hessian, gamma,
            solver_preference = control$solver_preference %||% c("CLARABEL", "ECOS", "SCS"),
            solver_opts = control$solver_opts %||% list(), coordinates = coordinates
          )
          record <- sqr_cv_precision_record_v2(precision, coordinates)
          if (!isTRUE(precision$feasible)) return_value <- list(record = record, loss = NA_real_) else {
            defect <- profile_dense_crossprod_v2(prepared$validation_hessian, precision$omega_symmetric) -
              diag(ncol(training$X))
            value <- sum(defect^2) / ncol(training$X)
            if (!is.finite(value)) stop("The heldout full-matrix SQR inverse defect is non-finite.")
            return_value <- list(record = record, loss = value)
          }
          return_value
        }, warning = function(condition) {
          native_warnings <<- c(native_warnings, conditionMessage(condition))
        }), error = function(condition) list(error = condition))
        if (!is.null(computed$error)) {
          failure_stage <- computed$error$failure_stage %||% "precision_solver"
          failure_message <- conditionMessage(computed$error)
          solver_status <- computed$error$diagnostics$status %||% "failed"
          precision_record <- list(status = "failed", condition = computed$error,
                                     failure_diagnostics = computed$error$diagnostics)
        } else {
          precision_record <- computed$record
          loss <- computed$loss
          if (isTRUE(precision_record$feasible) && is.finite(loss)) status <- "ok" else {
            failure_stage <- "precision_solver"
            failure_message <- "No feasible complete SQR precision matrix was obtained at the fold gamma."
          }
          solver_status <- if (status == "ok") "optimal_all_matrix_rows" else "infeasible_complete_matrix"
        }
      }
      diagnostic$precision_records[[row]] <- list(candidate = candidate, fold = fold,
                                                   multiplier = multiplier, gamma = gamma,
                                                   result = precision_record,
                                                   warning_messages = native_warnings)
      diagnostic$cells[[row]] <- data.frame(
        candidate = candidate, multiplier = multiplier, fold = fold, gamma = gamma,
        n_training = length(training_subjects), N_training = length(training$y),
        n_validation = length(validation_subjects), N_validation = length(validation$y),
        h_source = if (is.null(prepared$error)) prepared$fit$h else NA_real_,
        h_gamma_anchor = fold_tuning$h_inf,
        gamma_anchor = gamma / multiplier,
        lambda_multiplier = if (is.null(prepared$error)) prepared$selected$tuning$lambda_beta_multiplier else NA_real_,
        lambda_source = if (is.null(prepared$error)) prepared$fit$lambda else NA_real_,
        status = status, solver_status = solver_status, failure_stage = failure_stage,
        failure_message = failure_message, loss = loss,
        warning_count = length(native_warnings), stringsAsFactors = FALSE)
    }
  }
  diagnostic$cells <- do.call(rbind, diagnostic$cells)
  diagnostic$candidates <- do.call(rbind, lapply(seq_along(multipliers), function(candidate) {
    rows <- diagnostic$cells[diagnostic$cells$candidate == candidate, , drop = FALSE]
    eligible <- nrow(rows) == 5L && all(rows$status == "ok") && all(is.finite(rows$loss))
    data.frame(candidate = candidate, multiplier = multipliers[candidate], eligible = eligible,
                 successful_folds = sum(rows$status == "ok"),
                 loss = if (eligible) mean(rows$loss) else NA_real_, selected = FALSE)
  }))
  eligible <- which(diagnostic$candidates$eligible)
  if (!length(eligible)) sqr_source_tuning_stop_v2("No eligible complete five-fold SQR gamma candidate.", diagnostic)
  minimum <- min(diagnostic$candidates$loss[eligible])
  tolerance <- 1e-12 * max(1, abs(minimum))
  ties <- eligible[diagnostic$candidates$loss[eligible] <= minimum + tolerance]
  chosen <- ties[which.min(diagnostic$candidates$multiplier[ties])]
  selected_multiplier <- diagnostic$candidates$multiplier[chosen]
  selected_gamma <- tuning$mu_grid[match(selected_multiplier, configured)]
  diagnostic$selected_multiplier <- selected_multiplier
  diagnostic$selected_gamma <- selected_gamma
  diagnostic$selected_candidate <- chosen
  diagnostic$tie_tolerance <- tolerance
  diagnostic$candidates$selected[chosen] <- TRUE
  list(gamma = selected_gamma, multiplier = selected_multiplier, diagnostics = diagnostic)
}

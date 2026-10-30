submission_sqr_gamma_cv <- function(train, tau, target_coords, config, tuning, seed,
                                     lambda_selector, fold_tuning_builder, control) {
  multipliers <- sort(unique(as.numeric(config$dantzig_multipliers)))
  configured <- as.numeric(config$dantzig_multipliers)
  coordinates <- if (is.character(target_coords)) match(target_coords, colnames(train$X)) else as.integer(target_coords)
  if (!length(coordinates) || anyNA(coordinates) || !length(multipliers) ||
      length(tuning$mu_grid) != length(configured)) stop("Invalid SQR tuning definition.")
  role_seed <- seed_from_id_v2(seed, config$experiment_id, 1L, "sqr-gamma-subject-folds")
  folds <- make_subject_folds_v2(train$cluster_id, k = 5L, seed = role_seed)
  diagnostic <- list(method_version = "source_20261023_CV",
    selection_rule = "five_subject_folds_mean_full_matrix_inverse_defect",
    execution = if (identical(control$solver_opts$algorithm, "compiled_certified"))
      "certified_column_and_fold_loss_bounds" else "certified_nonnegative_fold_loss_pruning",
    loss_definition = "squared_Frobenius_H_validation_W_sym_minus_I_divided_by_p",
    configured_multipliers = configured, candidate_multipliers = multipliers, folds = folds, fold_seed = role_seed,
    roles = vector("list", 5L), inner_lambda_cv = vector("list", 5L),
    fold_initial_estimators = vector("list", 5L))
  prepared <- vector("list", 5L)
  preparation_started <- proc.time()[[3L]]
  for (fold in 1:5) {
    validation_subjects <- as.character(folds$cluster_id[folds$fold == fold])
    training_subjects <- as.character(folds$cluster_id[folds$fold != fold])
    diagnostic$roles[[fold]] <- list(fold = fold, training_subjects = training_subjects, validation_subjects = validation_subjects)
    training <- subset_clusters_v2(train, training_subjects)
    validation <- subset_clusters_v2(train, validation_subjects)
    fold_seed <- seed_from_id_v2(seed, config$experiment_id, fold, "sqr-gamma-inner-lambda")
    fold_tuning <- fold_tuning_builder(training)
    selected <- lambda_selector(training, fold_seed)
    N <- length(training$y)
    p <- ncol(training$X)
    h <- max(sqrt(tau * (1 - tau)) * sqrt(log(max(p, 2L))) * N^(-.3), .05)
    multiplier <- selected$tuning$lambda_beta_multiplier
    lambda <- multiplier * sqrt(tau * (1 - tau) * log(max(p, 2L)) / N)
    fitted <- fit_sqr_iid_lasso_v2(training$y, training$X, tau, h, lambda,
      penalty_factor = control$penalty_factor %||% rep(1, p), control = control$fit_control %||% list())
    diagnostic$inner_lambda_cv[[fold]] <- selected$diagnostics
    diagnostic$fold_initial_estimators[[fold]] <- list(beta = fitted$beta, score = fitted$score, h = fitted$h,
      lambda = fitted$lambda, lambda_multiplier = multiplier, converged = fitted$converged,
      iterations = fitted$iterations, kkt_residual = fitted$kkt_residual, history = fitted$history,
      n_training = length(training_subjects), N_training = N,
      n_validation = length(validation_subjects), N_validation = length(validation$y), warning_messages = character())
    if (!fitted$converged) sqr_source_tuning_stop_v2("An SQR fold initial estimator did not converge.", diagnostic)
    prepared[[fold]] <- list(H = fitted$hessian, path_cache = new.env(parent = emptyenv()),
      validation_H = sqr_iid_hessian_fn_v2(validation$y, validation$X, fitted$beta, tau, h), tuning = fold_tuning,
      n_training = length(training_subjects), N_training = N, n_validation = length(validation_subjects),
      N_validation = length(validation$y), h = h, lambda = lambda, lambda_multiplier = multiplier)
  }
  candidate_records <- vector("list", length(multipliers))
  diagnostic$initialisation_seconds <- proc.time()[[3L]] - preparation_started
  precision_started <- proc.time()[[3L]]
  cells <- precision_records <- vector("list", 5L * length(multipliers))
  incumbent <- Inf
  for (candidate in rev(seq_along(multipliers))) {
    accumulated <- 0
    successes <- 0L
    excluded <- failed <- pruned <- FALSE
    incumbent_at_pruning <- NA_real_
    for (fold in 1:5) {
      index <- (fold - 1L) * length(multipliers) + candidate
      item <- prepared[[fold]]
      gamma <- item$tuning$mu_grid[match(multipliers[candidate], configured)]
      status <- "failed"
      stage <- message <- solver_status <- ""
      loss <- NA_real_
      record <- NULL
      native_warnings <- character()
      if (gamma <= 0 || gamma >= 1 || !is.finite(gamma)) {
        excluded <- TRUE
        status <- "excluded_zero_direction"
        stage <- "gamma_domain"
      } else if (failed) {
        status <- "not_required_after_failed_fold"
      } else if (pruned) {
        status <- "pruned_by_loss_lower_bound"
      } else {
        # 单个gamma的数值失败使该候选不满足完整五折规则。
        computed <- tryCatch(withCallingHandlers({
          numerical_control <- control$solver_opts %||% list()
          numerical_control$path_cache <- item$path_cache
          if (identical(numerical_control$algorithm, "compiled_certified") &&
              isTRUE(numerical_control$column_loss_bounds)) {
            threshold <- incumbent + 1e-12 * max(1, abs(incumbent))
            budget <- if (is.finite(threshold)) max(0, 5 * threshold - accumulated) else Inf
            submission_clime_bounded_cv(item$H, gamma, item$validation_H, budget,
              numerical_control, coordinates)
          } else {
          precision <- solve_clime_matrix_iid_v2(item$H, gamma,
            solver_preference = control$solver_preference %||% c("CLARABEL", "ECOS", "SCS"),
            solver_opts = numerical_control, coordinates = coordinates)
          precision_record <- sqr_cv_precision_record_v2(precision, coordinates)
          if (!isTRUE(precision$feasible)) list(record = precision_record, feasible = FALSE) else {
            value <- if (identical(numerical_control$algorithm, "compiled_certified")) {
              submission_sqr_sparse_loss_cpp(item$validation_H, precision$omega_symmetric)
            } else {
              defect <- profile_dense_crossprod_v2(item$validation_H, precision$omega_symmetric) - diag(nrow(item$H))
              sum(defect^2) / nrow(item$H)
            }
            if (!is.finite(value) || value < 0) stop("Invalid SQR validation loss.")
            list(record = precision_record, feasible = TRUE, loss = value)
          }
          }
        }, warning = function(condition) {
          native_warnings <<- c(native_warnings, conditionMessage(condition))
        }), error = function(condition) list(condition = condition))
        if (!is.null(computed$condition)) {
          failed <- TRUE
          stage <- computed$condition$failure_stage %||% "precision_solver"
          message <- conditionMessage(computed$condition)
          solver_status <- computed$condition$diagnostics$status %||% "failed"
          record <- list(status = "failed", condition = computed$condition,
            failure_diagnostics = computed$condition$diagnostics)
        } else if (isTRUE(computed$pruned)) {
          record <- computed$record
          accumulated <- accumulated + computed$lower_bound
          pruned <- TRUE
          incumbent_at_pruning <- incumbent
          status <- "pruned_by_column_loss_lower_bound"
          solver_status <- "certified_nonnegative_column_loss_bound"
        } else {
          record <- computed$record
          if (!computed$feasible) {
            failed <- TRUE
            stage <- "precision_solver"
            solver_status <- "infeasible_complete_matrix"
          } else {
            status <- "ok"
            solver_status <- "optimal_all_matrix_rows"
            loss <- computed$loss
            accumulated <- accumulated + loss
            successes <- successes + 1L
            threshold <- incumbent + 1e-12 * max(1, abs(incumbent))
            if (fold < 5L && is.finite(incumbent) && accumulated / 5 > threshold) {
              pruned <- TRUE
              incumbent_at_pruning <- incumbent
            }
          }
        }
      }
      precision_records[[index]] <- list(candidate = candidate, fold = fold, multiplier = multipliers[candidate],
        gamma = gamma, result = record, warning_messages = native_warnings)
      cells[[index]] <- data.frame(candidate = candidate, multiplier = multipliers[candidate], fold = fold,
        gamma = gamma, n_training = item$n_training, N_training = item$N_training,
        n_validation = item$n_validation, N_validation = item$N_validation, h_source = item$h,
        h_gamma_anchor = item$tuning$h_inf, gamma_anchor = gamma / multipliers[candidate],
        lambda_multiplier = item$lambda_multiplier, lambda_source = item$lambda,
        status = status, solver_status = solver_status, failure_stage = stage, failure_message = message,
        loss = loss, warning_count = length(native_warnings), lower_bound = accumulated / 5)
    }
    eligible <- !excluded && !failed && !pruned && successes == 5L
    value <- if (eligible) accumulated / 5 else NA_real_
    if (eligible) incumbent <- min(incumbent, value)
    candidate_records[[candidate]] <- data.frame(candidate = candidate, multiplier = multipliers[candidate],
      eligible = eligible, successful_folds = successes, loss = value, selected = FALSE,
      certified_pruning = pruned, lower_bound = accumulated / 5, incumbent_at_pruning = incumbent_at_pruning)
  }
  diagnostic$cells <- do.call(rbind, cells)
  diagnostic$precision_records <- precision_records
  diagnostic$candidates <- do.call(rbind, candidate_records)
  eligible <- which(diagnostic$candidates$eligible)
  if (!length(eligible)) sqr_source_tuning_stop_v2("No eligible complete five-fold SQR gamma candidate.", diagnostic)
  minimum <- min(diagnostic$candidates$loss[eligible])
  tolerance <- 1e-12 * max(1, abs(minimum))
  ties <- eligible[diagnostic$candidates$loss[eligible] <= minimum + tolerance]
  chosen <- ties[which.min(diagnostic$candidates$multiplier[ties])]
  selected_multiplier <- diagnostic$candidates$multiplier[chosen]
  selected_gamma <- tuning$mu_grid[match(selected_multiplier, configured)]
  pruned <- diagnostic$candidates$certified_pruning
  if (any(diagnostic$candidates$lower_bound[pruned] <= minimum + tolerance))
    stop("An SQR pruning bound does not certify exclusion of its candidate.")
  diagnostic$candidates$selected[chosen] <- TRUE
  diagnostic$selected_multiplier <- selected_multiplier
  diagnostic$selected_gamma <- selected_gamma
  diagnostic$selected_candidate <- chosen
  diagnostic$tie_tolerance <- tolerance
  diagnostic$precision_seconds <- proc.time()[[3L]] - precision_started
  list(gamma = selected_gamma, multiplier = selected_multiplier, diagnostics = diagnostic)
}

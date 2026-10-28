# 冻结倍率使用训练 subject 的五折预测损失进行选择。
frozen_tuning_grid_v2 <- function(method_id, config, tuning) {
  paired <- identical(tuning$method_version, "profile_pair_20261022") &&
    method_id %in% c("PROFILE-DQR", "PROFILE-DQR-TRUE-SUPPORT")
  use_cv <- method_id %in% c("POOLED-QR-LASSO", "SQR-DEBIASED-IID") ||
    (paired && identical(config$experiment_id, "E07"))
  grid <- if (use_cv) config$lambda_beta_multipliers else
    tuning$lambda_multiplier %||% 1
  grid <- as.numeric(grid)
  if (!length(grid) || any(!is.finite(grid)) || any(grid <= 0)) {
    stop("Frozen penalty multipliers must be finite positive values.")
  }
  proposed <- identical(tuning$method_version, "profile_pair_20261022") &&
    method_id %in% c("PROFILE-DQR", "PROFILE-DQR-TRUE-SUPPORT", "PROFILE-DQR-TRUE-NUISANCE",
                     "PROFILE-DQR-POP-H", "POSTREFIT-EXACT-H", "PROFILE-DQR-SPLIT")
  list(use_cv = use_cv, proposed = proposed, configured = grid, values = sort(unique(grid)))
}

frozen_tuning_fit_summary_v2 <- function(fit) {
  if (is.null(fit)) return(NULL)
  fields <- c("beta", "beta_l1", "beta_refit", "selected_support", "refit_set",
              "converged", "kkt_residual", "beta_change", "iterations",
              "lambda_components", "h", "h_est", "h_inf", "lambda", "N",
              "post_refit_gradient_max", "post_refit_nuisance_gradient_max",
              "post_refit_hessian_min_eigenvalue", "post_refit_hessian_condition_number",
              "history", "score", "failure_stage", "opt", "component_convergence", "bootstrap_success", "B")
  out <- fit[intersect(fields, names(fit))]
  if (!is.null(fit$components)) {
    out$profile_identity_error <- fit$components$profile_identity_error
    out$max_nuisance_gradient <- fit$components$max_nuisance_gradient
  }
  out
}

fit_frozen_tuning_training_v2 <- function(method_id, train, tau, target_coords,
                                          tuning, multiplier, seed,
                                          context = list(), control = list()) {
  set.seed(seed)
  tuning$lambda_multiplier <- multiplier
  tuning$lambda_beta_multiplier <- multiplier
  if (method_id == "POOLED-QR-LASSO") {
    tuning$lambda_beta <- if (identical(context$config$method_version, "joint_profile_design_20261026")) {
      multiplier * sqrt(tau * (1 - tau) * length(train$y) * log(max(ncol(train$X), 2L)))
    } else multiplier * tuning$lambda_beta
    fit <- fit_benchmark_pooled_qr_lasso_v2(
      train, tau, target_coords, tuning, seed, control
    )
    return(list(status = fit$status, beta = fit$beta_hat,
                failure_stage = fit$failure_stage, failure_message = fit$failure_message,
                fit_summary = list(beta = fit$beta_hat, converged = fit$converged,
                                   lambda = tuning$lambda_beta,
                                   criterion = "pooled_check_loss_sum_plus_raw_l1")))
  }
  if (method_id == "SQR-DEBIASED-IID") {
    N <- length(train$y)
    p <- ncol(train$X)
    h <- max(sqrt(tau * (1 - tau)) * sqrt(log(max(p, 2))) * N^(-3 / 10), 0.05)
    lambda <- multiplier * sqrt(tau * (1 - tau) * log(max(p, 2)) / N)
    prediction_control <- control$fit_control %||% list()
    prediction_control$compute_hessian <- FALSE
    fit <- fit_sqr_iid_lasso_v2(
      train$y, train$X, tau, h, lambda,
      penalty_factor = control$penalty_factor %||% rep(1, p),
      control = prediction_control
    )
    return(list(status = if (isTRUE(fit$converged)) "ok" else "failed", beta = fit$beta,
                failure_stage = if (isTRUE(fit$converged)) "" else "penalised_fit",
                failure_message = if (isTRUE(fit$converged)) "" else
                  paste("SQR initial fit did not converge:", fit$failure_stage),
                fit_summary = frozen_tuning_fit_summary_v2(fit)))
  }
  if (!identical(tuning$method_version, "profile_pair_20261022") ||
      !method_id %in% c("PROFILE-DQR", "PROFILE-DQR-TRUE-SUPPORT")) {
    stop("Unsupported frozen prediction-CV method: ", method_id)
  }
  if (method_id == "PROFILE-DQR-TRUE-SUPPORT") {
    target_idx <- if (is.character(target_coords)) match(target_coords, colnames(train$X)) else
      as.integer(target_coords)
    if (anyNA(target_idx)) stop("Unknown oracle-support CV target coordinate.")
    keep <- sort(unique(c(context$active %||% stop("context$active is required"), target_idx)))
    target_coords <- colnames(train$X)[target_idx]
    train$X <- train$X[, keep, drop = FALSE]
    control$penalty_factor <- rep(0, length(keep))
  }
  control$shared_first_stage <- NULL
  fit <- profile_pair_fit_from_data_v2(train, tau, target_coords, tuning, control)
  list(status = "ok", beta = fit$beta, failure_stage = "", failure_message = "",
       fit_summary = frozen_tuning_fit_summary_v2(fit))
}

frozen_tuning_cluster_loss_v2 <- function(validation, beta, tau) {
  if (is.null(beta) || is.null(names(beta)) || any(!is.finite(beta)) ||
      !all(names(beta) %in% colnames(validation$X))) {
    stop("Prediction CV requires finite named training coefficients.")
  }
  residual <- validation$y - as.numeric(validation$X[, names(beta), drop = FALSE] %*% beta)
  loss <- residual * (tau - as.numeric(residual < 0))
  cluster_loss <- vapply(split(loss, as.character(validation$cluster_id)), mean, numeric(1))
  if (!length(cluster_loss) || any(!is.finite(cluster_loss))) {
    stop("Prediction CV produced invalid heldout cluster losses.")
  }
  cluster_loss
}

select_frozen_method_tuning_v2 <- function(method_id, train, tau, target_coords,
                                           config, tuning, seed,
                                           context = list(), control = list()) {
  grid <- frozen_tuning_grid_v2(method_id, config, tuning)
  diagnostic <- list(
    selection_rule = if (grid$use_cv) "five_fold_training_subject_mean_pinball" else
      if (grid$proposed) "configured_fixed_multiplier" else "source_method_tuning",
    configured_multipliers = grid$configured, candidate_multipliers = grid$values,
    selected_multiplier = if (grid$use_cv || !grid$proposed) NA_real_ else grid$values[[1]],
    method_id = method_id, training_subjects = unique(as.character(train$cluster_id)),
    validation_loss_weight = "equal_subject", folds = NULL, roles = NULL,
    fold_losses = NULL, candidates = NULL, training_fits = NULL
  )
  if (!grid$use_cv) return(list(tuning = tuning, diagnostics = diagnostic))
  role_seed <- seed_from_id_v2(seed, config$experiment_id, 1L, "penalty-cv-subject-folds")
  folds <- make_subject_folds_v2(train$cluster_id, k = 5L, seed = role_seed)
  diagnostic$fold_seed <- role_seed
  diagnostic$folds <- folds
  diagnostic$roles <- lapply(seq_len(5L), function(fold) list(
    fold = fold,
    training_subjects = folds$cluster_id[folds$fold != fold],
    validation_subjects = folds$cluster_id[folds$fold == fold]
  ))
  records <- vector("list", length(grid$values) * 5L)
  training_fits <- vector("list", length(records))
  candidate_records <- vector("list", length(grid$values))
  row <- 0L
  for (candidate in seq_along(grid$values)) {
    multiplier <- grid$values[[candidate]]
    for (fold in seq_len(5L)) {
      row <- row + 1L
      role <- diagnostic$roles[[fold]]
      training <- subset_clusters_v2(train, role$training_subjects)
      validation <- subset_clusters_v2(train, role$validation_subjects)
      fold_tuning <- build_tuning_for_analysis_v2(config, training)
      fold_seed <- seed_from_id_v2(seed, config$experiment_id, fold, "penalty-cv-training-fit")
      computed <- tryCatch({
        fitted <- fit_frozen_tuning_training_v2(
          method_id, training, tau, target_coords, fold_tuning, multiplier,
          fold_seed, context, control
        )
        if (!identical(fitted$status, "ok")) {
          list(fitted = fitted, losses = NULL)
        } else {
          list(fitted = fitted,
               losses = frozen_tuning_cluster_loss_v2(validation, fitted$beta, tau))
        }
      }, error = function(e) list(
        fitted = list(status = "failed", beta = NULL,
                      failure_stage = e$failure_stage %||% "tuning_fit",
                      failure_message = conditionMessage(e),
                      fit_summary = frozen_tuning_fit_summary_v2(e$fit),
                      failure_diagnostics = e$diagnostics), losses = NULL
      ))
      fitted <- computed$fitted
      training_fits[[row]] <- list(candidate = candidate, fold = fold,
                                    fit = fitted$fit_summary,
                                    failure_diagnostics = fitted$failure_diagnostics)
      losses <- computed$losses
      success <- identical(fitted$status, "ok") && !is.null(losses)
      records[[row]] <- data.frame(
        candidate = candidate, multiplier = multiplier, fold = fold, seed = fold_seed,
        n_training = training$n_clusters, N_training = length(training$y),
        n_validation = validation$n_clusters,
        h_est = fold_tuning$h_est, h_inf = fold_tuning$h_inf,
        raw_lambda_anchor = fold_tuning$lambda_beta,
        h_source = fitted$fit_summary[["h", exact = TRUE]] %||% NA_real_,
        lambda_source = fitted$fit_summary[["lambda", exact = TRUE]] %||% NA_real_,
        status = if (success) "ok" else "failed",
        failure_stage = fitted$failure_stage %||% "",
        failure_message = fitted$failure_message %||% "",
        loss_sum = if (success) sum(losses) else NA_real_,
        loss_mean = if (success) mean(losses) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    d <- do.call(rbind, records[(row - 4L):row])
    eligible <- all(d$status == "ok") && all(is.finite(d$loss_sum))
    candidate_records[[candidate]] <- data.frame(
      candidate = candidate, multiplier = multiplier, eligible = eligible,
      successful_folds = sum(d$status == "ok"),
      loss = if (eligible) sum(d$loss_sum) / sum(d$n_validation) else NA_real_,
      selected = FALSE,
      stringsAsFactors = FALSE
    )
  }
  diagnostic$fold_losses <- do.call(rbind, records)
  diagnostic$training_fits <- training_fits
  diagnostic$candidates <- do.call(rbind, candidate_records)
  eligible <- which(diagnostic$candidates$eligible)
  if (!length(eligible)) {
    stop(structure(list(
      message = paste("No eligible five-fold penalty candidate for", method_id),
      call = NULL, failure_stage = "tuning", diagnostics = list(tuning = diagnostic),
      method_version = if (method_id == "SQR-DEBIASED-IID") "source_20261023_CV" else
        if (method_id %in% c("PROFILE-DQR", "PROFILE-DQR-TRUE-SUPPORT")) tuning$method_version else NULL,
      score_definition = if (method_id == "SQR-DEBIASED-IID") "iid_Gaussian_smoothed_quantile_gradient" else
        if (method_id %in% c("PROFILE-DQR", "PROFILE-DQR-TRUE-SUPPORT")) "profile_bandwidth_pair_2_half_minus_full" else NULL
    ), class = c("frozen_tuning_error", "error", "condition")))
  }
  selected <- eligible[which.min(diagnostic$candidates$loss[eligible])]
  multiplier <- diagnostic$candidates$multiplier[[selected]]
  diagnostic$selected_multiplier <- multiplier
  diagnostic$selected_candidate <- selected
  diagnostic$candidates$selected <- seq_len(nrow(diagnostic$candidates)) == selected
  tuning$lambda_multiplier <- multiplier
  tuning$lambda_beta_multiplier <- multiplier
  if (method_id == "POOLED-QR-LASSO") tuning$lambda_beta <- if (identical(config$method_version, "joint_profile_design_20261026")) {
    multiplier * sqrt(tau * (1 - tau) * length(train$y) * log(max(ncol(train$X), 2L)))
  } else tuning$lambda_beta * multiplier
  list(tuning = tuning, diagnostics = diagnostic)
}

fit_benchmark_frozen_tuning_v2 <- function(method_id, train, tau, target_coords,
                                           tuning, seed, context = list(), control = list()) {
  if (startsWith(method_id, "JOINT-PROFILE")) {
    fitted <- submission_benchmark(method_id, train, tau, target_coords, tuning, seed, context, control)
    fitted$tuning_values <- list(h = tuning$h, h_inf = tuning$h, h_est = NA_real_,
      lambda_beta = NA_real_, lambda_gamma = tuning$lambda_gamma)
    fitted$tuning_diagnostics <- list(selection_rule = "independent_outcome_and_adjustment_learning",
      budget_multiplier = context$config$budget_multiplier, selected_multiplier = NA_real_,
      not_applicable = c("fixed_effect_lasso_penalty", "Dantzig_tolerance"))
    return(fitted)
  }
  working_intercept <- identical(context$config$method_version, "joint_profile_design_20261026") &&
    method_id %in% c("POOLED-QR-LASSO", "SQR-DEBIASED-IID", "QGEE-SCAD")
  original_names <- colnames(train$X)
  if (working_intercept) {
    train$X <- cbind(`(Intercept)` = 1, train$X)
    train$p <- ncol(train$X)
    control$penalty_factor <- c(0, rep(1, length(original_names)))
    control$intercept <- TRUE
  }
  tuning_started <- proc.time()[[3]]
  config <- context$config %||% stop("Frozen tuning requires context$config.")
  selected <- select_frozen_method_tuning_v2(
    method_id, train, tau, target_coords, config, tuning, seed, context, control
  )
  if (method_id == "SQR-DEBIASED-IID") {
    gamma <- tryCatch(
      select_sqr_gamma_cv_v2(
        train, tau, target_coords, config, selected$tuning, seed,
        lambda_selector = function(outer_training, outer_seed) {
          select_frozen_method_tuning_v2(
            method_id, outer_training, tau, target_coords, config,
            build_tuning_for_analysis_v2(config, outer_training), outer_seed, context, control
          )
        },
        fold_tuning_builder = function(outer_training) build_tuning_for_analysis_v2(config, outer_training),
        control = control
      ),
      error = function(condition) {
        condition$diagnostics <- modifyList(
          condition$diagnostics %||% list(),
          list(full_training_lambda_cv = selected$diagnostics)
        )
        stop(condition)
      }
    )
    selected$diagnostics$gamma_cv <- gamma$diagnostics
    selected$diagnostics$selected_gamma_multiplier <- gamma$multiplier
    selected$diagnostics$selected_gamma <- gamma$gamma
    selected$tuning$gamma_grid_requested <- selected$tuning$mu_grid
    selected$tuning$gamma_multipliers_requested <- config$dantzig_multipliers
    selected$tuning$mu_grid <- gamma$gamma
    selected$tuning$source_method_version <- "source_20261023_CV"
  }
  selected$diagnostics$elapsed_sec <- proc.time()[[3]] - tuning_started
  started <- proc.time()[[3]]
  fit <- tryCatch(
    fit_benchmark_by_id_v2(method_id, train, tau, target_coords, selected$tuning,
                           seed, context, control),
    error = function(e) {
      e$diagnostics <- modifyList(e$diagnostics %||% list(),
                                 list(tuning = selected$diagnostics))
      stop(e)
    }
  )
  fit$tuning_diagnostics <- selected$diagnostics
  if (method_id == "SQR-DEBIASED-IID") {
    fit$method_version <- "source_20261023_CV"
    fit$implementation_version <- "source_20261023_CV"
    fit$score_definition <- "iid_Gaussian_smoothed_quantile_gradient"
    fit$tuning_diagnostics$method_version <- "source_20261023_CV"
    fit$tuning_diagnostics$gamma_selection_rule <- gamma$diagnostics$selection_rule
  }
  fit$tuning_values <- selected$tuning
  fit$runtime_sec <- proc.time()[[3]] - tuning_started
  if (identical(method_id, "SQR-DEBIASED-IID") && !is.null(fit$fit_object)) {
    fit$tuning_values$lambda_beta <- fit$fit_object$lambda
    fit$tuning_values$h <- fit$fit_object$h
  }
  if (identical(method_id, "DOUBLE-PEN-QLMM")) {
    fit$tuning_diagnostics$configured_multipliers <- control$lambda_beta_grid
    fit$tuning_diagnostics$candidate_multipliers <- sort(unique(control$lambda_beta_grid))
    fit$tuning_diagnostics$lambda_alpha_multipliers <- control$lambda_alpha_grid
    source_fit <- fit$fit_object$fit
    if (!is.null(source_fit)) {
      fields <- c("lambda_beta", "lambda_alpha", "lambda_beta_grid", "lambda_alpha_grid",
                  "initialization_candidates", "tuning_candidates", "tuning_history")
      fit$tuning_diagnostics$source_block_tuning <- source_fit[fields]
      fit$tuning_values$lambda_beta <- source_fit$lambda_beta
      fit$tuning_diagnostics$selected_multiplier <- source_fit$lambda_beta /
        fit$fit_object$penalty_scales[["beta"]]
    }
  }
  fit$tuning_diagnostics$full_training <- list(
    subject_ids = unique(as.character(train$cluster_id)), seed = seed,
    multiplier = fit$tuning_diagnostics$selected_multiplier,
    raw_lambda_anchor = tuning$lambda_beta,
    lambda_beta = fit$tuning_values$lambda_beta,
    h = fit$tuning_values$h, beta_hat = fit$beta_hat,
    status = fit$status, failure_stage = fit$failure_stage,
    failure_message = fit$failure_message,
    fit = frozen_tuning_fit_summary_v2(fit$fit_object),
    elapsed_sec = proc.time()[[3]] - started
  )
  if (working_intercept && !is.null(fit$beta_hat)) {
    fit$prediction_intercept <- unname(fit$beta_hat["(Intercept)"])
    fit$working_dimension <- length(fit$beta_hat)
    fit$selected <- setdiff(names(fit$beta_hat)[fit$selected], "(Intercept)")
    fit$beta_hat <- fit$beta_hat[original_names]
  }
  fit
}

application_initial_fit_record_v2 <- function(fit) {
  fields <- c("beta_l1", "h_est", "method_version", "score_definition", "loadings",
              "lambda_stochastic", "lambda_mean", "lambda_base", "lambda_multiplier",
              "coordinate_penalty", "converged", "kkt_normalized", "kkt_residual",
              "iterations", "first_stage_beta_change", "profile_objective", "composite_objective",
              "initial_acceleration", "initial_optimizer_history")
  c(fit[fields], list(profile_identity_error = fit$components$profile_identity_error,
                      max_nuisance_gradient = fit$components$max_nuisance_gradient,
                      cluster_scores = fit$components$cluster_scores))
}

fit_application_prediction_candidate_v2 <- function(metadata, design, training_rows, test_rows,
                                                    tau, h_multiplier, lambda_gamma,
                                                    lambda_multiplier,
                                                    method = "PROFILE-INITIAL", seed = 1L,
                                                    initial_control = list()) {
  n <- length(unique(metadata$subject_id[training_rows]))
  h_est <- h_multiplier * (log(max(sum(design$penalty_factor > 0), 2L)) / n)^(1 / 4)
  started <- proc.time()[[3]]
  if (method == "PROFILE-INITIAL") {
    fit <- fit_application_profile_initial_v2(
      metadata$log_c3[training_rows], design$X[training_rows, , drop = FALSE],
      design$Z[training_rows, , drop = FALSE], metadata$subject_id[training_rows], tau,
      h_est = h_est, lambda_gamma = lambda_gamma, lambda_multiplier = lambda_multiplier,
      penalty_factor = design$penalty_factor, control = initial_control
    )
    record <- application_initial_fit_record_v2(fit)
    beta <- fit$beta_l1
    status <- "ok"
    converged <- fit$converged
    kkt <- fit$kkt_residual
    identity <- fit$components$profile_identity_error
    nuisance <- fit$components$max_nuisance_gradient
  } else {
    record <- fit_application_comparator_v2(method, metadata, design, training_rows,
                                             tau, lambda_multiplier, seed)
    beta <- record$beta_prediction
    status <- record$status
    converged <- record$converged
    kkt <- identity <- nuisance <- NA_real_
  }
  prediction <- as.numeric(design$X[test_rows, , drop = FALSE] %*% beta)
  y <- metadata$log_c3[test_rows]
  loss <- (y - prediction) * (tau - as.numeric(y < prediction))
  per_subject <- aggregate(
    cbind(loss, calibration = as.numeric(y <= prediction)),
    by = list(subject_id = as.character(metadata$subject_id[test_rows])), FUN = mean
  )
  train_size <- table(metadata$subject_id[training_rows])
  baseline <- quantreg::rq(metadata$log_c3[training_rows] ~ 1, tau = tau,
                            weights = 1 / as.numeric(train_size[as.character(metadata$subject_id[training_rows])]))
  baseline_prediction <- as.numeric(stats::coef(baseline))
  baseline_loss <- (y - baseline_prediction) * (tau - as.numeric(y < baseline_prediction))
  baseline_mean <- mean(tapply(baseline_loss, metadata$subject_id[test_rows], mean))
  list(
    record = record,
    metrics = data.frame(
      status = status, n_train_subjects = n, n_test_subjects = nrow(per_subject),
      pinball_loss = mean(per_subject$loss), intercept_pinball_loss = baseline_mean,
      intercept_skill = 1 - mean(per_subject$loss) / baseline_mean,
      quantile_calibration = mean(per_subject$calibration),
      selected_genes = sum(abs(beta[design$penalty_factor > 0]) > 1e-10),
      runtime_sec = proc.time()[[3]] - started, converged = converged,
      kkt_residual = kkt, profile_identity_error = identity,
      max_nuisance_gradient = nuisance, method = method,
      method_version = record$method_version,
      score_definition = if (method == "PROFILE-INITIAL") record$score_definition else record$target_type,
      stringsAsFactors = FALSE
    ),
    subject_metrics = per_subject,
    predictions = data.frame(sample_id = rownames(metadata)[test_rows],
                              subject_id = metadata$subject_id[test_rows],
                              observed_log_c3 = y, predicted_quantile = prediction, loss = loss,
                              stringsAsFactors = FALSE)
  )
}

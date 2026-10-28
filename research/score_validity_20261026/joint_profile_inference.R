ds_fit_profile <- function(y, X, Z, cluster_id, tau, coordinates, Lambda = diag(ncol(Z)),
                            penalty_factor = rep(1, ncol(X)),
                            h = length(unique(cluster_id))^(-0.3), confidence = 0.95,
                            keep_prepared = FALSE, budget_multiplier = 1) {
  started <- proc.time()[[3L]]
  validate_tau_h_v2(tau, h)
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  if (length(y) != nrow(X) || nrow(Z) != nrow(X) || length(cluster_id) != nrow(X) ||
      any(!is.finite(y)) || anyNA(cluster_id)) stop("Invalid observations.")
  if (is.null(colnames(X)) || anyNA(colnames(X)) || any(!nzchar(colnames(X))) || anyDuplicated(colnames(X))) {
    stop("X needs unique nonempty column names.")
  }
  if (is.character(coordinates)) coordinates <- match(coordinates, colnames(X))
  if (!length(coordinates) || anyNA(coordinates) || anyDuplicated(coordinates) ||
      any(coordinates < 1L | coordinates > ncol(X)) || any(coordinates != as.integer(coordinates))) stop("Invalid target coordinates.")
  if (length(penalty_factor) != ncol(X) || any(!penalty_factor %in% c(0, 1))) stop("Invalid inclusion indicators.")
  if (length(confidence) != 1L || !is.finite(confidence) || confidence <= 0 || confidence >= 1) stop("Invalid confidence level.")
  Lambda <- make_lambda_matrix_v2(Lambda, ncol(Z))
  prepared <- ds_prepare(y, X, Z, cluster_id, tau, coordinates, Lambda, penalty_factor, h,
                          include_design = TRUE, budget_multiplier = budget_multiplier)
  training_seconds <- proc.time()[[3L]] - started
  proportions <- vapply(prepared$folds, function(fold) length(fold$validation_subjects), integer(1)) / prepared$n
  details <- lapply(seq_along(coordinates), function(j) {
    fits <- lapply(prepared$folds, function(fold) {
      model <- fold$models[[j]]
      ds_pair_small(model$data, model$initial, tau, h, Lambda, gram = model$gram)
    })
    estimate <- sum(proportions * vapply(fits, `[[`, numeric(1), "estimate"))
    influence <- unlist(lapply(fits, function(fit) fit$influence - mean(fit$influence)), use.names = FALSE)
    se <- sqrt(sum(influence^2)) / prepared$n
    if (!is.finite(se) || se <= 0) stop("Nonpositive cluster variance.")
    critical <- qnorm((1 + confidence) / 2)
    list(estimate = estimate, se = se, lower = estimate - critical * se,
         upper = estimate + critical * se, fits = fits, influence = influence,
         max_kkt = max(vapply(fits, `[[`, numeric(1), "max_kkt")))
  })
  table <- data.frame(coordinate = colnames(X)[coordinates],
                       beta_tilde = vapply(details, `[[`, numeric(1), "estimate"),
                       se = vapply(details, `[[`, numeric(1), "se"),
                       lower = vapply(details, `[[`, numeric(1), "lower"),
                       upper = vapply(details, `[[`, numeric(1), "upper"),
                       max_kkt = vapply(details, `[[`, numeric(1), "max_kkt"))
  training <- lapply(prepared$folds, function(fold) fold[c("training_subjects", "validation_subjects", "initial", "graphs")])
  list(table = table, details = details, training = training, h = h, confidence = confidence,
       timing = c(training = training_seconds, inference = proc.time()[[3L]] - started - training_seconds,
                  total = proc.time()[[3L]] - started), method_version = "joint_profile_design_20261026",
       budget_multiplier = budget_multiplier, prepared = if (keep_prepared) prepared else NULL)
}

qif_see_stop_v2 <- function(stage, message, diagnostics, partial_fit = NULL) {
  condition <- simpleError(message)
  condition$failure_stage <- stage
  condition$diagnostics <- diagnostics
  condition$fit <- partial_fit
  class(condition) <- c("qif_see_source_error", class(condition))
  stop(condition)
}

qif_see_basis_v2 <- function(m, corstr = "combined_known_basis") {
  identity <- diag(m)
  exchangeable <- matrix(1, m, m) - identity
  adjacent <- matrix(0, m, m)
  if (m > 1L) {
    adjacent[cbind(seq_len(m - 1L), seq.int(2L, m))] <- 1
    adjacent[cbind(seq.int(2L, m), seq_len(m - 1L))] <- 1
  }
  switch(corstr, independence = list(identity),
         exchangeable = list(identity, exchangeable),
         combined_known_basis = list(identity, exchangeable, adjacent),
         stop("Unknown declared QIF basis specification."))
}

qif_see_extended_source_v2 <- function(y, X, cluster_id, beta, tau, h,
                                        density_weights, delta,
                                        corstr = "combined_known_basis") {
  clusters <- split(seq_along(y), factor(as.character(cluster_id), levels = unique(as.character(cluster_id))))
  p <- ncol(X)
  B <- length(qif_see_basis_v2(length(clusters[[1L]]), corstr))
  n <- length(clusters)
  G <- matrix(0, B * p, n)
  derivative <- matrix(0, B * p, p)
  lift <- kronecker(diag(p), matrix(1, B, 1L))
  weights <- rep(1 - delta, each = B)
  for (i in seq_along(clusters)) {
    rows <- clusters[[i]]
    Xi <- X[rows, , drop = FALSE]
    residual <- y[rows] - as.numeric(Xi %*% beta)
    psi <- stats::pnorm(residual / h) - (1 - tau)
    phi <- stats::dnorm(residual / h) / h
    basis <- qif_see_basis_v2(length(rows), corstr)
    raw <- matrix(0, p, B)
    dRaw <- matrix(0, B * p, p)
    for (b in seq_len(B)) {
      raw[, b] <- as.numeric(crossprod(Xi, density_weights[rows] * as.numeric(basis[[b]] %*% psi)))
      positions <- seq.int(b, B * p, by = B)
      dRaw[positions, ] <- -crossprod(Xi, density_weights[rows] * (basis[[b]] %*% (phi * Xi)))
    }
    # Eq3按coordinate排列各basis分量，Eq5逐coordinate重复加入delta乘beta。
    G[, i] <- weights * as.numeric(t(raw)) + as.numeric(lift %*% (delta * beta))
    derivative <- derivative + weights * dRaw + sweep(lift, 2L, delta, "*")
  }
  list(G = G, mean = rowMeans(G), covariance = tcrossprod(G) / n,
       derivative = derivative / n, B = B, n = n,
       score_definition = "source_Eq5_modified_extended_quasi_score")
}

qif_see_source_step_v2 <- function(component, beta, context = list()) {
  C <- component$covariance
  dimension <- nrow(C)
  rank <- qr(component$G)$rank
  diagnostics <- c(context, list(matrix_dimension = dimension, rank_upper_bound = component$n,
                                  actual_rank = rank, covariance = C,
                                  score_mean = component$mean, derivative = component$derivative))
  if (rank < dimension) qif_see_stop_v2(
    "qif_covariance_rank", "The source QIF covariance has no ordinary inverse.",
    c(diagnostics, list(status = "rank_deficient")), list(beta = beta, component = component)
  )
  withCallingHandlers({
    weighted_mean <- solve(C, component$mean)
    weighted_derivative <- solve(C, component$derivative)
    gradient <- 2 * as.numeric(crossprod(component$derivative, weighted_mean))
    hessian <- 2 * crossprod(component$derivative, weighted_derivative)
    step <- -as.numeric(solve(hessian, gradient))
    if (any(!is.finite(step))) stop("The source QIF Newton step is non-finite.")
    list(step = step, gradient = gradient, hessian = hessian,
         q_value = component$n * sum(component$mean * weighted_mean), rank = rank)
  }, error = function(condition) {
    if (inherits(condition, "qif_see_source_error")) return(invisible(NULL))
    condition$failure_stage <- "qif_ordinary_inverse"
    condition$diagnostics <- diagnostics
    condition$fit <- list(beta = beta, component = component)
    stop(condition)
  })
}

fit_qif_see_source_v2 <- function(y, X, cluster_id, tau, beta_init, density_weights,
                                   lambda, gamma, corstr = "combined_known_basis",
                                   max_iter = 50L, tol = 1e-5) {
  n <- length(unique(cluster_id))
  h <- n^(-1 / 2)
  beta <- as.numeric(beta_init)
  history <- list()
  converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    # Algorithm1第3步每轮更新delta，第11步在该轮固定delta计算导数。
    denominator <- abs(beta)^(1 + gamma)
    delta <- rep(1, length(beta))
    positive <- denominator > 0
    delta[positive] <- pmin(1, lambda / denominator[positive])
    component <- qif_see_extended_source_v2(y, X, cluster_id, beta, tau, h,
                                             density_weights, delta, corstr)
    calculation <- qif_see_source_step_v2(
      component, beta, list(lambda = lambda, gamma = gamma, iteration = iteration,
                            completed_iterations = history, delta = delta)
    )
    beta_new <- beta + calculation$step
    change <- sqrt(sum(calculation$step^2))
    history[[iteration]] <- list(beta = beta, beta_new = beta_new, delta = delta,
                                  change = change, q_value = calculation$q_value,
                                  covariance_rank = calculation$rank)
    beta <- beta_new
    if (change < tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) qif_see_stop_v2(
    "qif_iteration_limit", "The source QIF algorithm reached its numerical iteration limit.",
    list(lambda = lambda, gamma = gamma, history = history), list(beta = beta)
  )
  denominator <- abs(beta)^(1 + gamma)
  delta <- rep(1, length(beta))
  positive <- denominator > 0
  delta[positive] <- pmin(1, lambda / denominator[positive])
  final <- qif_see_extended_source_v2(y, X, cluster_id, beta, tau, h, density_weights, delta, corstr)
  calculation <- qif_see_source_step_v2(final, beta, list(lambda = lambda, gamma = gamma, iteration = iteration))
  df <- sum(abs(beta) > 1e-6)
  list(beta = beta, lambda = lambda, gamma = gamma, h = h, delta = delta,
       q_value = calculation$q_value, bic = calculation$q_value + df * log(n),
       df = df, converged = TRUE, iterations = iteration, history = history,
       final_component = final, final_step = calculation,
       source_algorithm = "Algorithm1_iterated_delta_Gauss_Newton_ordinary_inverse")
}

fit_benchmark_qif_see_v2 <- function(train, tau, target_coords, tuning, seed,
                                     context = list(), control = list()) {
  set.seed(seed)
  started <- proc.time()[[3L]]
  X <- assert_numeric_matrix_v2(train$X, "X")
  y <- as.numeric(train$y)
  cluster <- as.character(train$cluster_id)
  n <- length(unique(cluster))
  p <- ncol(X)
  corstr <- control$corstr %||% "combined_known_basis"
  B <- length(qif_see_basis_v2(2L, corstr))
  if (length(y) != nrow(X) || length(cluster) != length(y) ||
      any(!is.finite(y)) || is.null(colnames(X))) stop("Invalid QIF data dimensions or coordinates.")
  if (B * p > n) {
    ans <- benchmark_failure_v2("QIF-SEE", "qif_covariance_rank",
                                "The source Bp-dimensional covariance has rank at most the number of clusters.",
                                proc.time()[[3L]] - started)
    ans$fit_object <- list(matrix_dimension = B * p, n_clusters = n, rank_upper_bound = n,
                           actual_rank = NA_integer_, status = "rank_deficient_by_dimension",
                           covariance_definition = "n_inverse_sum_modified_score_outer_product",
                           inverse_definition = "ordinary_inverse", basis = corstr,
                           density_definition = "source_conditional_error_density_at_zero")
    return(benchmark_add_metadata_v2(ans,
      reference_identifier = "Bhattacharya-Bhuiyan-Chatla-2026-10.1111/sjos.70077",
      fidelity_status = "source_Algorithm1_domain_failure"))
  }
  density <- control$density_weights
  if (is.null(density)) qif_see_stop_v2(
    "source_density_unspecified", "Source conditional error-density weights must be supplied explicitly.",
    list(source = "Eq3 Gamma_i", source_estimator_specification = "not_supplied_in_article")
  )
  if (length(density) != length(y) || any(!is.finite(density)) || any(density < 0)) {
    stop("QIF conditional density weights must be nonnegative finite values for every observation.")
  }
  folds <- make_subject_folds_v2(cluster, k = 10L, seed = seed)
  fold_id <- folds$fold[match(cluster, folds$cluster_id)]
  init <- rqPen::rq.pen.cv(x = X, y = y, tau = tau, penalty = "LASSO", nfolds = 10L,
                           foldid = fold_id,
                           printProgress = FALSE)
  coefficients <- as.numeric(stats::coef(init))
  if (length(coefficients) != p + 1L || any(!is.finite(coefficients))) {
    qif_see_stop_v2("initial_fit", "The rqPen initial coefficient vector has an unexpected dimension.", list(initial = init))
  }
  beta_init <- coefficients[-1L]
  grid <- expand.grid(lambda = control$lambda_grid %||% exp(seq(log(.05), log(1), length.out = 8L)),
                       gamma = control$gamma_grid %||% control$gamma %||% 1,
                       stringsAsFactors = FALSE)
  if (any(!is.finite(as.matrix(grid))) || any(grid$lambda <= 0) || any(grid$gamma <= 0)) stop("Invalid QIF tuning grid.")
  candidates <- vector("list", nrow(grid))
  for (index in seq_len(nrow(grid))) {
    candidates[[index]] <- withCallingHandlers(
      fit_qif_see_source_v2(y, X, cluster, tau, beta_init, density,
                            grid$lambda[index], grid$gamma[index], corstr,
                            max_iter = control$max_iter %||% 50L, tol = control$tol %||% 1e-5),
      error = function(condition) {
        condition$diagnostics <- c(condition$diagnostics, list(candidate = index,
          tuning_grid = grid, completed_candidates = candidates[seq_len(index - 1L)], initial = init))
        stop(condition)
      }
    )
  }
  chosen <- which.min(vapply(candidates, `[[`, numeric(1), "bic"))
  best <- candidates[[chosen]]
  beta <- setNames(best$beta, colnames(X))
  ans <- list(method_id = "QIF-SEE", status = "ok", failure_stage = "", failure_message = "",
              beta_hat = beta, beta_tilde = NULL, se = NULL, selected = which(abs(beta) > 1e-6),
              fit_object = list(initial = init, initial_subject_folds = folds,
                                beta_tilde_init = beta_init,
                                tuning_grid = grid, candidates = candidates, selected_candidate = chosen,
                                density_weights = density, basis = corstr, best = best),
              inference_object = list(table = data.frame(), target_scope = "longitudinal_qif_selection"),
              runtime_sec = proc.time()[[3L]] - started, converged = TRUE, kkt_residual = NA_real_,
              warning_messages = character(), implementation_version = "qif_see_source_20261022",
              target_scope = "longitudinal_qif_selection",
              prediction_function = local({ coefficient <- beta; function(newX) as.numeric(as.matrix(newX) %*% coefficient) }))
  benchmark_add_metadata_v2(ans,
    reference_identifier = "Bhattacharya-Bhuiyan-Chatla-2026-10.1111/sjos.70077",
    fidelity_status = "source_Algorithm1_with_explicit_density")
}

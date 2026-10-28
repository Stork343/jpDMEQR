ds_profile <- function(data, beta, tau, h, Lambda, gamma = NULL) {
  if (is.null(gamma)) gamma <- matrix(0, data$n, ncol(data$Z))
  .Call("jpdmeqr_cluster_profiles", as.numeric(data$y - data$X %*% beta), data$Z,
        data$order, data$pointers, gamma, as.numeric(tau), as.numeric(h), Lambda,
        list(grad_tol = min(1e-10, data$n^(-3)), newton_maxit = 150L), PACKAGE = "jpDMEQR")
}

ds_fit_small <- function(data, initial, tau, h, Lambda, weights = rep(1, data$n),
                          ridge = data$n^(-2), gram = crossprod(data$X, data$X * data$weight) / data$n) {
  if (length(weights) != data$n || any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - data$n) > 1e-8) stop("Invalid cluster weights.")
  row_weight <- data$weight * weights[data$group]
  beta <- initial
  evaluate <- function(coefficient, gamma = NULL) {
    profile <- ds_profile(data, coefficient, tau, h, Lambda, gamma)
    gradient <- -drop(crossprod(data$X, profile$psi * row_weight)) / data$n + ridge * drop(gram %*% coefficient)
    objective <- sum(profile$objective * weights) / data$n + ridge * sum(coefficient * (gram %*% coefficient)) / 2
    list(profile = profile, gradient = gradient, objective = objective)
  }
  current <- evaluate(beta)
  for (iteration in seq_len(300L)) {
    H <- sb_geometry_cpp(data$X, data$Z, current$profile$phi, data$order, data$pointers,
                           Lambda, FALSE, weights) + ridge * gram
    factor <- chol(H)
    direction <- -drop(backsolve(factor, forwardsolve(t(factor), current$gradient)))
    decrement <- -sum(current$gradient * direction)
    if (max(abs(current$gradient)) < min(1e-8, data$n^(-2)) &&
        max(abs(direction)) < min(1e-7, data$n^(-1.5))) break
    if (!is.finite(decrement) || decrement <= 0) stop("Invalid Newton decrement.")
    accepted <- FALSE
    for (backtrack in 0:50) {
      step <- 2^(-backtrack)
      candidate <- evaluate(beta + step * direction, current$profile$gamma)
      allowance <- 64 * .Machine$double.eps * max(1, abs(current$objective))
      if (candidate$objective <= current$objective - 1e-4 * step * decrement + allowance) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) stop("Weighted profile linesearch did not converge.")
    beta <- beta + step * direction
    current <- candidate
    if (iteration == 300L) stop("Weighted profile Newton iteration limit.")
  }
  target <- c(1, rep(0, ncol(data$X) - 1L))
  inverse <- drop(backsolve(factor, forwardsolve(t(factor), target)))
  scores <- -as.numeric(rowsum(drop(data$X %*% inverse) * current$profile$psi * data$weight,
                                data$group, reorder = FALSE))
  list(beta = beta, projected_scores = scores, gamma = current$profile$gamma,
       iterations = iteration, kkt = max(abs(current$gradient)), objective = current$objective,
       hessian = H, inverse = inverse)
}

ds_pair_small <- function(data, initial, tau, h, Lambda, weights = rep(1, data$n),
                           ridge = data$n^(-2), gram = crossprod(data$X, data$X * data$weight) / data$n) {
  full <- ds_fit_small(data, initial, tau, h, Lambda, weights, ridge, gram)
  half <- ds_fit_small(data, full$beta, tau, h / 2, Lambda, weights, ridge, gram)
  influence <- -2 * half$projected_scores + full$projected_scores
  list(estimate = 2 * half$beta[1L] - full$beta[1L], full = full, half = half,
       influence = influence, max_kkt = max(full$kkt, half$kkt))
}

ds_prepare <- function(y, X, Z, cluster_id, tau, coordinates, Lambda, penalty_factor, h,
                       include_design = TRUE, budget_multiplier = 1) {
  data <- sb_data(y, X, Z, cluster_id)
  assignment <- pc_folds(data$rows, 2L)
  folds <- vector("list", 2L)
  for (f in 1:2) {
    training <- sb_subset(data, unlist(data$rows[assignment != f], use.names = FALSE))
    validation <- sb_subset(data, unlist(data$rows[assignment == f], use.names = FALSE))
    initial <- sb_initial(training, tau, Lambda, penalty_factor, budget_multiplier)
    graphs <- sb_graph(initial, coordinates, tau, h, Lambda, budget_multiplier)
    if (include_design) {
      budget <- min(vapply(initial$training, function(data) sb_budget(data$n, ncol(data$X), multiplier = budget_multiplier), integer(1)))
      factor <- training$X * sqrt(training$weight)
      for (j in seq_along(coordinates)) {
        k <- coordinates[j]
        direction <- pc_nodewise(factor, training$rows, k)
        other <- setdiff(seq_len(ncol(X)), k)
        selected <- which(abs(direction$coefficient) > 1e-10)
        selected <- head(selected[order(abs(direction$coefficient[selected]), decreasing = TRUE)], budget)
        design_support <- other[selected]
        graphs[[j]]$curvature_support <- graphs[[j]]$support
        graphs[[j]]$design_support <- design_support
        graphs[[j]]$design_tuning <- direction$tuning
        graphs[[j]]$design_kkt <- direction$kkt
        graphs[[j]]$support <- c(k, sort(setdiff(union(graphs[[j]]$support, design_support), k)))
      }
    }
    beta <- Reduce(`+`, initial$coefficients) / 3
    models <- lapply(graphs, function(graph) {
      small <- validation
      small$X <- validation$X[, graph$support, drop = FALSE]
      list(data = small, initial = beta[graph$support], support = graph$support,
           gram = crossprod(small$X, small$X * small$weight) / small$n)
    })
    initial$training <- initial$validation <- initial$paths <- NULL
    folds[[f]] <- list(training_subjects = names(training$rows),
                       validation_subjects = names(validation$rows), initial = initial,
                       graphs = graphs, models = models)
  }
  list(folds = folds, coordinates = coordinates, n = data$n)
}

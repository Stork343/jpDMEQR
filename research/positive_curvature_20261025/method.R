library(glmnet)

pc_compile <- function(cache_dir = "results/scratch/positive_curvature_cache") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp("research/positive_curvature_20261025/joint_pdal.cpp",
                  cacheDir = cache_dir, rebuild = FALSE, showOutput = FALSE,
                  env = .GlobalEnv)
}

pc_cluster_scores <- function(X, dual, rows) {
  t(vapply(rows, function(index) -drop(crossprod(X[index, , drop = FALSE], dual[index])),
           numeric(ncol(X))))
}

pc_loadings <- function(scores) {
  centred <- sweep(scores, 2, colMeans(scores))
  loading <- sqrt(colMeans(centred^2))
  if (any(!is.finite(loading))) stop("Nonfinite cluster score loadings.")
  loading
}

pc_ratio_path <- function(n, length_out) {
  exp(seq(log(1), log(min(0.05, 1 / sqrt(n))), length.out = length_out))
}

pc_exact_profile <- function(y, Z, rows, Lambda, tau) {
  exact_cluster_profile_cpp(as.numeric(y), Z,
                            as.integer(unlist(rows, use.names = FALSE) - 1L),
                            as.integer(c(0, cumsum(lengths(rows)))), Lambda, tau)
}

pc_folds <- function(rows, count = 5L) {
  if (length(rows) < 2L * count) stop("Insufficient clusters for the declared folds.")
  ordered <- sort(names(rows))
  assignment <- setNames(rep(seq_len(count), length.out = length(ordered)), ordered)
  as.integer(assignment[names(rows)])
}

pc_prepare_initial <- function(y, X, Z, cluster_id, tau, Lambda, penalty_factor,
                                max_iterations = 30000L, fit_always_included = TRUE) {
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  y <- as.numeric(y)
  if (length(y) != nrow(X) || nrow(Z) != nrow(X) || length(cluster_id) != length(y) ||
      any(!is.finite(y)) || tau <= 0 || tau >= 1 || length(penalty_factor) != ncol(X) ||
      any(!is.finite(penalty_factor)) || any(penalty_factor < 0)) stop("Invalid initial data.")
  clusters <- cluster_index_v2(cluster_id)
  n <- clusters$n
  p <- ncol(X)
  q <- ncol(Z)
  Lambda <- make_lambda_matrix_v2(Lambda, q)
  row_cluster <- integer(length(y))
  inverse_size <- numeric(length(y))
  for (i in seq_along(clusters$rows)) {
    index <- clusters$rows[[i]]
    row_cluster[index] <- i - 1L
    inverse_size[index] <- 1 / length(index)
  }
  scale <- sqrt(colSums(X^2 * inverse_size) / n)
  if (any(!is.finite(scale)) || any(scale <= 0)) stop("Every fixed effect column needs positive RMS.")
  Xs <- sweep(X, 2, scale, "/")
  null <- pc_exact_profile(y, Z, clusters$rows, Lambda, tau)
  theta0 <- numeric(p)
  U <- which(penalty_factor == 0)
  fitted_always <- length(U) > 0L && isTRUE(fit_always_included)
  if (fitted_always) {
    null <- joint_profile_pdal_cpp(
      Xs[, U, drop = FALSE], Z, y, row_cluster, inverse_size, Lambda,
      numeric(length(U)), tau, numeric(length(U)), null$gamma, null$dual,
      max_iterations, 1e-8
    )
    theta0[U] <- null$theta
  }
  scores0 <- pc_cluster_scores(Xs, null$dual, clusters$rows)
  loading <- pc_loadings(scores0)
  P <- which(penalty_factor > 0)
  if (length(P) && any(loading[P] <= 1e-10)) stop("Initial loadings are degenerate.")
  lambda_max <- if (length(P)) max(abs(colMeans(scores0)[P]) / (loading[P] * penalty_factor[P])) else 0
  list(y = y, Xs = Xs, Z = Z, rows = clusters$rows, n = n, p = p, q = q,
       scale = scale, Lambda = Lambda, row_cluster = row_cluster,
       inverse_size = inverse_size, theta0 = theta0, gamma0 = null$gamma,
       dual0 = null$dual, loadings = loading, penalty_factor = penalty_factor,
       null_gap = if (fitted_always) null$duality_gap else mean(null$gap),
       null_iterations = if (fitted_always) null$iterations else sum(null$iterations),
       lambda_max = lambda_max)
}

pc_solve_initial <- function(prepared, ratio, tau, gap_tolerance = 1e-8,
                              max_iterations = 30000L, warm = NULL) {
  if (is.null(warm)) warm <- list(theta = prepared$theta0, gamma = prepared$gamma0, dual = prepared$dual0)
  base <- ratio * prepared$lambda_max
  fitted <- joint_profile_pdal_cpp(
    prepared$Xs, prepared$Z, prepared$y, prepared$row_cluster,
    prepared$inverse_size, prepared$Lambda,
    base * prepared$loadings * prepared$penalty_factor,
    tau, warm$theta, warm$gamma, warm$dual,
    max_iterations, gap_tolerance
  )
  names(fitted$theta) <- colnames(prepared$Xs)
  rownames(fitted$gamma) <- names(prepared$rows)
  fitted$beta <- fitted$theta / prepared$scale
  fitted$scale <- prepared$scale
  fitted$Xs <- prepared$Xs
  fitted$rows <- prepared$rows
  fitted$Lambda <- prepared$Lambda
  fitted$n_clusters <- prepared$n
  fitted$lambda_base <- base
  fitted$lambda_ratio <- ratio
  fitted$lambda_max <- prepared$lambda_max
  fitted$coordinate_penalty <- base * prepared$loadings * prepared$penalty_factor
  fitted$loadings <- prepared$loadings
  fitted$initial_ridge <- 1 / prepared$n
  fitted
}

pc_fit_initial <- function(y, X, Z, cluster_id, tau, Lambda = diag(ncol(Z)),
                           penalty_factor = rep(1, ncol(X)),
                           gap_tolerance = 1e-8, max_iterations = 30000L) {
  started <- proc.time()[[3]]
  prepared <- pc_prepare_initial(y, X, Z, cluster_id, tau, Lambda, penalty_factor)
  ratios <- if (any(penalty_factor > 0)) pc_ratio_path(prepared$n, 8L) else 1
  folds <- pc_folds(prepared$rows)
  loss <- matrix(0, 5L, length(ratios))
  fold_coefficients <- vector("list", 5L)
  fold_diagnostics <- vector("list", 5L)
    for (f in seq_len(5L)) {
      training <- unlist(prepared$rows[folds != f], use.names = FALSE)
      validation <- unlist(prepared$rows[folds == f], use.names = FALSE)
      train <- pc_prepare_initial(y[training], X[training, , drop = FALSE],
                                   Z[training, , drop = FALSE], cluster_id[training],
                                   tau, Lambda, penalty_factor)
      validation_rows <- cluster_index_v2(cluster_id[validation])$rows
      warm <- NULL
      coefficient_path <- matrix(0, ncol(X), length(ratios))
      path_diagnostics <- vector("list", length(ratios))
      for (j in seq_along(ratios)) {
        warm <- pc_solve_initial(train, ratios[j], tau, gap_tolerance = 1e-7,
                                  max_iterations = max_iterations, warm = warm)
        evaluated <- pc_exact_profile(
          y[validation] - drop(X[validation, , drop = FALSE] %*% warm$beta),
          Z[validation, , drop = FALSE], validation_rows, Lambda, tau
        )
        loss[f, j] <- sum(evaluated$loss)
        coefficient_path[, j] <- warm$beta
        path_diagnostics[[j]] <- c(lambda_base = warm$lambda_base, gap = warm$duality_gap,
                                   iterations = warm$iterations)
      }
      fold_coefficients[[f]] <- coefficient_path
      fold_diagnostics[[f]] <- path_diagnostics
    }
  selected <- which.min(colSums(loss))
  fitted <- pc_solve_initial(prepared, ratios[selected], tau, gap_tolerance, max_iterations)
  fitted$tuning <- list(rule = "cluster_profile_risk_cv", ratios = ratios,
                        fold_loss_sum = loss, fold_subjects = lapply(seq_len(5L), function(f) names(prepared$rows)[folds == f]),
                        selected = selected, validation_loss = colSums(loss) / prepared$n)
  fitted$fold_coefficients <- lapply(fold_coefficients, function(path) path[, selected])
  fitted$fold_coefficient_paths <- fold_coefficients
  fitted$fold_diagnostics <- fold_diagnostics
  fitted$fold_assignment <- folds
  fitted$runtime_sec <- proc.time()[[3]] - started
  fitted
}

pc_profile_geometry <- function(y, Xs, Z, rows, theta, gamma, Lambda, tau, h) {
  residual <- as.numeric(y - Xs %*% theta)
  full <- compiled_cluster_profiles_v2(
    residual, Z, rows, tau, h, Lambda, gamma_init = gamma,
    control = list(grad_tol = 1e-9, newton_maxit = 150L)
  )
  half <- compiled_cluster_profiles_v2(
    residual, Z, rows, tau, h / 2, Lambda, gamma_init = full$gamma,
    control = list(grad_tol = 1e-9, newton_maxit = 150L)
  )
  n <- length(rows)
  p <- ncol(Xs)
  factor <- matrix(0, nrow(Xs), p)
  scores <- matrix(0, n, p, dimnames = list(names(rows), colnames(Xs)))
  for (i in seq_along(rows)) {
    index <- rows[[i]]
    Xi <- Xs[index, , drop = FALSE]
    Zi <- Z[index, , drop = FALSE]
    w <- full$phi[index] / length(index)
    WZ <- Zi * w
    B <- crossprod(Zi, WZ) + Lambda
    L <- diag(w, nrow = length(index)) - WZ %*% solve(B, t(WZ))
    eig <- eigen((L + t(L)) / 2, symmetric = TRUE)
    if (min(eig$values) < -1e-10 * max(1, max(eig$values))) stop("The profile factor is not positive semidefinite.")
    root <- eig$vectors %*% (sqrt(pmax(eig$values, 0)) * t(eig$vectors))
    factor[index, ] <- root %*% Xi
    score <- 2 * half$psi[index] - full$psi[index]
    scores[i, ] <- -drop(crossprod(Xi, score)) / length(index)
  }
  colnames(factor) <- colnames(Xs)
  list(factor = factor, cluster_scores = scores, full = full, half = half,
       h = h, rows = rows, n_clusters = n)
}

pc_nodewise <- function(factor, rows, coordinate) {
  n <- length(rows)
  p <- ncol(factor)
  N <- nrow(factor)
  diagonal <- sqrt(colSums(factor^2) / n)
  if (any(!is.finite(diagonal)) || any(diagonal <= 1e-10)) stop("A profile curvature coordinate is unidentified.")
  standardised <- sweep(factor, 2, diagonal, "/")
  other <- setdiff(seq_len(p), coordinate)
  response <- standardised[, coordinate]
  design <- standardised[, other, drop = FALSE]
  solve_path <- function(indices, clusters, ratios) {
    penalty_max <- max(abs(drop(crossprod(design[indices, , drop = FALSE], response[indices]))) / clusters)
    penalties <- penalty_max * ratios
    if (length(other) == 1L) {
      coefficient <- vapply(penalties, function(penalty) {
        soft_threshold_weighted_v2(sum(design[indices, 1] * response[indices]) / clusters, penalty) /
          (sum(design[indices, 1]^2) / clusters)
      }, numeric(1))
      return(list(coefficient = matrix(coefficient, 1), penalties = penalties))
    }
    fit <- glmnet::glmnet(
      x = design[indices, , drop = FALSE] * sqrt(length(indices) / clusters),
      y = response[indices] * sqrt(length(indices) / clusters),
      family = "gaussian", alpha = 1, lambda = penalties,
      intercept = FALSE, standardize = FALSE, thresh = 1e-13, maxit = 1000000L
    )
    if (fit$jerr != 0L) stop("The nodewise solver did not converge.")
    list(coefficient = as.matrix(fit$beta), penalties = fit$lambda)
  }
  ratios <- pc_ratio_path(n, 24L)
  folds <- pc_folds(rows)
  loss <- matrix(0, 5L, length(ratios))
  if (length(other)) {
    for (f in seq_len(5L)) {
      training <- unlist(rows[folds != f], use.names = FALSE)
      validation <- unlist(rows[folds == f], use.names = FALSE)
      path <- solve_path(training, sum(folds != f), ratios)
      residual <- response[validation] - design[validation, , drop = FALSE] %*% path$coefficient
      loss[f, ] <- colSums(residual^2)
    }
    selected <- which.min(colSums(loss))
    path <- solve_path(seq_len(N), n, ratios[selected])
    coefficient <- as.numeric(path$coefficient[, 1])
    penalty <- rep(path$penalties[1], length(other))
  } else {
    coefficient <- penalty <- numeric()
    selected <- 1L
  }
  residual <- response - drop(design %*% coefficient)
  gradient <- -drop(crossprod(design, residual)) / n
  defect <- if (length(other)) {
    ifelse(abs(coefficient) > 1e-10,
           abs(gradient + penalty * sign(coefficient)), pmax(abs(gradient) - penalty, 0))
  } else 0
  if (max(defect) > 1e-6) stop("The nodewise KKT check failed.")
  residual_tau2 <- sum(residual^2) / n + sum(penalty * abs(coefficient))
  tau2 <- sum(response * residual) / n
  kkt_inner_product <- sum(coefficient * (gradient + penalty * sign(coefficient)))
  identity_error <- abs(residual_tau2 - tau2 - kkt_inner_product)
  if (identity_error > 1e-8) stop("The nodewise curvature identity check failed.")
  if (!is.finite(tau2) || tau2 <= 1e-10) stop("The nodewise directional curvature is nonpositive.")
  contrast <- numeric(p)
  contrast[coordinate] <- 1
  contrast[other] <- -coefficient
  omega <- contrast / diagonal / (diagonal[coordinate] * tau2)
  product <- drop(crossprod(factor, factor %*% omega)) / n
  normalisation_error <- abs(product[coordinate] - 1)
  if (normalisation_error > 1e-6) stop("The nodewise direction normalisation check failed.")
  list(omega = omega, coefficient = coefficient, penalty = penalty,
       kkt = max(defect), normalisation_error = normalisation_error,
       inverse_defect = max(abs(product - replace(numeric(p), coordinate, 1))),
       tau2 = tau2, residual_tau2 = residual_tau2, identity_error = identity_error, diagonal = diagonal,
       tuning = list(rule = "cluster_curvature_cv", ratios = ratios,
                     fold_loss_sum = loss, selected = selected,
                     selected_ratio = ratios[selected]))
}

pc_fit_infer <- function(y, X, Z, cluster_id, tau, coordinates,
                         Lambda = diag(ncol(Z)), penalty_factor = rep(1, ncol(X)),
                         h = length(unique(cluster_id))^(-0.3),
                         gap_tolerance = 1e-8, keep_factor = FALSE, initial_fit = NULL) {
  started <- proc.time()[[3]]
  if (is.character(coordinates)) coordinates <- match(coordinates, colnames(X))
  if (!length(coordinates) || anyNA(coordinates) || anyDuplicated(coordinates) ||
      any(coordinates < 1 | coordinates > ncol(X))) stop("Invalid requested coordinates.")
  initial <- if (is.null(initial_fit)) {
    pc_fit_initial(y, X, Z, cluster_id, tau, Lambda, penalty_factor,
                    gap_tolerance = gap_tolerance)
  } else initial_fit
  geometry_started <- proc.time()[[3]]
  factor <- matrix(0, nrow(X), ncol(X), dimnames = dimnames(X))
  scores <- matrix(0, initial$n_clusters, ncol(X),
                    dimnames = list(names(initial$rows), colnames(X)))
  theta_center <- numeric(ncol(X))
  for (f in seq_len(5L)) {
    subjects <- which(initial$fold_assignment == f)
    index <- unlist(initial$rows[subjects], use.names = FALSE)
    local_rows <- cluster_index_v2(cluster_id[index])$rows
    theta <- initial$scale * initial$fold_coefficients[[f]]
    evaluated <- pc_profile_geometry(
      y[index], initial$Xs[index, , drop = FALSE], Z[index, , drop = FALSE], local_rows,
      theta, matrix(0, length(subjects), ncol(Z)), initial$Lambda, tau, h
    )
    factor[index, ] <- evaluated$factor
    scores[match(rownames(evaluated$cluster_scores), rownames(scores)), ] <- evaluated$cluster_scores
    theta_center <- theta_center + length(subjects) / initial$n_clusters * theta
  }
  geometry <- list(factor = factor, cluster_scores = scores)
  geometry_sec <- proc.time()[[3]] - geometry_started
  inference_started <- proc.time()[[3]]
  directions <- lapply(coordinates, function(k) pc_nodewise(geometry$factor, initial$rows, k))
  omega <- do.call(cbind, lapply(directions, `[[`, "omega"))
  projected <- geometry$cluster_scores %*% omega
  centred <- sweep(projected, 2, colMeans(projected))
  n <- initial$n_clusters
  estimate <- (theta_center[coordinates] - colMeans(projected)) / initial$scale[coordinates]
  se <- sqrt(colMeans(centred^2) / n) / initial$scale[coordinates]
  if (any(!is.finite(estimate)) || any(!is.finite(se)) || any(se <= 0)) stop("Invalid coordinate inference.")
  critical <- stats::qnorm(0.975)
  table <- data.frame(coordinate = colnames(X)[coordinates], beta_tilde = estimate,
                      se = se, lower = estimate - critical * se, upper = estimate + critical * se)
  inference_sec <- proc.time()[[3]] - inference_started
  initial$Xs <- NULL
  list(table = table, initial = initial, directions = directions,
       cluster_scores = geometry$cluster_scores,
       factor = if (keep_factor) geometry$factor else NULL,
       projected_scores = projected, coordinates = coordinates,
       beta_center = theta_center / initial$scale,
       n_clusters = n, h = h, method_version = "crossfit_positive_curvature_20261025",
       score_definition = "profile_bandwidth_pair_2_half_minus_full",
       curvature_definition = "single_bandwidth_positive_profile",
       timing = c(initial = initial$runtime_sec, geometry = geometry_sec,
                  inference = inference_sec, total = proc.time()[[3]] - started))
}

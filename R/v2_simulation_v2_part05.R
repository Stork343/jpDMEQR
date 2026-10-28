profile_pair_population_eigen_init_v2 <- function(root) {
  if (!is.loaded("jpDMEQR_population_crossprod", PACKAGE = "jpDMEQR")) {
    stop("The installed jpDMEQR native population crossproduct routine is required.")
  }
  options(jpDMEQR.population_matrix_backend = "Eigen_single_thread")
  invisible(TRUE)
}

profile_pair_population_crossprod_v2 <- function(X, Y = NULL) {
  backend <- getOption("jpDMEQR.population_matrix_backend", "R_BLAS")
  if (identical(backend, "R_BLAS")) return(crossprod(X, Y))
  if (!identical(backend, "Eigen_single_thread") || !is.loaded("jpDMEQR_population_crossprod", PACKAGE = "jpDMEQR")) {
    stop("The requested population matrix backend is unavailable.")
  }
  .Call("jpDMEQR_population_crossprod", X, Y, PACKAGE = "jpDMEQR")
}

profile_pair_population_layout_v2 <- function(dat) {
  ci <- cluster_index_v2(dat$cluster_id)
  sizes <- lengths(ci$rows)
  groups <- lapply(sort(unique(sizes)), function(m) {
    clusters <- which(sizes == m)
    rows <- matrix(unlist(ci$rows[clusters], use.names = FALSE),
                   nrow = length(clusters), ncol = m, byrow = TRUE)
    list(clusters = clusters, rows = rows, m = m)
  })
  list(rows = ci$rows, levels = ci$levels, sizes = sizes, groups = groups,
       n = ci$n, row_cluster = as.integer(factor(dat$cluster_id, levels = ci$levels)))
}

profile_pair_population_nuisance_v2 <- function(residual_base, Z, layout, tau, h,
                                                Lambda, tolerance = 1e-9,
                                                max_iter = 100L) {
  q <- ncol(Z)
  if (!q %in% c(1L, 2L)) stop("Streamed population profiles require q=1 or q=2.")
  backend <- getOption("jpDMEQR.population_nuisance_backend", "vectorised_R")
  if (identical(backend, "compiled")) {
    result <- compiled_cluster_profiles_v2(
      residual_base, Z, layout$rows, tau, h, Lambda,
      control = list(grad_tol = tolerance, newton_maxit = max_iter))
    weights <- result$phi / layout$sizes[layout$row_cluster]
    sums <- function(value) as.numeric(rowsum(value, layout$row_cluster, reorder = FALSE))
    return(list(
      gamma = result$gamma, psi = result$psi, weights = weights,
      objective = mean(result$objective), max_nuisance_gradient = max(result$gradient_max),
      nuisance_converged = TRUE, iterations = result$iterations,
      B11 = Lambda[1, 1] + sums(weights * Z[, 1]^2),
      B12 = if (q == 2L) Lambda[1, 2] + sums(weights * Z[, 1] * Z[, 2]) else numeric(layout$n),
      B22 = if (q == 2L) Lambda[2, 2] + sums(weights * Z[, 2]^2) else numeric(layout$n),
      nuisance_backend = "compiled_cluster_profiles"))
  }
  if (!identical(backend, "vectorised_R")) stop("Unknown population nuisance backend.")
  gamma <- matrix(0, layout$n, q)
  psi <- weights <- numeric(length(residual_base))
  objectives <- gradients <- numeric(layout$n)
  B11 <- B12 <- B22 <- numeric(layout$n)
  iterations <- integer(layout$n)
  for (group in layout$groups) {
    ids <- group$clusters
    b <- length(ids)
    m <- group$m
    rows <- group$rows
    base <- matrix(residual_base[rows], b, m)
    z1 <- matrix(Z[as.vector(rows), 1], b, m)
    z2 <- if (q == 2L) matrix(Z[as.vector(rows), 2], b, m) else matrix(0, b, m)
    value_at <- function(g) {
      r <- base - z1 * g[, 1]
      if (q == 2L) r <- r - z2 * g[, 2]
      objective <- rowMeans(matrix(rho_smooth_v2(r, tau, h), b, m)) +
        0.5 * rowSums(g * (g %*% Lambda))
      list(r = r, objective = objective)
    }
    g <- matrix(0, b, q)
    state <- value_at(g)
    converged <- rep(FALSE, b)
    for (iteration in seq_len(max_iter)) {
      score <- matrix(psi_smooth_v2(state$r, tau, h), b, m)
      curvature <- matrix(phi_smooth_v2(state$r, tau, h), b, m) / m
      gradient <- g %*% Lambda
      gradient[, 1] <- gradient[, 1] - rowMeans(score * z1)
      if (q == 2L) gradient[, 2] <- gradient[, 2] - rowMeans(score * z2)
      gradient_norm <- apply(abs(gradient), 1, max)
      active <- gradient_norm > tolerance
      converged <- !active
      if (!any(active)) break
      a11 <- Lambda[1, 1] + rowSums(curvature * z1^2)
      direction <- matrix(0, b, q)
      if (q == 1L) {
        direction[, 1] <- -gradient[, 1] / a11
      } else {
        a12 <- Lambda[1, 2] + rowSums(curvature * z1 * z2)
        a22 <- Lambda[2, 2] + rowSums(curvature * z2^2)
        determinant <- a11 * a22 - a12^2
        if (any(!is.finite(determinant)) || any(determinant <= 0)) {
          stop("Population nuisance Hessian is not positive definite.")
        }
        direction[, 1] <- (-a22 * gradient[, 1] + a12 * gradient[, 2]) / determinant
        direction[, 2] <- (a12 * gradient[, 1] - a11 * gradient[, 2]) / determinant
      }
      direction[!active, ] <- 0
      slope <- rowSums(direction * gradient)
      steps <- rep(1, b)
      accepted <- !active
      for (backtrack in 0:50L) {
        candidate <- g + direction * steps
        trial <- value_at(candidate)
        rounding <- 16 * .Machine$double.eps * pmax(1, abs(state$objective))
        sufficient <- is.finite(trial$objective) &
          trial$objective <= state$objective + 1e-4 * steps * slope + rounding
        accepted <- accepted | sufficient
        if (all(accepted)) break
        steps[!accepted] <- steps[!accepted] / 2
      }
      if (!all(accepted)) {
        stop("Population nuisance line search failed at cluster ",
             layout$levels[ids[which(!accepted)[1]]], ".")
      }
      g <- candidate
      state <- trial
      iterations[ids[active]] <- iteration
    }
    score <- matrix(psi_smooth_v2(state$r, tau, h), b, m)
    curvature <- matrix(phi_smooth_v2(state$r, tau, h), b, m) / m
    gradient <- g %*% Lambda
    gradient[, 1] <- gradient[, 1] - rowMeans(score * z1)
    if (q == 2L) gradient[, 2] <- gradient[, 2] - rowMeans(score * z2)
    gradient_norm <- apply(abs(gradient), 1, max)
    if (any(!is.finite(gradient_norm)) || any(gradient_norm > tolerance)) {
      stop("Population nuisance gradient exceeds tolerance: ", max(gradient_norm), ".")
    }
    gamma[ids, ] <- g
    psi[as.vector(rows)] <- as.vector(score)
    weights[as.vector(rows)] <- as.vector(curvature)
    objectives[ids] <- state$objective
    gradients[ids] <- gradient_norm
    B11[ids] <- Lambda[1, 1] + rowSums(curvature * z1^2)
    if (q == 2L) {
      B12[ids] <- Lambda[1, 2] + rowSums(curvature * z1 * z2)
      B22[ids] <- Lambda[2, 2] + rowSums(curvature * z2^2)
    }
  }
  list(gamma = gamma, psi = psi, weights = weights, objective = mean(objectives),
       max_nuisance_gradient = max(gradients), nuisance_converged = TRUE,
       iterations = iterations, B11 = B11, B12 = B12, B22 = B22)
}

profile_pair_population_batch_hessian_v2 <- function(X, Z, ids, row_cluster,
                                                     nuisance, Lambda,
                                                     identity_check = TRUE) {
  w <- nuisance$weights
  q <- ncol(Z)
  U1 <- rowsum(X * (w * Z[, 1]), row_cluster, reorder = FALSE)
  b11 <- nuisance$B11[ids]
  t1 <- U1 / sqrt(b11)
  correction <- profile_pair_population_crossprod_v2(t1)
  if (q == 1L) {
    A1 <- U1 / b11
  } else {
    U2 <- rowsum(X * (w * Z[, 2]), row_cluster, reorder = FALSE)
    b12 <- nuisance$B12[ids]
    b22 <- nuisance$B22[ids]
    l21 <- b12 / sqrt(b11)
    l22 <- sqrt(b22 - l21^2)
    t2 <- (U2 - t1 * l21) / l22
    correction <- correction + profile_pair_population_crossprod_v2(t2)
    determinant <- b11 * b22 - b12^2
    A1 <- (U1 * b22 - U2 * b12) / determinant
    A2 <- (U2 * b11 - U1 * b12) / determinant
  }
  H <- profile_pair_population_crossprod_v2(X * sqrt(w)) - correction
  identity <- NULL
  if (identity_check) {
    Xt <- X - A1[row_cluster, , drop = FALSE] * Z[, 1]
    if (q == 2L) Xt <- Xt - A2[row_cluster, , drop = FALSE] * Z[, 2]
    R <- chol(Lambda)
    AL1 <- R[1, 1] * A1
    if (q == 2L) AL1 <- AL1 + R[1, 2] * A2
    identity <- profile_pair_population_crossprod_v2(Xt * sqrt(w)) + profile_pair_population_crossprod_v2(AL1)
    if (q == 2L) identity <- identity + profile_pair_population_crossprod_v2(R[2, 2] * A2)
  }
  list(hessian_sum = H, identity_sum = identity)
}

profile_pair_population_streamed_v2 <- function(dat, beta, tau, h, lambda_gamma = 1,
                                               batch_clusters = 1000L,
                                               need_hessian = TRUE,
                                               need_moments = TRUE,
                                               nuisance_tol = 1e-9) {
  validate_tau_h_v2(tau, h)
  p <- ncol(dat$X)
  if (length(beta) != p || any(!is.finite(beta))) stop("Population beta has wrong dimension.")
  layout <- dat$population_layout %||% profile_pair_population_layout_v2(dat)
  n <- layout$n
  Lambda <- make_lambda_matrix_v2(lambda_gamma, ncol(dat$Z))
  residual_base <- dat$y - as.numeric(dat$X %*% beta)
  half <- profile_pair_population_nuisance_v2(residual_base, dat$Z, layout, tau, h / 2,
                                              Lambda, nuisance_tol)
  full <- profile_pair_population_nuisance_v2(residual_base, dat$Z, layout, tau, h,
                                              Lambda, nuisance_tol)
  inv_m <- 1 / layout$sizes
  score_half <- -as.numeric(crossprod(dat$X, half$psi * inv_m[layout$row_cluster])) / n
  score_full <- -as.numeric(crossprod(dat$X, full$psi * inv_m[layout$row_cluster])) / n
  names(score_half) <- names(score_full) <- colnames(dat$X)
  output <- list(
    objective = 2 * half$objective - full$objective,
    score = 2 * score_half - score_full, n_clusters = n,
    h = h, bandwidths = c(h / 2, h), bandwidth_weights = c(2, -1),
    score_definition = "profile_bandwidth_pair_2_half_minus_full",
    nuisance_converged = TRUE,
    max_nuisance_gradient = max(half$max_nuisance_gradient, full$max_nuisance_gradient),
    profile_identity_error = NA_real_,
    population_engine = "streamed_profile_crossproducts_v1",
    population_nuisance_backend = getOption("jpDMEQR.population_nuisance_backend", "vectorised_R"),
    population_matrix_backend = getOption("jpDMEQR.population_matrix_backend", "R_BLAS"),
    batch_clusters = as.integer(batch_clusters)
  )
  if (!need_hessian && !need_moments) return(output)
  empty <- matrix(0, p, p, dimnames = list(colnames(dat$X), colnames(dat$X)))
  Hh <- Hf <- Ih <- If <- empty
  Ch <- Cf <- Ccross <- Cpair <- empty
  batches <- split(seq_len(n), ceiling(seq_len(n) / batch_clusters))
  for (ids in batches) {
    rows <- unlist(layout$rows[ids], use.names = FALSE)
    local_cluster <- rep(seq_along(ids), layout$sizes[ids])
    X <- dat$X[rows, , drop = FALSE]
    Z <- dat$Z[rows, , drop = FALSE]
    if (need_hessian) {
      half_batch <- half
      full_batch <- full
      half_batch$weights <- half$weights[rows]
      full_batch$weights <- full$weights[rows]
      bh <- profile_pair_population_batch_hessian_v2(X, Z, ids, local_cluster,
                                                     half_batch, Lambda)
      bf <- profile_pair_population_batch_hessian_v2(X, Z, ids, local_cluster,
                                                     full_batch, Lambda)
      Hh <- Hh + bh$hessian_sum
      Hf <- Hf + bf$hessian_sum
      Ih <- Ih + bh$identity_sum
      If <- If + bf$identity_sum
    }
    if (need_moments) {
      Gh <- -rowsum(X * (half$psi[rows] * inv_m[ids][local_cluster]),
                     local_cluster, reorder = FALSE)
      Gf <- -rowsum(X * (full$psi[rows] * inv_m[ids][local_cluster]),
                     local_cluster, reorder = FALSE)
      Gp <- 2 * Gh - Gf
      Ch <- Ch + profile_pair_population_crossprod_v2(Gh)
      Cf <- Cf + profile_pair_population_crossprod_v2(Gf)
      Ccross <- Ccross + profile_pair_population_crossprod_v2(Gh, Gf)
      Cpair <- Cpair + profile_pair_population_crossprod_v2(Gp)
    }
  }
  if (need_hessian) {
    Hh <- (Hh + t(Hh)) / (2 * n)
    Hf <- (Hf + t(Hf)) / (2 * n)
    output$hessian <- 2 * Hh - Hf
    output$hessian_identity <- (2 * Ih - If) / n
    output$profile_identity_error <- max(abs(output$hessian - output$hessian_identity))
    output$components_half_h <- list(hessian = Hh)
    output$components_h <- list(hessian = Hf)
  }
  if (need_moments) {
    Sh <- Ch / n - tcrossprod(score_half)
    Sf <- Cf / n - tcrossprod(score_full)
    cross <- Ccross / n - tcrossprod(score_half, score_full)
    Sigma <- Cpair / n - tcrossprod(output$score)
    identity <- max(abs(Sigma - (4 * Sh + Sf - 2 * cross - 2 * t(cross))))
    if (!is.finite(identity) || identity > 1e-8) stop("Streamed covariance identity failed.")
    output$population_moments <- list(
      n_clusters = n, score_mean = output$score,
      score_mean_half = score_half, score_mean_full = score_full,
      Sigma_population = Sigma, Sigma_half = Sh, Sigma_full = Sf,
      cross_covariance_half_full = cross, covariance_identity_error = identity,
      H_population = output$hessian, H_half = Hh, H_full = Hf
    )
  }
  output
}

profile_pair_population_evaluate_v2 <- function(dat, beta, tau, h, lambda_gamma,
                                                control, level = 2L) {
  if (identical(control$population_engine, "conditional_gaussian")) {
    return(profile_pair_conditional_gaussian_components_v2(
      dat, beta, tau, h, lambda_gamma, nuisance_tol = min(control$nuisance_tol, 1e-9)))
  }
  if (identical(control$population_engine, "reference")) {
    return(profile_bandwidth_pair_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, h, lambda_gamma,
      identity_check = TRUE, need_fold_sums = FALSE,
      nuisance_control = list(grad_tol = control$nuisance_tol, reltol = 1e-12, maxit = 500L)
    ))
  }
  if (!identical(control$population_engine, "streamed")) stop("Unknown population engine.")
  profile_pair_population_streamed_v2(
    dat, beta, tau, h, lambda_gamma, batch_clusters = control$batch_clusters,
    need_hessian = level >= 1L, need_moments = level >= 2L,
    nuisance_tol = min(control$nuisance_tol, 1e-9)
  )
}

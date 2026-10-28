jp_joint_qp <- function(X, Z, outcome, cluster, inverse_size, Lambda, coordinate_penalty,
                         quantile, theta, gamma, dual, max_iterations = 200L, gap_tolerance = 1e-8) {
  N <- nrow(X)
  d <- ncol(X)
  q <- ncol(Z)
  n <- nrow(gamma)
  if (any(coordinate_penalty != 0) || length(cluster) != N || any(cluster < 0 | cluster >= n)) stop("Invalid selected quadratic programme.")
  block <- Matrix::sparseMatrix(i = rep(seq_len(N), each = q),
                                 j = as.vector(t(outer(cluster * q, seq_len(q), "+"))),
                                 x = as.vector(t(Z)), dims = c(N, n * q))
  equality <- cbind(Matrix::Matrix(X, sparse = TRUE), block, Matrix::Diagonal(N), -Matrix::Diagonal(N))
  inequality <- cbind(Matrix::Matrix(0, 2L * N, d + n * q, sparse = TRUE), -Matrix::Diagonal(2L * N))
  A <- as(as(rbind(equality, inequality), "generalMatrix"), "CsparseMatrix")
  P <- Matrix::bdiag(Matrix::Diagonal(d), kronecker(Matrix::Diagonal(n), Matrix::Matrix(Lambda, sparse = TRUE)),
                      Matrix::Matrix(0, 2L * N, 2L * N, sparse = TRUE))
  P <- as(as(P, "generalMatrix"), "CsparseMatrix")
  linear <- c(numeric(d + n * q), quantile * inverse_size, (1 - quantile) * inverse_size)
  control <- clarabel::clarabel_control(verbose = FALSE, max_iter = 200L,
                                        tol_gap_abs = 1e-10, tol_gap_rel = 1e-10, tol_feas = 1e-10)
  messages <- capture.output(solution <- clarabel::clarabel(A = A, b = c(outcome, numeric(2L * N)),
                                                            q = linear, P = P,
                                                            cones = list(z = as.integer(N), l = as.integer(2L * N)),
                                                            control = control))
  status <- names(clarabel::solver_status_descriptions())[solution$status]
  if (!status %in% c("Solved", "AlmostSolved")) stop(paste("Selected quadratic solver status", status))
  theta <- solution$x[seq_len(d)]
  gamma <- matrix(solution$x[d + seq_len(n * q)], nrow = n, ncol = q, byrow = TRUE)
  dual <- pmin(quantile * inverse_size, pmax((quantile - 1) * inverse_size, -solution$z[seq_len(N)]))
  residual <- outcome - drop(X %*% theta) - rowSums(Z * gamma[cluster + 1L, , drop = FALSE])
  primal <- sum(inverse_size * residual * (quantile - (residual < 0))) +
    sum(theta^2) / 2 + sum(gamma * (gamma %*% Lambda)) / 2
  beta_gradient <- drop(crossprod(X, dual))
  gamma_gradient <- rowsum(Z * dual, cluster + 1L, reorder = TRUE)[as.character(seq_len(n)), , drop = FALSE]
  gamma_conjugate <- t(solve(Lambda, t(gamma_gradient)))
  dual_value <- sum(outcome * dual) - sum(beta_gradient^2) / 2 - sum(gamma_gradient * gamma_conjugate) / 2
  gap <- (primal - dual_value) / n
  if (!is.finite(gap) || gap < -1e-10 * max(1, abs(primal / n)) ||
      gap > gap_tolerance * max(1, abs(primal / n))) stop(paste("Selected quadratic gap certificate failed", gap))
  list(theta = theta, gamma = gamma, dual = dual, primal = primal / n, dual_value = dual_value / n,
       duality_gap = gap, beta_kkt = max(abs(theta - beta_gradient)) / n,
       gamma_kkt = max(abs(gamma %*% Lambda - gamma_gradient)), iterations = solution$iterations,
       solver_status = status, solver = "clarabel_selected_qp", messages = messages)
}

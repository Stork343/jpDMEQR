submission_clime_parametric_basis <- function(H, coordinate, gamma) {
  invisible(loadNamespace("fastclime"))
  p <- nrow(H)
  canonical <- numeric(p)
  canonical[coordinate] <- 1
  steps <- 4L * p
  path <- .C("dantzig", as.double(H), as.double(canonical), beta = numeric(p * steps),
    as.integer(p), as.double(gamma), as.integer(steps), gamma = numeric(steps), PACKAGE = "fastclime")
  valid <- which(path$gamma > 0)
  if (!length(valid) || any(!is.finite(path$beta))) stop("The parametric CLIME path is nonfinite or empty.")
  knots <- path$gamma[valid]
  coefficients <- matrix(path$beta, p, steps)[, valid, drop = FALSE]
  upper <- max(which(knots >= gamma))
  lower <- min(upper + 1L, length(knots))
  if (lower > upper && knots[upper] > knots[lower]) {
    weight <- (gamma - knots[lower]) / (knots[upper] - knots[lower])
    omega <- weight * coefficients[, upper] + (1 - weight) * coefficients[, lower]
  } else omega <- coefficients[, upper]
  submission_clime_basis(H, coordinate, gamma, omega)
}

submission_clime_basis <- function(H, coordinate, gamma, omega) {
  p <- nrow(H)
  canonical <- numeric(p)
  canonical[coordinate] <- 1
  # 候选方向提供basis，HiGHS随后完成原始LP并核查primal和dual条件。
  support <- which(abs(omega) > 0)
  col_status <- rep(0L, 2L * p)
  row_status <- rep(1L, p)
  if (length(support)) {
    defect <- as.numeric(H %*% omega - canonical)
    ordered <- order(abs(abs(defect) - gamma))
    factorisation <- qr(t(H[ordered, support, drop = FALSE]), tol = 1e-10, LAPACK = FALSE)
    if (factorisation$rank != length(support)) stop("The parametric CLIME support is rank deficient.")
    active <- ordered[factorisation$pivot[seq_along(support)]]
    col_status[support + p * (omega[support] < 0)] <- 1L
    row_status[active] <- ifelse(defect[active] < 0, 0L, 2L)
  }
  list(col_status = col_status, row_status = row_status)
}

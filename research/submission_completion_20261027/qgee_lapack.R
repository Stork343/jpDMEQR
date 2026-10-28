submission_qgee_inverse_lapack <- function(matrix, rhs, tolerance = 1.490116e-8) {
  if (nrow(matrix) < 1L || nrow(matrix) != ncol(matrix) || nrow(rhs) != nrow(matrix) ||
      any(!is.finite(matrix)) || any(!is.finite(rhs))) stop("Invalid QGEE inverse dimensions or entries.")
  decomposition <- La.svd(matrix, nu = nrow(matrix), nv = ncol(matrix))
  keep <- which(decomposition$d > tolerance * decomposition$d[1L])
  if (!length(keep)) return(matrix(0, ncol(matrix), ncol(rhs)))
  projected <- crossprod(decomposition$u[, keep, drop = FALSE], rhs) / decomposition$d[keep]
  crossprod(decomposition$vt[keep, , drop = FALSE], projected)
}

profile_dense_crossprod_v2 <- function(X, Y = NULL) {
  backend <- getOption("jpDMEQR.matrix_backend",
                       Sys.getenv("JPDMEQR_MATRIX_BACKEND", unset = "R_BLAS"))
  if (identical(backend, "R_BLAS")) return(base::crossprod(X, Y))
  if (!identical(backend, "Eigen_single_thread") ||
      !is.loaded("jpDMEQR_population_crossprod", PACKAGE = "jpDMEQR")) {
    stop("The declared profile matrix backend is unavailable.")
  }
  X <- as.matrix(X)
  if (!is.null(Y)) Y <- as.matrix(Y)
  output <- .Call("jpDMEQR_population_crossprod", X, Y, PACKAGE = "jpDMEQR")
  labels <- list(colnames(X), if (is.null(Y)) colnames(X) else colnames(Y))
  if (any(!vapply(labels, is.null, logical(1)))) dimnames(output) <- labels
  output
}

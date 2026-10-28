compiled_cluster_profiles_v2 <- function(y_minus_x_beta, Z, cluster_rows, tau, h,
                                          Lambda, gamma_init = NULL,
                                          control = list()) {
  if (!requireNamespace("jpDMEQR", quietly = TRUE) ||
      !is.loaded("jpdmeqr_cluster_profiles", PACKAGE = "jpDMEQR")) {
    stop("The compiled jpDMEQR cluster-profile engine must be installed and loaded.")
  }
  q <- ncol(Z)
  n <- length(cluster_rows)
  levels <- names(cluster_rows)
  initial <- if (is.null(gamma_init)) matrix(0, n, q) else {
    matrix(vapply(seq_len(n), function(i) {
      resolve_gamma_init_v2(gamma_init, levels[i], i, q)
    }, numeric(q)), nrow = n, ncol = q, byrow = TRUE)
  }
  .Call("jpdmeqr_cluster_profiles", as.numeric(y_minus_x_beta), as.matrix(Z),
        as.integer(unlist(cluster_rows, use.names = FALSE) - 1L),
        as.integer(c(0L, cumsum(lengths(cluster_rows)))), initial,
        as.numeric(tau), as.numeric(h), as.matrix(Lambda), control,
        PACKAGE = "jpDMEQR")
}

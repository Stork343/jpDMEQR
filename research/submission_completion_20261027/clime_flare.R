submission_flare_matrix <- function(H, gamma, precision = 1e-5, max_iterations = 10000L) {
  invisible(loadNamespace("flare"))
  p <- nrow(H)
  gammas <- sort(unique(gamma), decreasing = TRUE)
  count <- length(gammas)
  gram <- profile_dense_crossprod_v2(H, H)
  spectral_bound <- max(colSums(abs(gram)))
  first_size <- if (p <= 3L) p else ceiling(p / log(p))
  second_size <- if (p <= 3L) p else p - 1L
  matrices <- lapply(gammas, function(x) matrix(0, p, p))
  iterations <- vector("list", p)
  for (j in seq_len(p)) {
    ordering <- order(H[, j], decreasing = TRUE)
    first <- head(ordering, first_size)
    second <- head(ordering, second_size)
    fitted <- .C("sugm_clime_ladm_scr", as.double(H), as.double(gram[ordering, ordering]),
      coefficients = numeric(p * count), sparse_values = numeric(p * count), as.integer(p),
      as.double(spectral_bound), as.double(gammas), as.integer(count), as.double(1),
      nonzero = integer(1L), sparse_rows = integer(p * count),
      full_iterations = integer(count), first_iterations = integer(count), second_iterations = integer(count),
      as.integer(first_size), as.integer(second_size), as.integer(ordering), as.integer(first), as.integer(second),
      as.integer(max_iterations), as.double(precision), as.integer(j), PACKAGE = "flare")
    columns <- matrix(fitted$coefficients, p, count)
    if (any(!is.finite(columns))) stop("The flare CLIME direction is nonfinite.")
    for (index in seq_len(count)) matrices[[index]][, j] <- columns[, index]
    iterations[[j]] <- cbind(first = fitted$first_iterations, second = fitted$second_iterations,
      full = fitted$full_iterations)
  }
  list(matrices = matrices, gamma = gammas, iterations = iterations, precision = precision,
       maximum_iterations = max_iterations, package_version = as.character(packageVersion("flare")))
}

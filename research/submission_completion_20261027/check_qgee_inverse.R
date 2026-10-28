source("research/submission_completion_20261027/load_comparators.R")
source("research/submission_completion_20261027/qgee_lapack.R")
records <- list()
for (p in c(20L, 80L, 200L)) for (kind in c("full", "singular", "threshold")) {
  set.seed(202610280L + p)
  first <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  second <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  values <- seq(1, .1, length.out = p)
  if (kind == "singular") values[seq.int(p %/% 2L, p)] <- 0
  if (kind == "threshold") values[(p - 1L):p] <- c(2, .5) * 1.490116e-8
  matrix <- first %*% diag(values) %*% t(second)
  rhs <- matrix(rnorm(p * 2L), p, 2L)
  reference <- get("geninv", asNamespace("geeVerse"))(matrix) %*% rhs
  candidate <- submission_qgee_inverse_lapack(matrix, rhs)
  relative_error <- max(abs(reference - candidate)) / max(1, max(abs(reference)))
  if (relative_error > 2e-8) stop("The QGEE generalised inverse action differs from its source tolerance rule.")
  records[[length(records) + 1L]] <- data.frame(p = p, kind = kind, relative_error = relative_error)
}
write.csv(do.call(rbind, records), "results/scratch/submission_completion_20261027/qgee_lapack_inverse_checks.csv", row.names = FALSE)
print(do.call(rbind, records))

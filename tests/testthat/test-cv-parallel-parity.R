# Parity: parallel CLARABEL mu-CV (jpDMEQR.cv_cores > 1) vs serial reference.
# The (fold x mu) grid is solved concurrently (forked mclapply off Windows,
# deterministic PSOCK on Windows) and aggregated in the original fold order;
# the selected mu, defect statistics, and every candidate row must be
# bit-identical to the serial path. This is the correctness gate for the
# dominant cost term of the simulation (direction-1 optimisation).

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
for (mod in c("profile_v2", "simulation_v2", "metrics_v2")) {
  source_v2_module(root, mod, envir = .GlobalEnv)
}

make_synthetic_fold_inputs <- function(p, n, n_folds = 4L, seed = 7L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  Hhat <- crossprod(X) / n + diag(0.05, p)
  Hhat <- (Hhat + t(Hhat)) / 2
  colnames(Hhat) <- sprintf("x%04d", seq_len(p))
  bucket <- ((seq_len(n) - 1L) %% n_folds) + 1L
  fold_sums <- lapply(seq_len(n_folds), function(f) {
    idx <- which(bucket == f)
    crossprod(X[idx, , drop = FALSE]) / n + diag(0.05, p) * length(idx) / n
  })
  fold_counts <- tabulate(bucket, nbins = n_folds)
  list(H = Hhat, fold_sums = fold_sums, fold_counts = fold_counts, n = n)
}

run_cv_with_cores <- function(inp, coordinate, mu_grid, cv_cores) {
  old <- getOption("jpDMEQR.cv_cores", 1L)
  options(jpDMEQR.cv_cores = cv_cores)
  on.exit(options(jpDMEQR.cv_cores = old), add = TRUE)
  select_mu_inverse_hessian_cv_v2(
    H = inp$H, coordinate = coordinate, mu_grid = mu_grid,
    fold_sums = inp$fold_sums, fold_counts = inp$fold_counts
  )
}

selection_identical <- function(a, b) {
  if (!identical(a$selected, b$selected)) return(FALSE)
  if (isTRUE(a$selected) && isTRUE(b$selected)) {
    if (!identical(a$mu, b$mu)) return(FALSE)
    if (!identical(a$min_defect, b$min_defect)) return(FALSE)
    if (!identical(a$defect_mean, b$defect_mean)) return(FALSE)
    if (!identical(a$quad_mean, b$quad_mean)) return(FALSE)
    ca <- a$candidates; cb <- b$candidates
    if (!identical(dim(ca), dim(cb))) return(FALSE)
    if (!isTRUE(all.equal(unname(as.matrix(ca)), unname(as.matrix(cb)),
                          tolerance = 0, check.attributes = FALSE))) return(FALSE)
  }
  TRUE
}

testthat::test_that("parallel mu-CV is bit-identical to serial CLARABEL reference", {
  p <- 60L
  n <- 240L
  n_folds <- 4L
  mu_grid <- 10^seq(-1.6, 0.5, length.out = 4)
  inp <- make_synthetic_fold_inputs(p, n, n_folds = n_folds)
  for (k in c(1L, 17L, 41L)) {
    for (cores in c(2L, 4L)) {
      serial <- run_cv_with_cores(inp, k, mu_grid, 1L)
      par <- run_cv_with_cores(inp, k, mu_grid, cores)
      testthat::expect_true(
        selection_identical(serial, par),
        info = sprintf("coord=%d cores=%d: parallel mu-CV differs from serial", k, cores)
      )
    }
  }
})
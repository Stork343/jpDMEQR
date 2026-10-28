# Parity: compact Dantzig LP solved by HiGHS with its default solver choice vs
# the frozen CVXR/CLARABEL reference and vs the compact CLARABEL cone form. All
# three solve the SAME L1 program with the same acceptance rule (primal residual
# <= mu + 1e-6), so decisions must agree:
#   - same feasibility pattern over the mu grid,
#   - same selected mu (smallest feasible),
#   - omega equal far below the acceptance tolerance.
#
# IMPORTANT: the size range deliberately includes p >= 200. The Lever-1 chained
# OSQP path passed its parity test at p = 20/60/150 and only failed at p >= 200
# (see COMPUTATIONAL_OPTIMISATION_LEVER1_ADJUDICATION.md section 6), so a test
# capped below p = 200 does not protect a solver substitution. Measured HiGHS
# agreement is max|domega| = 7.8e-08 at p = 300 and 7.4e-08 at p = 600, and the
# mu-CV selected the same mu at both sizes.
#
# The path itself is adjudicated inadvisable: it is faster on these
# well-conditioned constructions and about 11.5x slower on the real profile
# Hessian (COMPUTATIONAL_OPTIMISATION_HIGHS_PATH_REPORT.md section 5). This file
# protects the decision equivalence of the switch, not its speed.

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
for (mod in c("profile_v2", "simulation_v2", "metrics_v2")) {
  source_v2_module(root, mod, envir = .GlobalEnv)
}

make_h <- function(p, seed, n_scale = 5L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n_scale * p * p), nrow = n_scale * p, ncol = p)
  crossprod(X) / (n_scale * p) + diag(0.05, p)
}

with_solver <- function(mode, expr) {
  old <- Sys.getenv("JPDMEQR_DANTZIG_SOLVER", unset = NA_character_)
  if (is.na(mode)) Sys.unsetenv("JPDMEQR_DANTZIG_SOLVER") else
    Sys.setenv(JPDMEQR_DANTZIG_SOLVER = mode)
  on.exit({
    if (is.na(old)) Sys.unsetenv("JPDMEQR_DANTZIG_SOLVER") else
      Sys.setenv(JPDMEQR_DANTZIG_SOLVER = old)
  }, add = TRUE)
  force(expr)
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

testthat::test_that("HiGHS compact LP reproduces CVXR/CLARABEL Dantzig decisions", {
  testthat::skip_if_not_installed("highs")
  for (p in c(20L, 60L, 150L, 200L)) {
    H <- make_h(p, seed = 100L + p)
    colnames(H) <- sprintf("x%05d", seq_len(p))
    mu_grid <- 10^seq(-1.6, 0.5, length.out = 9)
    ks <- if (p <= 150L) 1:3 else 1:2

    for (k in ks) {
      ref <- with_solver(NA, solve_dantzig_row_v2(H, k, mu_grid))
      hig <- with_solver("highs", solve_dantzig_row_v2(H, k, mu_grid))

      testthat::expect_equal(hig$feasible, ref$feasible,
                             info = sprintf("p=%d k=%d feasibility pattern", p, k))
      if (isTRUE(ref$feasible) && isTRUE(hig$feasible)) {
        testthat::expect_equal(hig$mu, ref$mu,
                               info = sprintf("p=%d k=%d selected mu", p, k))
        testthat::expect_true(hig$residual <= hig$mu + 1e-6,
                              info = sprintf("p=%d k=%d residual within mu", p, k))
        rel <- max(abs(hig$omega - ref$omega)) / max(1e-6, sum(abs(ref$omega)))
        testthat::expect_true(rel < 1e-6,
                              info = sprintf("p=%d k=%d omega rel diff %.3g", p, k, rel))
        testthat::expect_equal(hig$l1_norm, ref$l1_norm, tolerance = 1e-6,
                               info = sprintf("p=%d k=%d l1 norm", p, k))
      }
    }
  }
})

testthat::test_that("HiGHS mu-CV selects the same mu and the same defect table", {
  testthat::skip_if_not_installed("highs")
  p <- 60L
  n <- 240L
  mu_grid <- 10^seq(-1.6, 0.5, length.out = 4)
  inp <- make_synthetic_fold_inputs(p, n, n_folds = 4L)
  for (k in c(1L, 17L, 41L)) {
    ref <- with_solver(NA, select_mu_inverse_hessian_cv_v2(
      H = inp$H, coordinate = k, mu_grid = mu_grid,
      fold_sums = inp$fold_sums, fold_counts = inp$fold_counts))
    hig <- with_solver("highs", select_mu_inverse_hessian_cv_v2(
      H = inp$H, coordinate = k, mu_grid = mu_grid,
      fold_sums = inp$fold_sums, fold_counts = inp$fold_counts))

    testthat::expect_identical(hig$selected, ref$selected,
                               info = sprintf("coord=%d selected flag", k))
    if (isTRUE(ref$selected)) {
      testthat::expect_equal(hig$mu, ref$mu,
                             info = sprintf("coord=%d selected mu", k))
      testthat::expect_equal(hig$defect_mean, ref$defect_mean, tolerance = 1e-5,
                             info = sprintf("coord=%d defect mean", k))
      ca <- ref$candidates; cb <- hig$candidates
      testthat::expect_identical(ca$folds_used, cb$folds_used,
                                 info = sprintf("coord=%d folds used", k))
      dev <- max(abs(ca$defect_mean - cb$defect_mean), na.rm = TRUE)
      testthat::expect_lt(dev, 1e-5,
                          label = sprintf("coord=%d max defect deviation %.3g", k, dev))
    }
  }
})

testthat::test_that("highs path is inert when the dispatcher is not selected", {
  testthat::skip_if_not_installed("highs")
  p <- 60L
  H <- make_h(p, seed = 999L)
  colnames(H) <- sprintf("x%05d", seq_len(p))
  ref <- with_solver(NA, solve_dantzig_row_v2(H, 1L, 10^seq(-1.6, 0.5, length.out = 5)))
  testthat::expect_false(identical(ref$solver, "HiGHS"))
})

testthat::test_that("HiGHS status labels map the documented model-status enum", {
  testthat::expect_identical(highs_model_status_label_v2(7), "optimal")
  testthat::expect_identical(highs_model_status_label_v2(8), "infeasible")
  testthat::expect_identical(highs_model_status_label_v2(14), "iteration_limit")
  testthat::expect_identical(highs_model_status_label_v2(17), "interrupt")
  testthat::expect_identical(highs_model_status_label_v2(99), "status_99")
  testthat::expect_identical(highs_model_status_label_v2(NA_integer_), "status_NA")
})

testthat::test_that("highs path is bit-identical between serial and parallel mu-CV", {
  withr::local_envvar(c(JPDMEQR_CV_CORES = NA))
  # The parallel CV exports the solver chain to the workers by name, so the
  # highs functions must be in that list. Solving the same (fold, mu) tasks on
  # workers returns the same numbers, so a serial and a parallel run under the
  # highs dispatcher must agree exactly: a silent fall-back to the reference
  # solver on the workers would break bit-identity.
  testthat::skip_if_not_installed("highs")
  inp <- make_synthetic_fold_inputs(p = 30L, n = 80L, n_folds = 4L, seed = 11L)
  mu_grid <- 10^seq(-1.6, 0.5, length.out = 3)
  run_cv <- function(k, cores) {
    old <- getOption("jpDMEQR.cv_cores", 1L)
    options(jpDMEQR.cv_cores = cores)
    on.exit(options(jpDMEQR.cv_cores = old), add = TRUE)
    with_solver("highs", select_mu_inverse_hessian_cv_v2(
      H = inp$H, coordinate = k, mu_grid = mu_grid,
      fold_sums = inp$fold_sums, fold_counts = inp$fold_counts))
  }
  for (k in c(1L, 17L)) {
    serial <- run_cv(k, 1L)
    par <- run_cv(k, 2L)
    testthat::expect_identical(par$selected, serial$selected,
                               info = sprintf("coord=%d selected flag", k))
    testthat::expect_identical(par$mu, serial$mu,
                               info = sprintf("coord=%d selected mu", k))
    testthat::expect_identical(unname(as.matrix(par$candidates)),
                               unname(as.matrix(serial$candidates)),
                               info = sprintf("coord=%d candidate table", k))
  }
})

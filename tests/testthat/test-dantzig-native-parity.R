# Parity: compact native CLARABEL Dantzig path (layer-1) vs the frozen
# CVXR/CLARABEL reference path. Both solve the SAME L1 program; decisions must
# agree:
#   - same feasibility pattern over the mu grid,
#   - same selected mu (smallest feasible),
#   - omega equal to numerical precision.
#
# IMPORTANT: the size range deliberately includes p >= 200. The Lever-1 chained
# OSQP path passed its parity test at p = 20/60/150 and only failed at p >= 200
# (see COMPUTATIONAL_OPTIMISATION_LEVER1_ADJUDICATION.md section 6), so a test
# capped below p = 200 does not actually protect this substitution. Measured
# agreement for the native path is max|domega| ~ 2.2e-16 at p = 200 and 500.

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
for (mod in c("profile_v2", "simulation_v2")) {
  source_v2_module(root, mod, envir = .GlobalEnv)
}

make_h <- function(p, seed, n_scale = 5L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n_scale * p * p), nrow = n_scale * p, ncol = p)
  H <- crossprod(X) / (n_scale * p) + diag(0.05, p)
  H
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

testthat::test_that("native compact Dantzig path matches CVXR/CLARABEL decisions", {
  # p = 200 is the smallest size where the rejected Lever-1 path diverged, so it
  # is mandatory here; p = 500 is excluded from the CI test only for runtime.
  for (p in c(20L, 60L, 150L, 200L)) {
    H <- make_h(p, seed = 100L + p)
    colnames(H) <- sprintf("x%05d", seq_len(p))
    mu_grid <- 10^seq(-1.6, 0.5, length.out = 9)
    ks <- if (p <= 150L) 1:3 else 1:2

    for (k in ks) {
      ref <- with_solver(NA, solve_dantzig_row_v2(H, k, mu_grid))
      nat <- with_solver("clarabel_native", solve_dantzig_row_v2(H, k, mu_grid))

      testthat::expect_equal(nat$feasible, ref$feasible,
                             info = sprintf("p=%d k=%d feasibility pattern", p, k))
      if (isTRUE(ref$feasible) && isTRUE(nat$feasible)) {
        testthat::expect_equal(nat$mu, ref$mu,
                               info = sprintf("p=%d k=%d selected mu", p, k))
        # Same acceptance tolerance, so residuals must both satisfy mu + 1e-6.
        testthat::expect_true(nat$residual <= nat$mu + 1e-6,
                              info = sprintf("p=%d k=%d native residual within mu", p, k))
        rel <- max(abs(nat$omega - ref$omega)) / max(1e-6, sum(abs(ref$omega)))
        testthat::expect_true(rel < 1e-8,
                              info = sprintf("p=%d k=%d omega rel diff %.3g", p, k, rel))
      }
    }
  }
})

testthat::test_that("native path is inert when the dispatcher is not selected", {
  # The default path must remain the frozen CVXR reference: no env var, no
  # option -> solver reported by the reference path.
  p <- 60L
  H <- make_h(p, seed = 999L)
  colnames(H) <- sprintf("x%05d", seq_len(p))
  mu_grid <- 10^seq(-1.6, 0.5, length.out = 5)
  ref <- with_solver(NA, solve_dantzig_row_v2(H, 1L, mu_grid))
  testthat::expect_false(identical(ref$solver, "CLARABEL-native"))
})

testthat::test_that("compact form reproduces the redundant-constraint form within tolerance", {
  # Form C drops the implied `-t <= 0` block. The two forms have the same
  # feasible set and objective, so they agree to solver tolerance. Measured
  # difference is ~3e-10 (the interior-point iteration path differs once the
  # redundant block is removed), i.e. NOT bit-identical; the tight agreement
  # with the CVXR reference (2.2e-16) is a separate claim asserted above.
  p <- 40L
  H <- make_h(p, seed = 4242L)
  H <- (H + t(H)) / 2
  e <- numeric(p); e[1L] <- 1
  mu <- 0.05

  A_compact <- dantzig_lp_compact_build_v2(H, p)
  Ip <- Matrix::Diagonal(p); Zp <- Matrix::Matrix(0, p, p, sparse = TRUE)
  A_full <- methods::as(methods::as(
    rbind(cbind(Ip, -Ip), cbind(-Ip, -Ip), cbind(Zp, -Ip),
          cbind(H, Zp), cbind(-H, Zp)), "CsparseMatrix"), "dgCMatrix")

  q <- c(numeric(p), rep(1, p))
  r_c <- clarabel::clarabel(A = A_compact, b = dantzig_lp_compact_rhs_v2(p, e, mu),
                            q = q, P = NULL, cones = list(l = nrow(A_compact)),
                            control = list(verbose = FALSE))
  b_full <- c(numeric(p), numeric(p), numeric(p), e + mu, -e + mu)
  r_f <- clarabel::clarabel(A = A_full, b = b_full, q = q, P = NULL,
                            cones = list(l = nrow(A_full)),
                            control = list(verbose = FALSE))

  om_c <- as.numeric(r_c$x[seq_len(p)])
  om_f <- as.numeric(r_f$x[seq_len(p)])
  testthat::expect_lt(max(abs(om_c - om_f)), 1e-8)
  testthat::expect_lte(max(abs(as.numeric(H %*% om_c) - e)), mu + 1e-6)
})
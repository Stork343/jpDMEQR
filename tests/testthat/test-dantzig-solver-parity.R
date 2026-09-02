# Parity: chained-OSQP Dantzig path (Lever 1) vs CVXR/CLARABEL reference.
# The frozen program is solved by both; decisions must agree:
#   - same feasibility pattern over the mu grid,
#   - same selected mu (smallest feasible),
#   - omega within 1e-4 (max abs).

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

testthat::test_that("chained-OSQP Dantzig path matches CVXR/CLARABEL decisions", {
  for (p in c(20L, 60L, 150L)) {
    H <- make_h(p, seed = 100L + p)
    colnames(H) <- sprintf("x%05d", seq_len(p))
    mu_grid <- 10^seq(-1.6, 0.5, length.out = 9)
    for (k in 1:3) {
      Sys.unsetenv("JPDMEQR_DANTZIG_SOLVER")
      ref <- solve_dantzig_row_v2(H, k, mu_grid)                  # CLARABEL path
      Sys.setenv(JPDMEQR_DANTZIG_SOLVER = "osqp")
      fast <- solve_dantzig_row_v2(H, k, mu_grid)                 # chained OSQP
      Sys.unsetenv("JPDMEQR_DANTZIG_SOLVER")
      testthat::expect_equal(fast$feasible, ref$feasible,
                             info = sprintf("p=%d k=%d feasibility pattern", p, k))
      if (ref$feasible && fast$feasible) {
        testthat::expect_equal(fast$mu, ref$mu,
                               info = sprintf("p=%d k=%d selected mu", p, k))
        rel <- max(abs(fast$omega - ref$omega)) / max(1e-6, sum(abs(ref$omega)))
        testthat::expect_true(rel < 1e-4,
                              info = sprintf("p=%d k=%d omega rel diff %.3g", p, k, rel))
      }
    }
  }
})
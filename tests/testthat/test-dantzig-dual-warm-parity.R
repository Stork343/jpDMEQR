# Parity: the dual form min e_k'w + mu||w||_1 s.t. ||H'w||_inf <= 1 solved by a
# persistent HiGHS simplex model, with the basis carried down the descending mu
# grid, against the frozen CVXR/CLARABEL row path. Both solve the same L1 program
# with the same acceptance rule (primal residual <= mu + 1e-6), so decisions must
# agree:
#   - same feasibility pattern over the mu grid,
#   - same selected mu (smallest feasible),
#   - omega equal to a relative 1e-4, which is the solver-specific difference
#     measured on the registered cell (max|domega| = 1.1e-05 at p = 500 on a
#     matrix with condition number 4.9e+18).
# The grid sweep also settles every candidate below an unbounded one by
# monotonicity, so the statuses must show that settlement explicitly.

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

# 秩亏的 Hessian：列数多于行数，H omega = e_k 无解，小 mu 上的行规划不可行
make_singular_h <- function(p, seed) {
  set.seed(seed)
  X <- matrix(stats::rnorm((p %/% 2L) * p), nrow = p %/% 2L, ncol = p)
  H <- crossprod(X) / nrow(X)
  (H + t(H)) / 2
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

with_cores <- function(cores, expr) {
  withr::local_envvar(c(JPDMEQR_CV_CORES = NA))
  old <- getOption("jpDMEQR.cv_cores", 1L)
  options(jpDMEQR.cv_cores = cores)
  on.exit(options(jpDMEQR.cv_cores = old), add = TRUE)
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

testthat::test_that("dual warm path reproduces the frozen Dantzig row decisions", {
  testthat::skip_if_not_installed("highs")
  for (p in c(20L, 60L, 150L, 200L)) {
    H <- make_h(p, seed = 100L + p)
    colnames(H) <- sprintf("x%05d", seq_len(p))
    mu_grid <- 10^seq(-1.6, 0.5, length.out = 9)
    ks <- if (p <= 150L) 1:3 else 1:2

    for (k in ks) {
      ref <- with_solver(NA, solve_dantzig_row_v2(H, k, mu_grid))
      dw <- with_solver("dual_warm", solve_dantzig_row_v2(H, k, mu_grid))

      testthat::expect_true(dw$feasible,
                            info = sprintf("p=%d k=%d row feasible", p, k))
      testthat::expect_equal(dw$feasible, ref$feasible,
                             info = sprintf("p=%d k=%d feasibility pattern", p, k))
      if (isTRUE(ref$feasible) && isTRUE(dw$feasible)) {
        testthat::expect_equal(dw$mu, ref$mu,
                               info = sprintf("p=%d k=%d selected mu", p, k))
        testthat::expect_true(dw$residual <= dw$mu + 1e-6,
                              info = sprintf("p=%d k=%d residual within mu", p, k))
        rel <- max(abs(dw$omega - ref$omega)) / max(1e-6, sum(abs(ref$omega)))
        testthat::expect_true(rel < 1e-4,
                              info = sprintf("p=%d k=%d omega rel diff %.3g", p, k, rel))
        testthat::expect_equal(dw$l1_norm, ref$l1_norm, tolerance = 1e-4,
                               info = sprintf("p=%d k=%d l1 norm", p, k))
      }
    }
  }
})

testthat::test_that("an unbounded dual candidate settles every smaller candidate", {
  testthat::skip_if_not_installed("highs")
  p <- 60L
  H <- make_singular_h(p, seed = 41L)
  mu_grid <- c(0.02, 0.05, 0.1, 0.25, 0.5, 1, 2, 4)
  grid <- dantzig_dual_solve_grid_v2(H, 7L, mu_grid)

  testthat::expect_identical(grid$mu, mu_grid)
  # 降序扫描在 mu=0.1 上拿到无界证书，比它小的两个候选由单调性判定
  testthat::expect_identical(grid$status[1:2],
                             rep("infeasible_by_monotonicity", 2L))
  testthat::expect_identical(grid$status[3], "unbounded")
  testthat::expect_false(any(grid$feasible[1:3]))
  testthat::expect_true(all(grid$feasible[4:8]))
  testthat::expect_true(all(grid$residual[4:6] <= mu_grid[4:6] + 1e-6))

  # 冻结路径在每个候选上的可行性判定必须一致
  for (ii in seq_along(mu_grid)) {
    ref <- tryCatch(solve_dantzig_row_v2(H, 7L, mu_grid[ii]),
                    error = function(e) list(feasible = FALSE))
    testthat::expect_identical(isTRUE(ref$feasible), grid$feasible[ii],
                               info = sprintf("frozen feasibility at mu=%.4g", mu_grid[ii]))
  }
})

testthat::test_that("the row wrapper returns the smallest feasible candidate", {
  testthat::skip_if_not_installed("highs")
  p <- 60L
  H <- make_singular_h(p, seed = 42L)
  mu_grid <- c(0.02, 0.1, 0.5, 1, 2, 4)
  row <- solve_dantzig_row_dual_warm_v2(H, 3L, mu_grid)
  ok <- mu_grid[sapply(seq_along(mu_grid), function(ii) {
    r <- tryCatch(solve_dantzig_row_v2(H, 3L, mu_grid[ii]),
                  error = function(e) list(feasible = FALSE))
    isTRUE(r$feasible)
  })]
  testthat::expect_true(length(ok) >= 1L)
  testthat::expect_equal(row$mu, min(ok))
  testthat::expect_identical(row$solver, "HiGHS-dual-simplex")
  testthat::expect_true(row$residual <= row$mu + 1e-6)
  testthat::expect_identical(nrow(row$attempts), length(mu_grid))
})

testthat::test_that("an entirely infeasible grid returns the failure schema", {
  testthat::skip_if_not_installed("highs")
  p <- 40L
  H <- make_singular_h(p, seed = 43L)
  mu_grid <- c(0.01, 0.02)
  row <- solve_dantzig_row_dual_warm_v2(H, 3L, mu_grid)
  testthat::expect_false(row$feasible)
  testthat::expect_true(all(is.na(row$omega)))
  testthat::expect_true(is.na(row$mu))
  testthat::expect_identical(row$solver_status, "infeasible_on_grid")
  testthat::expect_identical(nrow(row$attempts), 2L)
})

testthat::test_that("dual warm mu-CV selects the same mu and the same defect table", {
  testthat::skip_if_not_installed("highs")
  p <- 60L
  n <- 240L
  mu_grid <- 10^seq(-1.6, 0.5, length.out = 4)
  inp <- make_synthetic_fold_inputs(p, n, n_folds = 4L)
  for (k in c(1L, 17L, 41L)) {
    ref <- with_solver(NA, select_mu_inverse_hessian_cv_v2(
      H = inp$H, coordinate = k, mu_grid = mu_grid,
      fold_sums = inp$fold_sums, fold_counts = inp$fold_counts))
    dw <- with_solver("dual_warm", select_mu_inverse_hessian_cv_v2(
      H = inp$H, coordinate = k, mu_grid = mu_grid,
      fold_sums = inp$fold_sums, fold_counts = inp$fold_counts))

    testthat::expect_identical(dw$selected, ref$selected,
                               info = sprintf("coord=%d selected flag", k))
    if (isTRUE(ref$selected)) {
      testthat::expect_equal(dw$mu, ref$mu,
                             info = sprintf("coord=%d selected mu", k))
      testthat::expect_equal(dw$defect_mean, ref$defect_mean, tolerance = 1e-5,
                             info = sprintf("coord=%d defect mean", k))
      ca <- ref$candidates; cb <- dw$candidates
      testthat::expect_identical(ca$folds_used, cb$folds_used,
                                 info = sprintf("coord=%d folds used", k))
      testthat::expect_identical(ca$infeasible_folds, cb$infeasible_folds,
                                 info = sprintf("coord=%d infeasible folds", k))
      dev <- max(abs(ca$defect_mean - cb$defect_mean), na.rm = TRUE)
      testthat::expect_lt(dev, 1e-5,
                          label = sprintf("coord=%d max defect deviation %.3g", k, dev))
    }
  }
})

testthat::test_that("the dual warm path is inert when the dispatcher is not selected", {
  testthat::skip_if_not_installed("highs")
  p <- 60L
  H <- make_h(p, seed = 999L)
  colnames(H) <- sprintf("x%05d", seq_len(p))
  ref <- with_solver(NA, solve_dantzig_row_v2(H, 1L, 10^seq(-1.6, 0.5, length.out = 5)))
  testthat::expect_false(identical(ref$solver, "HiGHS-dual-simplex"))
})

# 折与折之间没有共享状态，每个折在自己的训练矩阵上走完网格，因此把折分给
# 多个进程与逐个折串行调用逐位相同。
testthat::test_that("the fold-parallel dual warm mu-CV matches the serial run", {
  testthat::skip_if_not_installed("highs")
  p <- 60L
  n <- 240L
  mu_grid <- 10^seq(-1.6, 0.5, length.out = 4)
  inp <- make_synthetic_fold_inputs(p, n, n_folds = 4L)
  for (k in c(1L, 17L)) {
    ser <- with_cores(1L, with_solver("dual_warm", select_mu_inverse_hessian_cv_v2(
      H = inp$H, coordinate = k, mu_grid = mu_grid,
      fold_sums = inp$fold_sums, fold_counts = inp$fold_counts)))
    par <- with_cores(2L, with_solver("dual_warm", select_mu_inverse_hessian_cv_v2(
      H = inp$H, coordinate = k, mu_grid = mu_grid,
      fold_sums = inp$fold_sums, fold_counts = inp$fold_counts)))

    info <- sprintf("coord=%d", k)
    testthat::expect_identical(par$selected, ser$selected, info = info)
    testthat::expect_identical(par$mu, ser$mu, info = info)
    testthat::expect_identical(par$defect_mean, ser$defect_mean, info = info)
    testthat::expect_identical(par$candidates$mu, ser$candidates$mu, info = info)
    testthat::expect_identical(par$candidates$defect_mean,
                               ser$candidates$defect_mean, info = info)
    testthat::expect_identical(par$candidates$folds_used,
                               ser$candidates$folds_used, info = info)
    testthat::expect_identical(par$candidates$infeasible_folds,
                               ser$candidates$infeasible_folds, info = info)
  }
})

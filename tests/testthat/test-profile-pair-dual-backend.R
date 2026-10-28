root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", .GlobalEnv)
source_v2_module(root, "simulation_v2", .GlobalEnv)

testthat::test_that("dual warm grids solve the same real profile precision programs", {
  for (q in c(1L, 2L)) {
    dat <- generate_profile_qr_data_v2(
      n_clusters = 35L, p = 12L, s = 2L, tau = 0.5, q = q, seed = 20261400L + q
    )
    component <- profile_bandwidth_pair_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, dat$beta0, 0.5, 0.35
    )
    grid <- c(0.05, 0.15, 0.4, 0.7)
    dual <- profile_pair_solve_grid_v2(
      component$hessian, 1L, grid, list(backend = "dual_warm")
    )
    reference <- lapply(grid, function(mu) profile_pair_solve_row_v2(
      component$hessian, 1L, mu, list(backend = "clarabel_native")
    ))
    for (j in seq_along(grid)) {
      testthat::expect_identical(dual[[j]]$feasible, reference[[j]]$feasible)
      if (dual[[j]]$feasible) {
        testthat::expect_lte(dual[[j]]$residual, grid[j] + 1e-6)
        testthat::expect_equal(dual[[j]]$l1_norm, reference[[j]]$l1_norm, tolerance = 1e-6)
        testthat::expect_equal(dual[[j]]$omega, reference[[j]]$omega, tolerance = 1e-5)
      }
    }
  }
})

testthat::test_that("dual and native inference preserve the complete registered grid", {
  dat <- generate_profile_qr_data_v2(
    n_clusters = 60L, p = 10L, s = 2L, tau = 0.5, q = 2L, seed = 20261410L
  )
  targets <- colnames(dat$X)[c(1L, 4L)]
  fit <- fit_profile_pair_v2(dat$y, dat$X, dat$Z, dat$cluster_id, 0.5,
                            60^(-3 / 10), target_coordinates = targets)
  grid <- c(0.02, 0.05, 0.10, 0.25, 0.50, 1, 2, 4) *
    (sqrt(log(10) / (60 * fit$h_inf)) + fit$h_inf^2)
  native <- infer_profile_pair_v2(fit, targets, grid,
                                  solver_opts = list(backend = "clarabel_native"))
  dual <- infer_profile_pair_v2(fit, targets, grid,
                                solver_opts = list(backend = "dual_warm"))
  testthat::expect_equal(dual$mu_grid_requested, native$mu_grid_requested)
  testthat::expect_equal(dual$excluded_zero_direction, native$excluded_zero_direction)
  testthat::expect_equal(dual$table$mu, native$table$mu, tolerance = 1e-12)
  testthat::expect_equal(dual$table$beta_tilde, native$table$beta_tilde, tolerance = 1e-5)
  testthat::expect_equal(dual$table$se, native$table$se, tolerance = 1e-5)
  testthat::expect_true(all(dual$table$solver == "HiGHS-dual-simplex"))
  for (coordinate in targets) {
    testthat::expect_equal(dual$cv[[coordinate]]$candidates$eligible,
                           native$cv[[coordinate]]$candidates$eligible)
    testthat::expect_equal(dual$cv[[coordinate]]$candidates$defect_mean,
                           native$cv[[coordinate]]$candidates$defect_mean, tolerance = 1e-5)
  }
})

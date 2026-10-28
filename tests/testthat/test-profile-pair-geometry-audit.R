pair_geometry_root <- v2_require_root()
source(file.path(pair_geometry_root, "scripts", "00_source_v2.R"), local = TRUE)
source_v2_module(pair_geometry_root, "profile_v2", environment())
source_v2_module(pair_geometry_root, "simulation_v2", environment())
source(file.path(pair_geometry_root, "scripts", "research",
                 "audit_profile_pair_geometry_20261022.R"), local = TRUE)

testthat::test_that("the deterministic geometry matrix covers the declared factors", {
  grid <- profile_pair_geometry_grid_20261022()
  testthat::expect_equal(nrow(grid), 80L)
  testthat::expect_length(unique(grid$seed), 20L)
  testthat::expect_equal(sort(unique(grid$tau)), c(0.1, 0.25, 0.5, 0.75))
  testthat::expect_equal(sort(unique(grid$q)), c(1, 2))
  testthat::expect_equal(sort(unique(grid$h_multiplier)), c(0.75, 1, 1.25))
  testthat::expect_true(all(table(grid$tau, grid$q, grid$h_multiplier) > 0))
  testthat::expect_true(any(grid$balanced) && any(!grid$balanced))
  testthat::expect_true(any(grid$non_diagonal_lambda))
})

testthat::test_that("geometry accepts a correctly differentiated negative paired Hessian", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  dat <- list(y = rep(c(0.75, 0.8, 0.85), each = 2L),
              X = matrix(1, 6L, 1L, dimnames = list(NULL, "intercept")),
              Z = matrix(1, 6L, 1L), cluster_id = rep(1:3, each = 2L))
  result <- profile_pair_geometry_check_20261022(dat, 0, 0.5, 0.4, matrix(1, 1L),
                                                check_assembly = TRUE)
  testthat::expect_true(result$metrics$pass)
  testthat::expect_lt(result$metrics$min_paired_eigenvalue, 0)
  testthat::expect_lt(result$metrics$sandwich_identity_error, 1e-10)
  testthat::expect_lte(result$metrics$dantzig_residual, result$metrics$dantzig_mu + 1e-6)
})

testthat::test_that("singular paired geometry preserves verified Dantzig infeasibility", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "newton"))
  dat <- list(y = rep(c(0.75, 0.8, 0.85), each = 2L),
              X = matrix(1, 6L, 2L, dimnames = list(NULL, c("first", "second"))),
              Z = matrix(1, 6L, 1L), cluster_id = rep(1:3, each = 2L))
  result <- profile_pair_geometry_check_20261022(dat, c(0, 0), 0.5, 0.4,
                                                matrix(1, 1L))
  testthat::expect_true(result$metrics$pass)
  testthat::expect_false(result$metrics$dantzig_feasible)
  testthat::expect_true(result$metrics$dantzig_status_consistent)
  testthat::expect_equal(result$metrics$numerical_rank, 1L)
  testthat::expect_equal(result$metrics$independent_minimum_defect, 0.5, tolerance = 1e-7)
  testthat::expect_identical(result$metrics$independent_feasibility_solver, "ECOS")
})

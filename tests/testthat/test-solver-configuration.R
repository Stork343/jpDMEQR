root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)

unset_solver_env <- function(code) {
  withr::with_envvar(
    c(JPDMEQR_DANTZIG_SOLVER = NA, JPDMEQR_CLUSTER_SOLVER = NA,
      JPDMEQR_ASSEMBLY = NA, JPDMEQR_CV_CORES = NA,
      JPDMEQR_HESSIAN_REFRESH = NA),
    code
  )
}

unset_solver_options <- function(code) {
  withr::with_options(
    list(jpDMEQR.dantzig_solver = NULL, jpDMEQR.cluster_solver = NULL,
         jpDMEQR.assembly = NULL, jpDMEQR.cv_cores = NULL,
         jpDMEQR.hessian_refresh = NULL),
    code
  )
}

testthat::test_that("the runtime record reports the frozen defaults", {
  rec <- unset_solver_env(unset_solver_options(solver_configuration_v2()))
  testthat::expect_identical(rec$dantzig_solver, "cvxr")
  testthat::expect_identical(rec$cluster_solver, "bfgs")
  testthat::expect_identical(rec$assembly, "sequential")
  testthat::expect_identical(rec$cv_workers, 1L)
  testthat::expect_identical(rec$hessian_refresh, 5L)
})

testthat::test_that("the environment selects the recorded solver paths", {
  rec <- unset_solver_options(withr::with_envvar(
    c(JPDMEQR_DANTZIG_SOLVER = "dual_warm", JPDMEQR_CLUSTER_SOLVER = "newton",
      JPDMEQR_ASSEMBLY = "batch"),
    solver_configuration_v2()
  ))
  testthat::expect_identical(rec$dantzig_solver, "dual_warm")
  testthat::expect_identical(rec$cluster_solver, "newton")
  testthat::expect_identical(rec$assembly, "batch")
})

testthat::test_that("the option is recorded when the environment is unset", {
  rec <- unset_solver_env(withr::with_options(
    list(jpDMEQR.dantzig_solver = "clarabel_native",
         jpDMEQR.cluster_solver = "newton", jpDMEQR.assembly = "batch"),
    solver_configuration_v2()
  ))
  testthat::expect_identical(rec$dantzig_solver, "clarabel_native")
  testthat::expect_identical(rec$cluster_solver, "newton")
  testthat::expect_identical(rec$assembly, "batch")
})

testthat::test_that("the environment takes precedence over the option", {
  rec <- unset_solver_env(withr::with_options(
    list(jpDMEQR.assembly = "sequential"),
    withr::with_envvar(c(JPDMEQR_ASSEMBLY = "batch"),
                       solver_configuration_v2())
  ))
  testthat::expect_identical(rec$assembly, "batch")
})

testthat::test_that("the recorded worker count follows the cv cores option", {
  rec <- unset_solver_env(withr::with_options(
    list(jpDMEQR.cv_cores = 6L), solver_configuration_v2()
  ))
  testthat::expect_identical(rec$cv_workers, 6L)
})

testthat::test_that("the worker count follows the environment over the option", {
  rec <- withr::with_options(
    list(jpDMEQR.cv_cores = 2L),
    withr::with_envvar(c(JPDMEQR_CV_CORES = "5"), solver_configuration_v2())
  )
  testthat::expect_identical(rec$cv_workers, 5L)
})

testthat::test_that("a recorded configuration is compared against the current one", {
  current <- unset_solver_env(unset_solver_options(solver_configuration_v2()))
  testthat::expect_true(configuration_matches_v2(current, current))
  testthat::expect_false(configuration_matches_v2(NULL, current))

  other_assembly <- current
  other_assembly$assembly <- "batch"
  testthat::expect_false(configuration_matches_v2(other_assembly, current))

  other_workers <- current
  other_workers$cv_workers <- current$cv_workers + 1L
  testthat::expect_false(configuration_matches_v2(other_workers, current))

  missing_workers <- current
  missing_workers$cv_workers <- NULL
  testthat::expect_false(configuration_matches_v2(missing_workers, current))

  other_refresh <- current
  other_refresh$hessian_refresh <- 20L
  testthat::expect_false(configuration_matches_v2(other_refresh, current))

  missing_refresh <- current
  missing_refresh$hessian_refresh <- NULL
  testthat::expect_false(configuration_matches_v2(missing_refresh, current))
})

testthat::test_that("the recorded refresh interval follows the environment and the option", {
  rec <- unset_solver_options(withr::with_envvar(
    c(JPDMEQR_HESSIAN_REFRESH = "17"), solver_configuration_v2()
  ))
  testthat::expect_identical(rec$hessian_refresh, 17L)

  rec <- unset_solver_env(withr::with_options(
    list(jpDMEQR.hessian_refresh = 9L), solver_configuration_v2()
  ))
  testthat::expect_identical(rec$hessian_refresh, 9L)
})

testthat::test_that("the recorded values survive a JSON round trip", {
  rec <- unset_solver_env(unset_solver_options(solver_configuration_v2()))
  back <- jsonlite::fromJSON(
    jsonlite::toJSON(rec, auto_unbox = TRUE, null = "null", digits = NA),
    simplifyVector = TRUE
  )
  testthat::expect_true(configuration_matches_v2(back, rec))
})
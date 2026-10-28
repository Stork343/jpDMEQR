root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
sys.source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"),
           envir = environment())

actual_geometry_record <- function() {
  path <- Sys.getenv("JPDMEQR_UNIT_GEOMETRY_MANIFEST", unset = file.path(
    root, "results", "scratch", "profile_pair_geometry_development_20261022",
    "production_validation_20261022_02", "manifest.json"
  ))
  if (!file.exists(path)) stop("The completed actual production geometry record is required.")
  read_json_v2(path)
}

testthat::test_that("the actual completed geometry record is accepted at its recorded identity", {
  actual <- actual_geometry_record()
  decision <- validate_geometry_record_v2(
    actual, actual$commit_sha, actual$solver_configuration, final = TRUE,
    method_version = "profile_pair_20261022"
  )
  testthat::expect_true(decision$pass)
  testthat::expect_identical(decision$reason, "accepted")
})

testthat::test_that("the geometry gate requires the actual assembly comparison", {
  actual <- actual_geometry_record()
  altered <- actual
  altered$assembly_modes_checked <- FALSE
  testthat::expect_false(validate_geometry_record_v2(
    altered, actual$commit_sha, actual$solver_configuration)$pass)
  altered <- actual
  altered$assembly_pass <- FALSE
  testthat::expect_false(validate_geometry_record_v2(
    altered, actual$commit_sha, actual$solver_configuration)$pass)
  testthat::expect_true(validate_geometry_record_v2(
    actual, actual$commit_sha, actual$solver_configuration)$pass)
})

testthat::test_that("the geometry record checks its Dantzig result and expected identity", {
  actual <- actual_geometry_record()
  altered <- actual
  altered$dantzig_pass <- FALSE
  testthat::expect_false(validate_geometry_record_v2(
    altered, actual$commit_sha, actual$solver_configuration)$pass)
  testthat::expect_false(validate_geometry_record_v2(
    actual, paste0(actual$commit_sha, "-changed"), actual$solver_configuration)$pass)
  changed_configuration <- actual$solver_configuration
  changed_configuration$assembly <- "sequential"
  testthat::expect_false(validate_geometry_record_v2(
    actual, actual$commit_sha, changed_configuration)$pass)
})

testthat::test_that("the file wrapper always checks the actual current commit", {
  path <- geometry_manifest_path_v2(root, "profile_pair_20261022")
  decision <- validate_geometry_gate_v2(root, final = TRUE, method_version = "profile_pair_20261022")
  if (file.exists(path)) {
    expected <- validate_geometry_record_v2(
      read_json_v2(path), current_commit_v2(root), solver_configuration_v2(), TRUE,
      "profile_pair_20261022"
    )
    testthat::expect_identical(decision, expected)
  } else {
    testthat::expect_false(decision$pass)
    testthat::expect_identical(decision$reason, "strict geometry manifest missing")
  }
})

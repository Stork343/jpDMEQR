root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
sys.source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"),
           envir = environment())

# 结构性审计资产验收规则的纯判定测试。对象字段与 build_target_for_config 在
# role = "structural_target_audit" 下写出的资产一致，只需构造字段，不进行拟合。
make_audit_fixture <- function(target_displacement = c(x00001 = 1e-3, x00002 = -2e-3),
                               target_mc_se_by_coordinate = c(x00001 = 4e-3, x00002 = 4e-3),
                               n_population = 100000L,
                               repeats = 4L,
                               all_converged = TRUE,
                               max_kkt_residual = 1e-7,
                               max_nuisance_gradient = 1e-9,
                               dependency_hash = "target-construction-v2-smoothing-distance",
                               implementation_commit = "fixture-commit") {
  list(
    experiment_id = "E03",
    target_displacement = target_displacement,
    target_mc_se_by_coordinate = target_mc_se_by_coordinate,
    n_population = n_population,
    repeats = repeats,
    all_converged = all_converged,
    max_kkt_residual = max_kkt_residual,
    max_nuisance_gradient = max_nuisance_gradient,
    smoothing_target_distance = 1e-3,
    dependency_hash = dependency_hash,
    implementation_commit = implementation_commit
  )
}

testthat::test_that("structural audit validator accepts an audit that meets the frozen criteria", {
  out <- validate_structural_target_audit_v2(
    make_audit_fixture(),
    expected_dependency_hash = "target-construction-v2-smoothing-distance"
  )
  testthat::expect_true(out$valid)
  testthat::expect_length(out$problems, 0L)
})

testthat::test_that("structural audit validator accepts the KKT plateau branch", {
  out <- validate_structural_target_audit_v2(
    make_audit_fixture(all_converged = FALSE, max_kkt_residual = 1e-6,
                       max_nuisance_gradient = 1e-8)
  )
  testthat::expect_true(out$valid)
})

testthat::test_that("structural audit validator rejects fits above the acceptance tolerances", {
  out <- validate_structural_target_audit_v2(
    make_audit_fixture(all_converged = FALSE, max_kkt_residual = 1e-4,
                       max_nuisance_gradient = 1e-8)
  )
  testthat::expect_false(out$valid)
  testthat::expect_true(any(grepl("did not converge", out$problems)))
})

testthat::test_that("structural audit validator rejects a short repeat count and a small population", {
  short <- validate_structural_target_audit_v2(make_audit_fixture(repeats = 3L))
  testthat::expect_false(short$valid)
  testthat::expect_true(any(grepl("repeat count is below four", short$problems)))

  small <- validate_structural_target_audit_v2(make_audit_fixture(n_population = 50000L))
  testthat::expect_false(small$valid)
  testthat::expect_true(any(grepl("below 100000", small$problems)))
})

testthat::test_that("structural audit validator rejects a stale dependency hash", {
  out <- validate_structural_target_audit_v2(
    make_audit_fixture(),
    expected_dependency_hash = "target-construction-v2-other"
  )
  testthat::expect_false(out$valid)
  testthat::expect_true(any(grepl("dependency hash is stale", out$problems)))
})

testthat::test_that("structural audit validator rejects a displacement outside the frozen tolerance", {
  out <- validate_structural_target_audit_v2(
    make_audit_fixture(target_displacement = c(x00001 = 0.05, x00002 = -0.06),
                       target_mc_se_by_coordinate = c(x00001 = 1e-3, x00002 = 1e-3))
  )
  testthat::expect_false(out$valid)
  testthat::expect_true(any(grepl("not compatible with the frozen structural target tolerance",
                                  out$problems)))
})
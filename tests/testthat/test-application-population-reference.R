root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "application_v2", envir = .GlobalEnv)

testthat::test_that("独立读取的真实总体参考具有相同统计定义", {
  path <- file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1", "matched_population_model_A_q1.rds")
  testthat::skip_if_not(file.exists(path), "固定总体参考未提供。")
  fields <- c("X", "Z", "rows", "beta0", "preprocessing")
  first <- readRDS(path)[fields]
  second <- readRDS(path)[fields]
  testthat::expect_true(application_population_reference_equal_v2(first, second))
  asset <- file.path(root, "results", "raw", "application_matched",
                      "matched_targets_20261022_v1_AM_A_q1_normal", "matched_target_asset.rds")
  testthat::skip_if_not(file.exists(asset), "完整100000 by four总体目标尚未提供。")
  testthat::expect_true(application_population_reference_equal_v2(readRDS(asset)$population_reference, first))
})

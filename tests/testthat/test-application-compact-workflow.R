compact_root <- v2_require_root()
source(file.path(compact_root, "scripts", "00_source_v2.R"), local = TRUE)
source_v2_module(compact_root, "profile_v2", environment())
source_v2_module(compact_root, "application_v2", environment())

testthat::test_that("真实129 subject workflow的compact文件重建相同总体单位表格", {
  registry <- read.csv(file.path(compact_root, "data", "derived", "GSE65391",
                                 "execution_20261022_v1", "matched_registry.csv"), stringsAsFactors = FALSE)
  cases <- data.frame(
    case_id = c("AM_A_q1_normal_tau25", "AM_B_q2_t3_tau75"),
    complete = c("workflow_newton_A_q1_normal_20261022_v3", "workflow_newton_B_q2_t3_20261022_v2"),
    compact = c("workflow_compact_A_q1_normal_20261022_v1", "workflow_compact_B_q2_t3_20261022_v1"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(cases))) {
    case <- registry[registry$case_id == cases$case_id[i], , drop = FALSE]
    full_path <- file.path(compact_root, "results", "scratch", cases$complete[i], "workflow_result.rds")
    compact_path <- file.path(compact_root, "results", "scratch", cases$compact[i], "workflow_result.rds")
    complete <- readRDS(full_path)
    replayed <- read_compact_application_workflow_v2(compact_path)
    population <- readRDS(file.path(compact_root, "data", "derived", "GSE65391", "execution_20261022_v1",
                                     sprintf("matched_population_model_%s_q%d.rds", case$model, case$q)))
    asset <- readRDS(file.path(compact_root, "results", "raw", "application_matched",
                               sprintf("matched_targets_20261022_v1_AM_%s_q%d_%s", case$model, case$q, case$error_family),
                               "matched_target_asset.rds"))
    target <- asset$targets[asset$targets$tau == case$tau, ]
    original <- application_matched_coordinate_metrics_v2(complete, population, target, case, 1L)
    rebuilt <- application_matched_coordinate_metrics_v2(replayed, population, target, case, 1L)
    testthat::expect_identical(replayed$table, complete$table)
    testthat::expect_identical(rebuilt$table, original$table)
    testthat::expect_identical(rebuilt$mapped_jackknife, original$mapped_jackknife)
    testthat::expect_identical(replayed$fit$components$cluster_scores, complete$fit$components$cluster_scores)
    testthat::expect_identical(replayed$fit$components$hessian, complete$fit$components$hessian)
    testthat::expect_identical(names(replayed$jackknife$leave_one), names(complete$jackknife$leave_one))
    testthat::expect_identical(replayed$jackknife$estimates, complete$jackknife$estimates)
    testthat::expect_equal(nrow(replayed$jackknife$diagnostics), 129L)
  }
})

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "application_v2", envir = .GlobalEnv)

testthat::test_that("预测用初始L1与paired方法的相同两遍loadings结果一致", {
  path <- file.path(root, "data", "interim", "GSE65391", "expression_processed_probe.rds")
  testthat::skip_if_not(file.exists(path), "本地GSE65391表达未提供。")
  rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
  members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
  signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
  annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
  expression <- readRDS(path)
  md <- read.csv(file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv"), stringsAsFactors = FALSE)
  audit <- audit_primary_outcome_v2(md)
  eligible <- md[audit$eligible_visit, ]
  subjects <- unlist(lapply(split(eligible$subject_id, eligible$batch), function(ids) {
    utils::head(sort(unique(ids)), 6L)
  }), use.names = FALSE)
  md <- eligible[eligible$subject_id %in% subjects, , drop = FALSE]
  rownames(md) <- md$sample_id
  design <- prepare_gse65391_lowdim_design_v2(
    md, expression[, rownames(md), drop = FALSE], annotation, signatures,
    rep(TRUE, nrow(md)), rule, module_score_members = members
  )
  X <- design$X[, c("(Intercept)", "module__M1.2", "module__M3.4"), drop = FALSE]
  initial <- fit_application_profile_initial_v2(log(md$c3), X, design$Z, md$subject_id, 0.5,
                                                penalty_factor = c(0, 1, 1))
  paired <- fit_profile_pair_v2(log(md$c3), X, design$Z, md$subject_id, 0.5,
                                 h_inf = 12^(-3 / 10), penalty_factor = c(0, 1, 1))
  testthat::expect_identical(initial$beta_l1, paired$beta_l1)
  testthat::expect_identical(initial$loadings, paired$loadings)
  testthat::expect_identical(initial$lambda_stochastic, paired$lambda_stochastic)
  testthat::expect_identical(initial$lambda_mean, paired$lambda_mean)
  testthat::expect_identical(initial$coordinate_penalty, paired$lambda_components$coordinate)
  testthat::expect_identical(initial$h_est, paired$h_est)
  testthat::expect_true(initial$converged)
  testthat::expect_true(initial$kkt_normalized <= 1e-3)
})

testthat::test_that("冻结的subject folds完整包含每名患者并保持角色分离", {
  path <- file.path(root, "data", "derived", "GSE65391", "subject_design_20261022_v1", "analysis_inputs.rds")
  testthat::skip_if_not(file.exists(path), "本地冻结subject设计未提供。")
  obj <- readRDS(path)
  summary <- validate_gse65391_subject_folds_v2(obj$eligible_metadata, obj$subject_folds)
  testthat::expect_equal(nrow(summary), 25L)
  testthat::expect_equal(nrow(obj$subject_folds), 645L)
  testthat::expect_true(all(summary$overlap_subjects == 0))
  testthat::expect_true(all(summary$training_subjects + summary$test_subjects == 129L))
  testthat::expect_true(all(summary$training_visits + summary$test_visits == 867L))
})

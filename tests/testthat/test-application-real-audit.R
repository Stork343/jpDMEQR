root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "application_v2", envir = .GlobalEnv)

testthat::test_that("真实GSE记录的标识检查保留观测字段", {
  metadata_path <- file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv")
  testthat::skip_if_not(file.exists(metadata_path), "本地GSE65391 metadata未提供。")
  md <- read.csv(metadata_path, stringsAsFactors = FALSE)
  before <- md
  identifiers <- audit_gse65391_identifiers_v2(md)
  testthat::expect_identical(md, before)
  testthat::expect_equal(identifiers$summary$critical_conflicts, 0)
  canonical <- identifiers$canonical_keys
  testthat::expect_false(anyDuplicated(canonical$canonical_visit_key[canonical$biological]) > 0)
  missing_visit <- canonical$biological & !is.finite(canonical$recorded_visit)
  testthat::expect_true(any(missing_visit))
  testthat::expect_true(all(canonical$key_source[missing_visit] == "subject_time"))
  testthat::expect_true(all(is.finite(canonical$recorded_cumulative_time[missing_visit])))
  relations <- audit_gse65391_replicate_relations_v2(md)
  flagged_sle <- relations$records$group_type == "SLE" & relations$records$batch_replicate_flag
  testthat::expect_true(any(flagged_sle))
  testthat::expect_true(all(relations$records$matching_subject_visit_records[flagged_sle] == 0))
  testthat::expect_true(all(relations$records$relation_status[flagged_sle] ==
                             "sole_observed_biological_record"))
})

testthat::test_that("真实C3审计生成诊断并阻止未冻结的曲率规则", {
  metadata_path <- file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv")
  testthat::skip_if_not(file.exists(metadata_path), "本地GSE65391 metadata未提供。")
  md <- read.csv(metadata_path, stringsAsFactors = FALSE)
  audit <- audit_primary_outcome_v2(md)
  testthat::expect_true(audit$summary$scalar_gates_pass)
  testthat::expect_false(audit$summary$all_gates_pass)
  testthat::expect_true(audit$summary$review_required)
  testthat::expect_true(all(audit$batch$summary$pass))
  testthat::expect_true(all(audit$curvature$diagnostics$status == "review_required"))
  testthat::expect_true(all(is.finite(audit$curvature$diagnostics$cluster_average_curvature)))
  testthat::expect_true(all(audit$curvature$diagnostics$local_subjects <= audit$summary$final_subjects))
  eligible <- md[audit$eligible_visit, , drop = FALSE]
  for (include_sledai in c(FALSE, TRUE)) {
    clinical <- prepare_clinical_design_split_v2(
      eligible, rep(TRUE, nrow(eligible)), include_sledai = include_sledai
    )
    testthat::expect_equal(qr(clinical$X)$rank, ncol(clinical$X))
    testthat::expect_lte(ncol(clinical$X), floor(audit$summary$final_subjects /
                                                 log(audit$summary$final_subjects)))
  }
  testthat::expect_error(validate_gse65391_fit_inputs_v2(
    audit, file.path(root, "config", "module_signatures.csv")
  ), "EMPIRICAL_GATE_INCOMPLETE")
  testthat::expect_error(read_frozen_module_definitions_v2(
    file.path(root, "config", "module_signatures_template.csv")
  ), "MODULE_GATE_FAILED.*nonempty")
  previous_path <- file.path(root, "data", "derived", "GSE65391", "outcome_audit_summary.csv")
  testthat::skip_if_not(file.exists(previous_path), "原队列摘要未提供。")
  previous <- read.csv(previous_path)
  testthat::expect_equal(audit$summary$final_subjects, previous$final_subjects)
  testthat::expect_equal(audit$summary$final_visits, previous$final_visits)
})

testthat::test_that("真实临床小样本完成所有subject的profile jackknife", {
  metadata_path <- file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv")
  testthat::skip_if_not(file.exists(metadata_path), "本地GSE65391 metadata未提供。")
  md <- read.csv(metadata_path, stringsAsFactors = FALSE)
  audit <- audit_primary_outcome_v2(md)
  subjects <- sort(unique(md$subject_id[audit$eligible_visit]))[seq_len(12L)]
  selected <- audit$eligible_visit & md$subject_id %in% subjects
  clinical <- md[selected, , drop = FALSE]
  X <- cbind(intercept = 1, centered_time = as.numeric(scale(clinical$cumulative_time)))
  Z <- matrix(1, nrow(X), 1L)
  result <- lowdim_profile_inference_v2(
    log(clinical$c3), X, Z, clinical$subject_id,
    tau = 0.50, h = 12^(-3 / 10), score_source = "smoothed"
  )
  testthat::expect_equal(result$table$jackknife_completed, rep(12L, ncol(X)))
  testthat::expect_setequal(result$jackknife$diagnostics$omitted_subject, subjects)
  testthat::expect_true(all(result$jackknife$diagnostics$converged))
  testthat::expect_true(all(result$table$sandwich_se > 0))
  testthat::expect_true(all(result$table$jackknife_se > 0))
  testthat::expect_true(all(result$table$target_type == "finite_h_regularised_profile"))
  testthat::expect_true(all(result$table$score_source == "smoothed"))
  testthat::expect_error(lowdim_profile_inference_v2(
    log(clinical$c3), X, Z, clinical$subject_id, 0.50, 12^(-3 / 10),
    score_source = "unsmoothed"
  ), "SANDWICH_RULE_REQUIRED")
})

testthat::test_that("完整subject bootstrap共同计算quantile与kernel均值", {
  metadata_path <- file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv")
  testthat::skip_if_not(file.exists(metadata_path), "本地GSE65391 metadata未提供。")
  md <- read.csv(metadata_path, stringsAsFactors = FALSE)
  rule <- yaml::read_yaml(file.path(root, "config", "application_curvature_20261022.yml"))
  audit <- audit_primary_outcome_v2(md, curvature_rule = rule)
  diagnostic <- audit$curvature$diagnostics
  testthat::expect_true(all(diagnostic$pass[diagnostic$primary]))
  testthat::expect_equal(dim(audit$curvature$bootstrap_means), c(10000L, 12L))
  testthat::expect_equal(unique(diagnostic$minimum_subjects), 27)
  testthat::expect_equal(unique(diagnostic$family_size[diagnostic$primary]), 9)
  testthat::expect_equal(unique(diagnostic$family_size[!diagnostic$primary]), 3)

  eligible <- md[audit$eligible_visit, , drop = FALSE]
  by_subject <- split(log(eligible$c3), factor(eligible$subject_id,
                                             levels = unique(eligible$subject_id)))
  set.seed(rule$bootstrap_seed)
  selected <- sample.int(length(by_subject), length(by_subject), replace = TRUE)
  replicated_subjects <- by_subject[selected]
  for (j in seq_len(nrow(diagnostic))) {
    location <- as.numeric(stats::quantile(unlist(replicated_subjects), diagnostic$tau[j],
                                           type = rule$quantile_type, names = FALSE))
    direct_mean <- mean(vapply(replicated_subjects, function(y) {
      mean(phi_smooth_v2(y - location, diagnostic$tau[j], diagnostic$h[j]))
    }, numeric(1)))
    testthat::expect_equal(audit$curvature$bootstrap_means[1, j], direct_mean, tolerance = 1e-12)
  }
})

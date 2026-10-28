root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "application_v2", envir = .GlobalEnv)

testthat::test_that("真实全流程jackknife重新计算每次clinical及全部module preprocessing", {
  expression_path <- file.path(root, "data", "interim", "GSE65391", "expression_processed_probe.rds")
  testthat::skip_if_not(file.exists(expression_path), "本地GSE65391表达未提供。")
  rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
  members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
  signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
  annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
  expression <- readRDS(expression_path)
  metadata <- read.csv(file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv"),
                       stringsAsFactors = FALSE)
  audit <- audit_primary_outcome_v2(metadata)
  eligible <- metadata[audit$eligible_visit, , drop = FALSE]
  subjects <- unlist(lapply(split(eligible$subject_id, eligible$batch), function(ids) {
    utils::head(sort(unique(ids)), 6L)
  }), use.names = FALSE)
  testthat::expect_equal(length(unique(subjects)), 12L)
  metadata <- metadata[audit$eligible_visit & metadata$subject_id %in% subjects, , drop = FALSE]
  rownames(metadata) <- metadata$sample_id
  expression <- expression[, rownames(metadata), drop = FALSE]
  coordinates <- c("(Intercept)", "module__M1.2")
  h <- 12^(-3 / 10)
  result <- lowdim_profile_pair_workflow_inference_v2(
    metadata, expression, annotation, signatures, members, rule,
    tau = 0.5, h = h, coordinate_subset = coordinates
  )
  testthat::expect_identical(colnames(result$jackknife$estimates), coordinates)
  testthat::expect_equal(nrow(result$jackknife$estimates), 12L)
  testthat::expect_true(all(result$table$jackknife_scope == "complete_sample_standardised_workflow"))
  testthat::expect_true(all(result$table$coordinate_scope == "declared_research_coordinate_subset"))
  testthat::expect_true(all(result$jackknife$diagnostics$converged))
  testthat::expect_true(all(result$jackknife$diagnostics$h == h))
  testthat::expect_true(all(result$jackknife$diagnostics$h_half == h / 2))
  testthat::expect_true(all(result$table$jackknife_se > 0))
  for (subject in names(result$jackknife$leave_one)) {
    stage <- result$jackknife$leave_one[[subject]]
    testthat::expect_false(subject %in% stage$preprocessing$subject_ids)
    testthat::expect_setequal(stage$preprocessing$module_audit$module_id, rule$module_order)
    testthat::expect_true(all(stage$preprocessing$module_audit$eligible))
    testthat::expect_identical(stage$preprocessing$gene_preprocessing$training_sample_ids,
                               stage$preprocessing$sample_ids)
  }
  subject <- names(result$jackknife$leave_one)[1]
  keep <- as.character(metadata$subject_id) != subject
  direct_design <- prepare_gse65391_lowdim_design_v2(
    metadata[keep, , drop = FALSE], expression[, keep, drop = FALSE], annotation,
    signatures, rep(TRUE, sum(keep)), rule, module_score_members = members
  )
  saved <- result$jackknife$leave_one[[subject]]
  testthat::expect_identical(saved$preprocessing$gene_preprocessing, direct_design$gene_preprocessing)
  testthat::expect_identical(saved$preprocessing$clinical$time_center, direct_design$clinical$time_center)
  testthat::expect_true(any(abs(saved$preprocessing$gene_preprocessing$center -
                                  result$preprocessing$gene_preprocessing$center) > 1e-12))
  testthat::expect_true(any(abs(saved$preprocessing$gene_preprocessing$scale -
                                  result$preprocessing$gene_preprocessing$scale) > 1e-12))
  direct_fit <- fit_lowdim_profile_pair_v2(
    log(metadata$c3[keep]), direct_design$X[, coordinates, drop = FALSE], direct_design$Z,
    metadata$subject_id[keep], tau = 0.5, h = h
  )
  testthat::expect_equal(unname(result$jackknife$estimates[subject, ]),
                         unname(direct_fit$fit$beta_pair), tolerance = 1e-12)
  centered <- sweep(result$jackknife$estimates, 2, colMeans(result$jackknife$estimates))
  testthat::expect_equal(unname(result$table$jackknife_se),
                         unname(sqrt(11 / 12 * colSums(centered^2))), tolerance = 1e-12)
})

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "application_v2", envir = .GlobalEnv)

testthat::test_that("matched模拟保留真实源记录并明确生成结果与总体单位", {
  directory <- file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1")
  path <- file.path(directory, "matched_population_model_B_q2.rds")
  testthat::skip_if_not(file.exists(path), "固定总体设计未提供。")
  population <- readRDS(path)
  before <- population$metadata
  config <- yaml::read_yaml(file.path(root, "config", "application_matched_prospective_20261022.yml"))
  expression <- readRDS(file.path(root, "data", "derived", "GSE65391", "module_measurements_20261022_v1",
                                    "observed_module_probe_expression.rds"))
  draw <- resample_application_matched_workflow_v2(population, expression, config, "normal", 314159L)
  repeat_draw <- resample_application_matched_workflow_v2(population, expression, config, "normal", 314159L)
  testthat::expect_identical(draw$generated$y, repeat_draw$generated$y)
  testthat::expect_identical(population$metadata, before)
  testthat::expect_equal(length(unique(draw$metadata$subject_id)), 129L)
  testthat::expect_false(anyDuplicated(rownames(draw$metadata)) > 0)
  testthat::expect_identical(colnames(draw$expression), rownames(draw$metadata))
  testthat::expect_identical(draw$metadata$observed_c3, before$c3[draw$generated$source_rows])
  testthat::expect_true(all(draw$metadata$outcome_source == "declared_matched_simulation"))
  dat <- draw$generated
  reconstructed <- as.numeric(dat$X %*% population$beta0) +
    rowSums(dat$Z * dat$random_effects[dat$cluster_id, , drop = FALSE]) + dat$errors
  testthat::expect_equal(dat$y, reconstructed, tolerance = 1e-12)
  identity <- application_population_unit_map_v2(population$preprocessing, population$preprocessing,
                                                 names(population$beta0))
  testthat::expect_equal(unname(identity$matrix), diag(17L))
  testthat::expect_true(identity$exact_affine_measurement_map)
  rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
  members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
  signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
  annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
  design <- prepare_gse65391_lowdim_design_v2(
    draw$metadata, draw$expression, annotation, signatures, rep(TRUE, nrow(draw$metadata)),
    rule, include_sledai = TRUE, module_score_members = members
  )
  preprocessing <- list(clinical = design$clinical, gene_preprocessing = design$gene_preprocessing)
  map <- application_population_unit_map_v2(preprocessing, population$preprocessing, names(population$beta0))
  population_units <- design$X
  for (coordinate in grep("^module__", colnames(design$X), value = TRUE)) {
    id <- sub("^module__", "", coordinate)
    raw <- design$X[, coordinate] * preprocessing$gene_preprocessing$scale[id] +
      preprocessing$gene_preprocessing$center[id]
    population_units[, coordinate] <- (raw - population$preprocessing$gene_preprocessing$center[id]) /
      population$preprocessing$gene_preprocessing$scale[id]
  }
  population_units[, "centered_time"] <- draw$metadata$cumulative_time -
    population$preprocessing$clinical$time_center
  beta_population_units <- as.numeric(map$matrix %*% population$beta0)
  testthat::expect_equal(as.numeric(design$X %*% population$beta0),
                         as.numeric(population_units %*% beta_population_units), tolerance = 1e-12)
})

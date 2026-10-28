root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "application_v2", envir = .GlobalEnv)

testthat::test_that("真实固定研究坐标的comparator调用保持官方参数和临床零penalty", {
  path <- file.path(root, "data", "interim", "GSE65391", "expression_processed_probe.rds")
  testthat::skip_if_not(file.exists(path), "本地GSE65391表达未提供。")
  md <- read.csv(file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv"), stringsAsFactors = FALSE)
  audit <- audit_primary_outcome_v2(md)
  eligible <- md[audit$eligible_visit, ]
  subjects <- unlist(lapply(split(eligible$subject_id, eligible$batch), function(ids) {
    utils::head(sort(unique(ids)), 6L)
  }), use.names = FALSE)
  md <- eligible[eligible$subject_id %in% subjects, , drop = FALSE]
  rownames(md) <- md$sample_id
  md$log_c3 <- log(md$c3)
  annotation <- read.csv(file.path(root, "data", "interim", "GSE65391", "probe_annotation_20261022.csv"),
                         stringsAsFactors = FALSE)
  expression <- readRDS(path)
  rows <- which(annotation$gene_symbol == "DDX58" & annotation$mapping_status == "unique_gene")
  rows <- rows[!is.na(rows)]
  probes <- expression[annotation$probe_id[rows], rownames(md), drop = FALSE]
  selected <- which.max(apply(probes, 1, stats::mad))
  gene <- as.numeric(scale(probes[selected, ]))
  design <- list(X = cbind("(Intercept)" = 1, centered_time = md$cumulative_time - mean(md$cumulative_time),
                            gene__DDX58 = gene), penalty_factor = c(0, 0, 1))
  train <- rep(TRUE, nrow(md))
  order <- order(as.character(md$subject_id), md$cumulative_time)
  X <- design$X[order, , drop = FALSE]
  y <- md$log_c3[order]
  ids <- as.character(md$subject_id[order])
  pooled <- fit_application_comparator_v2("POOLED-QR-LASSO", md, design, train, .5, seed = 8128L)
  set.seed(8128L)
  penalty <- quantreg::LassoLambdaHat(X, tau = .5)
  penalty[1:2] <- 0
  direct <- quantreg::rq.fit.lasso(X, y, tau = .5, lambda = penalty)
  testthat::expect_equal(pooled$beta_prediction, direct$coefficients, tolerance = 1e-12)
  testthat::expect_identical(pooled$lambda, penalty)
  testthat::expect_equal(unname(pooled$lambda[1:2]), c(0, 0))
  qgee <- fit_application_comparator_v2("QGEE-SCAD", md, design, train, .5, seed = 8128L)
  set.seed(8128L)
  attr(X, "assign") <- c(0L, 1L, 2L)
  kernel <- get(".fit_single_lambda", envir = asNamespace("geeVerse"))
  direct_qgee <- kernel(x = X, y = y, nobs = as.integer(table(factor(ids, levels = unique(ids)))),
                         tau = .5, corstr = "exchangeable", lambda = qgee$lambda,
                         penalty_factor = c(0, 1), intercept = TRUE, betaint = NULL,
                         control = geeVerse::qpgeeControl())
  testthat::expect_equal(qgee$beta_prediction, direct_qgee$coefficients, tolerance = 1e-12)
  testthat::expect_identical(qgee$converged, direct_qgee$converged)
  testthat::expect_equal(nrow(qgee$lambda_path), 30L)
  testthat::expect_identical(qgee$lambda, qgee$lambda_path$lambda[which.min(qgee$lambda_path$hbic)])
  testthat::expect_identical(qgee$penalty_factor_excluding_intercept, c(0, 1))
  lqmm <- fit_application_comparator_v2("LQMM", md, design, train, .5, seed = 8128L)
  fixed_names <- sprintf("covariate_%03d", seq_len(ncol(X)))
  data <- data.frame(response = y, subject = factor(ids), X, check.names = FALSE)
  names(data)[seq_len(ncol(X)) + 2L] <- fixed_names
  set.seed(8128L)
  direct_lqmm <- lqmm::lqmm(fixed = stats::reformulate(fixed_names, response = "response", intercept = FALSE),
                            random = ~1, group = "subject", tau = .5, data = data,
                            covariance = "pdSymm", nK = 7L,
                            control = lqmm::lqmmControl(LP_tol_ll = 1e-6))
  testthat::expect_equal(unname(lqmm$beta_prediction), as.numeric(stats::coef(direct_lqmm)), tolerance = 1e-10)
  testthat::expect_true(all(is.finite(lqmm$beta_prediction)))
})

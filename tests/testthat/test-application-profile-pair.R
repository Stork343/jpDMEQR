root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "application_v2", envir = .GlobalEnv)

testthat::test_that("真实C3数据的paired区间使用同一cluster score并完成subject jackknife", {
  path <- file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv")
  testthat::skip_if_not(file.exists(path), "本地GSE65391 metadata未提供。")
  metadata <- read.csv(path, stringsAsFactors = FALSE)
  audit <- audit_primary_outcome_v2(metadata)
  subjects <- sort(unique(metadata$subject_id[audit$eligible_visit]))[seq_len(12L)]
  selected <- audit$eligible_visit & metadata$subject_id %in% subjects
  data <- metadata[selected, , drop = FALSE]
  X <- cbind(intercept = 1, centered_time = as.numeric(scale(data$cumulative_time)))
  Z <- matrix(1, nrow(X), 1L)
  h <- 12^(-3 / 10)
  result <- lowdim_profile_pair_inference_v2(
    log(data$c3), X, Z, data$subject_id, tau = 0.5, h = h
  )
  testthat::expect_true(all(result$table$target_type == "unsmoothed_regularised_profile"))
  testthat::expect_true(all(result$table$score_definition ==
                             "profile_bandwidth_pair_2_half_minus_full"))
  testthat::expect_equal(result$table$jackknife_completed, rep(12L, ncol(X)))
  testthat::expect_setequal(result$jackknife$diagnostics$omitted_subject, subjects)
  testthat::expect_true(all(result$jackknife$diagnostics$converged))
  testthat::expect_true(all(result$table$sandwich_se > 0))
  testthat::expect_true(all(result$table$jackknife_se > 0))
  G <- result$fit$components$cluster_scores
  centered <- sweep(G, 2, colMeans(G))
  projected <- centered %*% result$precision
  testthat::expect_equal(result$table$sandwich_se,
                         unname(sqrt(colSums(projected^2) / 12^2)), tolerance = 1e-12)
  direct <- result$fit$beta_refit -
    as.numeric(crossprod(result$precision, colMeans(G)))
  testthat::expect_equal(result$table$beta_hat, unname(direct), tolerance = 1e-12)
  first <- rownames(result$jackknife$estimates)[1]
  keep <- data$subject_id != first
  leave <- fit_lowdim_profile_pair_v2(
    log(data$c3[keep]), X[keep, , drop = FALSE], Z[keep, , drop = FALSE],
    data$subject_id[keep], tau = 0.5, h = h
  )
  testthat::expect_equal(unname(result$jackknife$estimates[1, ]),
                         unname(leave$fit$beta_pair), tolerance = 1e-12)
})

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "application_v2", envir = .GlobalEnv)

make_exploratory_dgp_v2 <- function() {
  n_clusters <- 400L
  dat <- generate_profile_qr_data_v2(n_clusters = n_clusters, p = 8L, s = 2L, q = 1L,
                                      tau = .5, m_values = 3L, seed = 22610L)
  subject <- as.character(dat$cluster_id)
  index <- match(subject, unique(subject))
  set.seed(22611L)
  metadata <- data.frame(subject_id = subject,
                         cumulative_time = ave(seq_along(subject), subject, FUN = seq_along),
                         age = stats::rnorm(length(subject)),
                         gender = rep(c("F", "M"), length.out = n_clusters)[index],
                         race = rep(c("R1", "R2"), each = 3L, length.out = n_clusters)[index],
                         treatment = rep(c("T1", "T2"), each = 5L, length.out = n_clusters)[index],
                         batch = rep(c("B1", "B2"), each = 7L, length.out = n_clusters)[index],
                         sledai = stats::rnorm(length(subject)), log_c3 = dat$y,
                         stringsAsFactors = FALSE)
  metadata$c3 <- exp(metadata$log_c3)
  rownames(metadata) <- paste0("sample", seq_len(nrow(metadata)))
  expression <- t(dat$X)
  rownames(expression) <- paste0("probe", seq_len(nrow(expression)))
  colnames(expression) <- rownames(metadata)
  annotation <- data.frame(probe_id = rownames(expression), gene_symbol = paste0("G", seq_len(nrow(expression))),
                           is_control = FALSE, is_multimapping = FALSE, stringsAsFactors = FALSE)
  original <- data.frame(subject_id = unique(subject),
                          role = ifelse(seq_len(n_clusters) <= n_clusters / 2L, "selection", "inference"),
                          stringsAsFactors = FALSE)
  list(metadata = metadata, expression = expression, annotation = annotation, assignment = original)
}

testthat::test_that("exploratory preprocessing uses the frozen selection subjects", {
  dat <- make_exploratory_dgp_v2()
  a <- prepare_exploratory_gene_split_v2(dat$metadata, dat$expression, dat$annotation, dat$assignment, 8L)
  changed <- dat$metadata
  changed$log_c3[a$inference_rows] <- changed$log_c3[a$inference_rows] + 2
  changed$c3 <- exp(changed$log_c3)
  b <- prepare_exploratory_gene_split_v2(changed, dat$expression, dat$annotation, dat$assignment, 8L)
  testthat::expect_identical(a$genes, b$genes)
  testthat::expect_identical(a$preprocessing, b$preprocessing)
  testthat::expect_identical(a$clinical$A$X, b$clinical$A$X)
  testthat::expect_identical(a$clinical$B$X, b$clinical$B$X)
  altered <- dat$expression
  altered[, a$inference_rows] <- altered[, a$inference_rows] + 100
  c <- prepare_exploratory_gene_split_v2(dat$metadata, altered, dat$annotation, dat$assignment, 8L)
  testthat::expect_identical(a$preprocessing, c$preprocessing)
  testthat::expect_identical(a$genes[a$selection_rows, ], c$genes[c$selection_rows, ])
  testthat::expect_length(intersect(a$selection_subjects, a$inference_subjects), 0L)
  testthat::expect_error(prepare_exploratory_gene_split_v2(dat$metadata, dat$expression[, -1L],
                                                         dat$annotation, dat$assignment), "must agree")
})

testthat::test_that("exploratory real DGP fit preserves independent paired intervals and ranking", {
  dat <- make_exploratory_dgp_v2()
  prepared <- prepare_exploratory_gene_split_v2(dat$metadata, dat$expression, dat$annotation, dat$assignment, 8L)
  output <- file.path(root, "results", "scratch", "application_exploratory_tests_20261022")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  checkpoint <- tempfile(pattern = "actual_dgp_", tmpdir = output)
  result <- fit_exploratory_gene_role_v2(dat$metadata, prepared, "A", .5, checkpoint,
                                         source_inputs = list(source = "generated_profile_qr_DGP_seed22610"))
  testthat::expect_identical(result$status, "ok")
  testthat::expect_equal(nrow(result$table), 3L)
  testthat::expect_equal(sum(result$ranking$selected), 3L)
  testthat::expect_true(all(is.finite(result$table$se) & result$table$se > 0))
  selected <- readRDS(file.path(checkpoint, "selected_design.rds"))
  fitted <- readRDS(file.path(checkpoint, "inference_fit.rds"))
  testthat::expect_identical(selected$cluster_inference, dat$metadata$subject_id[prepared$inference_rows])
  testthat::expect_identical(colnames(selected$X_inference), names(fitted$fit$beta_refit))
  direct <- profile_bandwidth_pair_wald_v2(fitted$fit$components, fitted$fit$beta_refit,
                                            solve(fitted$fit$components$hessian), seq_len(ncol(selected$X_inference)))
  index <- match(result$table$coordinate, colnames(selected$X_inference))
  testthat::expect_equal(result$table$beta_tilde, direct$beta_tilde[index], tolerance = 1e-12)
  testthat::expect_equal(result$table$se, direct$se[index], tolerance = 1e-12)
  testthat::expect_identical(result$ranking$selection_rank, seq_len(8L))
  testthat::expect_true(all(diff(abs(result$ranking$beta_l1)) <= 0))
  jobs <- data.frame(job_id = 1:2, split = 1L, swapped = c(FALSE, TRUE), model = "A", tau = .5,
                     role_id = c("split01_swap0", "split01_swap1"))
  summary <- summarise_exploratory_gene_roles_v2(jobs, list(result, NULL))
  testthat::expect_equal(sum(summary$selection_frequency$selected_roles), 3L)
  testthat::expect_true(all(summary$selection_frequency$scheduled_roles == 2L))
  testthat::expect_true(all(summary$selection_frequency$selection_completed_roles == 1L))
  testthat::expect_true(all(summary$coefficient_median_iqr$n_inference_ok_roles == 1L))
  testthat::expect_equal(nrow(summary$intervals), 3L)
  testthat::expect_identical(summary$role_metrics$status, c("ok", "not_attempted"))
})

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)

# 四折沿用生成数据的 cluster 分组，空折对应移除该组全部 clusters。
cv_fold_fixture <- function(counts = NULL, seed = 20260928L, n_clusters = 14L) {
  dat <- generate_profile_qr_data_v2(n_clusters, 6L, 2L, 0.5, q = 1L,
                                     m_values = 3:5, seed = seed)
  beta <- dat$beta0 + c(0.06, -0.03, 0.02, -0.01, 0, 0)
  comp <- profile_components_v2(dat$y, dat$X, dat$Z, dat$cluster_id, beta,
                                0.5, 0.4, need_hessian = TRUE)
  real <- as.integer(comp$hessian_fold_counts)
  if (is.null(counts)) counts <- real
  if (length(counts) != length(real) || any(!counts %in% c(0L, real)) ||
      any(counts != 0L & counts != real)) stop("Only complete physical fold removal is supported.")
  sums <- comp$hessian_fold_sums
  if (any(counts == 0L)) {
    ids <- unique(as.character(dat$cluster_id))
    bucket <- rep(seq_along(real), length.out = length(ids))
    dat <- subset_clusters_v2(dat, ids[counts[bucket] > 0L])
    for (j in which(counts == 0L)) sums[[j]][] <- 0
    comp <- profile_components_v2(dat$y, dat$X, dat$Z, dat$cluster_id, beta,
                                  0.5, 0.4, need_hessian = TRUE)
  }
  if (max(abs(sum(counts) * comp$hessian - Reduce(`+`, sums))) > 1e-12) {
    stop("Physical fold sums do not reproduce the profile Hessian.")
  }
  list(dat = dat, beta = beta, comp = comp,
       H = comp$hessian,
       fold_sums = sums,
       fold_counts = as.integer(counts),
       real_counts = real)
}

# 独立重算：对给定 mu，在参与折上用同一个公开求解器复算留出缺陷。
held_out_defects <- function(fx, mu, coordinate = 1L, folds = NULL) {
  H <- fx$H
  counts <- fx$fold_counts
  n <- sum(counts)
  if (is.null(folds)) folds <- which(counts > 0L & counts < n)
  e <- numeric(ncol(H)); e[coordinate] <- 1
  vapply(folds, function(f) {
    cnt <- counts[f]
    s <- fx$fold_sums[[f]]
    H_val <- (s + t(s)) / 2 / cnt
    H_train <- (n * H - cnt * H_val) / (n - cnt)
    H_train <- (H_train + t(H_train)) / 2
    sol <- solve_dantzig_row_v2(H_train, coordinate = coordinate, mu_grid = mu)
    testthat::expect_true(isTRUE(sol$feasible))
    max(abs(as.numeric(H_val %*% as.numeric(sol$omega)) - e))
  }, numeric(1))
}

testthat::test_that("a fold without held-out clusters is not averaged in as a zero defect", {
  fx0 <- cv_fold_fixture()
  real <- fx0$real_counts
  testthat::expect_true(all(real > 0L))
  counts <- real
  counts[3L] <- 0L                      # 第三个 fold 没有留出样本
  fx <- cv_fold_fixture(counts = counts)
  testthat::expect_equal(sum(fx$fold_counts) * fx$H, Reduce(`+`, fx$fold_sums),
                         tolerance = 1e-12)
  testthat::expect_equal(max(abs(fx$fold_sums[[3L]])), 0)
  testthat::expect_equal(fx$dat$n_clusters, sum(counts))
  used <- which(counts > 0L & counts < sum(counts))
  testthat::expect_identical(used, c(1L, 2L, 4L))

  cv <- select_mu_inverse_hessian_cv_v2(
    H = fx$H, coordinate = 1L, mu_grid = c(0.05, 0.25), fold_sums = fx$fold_sums,
    fold_counts = fx$fold_counts, n_folds = 4L
  )

  testthat::expect_true(cv$selected)
  testthat::expect_identical(cv$folds_used, 3L)
  testthat::expect_identical(cv$candidates$folds_used, rep(3L, nrow(cv$candidates)))
  testthat::expect_identical(cv$candidates$skipped_folds, rep(1L, nrow(cv$candidates)))
  testthat::expect_identical(cv$candidates$infeasible_folds, rep(0L, nrow(cv$candidates)))

  for (i in seq_len(nrow(cv$candidates))) {
    mu <- cv$candidates$mu[i]
    d <- held_out_defects(fx, mu, folds = used)
    testthat::expect_equal(cv$candidates$defect_mean[i], mean(d), tolerance = 1e-12)
    # 把空折当作零的旧算法会得到 sum(d)/4，必须与当前记录不同
    testthat::expect_false(isTRUE(all.equal(cv$candidates$defect_mean[i], sum(d) / 4)))
  }
  testthat::expect_equal(cv$defect_mean, min(cv$candidates$defect_mean), tolerance = 1e-12)
})

testthat::test_that("a candidate requires at least half of the folds", {
  fx0 <- cv_fold_fixture()
  counts <- c(fx0$real_counts[1L], 0L, 0L, 0L)
  fx <- cv_fold_fixture(counts = counts)

  cv <- select_mu_inverse_hessian_cv_v2(
    H = fx$H, coordinate = 1L, mu_grid = c(0.05, 0.25), fold_sums = fx$fold_sums,
    fold_counts = fx$fold_counts, n_folds = 4L
  )
  testthat::expect_false(cv$selected)
  testthat::expect_true(is.na(cv$mu))
  testthat::expect_true(is.na(cv$folds_used))
  testthat::expect_identical(nrow(cv$candidates), 0L)
})

testthat::test_that("exactly half of the folds keeps the candidate", {
  fx0 <- cv_fold_fixture(n_clusters = 40L)
  counts <- c(fx0$real_counts[1L], fx0$real_counts[2L], 0L, 0L)
  fx <- cv_fold_fixture(counts = counts, n_clusters = 40L)
  for (j in c(1L, 2L)) {
    testthat::expect_gt(min(eigen(fx$fold_sums[[j]] / counts[[j]],
                                 symmetric = TRUE, only.values = TRUE)$values), 1e-8)
  }

  cv <- select_mu_inverse_hessian_cv_v2(
    H = fx$H, coordinate = 1L, mu_grid = c(0.05), fold_sums = fx$fold_sums,
    fold_counts = fx$fold_counts, n_folds = 4L
  )
  testthat::expect_true(cv$selected)
  testthat::expect_identical(cv$folds_used, 2L)
  d <- held_out_defects(fx, cv$mu, folds = c(1L, 2L))
  testthat::expect_equal(cv$defect_mean, mean(d), tolerance = 1e-12)
})

testthat::test_that("the parallel aggregation keeps the same fold semantics", {
  withr::local_envvar(c(JPDMEQR_CV_CORES = NA))
  fx0 <- cv_fold_fixture()
  counts <- fx0$real_counts
  counts[3L] <- 0L
  fx <- cv_fold_fixture(counts = counts)

  serial <- select_mu_inverse_hessian_cv_v2(
    H = fx$H, coordinate = 1L, mu_grid = c(0.05, 0.25), fold_sums = fx$fold_sums,
    fold_counts = fx$fold_counts, n_folds = 4L
  )
  old <- getOption("jpDMEQR.cv_cores", 1L)
  on.exit(options(jpDMEQR.cv_cores = old), add = TRUE)
  options(jpDMEQR.cv_cores = 2L)
  par <- select_mu_inverse_hessian_cv_v2(
    H = fx$H, coordinate = 1L, mu_grid = c(0.05, 0.25), fold_sums = fx$fold_sums,
    fold_counts = fx$fold_counts, n_folds = 4L
  )

  testthat::expect_identical(par$folds_used, serial$folds_used)
  testthat::expect_identical(par$mu, serial$mu)
  testthat::expect_identical(par$candidates$defect_mean, serial$candidates$defect_mean)
  testthat::expect_identical(par$candidates$folds_used, serial$candidates$folds_used)
})

testthat::test_that("OSQP rejects actual unconverged empty-fold solutions at the declared residual", {
  fx0 <- cv_fold_fixture()
  counts <- fx0$real_counts
  counts[3L] <- 0L
  fx <- cv_fold_fixture(counts = counts)

  old <- Sys.getenv("JPDMEQR_DANTZIG_SOLVER", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("JPDMEQR_DANTZIG_SOLVER") else
      Sys.setenv(JPDMEQR_DANTZIG_SOLVER = old)
  }, add = TRUE)
  Sys.setenv(JPDMEQR_DANTZIG_SOLVER = "osqp")

  cv <- select_mu_inverse_hessian_cv_v2(
    H = fx$H, coordinate = 1L, mu_grid = c(0.05, 0.25), fold_sums = fx$fold_sums,
    fold_counts = fx$fold_counts, n_folds = 4L
  )
  testthat::expect_false(cv$selected)
  testthat::expect_identical(nrow(cv$candidates), 0L)
  testthat::expect_true(is.na(cv$mu))
  testthat::expect_true(is.na(cv$folds_used))
  used <- which(counts > 0L & counts < sum(counts))
  testthat::expect_identical(used, c(1L, 2L, 4L))
  accepted <- matrix(FALSE, length(used), 2L)
  solutions <- vector("list", length(used))
  e <- c(1, rep(0, ncol(fx$H) - 1L))
  for (j in seq_along(used)) {
    fold <- used[[j]]
    n_f <- counts[[fold]]
    H_validation <- (fx$fold_sums[[fold]] + t(fx$fold_sums[[fold]])) / 2 / n_f
    H_all <- (fx$H + t(fx$H)) / 2
    H_training <- (sum(counts) * H_all - n_f * H_validation) / (sum(counts) - n_f)
    H_training <- (H_training + t(H_training)) / 2
    parts <- dantzig_build_osqp_model_v2(H_training, ncol(fx$H))
    solutions[[j]] <- dantzig_osqp_solve_chain_v2(parts, e, c(0.05, 0.25),
                                                feas_tol = 1e-6, max_iter = 20000L)
    for (k in seq_len(2L)) {
      sol <- solutions[[j]][[k]]
      residual <- max(abs(as.numeric(H_training %*% sol$omega) - e))
      testthat::expect_equal(sol$residual, residual, tolerance = 1e-12)
      testthat::expect_identical(sol$feasible, residual <= sol$mu + 1e-6)
      testthat::expect_true(nzchar(sol$status))
      accepted[j, k] <- sol$feasible
    }
  }
  testthat::expect_true(all(colSums(accepted) < ceiling(4L / 2)))
  artifact_dir <- getOption("jpDMEQR.unit_artifact_dir")
  if (!is.null(artifact_dir)) saveRDS(list(input = fx, cv = cv, solutions = solutions,
                                          accepted = accepted, max_iter = 20000L, residual_tolerance = 1e-6),
                                      file.path(artifact_dir, "osqp_empty_fold_actual.rds"))
})

debias_fixture <- function(with_folds = TRUE, seed = 20260928L) {
  dat <- generate_profile_qr_data_v2(14L, 6L, 2L, 0.5, q = 1L,
                                     m_values = 3:5, seed = seed)
  beta <- setNames(dat$beta0 + c(0.06, -0.03, 0.02, -0.01, 0, 0), colnames(dat$X))
  comp <- profile_components_v2(dat$y, dat$X, dat$Z, dat$cluster_id,
                                as.numeric(beta), 0.5, 0.4, need_hessian = TRUE)
  if (!with_folds) {
    comp$hessian_fold_sums <- NULL
    comp$hessian_fold_counts <- NULL
  }
  list(fit = list(beta = beta, components = comp), comp = comp)
}

testthat::test_that("the adjacent direction diagnostic reads the cross-validation column", {
  fx <- debias_fixture(with_folds = TRUE)
  inf <- debias_profile_coordinates_v2(fx$fit, coordinates = "x00001",
                                       mu_grid = c(0.05, 0.25))

  # 14 个 cluster 落在 n < 200 的分支，交叉验证使用 2 折
  testthat::expect_identical(inf$table$mu_source, "fold_cv")
  testthat::expect_identical(inf$table$mu_cv_folds_used, 2L)
  testthat::expect_true(is.finite(inf$table$adjacent_stability))
  testthat::expect_true(is.finite(inf$table$mu_cv_defect))
  testthat::expect_true(inf$table$mu %in% c(0.05, 0.25))
  testthat::expect_equal(inf$table$mu_cv_defect, inf$table$mu_cv_min_defect,
                         tolerance = 1e-12)
})

testthat::test_that("a fit without fold information is labelled as the smallest feasible value", {
  fx <- debias_fixture(with_folds = FALSE)
  inf <- debias_profile_coordinates_v2(fx$fit, coordinates = "x00001",
                                       mu_grid = c(0.05, 0.25))

  testthat::expect_identical(inf$table$mu_source, "smallest_feasible")
  testthat::expect_true(is.na(inf$table$mu_cv_folds_used))
  testthat::expect_true(is.na(inf$table$mu_cv_defect))
  testthat::expect_true(is.na(inf$table$adjacent_stability))
  testthat::expect_true(inf$table$feasible)
  testthat::expect_equal(inf$table$mu, min(c(0.05, 0.25)), tolerance = 1e-12)
})

small_fit_data <- function(seed = 20260928L) {
  generate_profile_qr_data_v2(8L, 5L, 2L, 0.5, q = 1L, m_values = 3:5, seed = seed)
}

testthat::test_that("an unknown fit control field raises an error", {
  dat <- small_fit_data()
  testthat::expect_error(
    fit_profile_lasso_v2(dat$y, dat$X, dat$Z, dat$cluster_id, 0.5, 0.4,
                         lambda_beta = 0.05, lambda_gamma = 1,
                         penalty_factor = rep(1, ncol(dat$X)),
                         control = list(kkt_tol = 1e-5)),
    "Unknown fit control field"
  )
})

testthat::test_that("the fit records the requested normalised KKT tolerance", {
  dat <- small_fit_data()
  fit <- fit_profile_lasso_v2(dat$y, dat$X, dat$Z, dat$cluster_id, 0.5, 0.4,
                              lambda_beta = 0.05, lambda_gamma = 1,
                              penalty_factor = rep(1, ncol(dat$X)),
                              control = list(max_iter = 200L,
                                             kkt_normalized_tol = 1e-3))
  testthat::expect_identical(fit$kkt_normalized_tol, 1e-3)
  testthat::expect_true(is.finite(fit$kkt_normalized))
})

testthat::test_that("the target call sites use the field the fit reads", {
  dat <- small_fit_data()
  target_control <- eval(formals(approximate_profile_target_v2)$fit_control)
  fit <- fit_profile_lasso_v2(dat$y, dat$X, dat$Z, dat$cluster_id, 0.5, 0.4,
                              lambda_beta = 0, lambda_gamma = 1,
                              penalty_factor = rep(0, ncol(dat$X)),
                              control = target_control)
  testthat::expect_identical(fit$kkt_normalized_tol, target_control$kkt_normalized_tol)

  for (f in c("R/v2_simulation_v2_part03.R", "R/v2_profile_v2_part05.R",
              "scripts/simulation/04_build_profile_targets.R",
              "scripts/application/03_fit_GSE65391.R")) {
    txt <- readLines(file.path(root, f), warn = FALSE)
    testthat::expect_false(any(grepl("kkt_tol\\s*=", txt)), info = f)
  }
})

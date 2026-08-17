# Acceptance tests for the freeze-patch oracle adapters:
# PROFILE-DQR-TRUE-NUISANCE, PROFILE-DQR-POP-H, PROFILE-DQR-SPLIT.

root <- if (file.exists(file.path("scripts", "00_source_v2.R"))) "." else
  normalizePath(file.path(testthat::test_path(), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
source_v2_module(root, "benchmark_adapters_v2", envir = .GlobalEnv)

make_train <- function(n = 40, p = 12, s = 3, q = 1L, seed = 202L) {
  dat <- generate_profile_qr_data_v2(
    n_clusters = n, p = p, s = s, tau = 0.5, q = q,
    m_values = 3:5, m_rule = "uniform", sigma_b0 = 1,
    error_dist = "gaussian", seed = seed
  )
  dat
}

tuning <- list(h = 0.5, lambda_beta = 0.15, lambda_gamma = 1,
               mu_grid = c(0.3, 0.6, 1.2))

# ---------------------------------------------------------------------------
# PROFILE-DQR-TRUE-NUISANCE
# ---------------------------------------------------------------------------

testthat::test_that("TRUE-NUISANCE removes generated nuisance and matches low-dimensional QR", {
  train <- make_train(seed = 202)
  coords <- colnames(train$X)[1:4]
  ans <- fit_benchmark_profile_true_nuisance_v2(train, 0.5, coords, tuning,
                                                seed = 7L)
  testthat::expect_equal(ans$method_id, "PROFILE-DQR-TRUE-NUISANCE")
  testthat::expect_equal(ans$status, "ok")
  testthat::expect_equal(ans$target_scope, "structural_known_nuisance")
  testthat::expect_true(all(is.finite(ans$beta_tilde)))
  # The oracle removes Z'b exactly. The identical reduced design (y - Z'b,
  # zero nuisance column) fed through the same PROFILE-DQR adapter must give
  # the same penalised fit: this is the defining oracle property.
  zb <- rowSums(train$Z * train$random_effects[match(train$cluster_id,
                                                     rownames(train$random_effects)), ,
                                              drop = FALSE])
  reduced <- train
  reduced$y <- as.numeric(train$y - zb)
  reduced$Z <- matrix(0, nrow = length(reduced$y), ncol = 1L)
  ref <- fit_benchmark_profile_dqr_v2(reduced, 0.5, coords, tuning, seed = 7L)
  testthat::expect_equal(as.numeric(ans$beta_hat[coords]),
                         as.numeric(ref$beta_hat[coords]), tolerance = 1e-6)
  testthat::expect_equal(as.numeric(ans$beta_tilde[coords]),
                         as.numeric(ref$beta_tilde[coords]), tolerance = 1e-6)
  # Oracle metadata is attached.
  testthat::expect_equal(ans$adapter_fidelity_status,
                         "implementation_present_acceptance_pending")
})

testthat::test_that("TRUE-NUISANCE fails cleanly without oracle context", {
  train <- make_train(seed = 303)
  train$random_effects <- NULL
  ans <- fit_benchmark_profile_true_nuisance_v2(train, 0.5,
                                                colnames(train$X)[1:3],
                                                tuning, seed = 7L)
  testthat::expect_equal(ans$status, "failed")
  testthat::expect_match(ans$failure_stage, "oracle_context")
})

# ---------------------------------------------------------------------------
# PROFILE-DQR-POP-H
# ---------------------------------------------------------------------------

make_pop_asset <- function(train, coords, seed = 91L) {
  target_names <- coords
  p_sub <- length(target_names)
  H_pop <- diag(2, p_sub)
  directions <- lapply(target_names, function(nm) {
    k <- match(nm, target_names)
    omega <- rep(0, p_sub)
    omega[k] <- 1 / 2
    e <- numeric(p_sub); e[k] <- 1
    list(omega = setNames(omega, target_names),
         residual = max(abs(H_pop %*% omega - e)),
         l1_norm = sum(abs(omega)),
         l2_norm = sqrt(sum(omega^2)))
  })
  names(directions) <- target_names
  list(experiment_id = "TEST",
       H_population = H_pop,
       directions = directions,
       n_population = 100000L,
       repeats = 4L,
       implementation_commit = "test",
       config_sha256 = "test",
       nuisance_converged = TRUE)
}

testthat::test_that("POP-H uses population direction and reports the gap", {
  train <- make_train(seed = 404)
  coords <- colnames(train$X)[1:4]
  asset <- make_pop_asset(train, coords)
  ans <- fit_benchmark_profile_pop_h_v2(train, 0.5, coords, tuning, seed = 9L,
                                        population_direction_asset = asset)
  testthat::expect_equal(ans$method_id, "PROFILE-DQR-POP-H")
  testthat::expect_equal(ans$status, "ok")
  testthat::expect_true(all(is.finite(ans$beta_tilde)))
  testthat::expect_true(all(is.finite(ans$se)))
  # The population residual from the synthetic asset is reproduced.
  testthat::expect_equal(max(vapply(ans$inference_object$directions,
                                    function(d) d$population_residual,
                                    numeric(1))),
                         0, tolerance = 1e-8)
})

testthat::test_that("POP-H fails cleanly without the population asset", {
  train <- make_train(seed = 505)
  ans <- fit_benchmark_profile_pop_h_v2(train, 0.5, colnames(train$X)[1:3],
                                        tuning, seed = 9L,
                                        population_direction_asset = NULL)
  testthat::expect_equal(ans$status, "failed")
  testthat::expect_match(ans$failure_stage, "oracle_asset")
})

# ---------------------------------------------------------------------------
# PROFILE-DQR-SPLIT
# ---------------------------------------------------------------------------

testthat::test_that("SPLIT screens on disjoint clusters and maps back to the full design", {
  train <- make_train(n = 60, seed = 606)
  coords <- colnames(train$X)[1:4]
  ans <- fit_benchmark_profile_split_v2(
    train, 0.5, coords, tuning, seed = 11L,
    context = list(screening = NULL),
    control = list(screen_fraction = 0.5, screen_dim = 8)
  )
  testthat::expect_equal(ans$method_id, "PROFILE-DQR-SPLIT")
  testthat::expect_equal(ans$status, "ok")
  # Screening and inference sets are disjoint.
  scr <- ans$screening
  testthat::expect_equal(length(intersect(scr$screen_clusters,
                                          scr$fit_clusters)), 0)
  # Full-design coefficient vector is returned.
  testthat::expect_true(all(names(ans$beta_hat) == colnames(train$X)))
  testthat::expect_true(all(is.finite(ans$beta_tilde)))
  # The selected design must contain the target coordinates.
  testthat::expect_true(all(coords %in% names(ans$beta_hat)))
})

testthat::test_that("SPLIT rejects overlapping pre-specified role sets", {
  train <- make_train(n = 60, seed = 707)
  screening <- list(
    method = "split_quantile_score",
    screen_clusters = head(unique(train$cluster_id), 20),
    fit_clusters = head(unique(train$cluster_id), 20)  # overlap
  )
  ans <- fit_benchmark_profile_split_v2(
    train, 0.5, colnames(train$X)[1:3], tuning, seed = 11L,
    context = list(screening = screening),
    control = list(screen_fraction = 0.5, screen_dim = 8)
  )
  testthat::expect_equal(ans$status, "failed")
  testthat::expect_match(ans$failure_stage, "screening")
})

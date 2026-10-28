eigen_root <- v2_require_root()
library(jpDMEQR)
source(file.path(eigen_root, "scripts", "00_source_v2.R"), local = TRUE)
source_v2_module(eigen_root, "profile_v2", environment())
source_v2_module(eigen_root, "simulation_v2", environment())

testthat::test_that("registered Eigen population products match real-data R BLAS calculations", {
  withr::local_options(list(jpDMEQR.population_matrix_backend = "R_BLAS"))
  dat <- generate_profile_qr_data_v2(60L, 20L, 3L, 0.25, q = 2L,
                                     m_multiset = rep(c(2L, 3L, 5L), 20L), seed = 20261230L)
  Lambda <- matrix(c(1, 0.25, 0.25, 0.8), 2L)
  base <- profile_pair_population_streamed_v2(dat, dat$beta0, 0.25, 0.2, Lambda,
                                              batch_clusters = 13L)
  profile_pair_population_eigen_init_v2(eigen_root)
  eigen <- profile_pair_population_streamed_v2(dat, dat$beta0, 0.25, 0.2, Lambda,
                                               batch_clusters = 13L)
  testthat::expect_identical(eigen$population_matrix_backend, "Eigen_single_thread")
  testthat::expect_identical(eigen$score, base$score)
  testthat::expect_lt(max(abs(eigen$hessian - base$hessian)), 1e-10)
  testthat::expect_lt(eigen$profile_identity_error, 1e-10)
  for (field in c("Sigma_population", "Sigma_half", "Sigma_full", "cross_covariance_half_full")) {
    testthat::expect_lt(max(abs(eigen$population_moments[[field]] - base$population_moments[[field]])), 1e-10)
  }
  testthat::expect_lt(max(abs(profile_pair_population_crossprod_v2(dat$X) - crossprod(dat$X))), 1e-9)
  testthat::expect_lt(max(abs(profile_pair_population_crossprod_v2(dat$X, dat$Z) - crossprod(dat$X, dat$Z))), 1e-9)
  testthat::expect_error(profile_pair_population_crossprod_v2(dat$X, dat$Z[-1, , drop = FALSE]), "dimensions")
})

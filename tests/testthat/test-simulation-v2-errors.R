root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)

testthat::test_that("quantile-centred error generators have the requested quantile", {
  set.seed(20260817)
  distributions <- c("gaussian", "t3", "laplace", "skew_chisq3",
                     "contaminated_normal", "asymmetric_laplace")
  for (tau in c(0.25, 0.5, 0.75)) {
    for (dist in distributions) {
      e <- r_quantile_centered_error_v2(120000L, tau, dist)
      testthat::expect_true(
        abs(as.numeric(stats::quantile(e, tau, names = FALSE))) < 0.035,
        info = paste(dist, tau)
      )
    }
  }
})

testthat::test_that("positive heteroskedastic scaling preserves zero conditional quantile", {
  set.seed(280817)
  n <- 100000L
  x <- stats::rnorm(n)
  e <- r_quantile_centered_error_v2(n, 0.25, "t3", x1 = x, heteroskedastic = TRUE)
  groups <- cut(x, breaks = stats::quantile(x, seq(0, 1, length.out = 6)), include.lowest = TRUE)
  q <- tapply(e, groups, stats::quantile, probs = 0.25)
  testthat::expect_lt(max(abs(q)), 0.06)
})

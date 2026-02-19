test_that("default_lambda_clime validates h and returns positive value", {
  expect_error(default_lambda_clime(p = 20, N = 100, h = 0), "h must be positive")
  expect_error(default_lambda_clime(p = 20, N = 100, h = -0.1), "h must be positive")
  expect_gt(default_lambda_clime(p = 20, N = 100, h = 0.2), 0)
})

test_that("data generation switches random-effects dimension by scenario", {
  set.seed(11)
  s1 <- generate_mixed_qr_data(
    n = 15,
    p = 30,
    s = 4,
    tau = 0.5,
    scenario = "S1",
    m_choices = c(2, 3)
  )
  s2 <- generate_mixed_qr_data(
    n = 15,
    p = 30,
    s = 4,
    tau = 0.5,
    scenario = "S2",
    m_choices = c(2, 3)
  )

  expect_equal(ncol(s1$Z), 1L)
  expect_equal(ncol(s1$gamma0), 1L)
  expect_equal(ncol(s2$Z), 2L)
  expect_equal(ncol(s2$gamma0), 2L)
})

test_that("screening and fitting input checks fail on malformed inputs", {
  set.seed(12)
  dat <- generate_mixed_qr_data(
    n = 10,
    p = 20,
    s = 4,
    tau = 0.5,
    scenario = "S1",
    m_choices = c(2, 3)
  )

  expect_error(
    CA_IQR_SIS(
      y = dat$y,
      X = dat$X,
      cluster_id = dat$cluster_id[-1],
      d_target = 5
    ),
    "cluster_id length mismatch"
  )

  expect_error(
    CA_IQR_SIS(
      y = dat$y,
      X = dat$X,
      cluster_id = dat$cluster_id,
      tau_grid = c(-0.1, 0.5),
      d_target = 5
    ),
    "tau_grid must be in"
  )

  expect_error(
    JP_DME_QR(
      y = dat$y,
      X = dat$X,
      Z = dat$Z,
      cluster_id = dat$cluster_id[-1],
      tau = 0.5,
      h = dat$N^(-1 / 3),
      lambda1 = 0.1,
      lambda2 = 1
    ),
    "Dimension mismatch"
  )
})

test_that("core pipeline runs on a small instance", {
  set.seed(101)
  dat <- generate_mixed_qr_data(
    n = 20,
    p = 40,
    s = 5,
    tau = 0.5,
    scenario = "S1",
    m_choices = c(3, 4)
  )

  screened <- CA_IQR_SIS(
    y = dat$y,
    X = dat$X,
    cluster_id = dat$cluster_id,
    d_target = 15,
    T_iter = 2
  )
  expect_true(length(screened) >= 1)

  X_red <- dat$X[, screened, drop = FALSE]
  fit <- JP_DME_QR(
    y = dat$y,
    X = X_red,
    Z = dat$Z,
    cluster_id = dat$cluster_id,
    tau = 0.5,
    h = dat$N^(-1 / 3),
    lambda1 = sqrt(log(max(2, ncol(X_red))) / dat$N),
    lambda2 = 1,
    max_iter = 20,
    verbose = FALSE
  )

  expect_length(fit$beta, ncol(X_red))
  expect_true(all(is.finite(fit$beta)))

  inf <- Debias_Inference(
    y = dat$y,
    X = X_red,
    Z = dat$Z,
    cluster_id = dat$cluster_id,
    beta_hat = fit$beta,
    gamma_hat = fit$gamma,
    tau = 0.5,
    h = dat$N^(-1 / 3),
    lambda2 = 1,
    method = "RIDGE"
  )

  expect_length(inf$beta_tilde, ncol(X_red))
  expect_length(inf$se, ncol(X_red))
  expect_true(all(is.finite(inf$beta_tilde)))
})

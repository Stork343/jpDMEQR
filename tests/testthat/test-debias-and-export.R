test_that("debiased inference works with RIDGE and CLIME-ADMM backends", {
  set.seed(21)
  dat <- generate_mixed_qr_data(
    n = 12,
    p = 24,
    s = 5,
    tau = 0.5,
    scenario = "S1",
    m_choices = c(2, 3)
  )

  screened <- CA_IQR_SIS(
    y = dat$y,
    X = dat$X,
    cluster_id = dat$cluster_id,
    d_target = 6,
    T_iter = 1
  )
  X_red <- dat$X[, screened, drop = FALSE]

  fit <- JP_DME_QR(
    y = dat$y,
    X = X_red,
    Z = dat$Z,
    cluster_id = dat$cluster_id,
    tau = 0.5,
    h = dat$N^(-1 / 3),
    lambda1 = 0.1,
    lambda2 = 1,
    max_iter = 15,
    tol = 1e-3,
    obj_tol = 1e-4,
    verbose = FALSE
  )

  inf_ridge <- Debias_Inference(
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
  expect_length(inf_ridge$beta_tilde, ncol(X_red))
  expect_true(all(is.finite(inf_ridge$beta_tilde)))
  expect_true(all(is.finite(inf_ridge$se)))

  inf_clime <- Debias_Inference(
    y = dat$y,
    X = X_red,
    Z = dat$Z,
    cluster_id = dat$cluster_id,
    beta_hat = fit$beta,
    gamma_hat = fit$gamma,
    tau = 0.5,
    h = dat$N^(-1 / 3),
    lambda2 = 1,
    lambda_clime = 0.5,
    method = "CLIME",
    clime_method = "admm",
    clime_admm_control = list(max_iter = 200, tol = 1e-3, ista_max_iter = 100, ista_tol = 1e-5)
  )
  expect_length(inf_clime$beta_tilde, ncol(X_red))
  expect_true(all(is.finite(inf_clime$beta_tilde)))
  expect_true(all(is.finite(inf_clime$se)))
})

test_that("simulation summary and CSV export interfaces work", {
  set.seed(22)
  out <- Simulation_Study(
    B = 1,
    n = 20,
    m_choices = c(3, 4),
    n_test = 150,
    p_grid = c(30),
    tau_grid = c(0.5),
    scenarios = c("S1"),
    methods = c("JP-DME-QR", "Naive QR-Lasso"),
    inference_mode = "RIDGE",
    jp_target = "beta0",
    verbose = FALSE
  )

  sm <- summarise_simulation_results(out$raw_results)
  expect_gt(nrow(sm), 0)
  expect_true(all(c("n_runs", "n_ok") %in% names(sm)))

  out_dir <- file.path(tempdir(), "jpDMEQR_export_test")
  paths <- export_simulation_outputs(out, out_dir = out_dir, prefix = "unit")
  expect_true(file.exists(paths$raw))
  expect_true(file.exists(paths$summary))
})

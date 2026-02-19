test_that("simulation API returns expected structure", {
  set.seed(202)
  out <- Simulation_Study(
    B = 1,
    n = 20,
    m_choices = c(3, 4),
    n_test = 200,
    p_grid = c(30),
    tau_grid = c(0.5),
    scenarios = c("S1"),
    methods = c("JP-DME-QR", "Naive QR-Lasso"),
    inference_mode = "RIDGE",
    jp_target = "beta0",
    verbose = FALSE
  )

  expect_true(is.list(out))
  expect_true(all(c("config_grid", "raw_results", "summary") %in% names(out)))
  expect_gt(nrow(out$raw_results), 0)
  expect_gt(nrow(out$summary), 0)
})

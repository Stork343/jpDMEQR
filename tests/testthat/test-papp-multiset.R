# Round-6 PAPP wiring: frozen empirical cluster-size multiset in the DGP.

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)

testthat::test_that("m_multiset is used verbatim by the DGP", {
  ms <- c(5, 3, 8, 4, 6)   # small frozen multiset
  dat <- generate_profile_qr_data_v2(
    n_clusters = 5, p = 20, s = 3, tau = 0.5, q = 1,
    m_values = 2:9, m_multiset = ms, seed = 1L
  )
  testthat::expect_equal(as.integer(dat$m_i), ms)
  testthat::expect_equal(dat$m_summary$m_max, max(ms))
  # Without the multiset the sizes are drawn (support respected).
  dat2 <- generate_profile_qr_data_v2(
    n_clusters = 5, p = 20, s = 3, tau = 0.5, q = 1,
    m_values = 2:9, seed = 2L
  )
  testthat::expect_true(all(dat2$m_i %in% 2:9))
  # Length mismatch errors.
  testthat::expect_error(
    generate_profile_qr_data_v2(n_clusters = 5, p = 20, s = 3, tau = 0.5,
                                m_multiset = c(3, 4), seed = 1L),
    "exactly n_clusters"
  )
})

testthat::test_that("config row m_multiset parses to a numeric vector", {
  cfg <- read.csv(file.path(root, "config", "simulation_papp.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
  lst <- config_row_to_list_v2(cfg[1, , drop = FALSE])
  testthat::expect_equal(length(lst$m_multiset), 129L)
  testthat::expect_equal(sum(lst$m_multiset), 867L)
})
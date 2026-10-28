root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
sys.source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"),
           envir = environment())

testthat::test_that("frozen main registry expects every manuscript fragment and figure", {
  registry <- read.csv(file.path(root, "config", "simulation_main_b200.csv"),
                       stringsAsFactors = FALSE)
  expected <- expected_manuscript_assets_v2(registry)
  testthat::expect_equal(sort(expected$fragments), sort(c(
    "tab_simulation_main.tex", "tab_simulation_methods.tex",
    "tab_simulation_inference.tex", "tab_simulation_theory.tex",
    "tab_simulation_failures.tex", "tab_simulation_rates.tex",
    "tab_simulation_computation.tex", "tab_simulation_selection.tex"
  )))
  testthat::expect_equal(sort(expected$figures), sort(c(
    "fig_type1_power.pdf", "fig_null_pvalue_uniformity.pdf",
    "fig_studentised_quantiles.pdf", "fig_rate_slopes.pdf",
    "fig_coverage_clusters.pdf", "fig_se_sd_ratio.pdf",
    "fig_computation_time.pdf"
  )))
})

testthat::test_that("preflight registry expects the fragments its modules support", {
  registry <- read.csv(file.path(root, "config", "simulation_preflight.csv"),
                       stringsAsFactors = FALSE)
  expected <- expected_manuscript_assets_v2(registry)
  testthat::expect_equal(sort(expected$fragments), sort(c(
    "tab_simulation_main.tex", "tab_simulation_methods.tex",
    "tab_simulation_inference.tex", "tab_simulation_theory.tex",
    "tab_simulation_failures.tex"
  )))
  testthat::expect_equal(sort(expected$figures), sort(c(
    "fig_type1_power.pdf", "fig_null_pvalue_uniformity.pdf",
    "fig_studentised_quantiles.pdf"
  )))
})

testthat::test_that("module contents drive the expected asset set", {
  module_a <- data.frame(module = "A_theory_scaling", target_coordinate_count = 0L)
  expected_a <- expected_manuscript_assets_v2(module_a)
  testthat::expect_true("tab_simulation_rates.tex" %in% expected_a$fragments)
  testthat::expect_false("tab_simulation_computation.tex" %in% expected_a$fragments)
  testthat::expect_false("tab_simulation_selection.tex" %in% expected_a$fragments)
  testthat::expect_equal(sort(expected_a$figures), sort(c(
    "fig_rate_slopes.pdf", "fig_coverage_clusters.pdf", "fig_se_sd_ratio.pdf"
  )))

  module_g <- data.frame(module = "G_computation", target_coordinate_count = 12L)
  expected_g <- expected_manuscript_assets_v2(module_g)
  testthat::expect_true("tab_simulation_computation.tex" %in% expected_g$fragments)
  testthat::expect_false("tab_simulation_rates.tex" %in% expected_g$fragments)
  testthat::expect_equal(sort(expected_g$figures), sort(c(
    "fig_type1_power.pdf", "fig_null_pvalue_uniformity.pdf",
    "fig_studentised_quantiles.pdf", "fig_computation_time.pdf"
  )))

  module_f <- data.frame(module = "F_screening", target_coordinate_count = 0L)
  expected_f <- expected_manuscript_assets_v2(module_f)
  testthat::expect_true("tab_simulation_selection.tex" %in% expected_f$fragments)
  testthat::expect_equal(expected_f$figures, character(0))
})
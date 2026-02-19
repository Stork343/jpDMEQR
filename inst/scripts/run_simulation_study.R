#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "simulation_outputs"
B <- if (length(args) >= 2) as.integer(args[[2]]) else 100L

if (!requireNamespace("jpDMEQR", quietly = TRUE)) {
  stop("Package 'jpDMEQR' is not installed. Install it first.")
}

set.seed(2026)
sim <- jpDMEQR::Simulation_Study(
  B = B,
  n = 200,
  m_choices = c(4, 5, 6),
  p_grid = c(200, 500, 1000, 2000),
  tau_grid = c(0.25, 0.5, 0.75),
  scenarios = c("S1", "S2", "S3", "S4", "S5", "S6"),
  methods = c("JP-DME-QR", "Naive QR-Lasso", "Bayesian Mixed-Lasso"),
  c_screen = 3,
  T_iter = 2,
  inference_mode = "auto",
  jp_target = "beta_star_ref",
  verbose = TRUE
)

jpDMEQR::export_simulation_outputs(sim, out_dir = out_dir, prefix = "jpdmeqr_main")
message("Simulation finished. Outputs written to: ", normalizePath(out_dir, mustWork = FALSE))

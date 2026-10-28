.libPaths(c(normalizePath("results/scratch/profile_pair_package_20261022/library", mustWork = TRUE), .libPaths()))
source("research/score_validity_20261026/load_joint.R")
library(highs)
for (module in c("metrics_v2", "benchmark_adapters_v2")) source_v2_module(".", module, .GlobalEnv)
root <- normalizePath(".")
source("scripts/simulation/_run_registry_helpers.R")
compile_submission <- function() {
  dir.create("results/scratch/submission_completion_20261027", recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock("results/scratch/submission_completion_20261027/compile.lock", timeout = 60000)
  if (is.null(lock)) stop("Unable to acquire the submission compilation lock.")
  on.exit(filelock::unlock(lock), add = TRUE)
  Rcpp::sourceCpp("research/submission_completion_20261027/sqr_fista.cpp",
    cacheDir = "results/scratch/submission_completion_20261027/cache", rebuild = FALSE, showOutput = FALSE, env = .GlobalEnv)
  Rcpp::sourceCpp("research/submission_completion_20261027/clime_anchor.cpp",
    cacheDir = "results/scratch/submission_completion_20261027/cache", rebuild = FALSE, showOutput = FALSE, env = .GlobalEnv)
}
compile_submission()
rm(compile_submission)
source("research/joint_profile_manuscript_20261026/fixed_inference.R")
source("research/submission_completion_20261027/simulation.R")
source("research/submission_completion_20261027/clime_column_generation.R")

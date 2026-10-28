library(filelock)
load_joint_profile <- function() {
  dir.create("results/scratch/score_validity_20261026", recursive = TRUE, showWarnings = FALSE)
  handle <- filelock::lock("results/scratch/score_validity_20261026/compile.lock", timeout = 60000)
  if (is.null(handle)) stop("Unable to acquire the compilation cache lock.")
  on.exit(filelock::unlock(handle), add = TRUE)
  library_path <- "results/scratch/score_validity_20261026/library"
  if (!dir.exists(file.path(library_path, "jpDMEQR"))) stop("Run research/score_validity_20261026/setup.R before loading the estimator.")
  .libPaths(c(normalizePath(library_path, mustWork = TRUE), .libPaths()))
  library(jpDMEQR)
  library(Rcpp)
  library(RcppEigen)
  sys.source("scripts/00_source_v2.R", envir = .GlobalEnv)
  source_v2_module(".", "profile_v2", .GlobalEnv)
  source_v2_module(".", "simulation_v2", .GlobalEnv)
  sys.source("research/positive_curvature_20261025/method.R", envir = .GlobalEnv)
  pc_compile()
  Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled", JPDMEQR_ASSEMBLY = "batch",
             JPDMEQR_MATRIX_BACKEND = "Eigen_single_thread", JPDMEQR_IDENTITY_DETAIL = "aggregate")
  sys.source("research/positive_curvature_20261025/orthogonal_greedy.R", envir = .GlobalEnv)
  Rcpp::sourceCpp("research/score_validity_20261026/geometry.cpp",
                  cacheDir = "results/scratch/score_validity_20261026/cache",
                  rebuild = FALSE, showOutput = FALSE, env = .GlobalEnv)
  sys.source("research/score_validity_20261026/method.R", envir = .GlobalEnv)
  sys.source("research/score_validity_20261026/small_profile.R", envir = .GlobalEnv)
  sys.source("research/score_validity_20261026/joint_profile_inference.R", envir = .GlobalEnv)
}
load_joint_profile()
rm(load_joint_profile)

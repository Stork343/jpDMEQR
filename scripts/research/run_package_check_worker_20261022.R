args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Expected repository, check output directory and scratch alias.")
root <- normalizePath(args[[1]], mustWork = TRUE)
setwd(root)
output <- normalizePath(args[[2]], mustWork = TRUE)
alias <- args[[3]]
scratch <- normalizePath(file.path(root, "results", "scratch"), mustWork = TRUE)
if (!identical(normalizePath(alias, mustWork = TRUE), scratch) ||
    !identical(Sys.getenv("TMPDIR"), alias)) stop("Package check scratch identity differs.")
library_path <- normalizePath("results/scratch/profile_pair_package_20261022/library", mustWork = TRUE)
.libPaths(c(library_path, .libPaths()))
check_env <- c(
  TMPDIR = alias, JPDMEQR_CLUSTER_SOLVER = "compiled", JPDMEQR_ASSEMBLY = "batch",
  JPDMEQR_MATRIX_BACKEND = "Eigen_single_thread", JPDMEQR_IDENTITY_DETAIL = "aggregate",
  JPDMEQR_UNPENALISED_OPTIMIZER = "newton_armijo", JPDMEQR_PROXIMAL_SCALING = "cluster_rms",
  JPDMEQR_DANTZIG_SOLVER = "dual_warm", JPDMEQR_HESSIAN_REFRESH = "20",
  JPDMEQR_SHARED_FIRST_STAGE = "TRUE", JPDMEQR_CV_CORES = "1",
  JPDMEQR_INITIAL_ACCELERATION = "clinical_kkt",
  `_R_CHECK_FORCE_SUGGESTS_` = "false", NOT_CRAN = "true"
)
do.call(Sys.setenv, as.list(check_env))
source(file.path(root, "scripts", "research", "prepare_package_source_20261022.R"))
staging <- prepare_package_source_20261022(root, file.path(output, "source", "jpDMEQR"))
writeLines(staging$files, file.path(output, "package_source_inventory.txt"))
result <- rcmdcheck::rcmdcheck(
  path = staging$path, check_dir = file.path(output, "check"),
  libpath = .libPaths(), args = character(), build_args = character(),
  env = check_env, error_on = "never", quiet = FALSE, timeout = Inf
)
saveRDS(result, file.path(output, "rcmdcheck.rds"))
writeLines(result$stdout, file.path(output, "check_stdout.txt"))
writeLines(result$stderr, file.path(output, "check_stderr.txt"))
manifest <- list(
  pass = length(result$errors) == 0L && length(result$warnings) == 0L,
  errors = result$errors, warnings = result$warnings, notes = result$notes,
  scratch_alias = alias, scratch_realpath = scratch, temp_directory = tempdir(),
  configuration = as.list(check_env), tests_enabled = TRUE, examples_enabled = TRUE,
  staging_test_file_count = staging$test_file_count, staging_rules = staging$rules,
  vignettes_enabled = TRUE, manual_enabled = TRUE,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
jsonlite::write_json(manifest, file.path(output, "manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null")
if (!manifest$pass) quit(status = 1L)

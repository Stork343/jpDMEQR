args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
setwd(root)
run_id <- if (length(args) >= 2L) args[[2L]] else format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
output <- file.path(root, "results", "scratch", "native_warning_records_20261022", run_id)
if (dir.exists(output)) stop("The source warning check directory already exists.")
dir.create(output, recursive = TRUE)
library_path <- normalizePath("results/scratch/profile_pair_package_20261022/library", mustWork = TRUE)
.libPaths(c(library_path, .libPaths()))
alias <- file.path(path.expand("~"), ".codex", "jpdmeqr-scratch-20261022")
if (!identical(normalizePath(alias, mustWork = TRUE), normalizePath("results/scratch"))) {
  stop("The check scratch alias differs from the project directory.")
}
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
dirty <- system2("git", c("status", "--porcelain"), stdout = TRUE)
if (length(args) >= 3L && (!identical(commit, args[[3L]]) || length(dirty))) {
  stop("Final source warning checks require the requested commit and clean source files.")
}
jsonlite::write_json(list(implementation_commit = commit, run_id = run_id,
                          dirty_source = dirty,
                          scratch_alias = alias, scratch_realpath = normalizePath(alias)),
                     file.path(output, "request.json"), auto_unbox = TRUE, pretty = TRUE)
child_env <- Sys.getenv()
child_env[c("TMPDIR", "R_LIBS_USER", "LC_CTYPE")] <- c(alias, library_path, "en_US.UTF-8")
invocation <- processx::run(
  file.path(R.home("bin"), "Rscript"),
  c("scripts/research/run_unit_file_worker_20261022.R", root,
    "tests/testthat/test-native-warning-records.R", output),
  wd = root, env = child_env, stdout = file.path(output, "stdout.log"),
  stderr = file.path(output, "stderr.log"), error_on_status = FALSE, timeout = Inf
)
manifest <- jsonlite::read_json(file.path(output, "manifest.json"), simplifyVector = TRUE)
if (invocation$status != 0L || !isTRUE(manifest$assertion_pass)) {
  stop("The real source warning comparison failed; inspect ", output)
}
message("Source warning comparison assertions: ", manifest$passed,
        "; recorded test warning events: ", manifest$warning_count)

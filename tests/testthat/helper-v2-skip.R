# Shared helper: locate the repository root for v2 source-dependent tests.
#
# These tests source the repository's modular R implementation directly
# (scripts/00_source_v2.R). Under `R CMD check` the tests run against the
# installed package, where the repository files are not available; the tests
# are then skipped rather than failing. Under local `test_dir()` runs from
# the repository root they execute fully.
v2_repo_root <- function() {
  if (file.exists(file.path("scripts", "00_source_v2.R"))) {
    return(".")
  }
  # Fallback: parent of the tests directory (local test_dir from repo root).
  cand <- normalizePath(file.path(testthat::test_path(), "../.."), mustWork = FALSE)
  if (file.exists(file.path(cand, "scripts", "00_source_v2.R"))) {
    return(cand)
  }
  NULL
}

v2_require_root <- function() {
  root <- v2_repo_root()
  if (is.null(root)) {
    testthat::skip("Repository v2 source files are not available in this environment.")
  }
  invisible(root)
}

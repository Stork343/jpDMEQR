args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Expected repository and a new run identifier.")
root <- normalizePath(args[[1]], mustWork = TRUE)
setwd(root)
run_id <- args[[2]]
phase <- if (length(args) >= 3L) args[[3]] else "all"
if (!phase %in% c("unit", "check", "all")) stop("Unknown check phase.")
output <- file.path(root, "results", "scratch", "repository_package_checks_20261022", run_id)
if (dir.exists(output)) stop("The check output directory already exists.")
dir.create(output, recursive = TRUE)
library_path <- normalizePath("results/scratch/profile_pair_package_20261022/library", mustWork = TRUE)
.libPaths(c(library_path, .libPaths()))
scratch <- normalizePath("results/scratch", mustWork = TRUE)
alias <- file.path(path.expand("~"), ".codex", "jpdmeqr-scratch-20261022")
link_target <- Sys.readlink(alias)
if (!file.exists(alias) && (is.na(link_target) || !nzchar(link_target))) {
  if (!file.symlink(scratch, alias)) stop("Cannot create the scratch path alias.")
}
if (!identical(normalizePath(alias, mustWork = TRUE), scratch)) {
  stop("The scratch alias does not point to the expected project directory.")
}
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
dirty <- system2("git", c("status", "--porcelain"), stdout = TRUE)
if (length(args) >= 4L && (!identical(commit, args[[4]]) || length(dirty))) {
  stop("Final checks require the requested commit and a clean source checkout.")
}
request <- list(run_id = run_id, phase = phase, implementation_commit = commit,
                dirty_source = dirty, package_version = as.character(packageVersion("jpDMEQR")),
                scratch_alias = alias, scratch_realpath = scratch,
                started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
jsonlite::write_json(request, file.path(output, "request.json"), auto_unbox = TRUE, pretty = TRUE)
unit_pass <- NA
check_pass <- NA
if (phase %in% c("unit", "all")) {
  files <- sort(list.files("tests/testthat", pattern = "^test-.*\\.R$", full.names = TRUE))
  if (!length(files)) stop("No repository test files are available.")
  writeLines(files, file.path(output, "test_file_inventory.txt"))
  summaries <- vector("list", length(files))
  warning_records <- vector("list", length(files))
  for (i in seq_along(files)) {
    case_output <- file.path(output, "unit", sub("\\.R$", "", basename(files[[i]])))
    dir.create(case_output, recursive = TRUE)
    child_env <- Sys.getenv()
    child_env[c("TMPDIR", "R_LIBS_USER", "LC_CTYPE")] <-
      c(alias, library_path, "en_US.UTF-8")
    command <- processx::run(
      file.path(R.home("bin"), "Rscript"),
      c("scripts/research/run_unit_file_worker_20261022.R", root, files[[i]], case_output),
      wd = root, env = child_env,
      stdout = file.path(case_output, "stdout.log"),
      stderr = file.path(case_output, "stderr.log"), error_on_status = FALSE,
      echo = FALSE, timeout = Inf
    )
    manifest_path <- file.path(case_output, "manifest.json")
    if (file.exists(manifest_path)) {
      observed <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
      summaries[[i]] <- data.frame(
        file = basename(files[[i]]), process_status = command$status,
        pass = isTRUE(observed$assertion_pass) && command$status == 0L,
        assertion_pass = isTRUE(observed$assertion_pass) && command$status == 0L,
        strict_pass = isTRUE(observed$strict_pass) && command$status == 0L,
        passed = observed$passed, failed = observed$failed, errors = observed$errors,
        warnings = observed$warnings, skipped = observed$skipped,
        elapsed_sec = observed$elapsed_sec, result_available = TRUE,
        stringsAsFactors = FALSE
      )
      warning_records[[i]] <- utils::read.csv(file.path(case_output, "warning_details.csv"),
                                               stringsAsFactors = FALSE)
    } else {
      summaries[[i]] <- data.frame(
        file = basename(files[[i]]), process_status = command$status, pass = FALSE,
        assertion_pass = FALSE, strict_pass = FALSE,
        passed = NA_integer_, failed = NA_integer_, errors = NA_integer_,
        warnings = NA_integer_, skipped = NA_integer_, elapsed_sec = command$timeout,
        result_available = FALSE, stringsAsFactors = FALSE
      )
    }
    current <- do.call(rbind, summaries[seq_len(i)])
    utils::write.csv(current, file.path(output, "unit_file_summary.csv"), row.names = FALSE)
    message(sprintf("[%d/%d] %s: pass=%s, assertions=%s, failures=%s, warnings=%s, skips=%s",
                    i, length(files), basename(files[[i]]), current$pass[[i]],
                    current$passed[[i]], current$failed[[i]], current$warnings[[i]], current$skipped[[i]]))
  }
  unit_pass <- all(current$pass) && nrow(current) == length(files)
  unit_manifest <- list(
    pass = unit_pass, assertion_pass = unit_pass, strict_pass = all(current$strict_pass),
    files_scheduled = length(files), files_attempted = nrow(current),
    files_with_results = sum(current$result_available), observed_assertions = sum(current$passed, na.rm = TRUE),
    failures = sum(current$failed, na.rm = TRUE), errors = sum(current$errors, na.rm = TRUE),
    warnings = sum(current$warnings, na.rm = TRUE), warning_count = sum(current$warnings, na.rm = TRUE),
    skipped = sum(current$skipped, na.rm = TRUE)
  )
  utils::write.csv(do.call(rbind, warning_records), file.path(output, "unit_warning_details.csv"),
                   row.names = FALSE)
  jsonlite::write_json(unit_manifest, file.path(output, "unit_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
}
if (phase %in% c("check", "all")) {
  check_output <- file.path(output, "package")
  dir.create(check_output)
  child_env <- Sys.getenv()
  child_env[c("TMPDIR", "R_LIBS_USER", "LC_CTYPE")] <-
    c(alias, library_path, "en_US.UTF-8")
  invocation <- processx::run(
    file.path(R.home("bin"), "Rscript"),
    c("scripts/research/run_package_check_worker_20261022.R", root, check_output, alias),
    wd = root, env = child_env,
    stdout = file.path(check_output, "stdout.log"),
    stderr = file.path(check_output, "stderr.log"), error_on_status = FALSE,
    echo = FALSE, timeout = Inf
  )
  check_manifest_path <- file.path(check_output, "manifest.json")
  check_pass <- invocation$status == 0L && file.exists(check_manifest_path) &&
    isTRUE(jsonlite::read_json(check_manifest_path)$pass)
}
jsonlite::write_json(list(
  pass = (is.na(unit_pass) || unit_pass) && (is.na(check_pass) || check_pass),
  unit_pass = unit_pass, package_check_pass = check_pass,
  assertion_pass = if (phase %in% c("unit", "all")) unit_pass else NA,
  implementation_commit = commit, phase = phase,
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")

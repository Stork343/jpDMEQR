args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("需要提供 unit 目录、package 目录、理论审阅文件、initial 检查文件和输出文件。")
root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
source(file.path(root, "scripts", "application", "_paired_method_gate_20261022.R"))
read_record <- function(path) jsonlite::read_json(normalizePath(path, mustWork = TRUE), simplifyVector = TRUE)
unit_dir <- normalizePath(args[[1L]], mustWork = TRUE)
package_dir <- normalizePath(args[[2L]], mustWork = TRUE)
review_path <- normalizePath(args[[3L]], mustWork = TRUE)
initial_path <- normalizePath(args[[4L]], mustWork = TRUE)
output <- file.path(root, args[[5L]])
if (file.exists(output)) stop("应用检查文件已经存在，需要使用新的输出文件。")
unit <- read_record(file.path(unit_dir, "unit_manifest.json"))
unit_request <- read_record(file.path(unit_dir, "request.json"))
package <- read_record(file.path(package_dir, "package", "manifest.json"))
package_request <- read_record(file.path(package_dir, "request.json"))
review <- read_record(review_path)
initial <- read_record(initial_path)
commit <- current_commit_v2(root)
if (!isTRUE(unit$assertion_pass) || unit$failures != 0L || unit$errors != 0L ||
    unit$skipped != 0L || unit$files_scheduled != unit$files_attempted ||
    unit$files_scheduled != unit$files_with_results || unit$observed_assertions < 1L ||
    !isTRUE(package$pass) || length(package$errors) || length(package$warnings) ||
    !identical(unit_request$implementation_commit, commit) ||
    !identical(package_request$implementation_commit, commit)) {
  stop("当前 commit 的真实 unit 或 package 检查未通过。")
}
numeric_paths <- sort(unique(c(
  application_paired_gate_dependency_paths_20261022(root),
  file.path("R", list.files(file.path(root, "R"), "^v2_benchmark_adapters_v2_part[0-9]+[.]R$"))
)))
dirty <- system2("git", c("diff", "--quiet", "--ignore-cr-at-eol", "HEAD", "--", shQuote(numeric_paths)))
tracked <- system2("git", c("ls-files", "--error-unmatch", shQuote(numeric_paths)),
                     stdout = FALSE, stderr = FALSE)
if (dirty != 0L || tracked != 0L) stop("应用数值源文件与当前 commit 不一致。")
if (!isTRUE(review$pass) || !identical(review$method_version, "profile_pair_20261022") ||
    !identical(review$score_definition, "profile_bandwidth_pair_2_half_minus_full") ||
    !length(review$source_paths) || !length(review$checks) || is.null(review$source_commit)) {
  stop("需要完整的理论审阅记录。")
}
review_paths <- as.character(review$source_paths)
if (any(grepl("^/|(^|/)[.][.](/|$)", review_paths)) ||
    any(!file.exists(file.path(root, review_paths)))) stop("理论审阅源文件路径无效。")
ancestor <- system2("git", c("merge-base", "--is-ancestor", review$source_commit, commit))
changed <- system2("git", c("diff", "--quiet", "--ignore-cr-at-eol", review$source_commit,
                            commit, "--", shQuote(review_paths)))
review_dirty <- system2("git", c("diff", "--quiet", "--ignore-cr-at-eol", "HEAD", "--", shQuote(review_paths)))
if (ancestor != 0L || changed != 0L || review_dirty != 0L) stop("理论审阅内容与当前源文件不一致。")
relative <- function(path) {
  normalized <- normalizePath(path, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(normalized, prefix)) stop("检查记录必须保存在当前项目。")
  substring(normalized, nchar(prefix) + 1L)
}
geometry <- geometry_manifest_path_v2(root, "profile_pair_20261022")
preflight <- file.path(profile_pair_gate_directory_v2(root), "application_method_preflight_manifest.json")
packages <- c("jpDMEQR", "CVXR", "quantreg", "lqmm", "geeVerse", "highs")
checks <- list(
  unit_checks_pass = unit$assertion_pass, package_checks_pass = package$pass,
  initial_checks_pass = initial$pass, initial_acceptance_manifest = relative(initial_path),
  generic_target_review_pass = review$pass, commit_sha = commit,
  method_version = "profile_pair_20261022", score_definition = "profile_bandwidth_pair_2_half_minus_full",
  geometry_manifest = relative(geometry), application_method_preflight_manifest = relative(preflight),
  numeric_dependency_paths = numeric_paths,
  benchmark_dependency_paths = numeric_paths[grepl("v2_benchmark_adapters", numeric_paths)],
  solver_configuration = solver_configuration_v2(),
  initial_configuration = list(clinical_acceleration = Sys.getenv("JPDMEQR_INITIAL_ACCELERATION", unset = "none")),
  package_versions = setNames(lapply(packages, function(name) as.character(packageVersion(name))), packages),
  unit_request = relative(file.path(unit_dir, "request.json")),
  unit_manifest = relative(file.path(unit_dir, "unit_manifest.json")),
  unit_warning_count = unit$warning_count,
  unit_warning_details = relative(file.path(unit_dir, "unit_warning_details.csv")),
  package_manifest = relative(file.path(package_dir, "package", "manifest.json")),
  theory_review = relative(review_path), completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write_json_v2(checks, output)
validated <- validate_application_paired_gate_20261022(root, output, inference = TRUE)
saveRDS(validated, paste0(output, ".validated.rds"))
message(output)

root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "application", "_paired_method_gate_20261022.R"))
args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args)) args[1] else "application_certificate_sources_20261023_v1"
output <- file.path(root, "results", "scratch", run_id)
if (dir.exists(output)) stop("The application source inventory directory already exists.")
directory <- file.path(root, "results", "scratch", "repository_package_checks_20261022", "development_full_03")
inventory <- readLines(file.path(directory, "test_file_inventory.txt"), warn = FALSE)
inventory <- inventory[nzchar(inventory)]
results <- read.csv(file.path(directory, "unit_file_summary.csv"), stringsAsFactors = FALSE)
record <- jsonlite::read_json(file.path(directory, "manifest.json"), simplifyVector = TRUE)
unit <- jsonlite::read_json(file.path(directory, "unit_manifest.json"), simplifyVector = TRUE)
if (length(inventory) != 65L || anyDuplicated(inventory) || nrow(results) != 65L ||
    !setequal(basename(inventory), results$file)) {
  stop("The actual 65-file repository unit inventory is incomplete.")
}
numerical <- application_paired_gate_dependency_paths_20261022(root)
paths <- sort(unique(c(inventory, numerical)))
if (any(!file.exists(file.path(root, paths)))) stop("An actual recorded test or application source is missing.")
git_lines <- function(arguments) {
  value <- system2("git", c("-C", shQuote(root), arguments), stdout = TRUE)
  if (!is.null(attr(value, "status")) && attr(value, "status") != 0L) stop("Cannot read the actual repository source state.")
  value
}
current <- git_lines(c("rev-parse", "HEAD"))
tracked <- git_lines(c("ls-tree", "-r", "--name-only", current, "--", shQuote(paths)))
recorded <- git_lines(c("ls-tree", "-r", "--name-only", record$implementation_commit, "--", shQuote(paths)))
changed <- git_lines(c("diff", "--name-only", "--ignore-cr-at-eol", record$implementation_commit, current,
                       "--", shQuote(paths)))
working <- git_lines(c("diff", "--name-only", "--ignore-cr-at-eol", "HEAD", "--", shQuote(paths)))
material_changes <- function(candidates, from, to = NULL) {
  candidates[vapply(candidates, function(path) {
    system2("git", c("-C", shQuote(root), "diff", "--quiet", "--ignore-cr-at-eol", from,
                      if (is.null(to)) character() else to, "--", shQuote(path))) != 0L
  }, logical(1))]
}
changed <- material_changes(changed, record$implementation_commit, current)
working <- material_changes(working, "HEAD")
info <- file.info(file.path(root, paths))
table <- data.frame(
  path = paths, checked_test_file = paths %in% inventory, application_numerical_source = paths %in% numerical,
  present_in_unit_record_commit = paths %in% recorded, tracked_in_current_HEAD = paths %in% tracked,
  committed_content_differs_from_unit_commit = paths %in% changed,
  working_content_differs_from_current_HEAD = paths %in% working,
  bytes = info$size, modified_utc = format(info$mtime, tz = "UTC", usetz = TRUE),
  unit_record_commit = record$implementation_commit, actual_source_commit = current,
  stringsAsFactors = FALSE
)
table$unchanged_committed_source <- table$present_in_unit_record_commit & table$tracked_in_current_HEAD &
  !table$committed_content_differs_from_unit_commit & !table$working_content_differs_from_current_HEAD
dir.create(output, recursive = TRUE)
write.csv(table, file.path(output, "all_observed_source_files.csv"), row.names = FALSE)
write.csv(table[table$checked_test_file, ], file.path(output, "repository_65_checked_test_files.csv"), row.names = FALSE)
write.csv(table[table$application_numerical_source, ], file.path(output, "application_numeric_source_files.csv"), row.names = FALSE)
write.csv(table[table$unchanged_committed_source, ], file.path(output, "unchanged_committed_source_files.csv"), row.names = FALSE)
jsonlite::write_json(list(
  actual_source_commit = current, unit_record_commit = record$implementation_commit,
  unit_record_directory = directory, test_files = length(inventory),
  application_numerical_source_files = length(numerical),
  observed_unit_record = unit, application_comparator_source = "R/v2_application_v2_part14.R",
  application_comparator_ids = names(yaml::read_yaml(file.path(root, "config", "application_execution_20261022.yml"))$comparators),
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "source_inventory_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

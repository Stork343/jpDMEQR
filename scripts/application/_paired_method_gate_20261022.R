application_paired_gate_dependency_paths_20261022 <- function(root) {
  profile <- list.files(file.path(root, "R"), "^v2_profile_v2_part[0-9]+[.]R$")
  application <- list.files(file.path(root, "R"), "^v2_application_v2_part[0-9]+[.]R$")
  native <- list.files(file.path(root, "src"), "[.](cpp|h|hpp)$", recursive = TRUE)
  sort(c(
    file.path("R", profile), file.path("R", application), file.path("src", native),
    "src/Makevars", "src/Makevars.win", "DESCRIPTION", "NAMESPACE",
    "scripts/00_source_v2.R", "scripts/application/_paired_method_gate_20261022.R",
    "scripts/application/03_fit_GSE65391.R", "scripts/application/03_fit_lowdim_GSE65391.R",
    "scripts/application/12_validate_GSE65391.R",
    "scripts/application/09_fit_exploratory_pair_GSE65391.R",
    "config/application.yml", "config/application_execution_20261022.yml",
    "config/application_curvature_20261022.yml", "config/application_matched_prospective_20261022.yml",
    "config/module_scoring_20261022.yml", "config/module_signatures.csv", "config/module_signatures.csv.freeze.csv",
    "config/module_score_members_20261022.csv", "config/module_score_members_20261022.csv.freeze.csv"
  ))
}

validate_application_paired_gate_20261022 <- function(
    root, checks_path, inference = TRUE, method_version = "profile_pair_20261022") {
  if (!file.exists(checks_path)) stop("APPLICATION_PAIRED_GATE_FAILED: actual application checks are missing.")
  checks <- jsonlite::read_json(checks_path, simplifyVector = TRUE)
  required <- c("unit_checks_pass", "package_checks_pass", "commit_sha", "method_version",
                "score_definition", "geometry_manifest", "application_method_preflight_manifest",
                "numeric_dependency_paths", "solver_configuration", "initial_configuration",
                "initial_acceptance_manifest")
  score_definition <- "profile_bandwidth_pair_2_half_minus_full"
  if (!all(required %in% names(checks)) ||
      !isTRUE(checks$unit_checks_pass) || !isTRUE(checks$package_checks_pass) ||
      !identical(checks$method_version, method_version) ||
      !identical(checks$score_definition, score_definition) ||
      (isTRUE(inference) && !isTRUE(checks$generic_target_review_pass))) {
    stop("APPLICATION_PAIRED_GATE_FAILED: matching actual package, unit and target checks are required.")
  }
  geometry_path <- file.path(root, checks$geometry_manifest)
  preflight_path <- file.path(root, checks$application_method_preflight_manifest)
  geometry <- jsonlite::read_json(geometry_path, simplifyVector = TRUE)
  preflight <- jsonlite::read_json(preflight_path, simplifyVector = TRUE)
  initial_configuration <- list(clinical_acceleration = Sys.getenv("JPDMEQR_INITIAL_ACCELERATION", unset = "none"))
  if (!initial_configuration$clinical_acceleration %in% c("none", "clinical_kkt") ||
      !identical(checks$initial_configuration, initial_configuration) ||
      !identical(preflight$initial_configuration, initial_configuration)) {
    stop("APPLICATION_PAIRED_GATE_FAILED: the actual initial-estimator configuration must match its certificate.")
  }
  for (record in list(geometry, preflight)) {
    if (!isTRUE(record$pass) || !identical(record$commit_sha, checks$commit_sha) ||
        !identical(record$method_version, method_version) ||
        !identical(record$score_definition, score_definition) ||
        !configuration_matches_v2(record$solver_configuration)) {
      stop("APPLICATION_PAIRED_GATE_FAILED: method records must certify the same source and production configuration.")
    }
  }
  if (!isTRUE(geometry$strict) || !isTRUE(geometry$dantzig_checked) ||
      !isTRUE(geometry$dantzig_pass) || !configuration_matches_v2(checks$solver_configuration)) {
    stop("APPLICATION_PAIRED_GATE_FAILED: strict paired geometry and the matching numerical configuration are required.")
  }
  paths <- as.character(checks$numeric_dependency_paths)
  mandatory <- application_paired_gate_dependency_paths_20261022(root)
  if (!length(paths) || anyNA(paths) || anyDuplicated(paths) ||
      any(grepl("^/|(^|/)[.][.](/|$)", paths)) || !all(mandatory %in% paths) ||
      any(!file.exists(file.path(root, paths)))) {
    stop("APPLICATION_PAIRED_GATE_FAILED: the complete declared numerical source paths are required.")
  }
  if (!is.null(checks$benchmark_dependency_paths) &&
      !all(as.character(checks$benchmark_dependency_paths) %in% paths)) {
    stop("APPLICATION_PAIRED_GATE_FAILED: a declared application comparator source was omitted.")
  }
  current <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
  certifying <- system2("git", c("-C", shQuote(root), "rev-parse", shQuote(paste0(checks$commit_sha, "^{commit}"))),
                         stdout = TRUE)
  tracked <- system2("git", c("-C", shQuote(root), "ls-files", "--error-unmatch", shQuote(paths)),
                     stdout = FALSE, stderr = FALSE)
  ancestry <- system2("git", c("-C", shQuote(root), "merge-base", "--is-ancestor", certifying, current))
  changed <- system2("git", c("-C", shQuote(root), "diff", "--quiet", "--ignore-cr-at-eol", certifying,
                               current, "--", shQuote(paths)))
  dirty <- system2("git", c("-C", shQuote(root), "diff", "--quiet", "--ignore-cr-at-eol", "HEAD", "--", shQuote(paths)))
  if (length(current) != 1L || length(certifying) != 1L || tracked != 0L ||
      ancestry != 0L || changed != 0L || dirty != 0L) {
    stop("APPLICATION_PAIRED_GATE_FAILED: the evaluated numerical source must remain unchanged in the current descendant.")
  }
  initial_path <- file.path(root, checks$initial_acceptance_manifest)
  initial <- jsonlite::read_json(initial_path, simplifyVector = TRUE)
  initial_paths <- as.character(initial$source_paths)
  required_initial <- c(
    file.path("R", list.files(file.path(root, "R"), "^v2_profile_v2_part[0-9]+[.]R$")),
    file.path("src", list.files(file.path(root, "src"), "[.](cpp|h|hpp)$", recursive = TRUE)),
    "src/Makevars", "src/Makevars.win", "DESCRIPTION", "NAMESPACE",
    file.path("R", sprintf("v2_application_v2_part%02d.R", c(3L, 5L, 6L, 9L, 10L, 13L)))
  )
  if (!isTRUE(initial$pass) || !identical(initial$method_version, "profile_pair_20261022_initial_l1") ||
      !identical(initial$initial_configuration, initial_configuration) ||
      !configuration_matches_v2(initial$solver_configuration) ||
      !length(initial_paths) || anyNA(initial_paths) || anyDuplicated(initial_paths) ||
      any(grepl("^/|(^|/)[.][.](/|$)", initial_paths)) || !all(required_initial %in% initial_paths)) {
    stop("APPLICATION_PAIRED_GATE_FAILED: actual initial-estimator numerical acceptance is incomplete.")
  }
  required_checks <- c("cache_and_full_design", "fixed_penalty_stationarity_and_gene_continuation",
                       "frozen_grid_global_kkt_and_timing", "matched_source_scope")
  if (!all(required_checks %in% names(initial$checks))) {
    stop("APPLICATION_PAIRED_GATE_FAILED: required initial numerical check records are missing.")
  }
  for (name in required_checks) {
    record <- initial$checks[[name]]
    table <- utils::read.csv(record$path, stringsAsFactors = FALSE)
    fields <- c("nb", "passed", "failed", "error", "skipped", "warning")
    if (!isTRUE(record$pass) || !nrow(table) || !all(fields %in% names(table)) ||
        anyNA(table[fields]) || any(table$failed != 0 | table$error | table$skipped | table$warning != 0) ||
        sum(table$nb) != sum(table$passed) || sum(table$passed) != record$assertions) {
      stop("APPLICATION_PAIRED_GATE_FAILED: an actual initial numerical check is incomplete: ", name)
    }
  }
  initial_ancestor <- system2("git", c("-C", shQuote(root), "merge-base", "--is-ancestor",
                                        initial$actual_source_commit, certifying))
  initial_changed <- system2("git", c("-C", shQuote(root), "diff", "--quiet", "--ignore-cr-at-eol",
                                       initial$actual_source_commit, "--", shQuote(initial_paths)))
  initial_tracked <- system2("git", c("-C", shQuote(root), "ls-files", "--error-unmatch", shQuote(initial_paths)),
                              stdout = FALSE, stderr = FALSE)
  if (initial_ancestor != 0L || initial_changed != 0L || initial_tracked != 0L) {
    stop("APPLICATION_PAIRED_GATE_FAILED: the actual initial-estimator source differs from its accepted ancestor.")
  }
  if (!is.null(checks$package_versions)) {
    for (name in names(checks$package_versions)) {
      if (!identical(as.character(utils::packageVersion(name)), as.character(checks$package_versions[[name]]))) {
        stop("APPLICATION_PAIRED_GATE_FAILED: the certified package version changed: ", name)
      }
    }
  }
  if (identical(solver_configuration_v2()$cluster_solver, "compiled") &&
      !is.loaded("jpdmeqr_cluster_profiles", PACKAGE = "jpDMEQR")) {
    stop("APPLICATION_PAIRED_GATE_FAILED: the certified compiled cluster routine is unavailable.")
  }
  list(
    commit_sha = current, actual_source_commit = current,
    numeric_certificate_commit = certifying, certifying_commit = certifying,
    method_version = method_version, score_definition = score_definition,
    solver_configuration = solver_configuration_v2(), numeric_dependency_paths = paths,
    initial_configuration = initial_configuration,
    initial_acceptance_manifest = normalizePath(initial_path),
    checks_path = normalizePath(checks_path), geometry_manifest = normalizePath(geometry_path),
    application_method_preflight_manifest = normalizePath(preflight_path)
  )
}

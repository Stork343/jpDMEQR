application_population_reference_equal_v2 <- function(left, right) {
  isTRUE(all.equal(left, right, tolerance = 0, check.attributes = TRUE, check.environment = FALSE))
}

application_numeric_dependency_paths_v2 <- function(root) {
  profile <- list.files(file.path(root, "R"), "^v2_profile_v2_part[0-9]+[.]R$")
  application <- setdiff(list.files(file.path(root, "R"), "^v2_application_v2_part[0-9]+[.]R$"),
                         c("v2_application_v2_part12.R", "v2_application_v2_part15.R"))
  native <- list.files(file.path(root, "src"), "[.](cpp|h|hpp)$", recursive = TRUE)
  c(file.path("R", profile), file.path("R", application), file.path("src", native),
    "scripts/00_source_v2.R", "scripts/application/13_run_matched_GSE65391.R",
    "scripts/application/17_launch_matched_study_GSE65391.R",
    "scripts/research/check_application_numeric_acceptance_20261022.R",
    "tests/testthat/test-profile-unpenalised-newton.R", "tests/testthat/test-profile-proximal-scaling.R",
    "tests/testthat/test-application-initial-estimator.R", "tests/testthat/test-application-compact-workflow.R",
    "config/application_matched_prospective_20261022.yml", "config/module_scoring_20261022.yml",
    "config/module_signatures.csv", "config/module_signatures.csv.freeze.csv",
    "config/module_score_members_20261022.csv", "config/module_score_members_20261022.csv.freeze.csv")
}

application_numeric_ast_exclusions_v2 <- function() {
  list(
    "R/v2_profile_v2_part03.R" = "solve_dantzig_row_v2",
    "R/v2_profile_v2_part04.R" = "dantzig_osqp_solve_chain_v2",
    "R/v2_application_v2_part10.R" = "fit_application_profile_initial_v2",
    "R/v2_application_v2_part13.R" = c("application_initial_fit_record_v2", "fit_application_prediction_candidate_v2")
  )
}

application_numeric_git_expressions_v2 <- function(root, commit, path) {
  lines <- system2("git", c("-C", shQuote(root), "show", shQuote(paste0(commit, ":", path))),
                    stdout = TRUE)
  if (!is.null(attr(lines, "status")) && attr(lines, "status") != 0L) {
    stop("The declared numerical source cannot be read from its actual commit.")
  }
  parse(text = sub("\r$", "", lines), keep.source = FALSE)
}

application_numeric_ast_equal_v2 <- function(reference, candidate, excluded_functions) {
  retained <- function(expressions) {
    declarations <- vapply(as.list(expressions), function(item) {
      if (is.call(item) && length(item) == 3L &&
          identical(item[[1]], as.name("<-")) && is.symbol(item[[2]]) &&
          is.call(item[[3]]) && identical(item[[3]][[1]], as.name("function"))) {
        as.character(item[[2]])
      } else ""
    }, character(1))
    if (any(vapply(excluded_functions, function(name) sum(declarations == name) != 1L, logical(1)))) {
      stop("Every excluded reference API must have exactly one top-level function declaration.")
    }
    as.expression(as.list(expressions)[!declarations %in% excluded_functions])
  }
  identical(retained(reference), retained(candidate))
}

application_numeric_sources_match_v2 <- function(root, from_commit, to_commit = NULL, paths) {
  exclusions <- application_numeric_ast_exclusions_v2()
  whole_paths <- setdiff(paths, names(exclusions))
  target <- if (is.null(to_commit)) character() else to_commit
  changed <- system2("git", c("-C", shQuote(root), "diff", "--quiet", "--ignore-cr-at-eol",
                               from_commit, target, "--", shQuote(whole_paths)))
  if (changed != 0L) return(FALSE)
  for (path in intersect(paths, names(exclusions))) {
    reference <- application_numeric_git_expressions_v2(root, from_commit, path)
    candidate <- if (is.null(to_commit)) parse(file.path(root, path), keep.source = FALSE) else
      application_numeric_git_expressions_v2(root, to_commit, path)
    if (!application_numeric_ast_equal_v2(reference, candidate, exclusions[[path]])) return(FALSE)
  }
  TRUE
}

validate_application_numeric_source_identity_v2 <- function(root, build_commit = NULL) {
  current <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
  if (is.null(build_commit) || !nzchar(build_commit)) build_commit <- current
  anchor <- system2("git", c("-C", shQuote(root), "rev-parse", shQuote(paste0(build_commit, "^{commit}"))),
                     stdout = TRUE)
  paths <- application_numeric_dependency_paths_v2(root)
  tracked <- system2("git", c("-C", shQuote(root), "ls-files", "--error-unmatch", shQuote(paths)),
                     stdout = FALSE, stderr = FALSE)
  ancestry <- system2("git", c("-C", shQuote(root), "merge-base", "--is-ancestor", anchor, current))
  unchanged <- application_numeric_sources_match_v2(root, anchor, current, paths)
  clean <- application_numeric_sources_match_v2(root, current, paths = paths)
  if (length(current) != 1L || length(anchor) != 1L ||
      tracked != 0L || ancestry != 0L || !unchanged || !clean) {
    stop("APPLICATION_NUMERIC_SOURCE_FAILED: committed numerical dependencies must match the declared ancestor build.")
  }
  path <- file.path(root, "results", "scratch", "numerical_solver_acceptance_20261022_v3", "test_summary.csv")
  checks <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("nb", "passed", "failed", "error", "skipped", "warning")
  if (!all(required %in% names(checks)) || anyNA(checks[required]) ||
      sum(checks$nb) != 148L || sum(checks$passed) != 148L ||
      any(checks$failed != 0 | checks$error | checks$skipped | checks$warning != 0)) {
    stop("APPLICATION_NUMERIC_SOURCE_FAILED: all actual production numerical checks are required.")
  }
  configuration_path <- file.path(dirname(path), "solver_configuration.json")
  recorded <- jsonlite::read_json(configuration_path, simplifyVector = TRUE)
  if (!configuration_matches_v2(recorded)) {
    stop("APPLICATION_NUMERIC_SOURCE_FAILED: the production numerical configuration changed.")
  }
  manifest <- jsonlite::read_json(file.path(dirname(path), "validation_manifest.json"), simplifyVector = TRUE)
  validation_ancestor <- system2("git", c("-C", shQuote(root), "merge-base", "--is-ancestor",
                                          manifest$implementation_commit, anchor))
  validation_unchanged <- application_numeric_sources_match_v2(root, manifest$implementation_commit, anchor, paths)
  metadata_paths <- "R/v2_application_v2_part15.R"
  validated_paths <- setdiff(as.character(manifest$dependency_paths), metadata_paths)
  if (validation_ancestor != 0L || !validation_unchanged || !identical(validated_paths, paths)) {
    stop("APPLICATION_NUMERIC_SOURCE_FAILED: evaluated source dependencies differ from the declared build.")
  }
  list(application_build_commit = anchor, actual_source_commit = current,
       dependency_paths = paths, numerical_checks_source = path, numerical_assertions = sum(checks$passed),
       numerical_configuration_source = configuration_path,
       numerical_validation_commit = manifest$implementation_commit,
       comparison_rule = "complete_top_level_AST_except_declared_unused_reference_API",
       excluded_reference_functions = application_numeric_ast_exclusions_v2(),
       exclusion_reason = list(
         solve_dantzig_row_v2 = "matched lowdim inference uses the direct paired Hessian inverse",
         dantzig_osqp_solve_chain_v2 = "matched lowdim inference uses the direct paired Hessian inverse",
         fit_application_profile_initial_v2 = "matched lowdim workflow has zero L1 penalties and never calls the prediction or exploratory initial estimator",
         application_initial_fit_record_v2 = "matched lowdim workflow never calls prediction candidate recording",
         fit_application_prediction_candidate_v2 = "matched lowdim workflow never calls prediction candidate fitting"
       ),
       metadata_source_paths = metadata_paths,
       metadata_exclusion_reason = "source identity and population validation functions do not compute beta or profile score")
}

validate_application_population_implementation_v2 <- function(root, core_build_commit = "5c47a62") {
  summaries <- c(
    compiled = file.path(root, "results", "scratch", "compiled_profiles_20261022", "test_summary.csv"),
    population = file.path(root, "results", "scratch", "profile_pair_population_performance_20261022", "Eigen", "equivalence_tests.csv")
  )
  expected <- c(compiled = 112L, population = 319L)
  counts <- lapply(names(summaries), function(name) {
    table <- utils::read.csv(summaries[[name]], stringsAsFactors = FALSE)
    required <- c("nb", "failed", "error", "skipped", "warning", "passed")
    if (!nrow(table) || !all(required %in% names(table)) || anyNA(table[required]) ||
        any(table$failed != 0) || any(table$error) || any(table$skipped) || any(table$warning != 0) ||
        sum(table$passed) != expected[[name]] || sum(table$nb) != expected[[name]]) {
      stop("POPULATION_CORE_GATE_FAILED: actual native and population checks must be complete.")
    }
    data.frame(check = name, assertions = sum(table$passed), source = summaries[[name]],
               source_modified_utc = format(file.info(summaries[[name]])$mtime, tz = "UTC", usetz = TRUE))
  })
  symbols <- c("jpdmeqr_cluster_profiles", "jpDMEQR_population_crossprod")
  if (!all(vapply(symbols, is.loaded, logical(1), PACKAGE = "jpDMEQR"))) {
    stop("POPULATION_CORE_GATE_FAILED: both registered installed native entry points are required.")
  }
  committed <- c("config/application_matched_prospective_20261022.yml",
                  "config/module_signatures.csv", "config/module_signatures.csv.freeze.csv",
                  "config/module_score_members_20261022.csv", "config/module_score_members_20261022.csv.freeze.csv",
                  "config/module_scoring_20261022.yml", "R/v2_application_v2_part11.R")
  for (path in committed) {
    tracked <- system2("git", c("-C", shQuote(root), "ls-files", "--error-unmatch", shQuote(path)),
                       stdout = FALSE, stderr = FALSE)
    unchanged <- system2("git", c("-C", shQuote(root), "diff", "--quiet", "HEAD", "--", shQuote(path)))
    if (tracked != 0L || unchanged != 0L) stop("Matched population definitions must be committed and unchanged.")
  }
  app_commit <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
  core_commit <- system2("git", c("-C", shQuote(root), "rev-parse", shQuote(paste0(core_build_commit, "^{commit}"))),
                         stdout = TRUE)
  ancestry <- system2("git", c("-C", shQuote(root), "merge-base", "--is-ancestor", core_commit, app_commit))
  if (length(app_commit) != 1L || length(core_commit) != 1L || ancestry != 0L) {
    stop("The actual application construction must descend from the declared core commit.")
  }
  config <- yaml::read_yaml(file.path(root, "config", "application_matched_prospective_20261022.yml"))
  if (config$n_subjects != 129L || config$replications != 200L ||
      config$population_targets$clusters_per_repeat != 100000L ||
      config$population_targets$independent_repeats != 4L ||
      !setequal(config$models, c("A", "B")) || !setequal(config$quantiles, c(.25, .50, .75)) ||
      !setequal(config$error_families, c("normal", "t3")) ||
      !setequal(config$random_effect_dimensions, c(1L, 2L))) {
    stop("The complete reviewed matched-study matrix is required.")
  }
  directory <- file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1")
  for (model in config$models) for (q in config$random_effect_dimensions) {
    population <- readRDS(file.path(directory, sprintf("matched_population_model_%s_q%d.rds", model, q)))
    d <- if (model == "A") 16L else 17L
    if (!identical(population$model, model) || population$q != q ||
        nrow(population$X) != 867L || ncol(population$X) != d ||
        nrow(population$Z) != 867L || ncol(population$Z) != q ||
        length(population$rows) != 129L || sum(lengths(population$rows)) != 867L ||
        qr(population$X)$rank != d || any(!is.finite(population$X)) || any(!is.finite(population$Z)) ||
        !identical(names(population$beta0), colnames(population$X)) ||
        any(!is.finite(population$beta0)) || nrow(population$preprocessing$module_audit) != 5L ||
        !all(population$preprocessing$module_audit$eligible)) {
      stop("The fixed complete population reference is invalid.")
    }
  }
  dll <- getLoadedDLLs()[["jpDMEQR"]][["path"]]
  list(commit_sha = app_commit, core_build_commit = core_commit,
       method_version = "profile_pair_20261022", scope = "population_construction_and_simulation_only",
       checks = do.call(rbind, counts), installed_dll = dll,
       installed_dll_modified_utc = format(file.info(dll)$mtime, tz = "UTC", usetz = TRUE),
       registered_native_symbols = symbols,
       installed_package_version = as.character(utils::packageVersion("jpDMEQR")))
}

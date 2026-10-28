root <- normalizePath(".", mustWork = TRUE)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled", JPDMEQR_UNPENALISED_OPTIMIZER = "newton_armijo",
           JPDMEQR_PROXIMAL_SCALING = "cluster_rms")
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args)) args[1] else "matched_source_identity_20261022_v1"
output <- file.path(root, "results", "scratch", run_id)
if (dir.exists(output)) stop("The source identity validation directory already exists.")
dir.create(output, recursive = TRUE)
change_body <- function(expressions, function_name) {
  result <- as.list(expressions)
  declarations <- vapply(result, function(item) {
    if (is.call(item) && length(item) == 3L && identical(item[[1]], as.name("<-")) &&
        is.symbol(item[[2]]) && is.call(item[[3]]) && identical(item[[3]][[1]], as.name("function"))) {
      as.character(item[[2]])
    } else ""
  }, character(1))
  index <- which(declarations == function_name)
  if (length(index) != 1L) stop("The actual source function must have one declaration.")
  # 仅修改内存中的语法树副本，检验来源比较；这些函数副本不会执行。
  result[[index]][[3]][[3]] <- quote(stop("source_identity_required_expression_change"))
  as.expression(result)
}
reporter <- testthat::ListReporter$new()
testthat::with_reporter(reporter, {
testthat::test_that("真实来源AST逐项比较并拒绝required expression变化", {
  exclusions <- application_numeric_ast_exclusions_v2()
  required <- c("R/v2_profile_v2_part03.R" = "fit_profile_lasso_v2",
                "R/v2_profile_v2_part04.R" = "cluster_sandwich_coordinate_v2")
  for (path in names(exclusions)) {
    original <- application_numeric_git_expressions_v2(root, "078559a", path)
    current <- parse(file.path(root, path), keep.source = FALSE)
    testthat::expect_true(application_numeric_ast_equal_v2(original, current, exclusions[[path]]))
    for (function_name in exclusions[[path]]) {
      testthat::expect_true(application_numeric_ast_equal_v2(
        original, change_body(current, function_name), exclusions[[path]]
      ))
    }
    altered_required <- if (path %in% names(required)) {
      change_body(current, required[[path]])
    } else {
      as.expression(c(as.list(current), as.list(parse(
        file.path(root, "R", "v2_profile_v2_part03.R"), keep.source = FALSE
      ))[1L]))
    }
    testthat::expect_false(application_numeric_ast_equal_v2(
      original, altered_required, exclusions[[path]]
    ))
  }
  identity <- validate_application_numeric_source_identity_v2(root, "078559a")
  testthat::expect_identical(identity$application_build_commit, "078559aa453622373dc765b4209510e519318f85")
  testthat::expect_identical(identity$actual_source_commit, system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
  testthat::expect_false("R/v2_application_v2_part15.R" %in% identity$dependency_paths)
  pending <- "lowdim_profile_pair_workflow_inference_v2"
  visited <- character()
  while (length(pending)) {
    name <- pending[1L]
    pending <- pending[-1L]
    if (name %in% visited) next
    visited <- c(visited, name)
    calls <- codetools::findGlobals(get(name, envir = .GlobalEnv), merge = FALSE)$functions
    available <- calls[grepl("_v2$", calls) & vapply(calls, exists, logical(1),
                                                          envir = .GlobalEnv, mode = "function")]
    pending <- unique(c(pending, setdiff(available, visited)))
  }
  unused <- c("fit_application_profile_initial_v2", "application_initial_fit_record_v2",
              "fit_application_prediction_candidate_v2")
  testthat::expect_length(intersect(visited, unused), 0L)
  saveRDS(list(entry = "lowdim_profile_pair_workflow_inference_v2", reachable_functions = visited,
               excluded_application_functions = unused), file.path(output, "actual_matched_call_graph.rds"))
  saveRDS(identity, file.path(output, "actual_source_identity.rds"))
})
})
results <- as.data.frame(reporter$get_results())
write.csv(results[c("test", "nb", "failed", "skipped", "error", "warning", "passed")],
           file.path(output, "test_summary.csv"), row.names = FALSE)
if (any(results$failed > 0 | results$error | results$skipped | results$warning > 0)) {
  stop("The actual numerical source identity checks failed.")
}

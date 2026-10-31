source("research/score_validity_20261026/load_joint.R")
source("research/milk_application_20261031/nodewise.R")
source("research/milk_application_20261031/weighted_solver.R")
source("research/soynam_application_20261031/refine_basis.R")
milk_reference_fit_small <- soy_refined_basis_fit
source("research/soynam_application_20261031/data.R")
load_soy_solver <- function() {
  lock <- filelock::lock("results/scratch/soynam_20261031/pair_qp_compile.lock", timeout = 60000)
  if (is.null(lock)) stop("The paired solver compilation lock was unavailable.")
  on.exit(filelock::unlock(lock), add = TRUE)
  Rcpp::sourceCpp("research/soynam_application_20261031/pair_qp.cpp",
    cacheDir = "results/scratch/soynam_20261031/pair_qp_cache", rebuild = FALSE, env = .GlobalEnv)
}
load_soy_solver()
rm(load_soy_solver)
soy_original_initial_preparation <- pc_prepare_initial
soy_original_path <- pc_orthogonal_path
pc_prepare_initial <- function(y, X, Z, cluster_id, tau, Lambda, penalty_factor,
                               max_iterations = 30000L, fit_always_included = FALSE) {
  soy_original_initial_preparation(y, X, Z, cluster_id, tau, Lambda, penalty_factor,
    max_iterations = max_iterations, fit_always_included = FALSE)
}
pc_orthogonal_path <- function(prepared, tau, sizes, gap_tolerance = 1e-8,
                               max_iterations = 30000L, solver = soy_pair_qp_cpp) {
  soy_original_path(prepared, tau, sizes, gap_tolerance, max_iterations, solver = solver)
}
soy_original_initial <- sb_initial
soy_original_prepare <- ds_prepare
ds_fit_small <- milk_fit_small
sb_lasso_supports <- milk_curvature_supports

soy_progress <- function(stage, detail = list()) {
  record <- c(list(stage = stage, time = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)), detail)
  jsonlite::write_json(record, file.path(soy_context$directory, "progress.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA)
  cat(record$time, soy_context$case_id, stage, "\n")
  flush.console()
}

soy_fold_number <- function(subjects) {
  fold <- which(vapply(soy_context$selection_subjects, identical, logical(1), y = subjects))
  if (length(fold) != 1L) stop("An unregistered selection sample was requested.")
  fold
}

sb_initial <- function(data, tau, Lambda, penalty_factor, budget_multiplier = 1) {
  fold <- soy_fold_number(names(data$rows))
  key <- sprintf("%s_tau%02d_lambda%03d_fold%d", soy_context$cohort,
    as.integer(round(100 * tau)), as.integer(round(100 * Lambda[1L])), fold)
  path <- file.path(soy_context$cache, paste0("outcome_", key, ".rds"))
  handle <- filelock::lock(paste0(path, ".lock"), timeout = 600000)
  if (is.null(handle)) stop("The outcome-search cache lock was unavailable.")
  on.exit(filelock::unlock(handle), add = TRUE)
  request <- list(subjects = names(data$rows), response = data$y, columns = colnames(data$X),
    penalty_factor = penalty_factor, Lambda = Lambda, budget_multiplier = budget_multiplier)
  if (file.exists(path)) {
    saved <- readRDS(path)
    if (!identical(saved$request, request)) stop("The outcome-search cache has a different input.")
    initial <- saved$initial
    inner <- pc_folds(data$rows, 3L)
    initial$training <- lapply(1:3, function(f) sb_subset(data, unlist(data$rows[inner != f], use.names = FALSE)))
    initial$validation <- lapply(1:3, function(f) sb_subset(data, unlist(data$rows[inner == f], use.names = FALSE)))
    soy_progress("outcome_search_cached", list(fold = fold, additions = initial$size))
    return(initial)
  }
  soy_progress("outcome_search", list(fold = fold, clusters = data$n, variables = ncol(data$X)))
  initial <- soy_original_initial(data, tau, Lambda, penalty_factor, budget_multiplier)
  compact <- initial
  compact$training <- compact$validation <- compact$paths <- NULL
  certificates <- lapply(initial$paths, function(path) path$diagnostics)
  saveRDS(list(request = request, initial = compact, certificates = certificates), path)
  soy_progress("outcome_search_completed", list(fold = fold, additions = initial$size))
  initial
}

pc_nodewise <- function(factor, rows, coordinate) {
  fold <- soy_fold_number(names(rows))
  key <- sprintf("%s_coordinate%d_fold%d", soy_context$cohort, coordinate, fold)
  path <- file.path(soy_context$cache, paste0("design_", key, ".rds"))
  handle <- filelock::lock(paste0(path, ".lock"), timeout = 600000)
  if (is.null(handle)) stop("The design-regression cache lock was unavailable.")
  on.exit(filelock::unlock(handle), add = TRUE)
  request <- list(subjects = names(rows), coordinate = coordinate, columns = colnames(factor),
    dimensions = dim(factor), rms = sqrt(colSums(factor^2) / length(rows)))
  if (file.exists(path)) {
    saved <- readRDS(path)
    if (!identical(saved$request, request)) stop("The conditional-design cache has a different input.")
    soy_progress("conditional_design_cached", list(fold = fold, coordinate = coordinate))
    return(saved$fit)
  }
  soy_progress("conditional_design", list(fold = fold, coordinate = coordinate))
  fit <- milk_nodewise_compute(factor, rows, coordinate)
  saveRDS(list(request = request, fit = fit), path)
  soy_progress("conditional_design_completed", list(fold = fold, coordinate = coordinate, kkt = fit$kkt))
  fit
}

ds_prepare <- function(y, X, Z, cluster_id, tau, coordinates, Lambda, penalty_factor, h,
                       include_design = TRUE, budget_multiplier = 1) {
  prepared <- soy_original_prepare(y, X, Z, cluster_id, tau, coordinates, Lambda,
    penalty_factor, h, include_design, budget_multiplier)
  saveRDS(prepared, file.path(soy_context$directory, "prepared.rds"), compress = FALSE)
  soy_progress("paired_coefficient_fits")
  prepared
}

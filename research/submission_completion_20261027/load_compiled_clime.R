compile_clime_driver <- function() {
  prefix <- normalizePath("results/scratch/dependencies/highs-native", mustWork = TRUE)
  cache <- "results/scratch/sqr_acceleration_20261030/compiled"
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(file.path(cache, "compile.lock"), timeout = 60000)
  if (is.null(lock)) stop("The compiled CLIME driver is being prepared by another process.")
  on.exit(filelock::unlock(lock), add = TRUE)
  previous <- Sys.getenv(c("PKG_CPPFLAGS", "PKG_LIBS"), unset = NA_character_)
  on.exit({
    for (name in names(previous)) {
      if (is.na(previous[[name]])) Sys.unsetenv(name) else
        do.call(Sys.setenv, setNames(list(previous[[name]]), name))
    }
  }, add = TRUE)
  Sys.setenv(PKG_CPPFLAGS = paste0("-I", shQuote(file.path(prefix, "include/highs"))),
    PKG_LIBS = paste0("-L", shQuote(file.path(prefix, "lib")), " -lhighs -Wl,-rpath,", shQuote(file.path(prefix, "lib"))))
  Rcpp::sourceCpp("research/submission_completion_20261027/clime_highs.cpp",
    cacheDir = cache, rebuild = FALSE, showOutput = FALSE, env = .GlobalEnv)
}
compile_clime_driver()
rm(compile_clime_driver)
compiled_identity <- submission_clime_highs_cpp(diag(2), .5, FALSE, NULL, Inf)
if (!isTRUE(compiled_identity$feasible) || compiled_identity$highs_version != "1.15.1")
  stop("The compiled CLIME dependency does not match its declared version.")
rm(compiled_identity)
source("research/submission_completion_20261027/clime_path_cache.R")

submission_clime_compiled_grid <- function(H, gamma_grid, solver_opts = list()) {
  p <- nrow(H)
  gamma_grid <- sort(unique(gamma_grid))
  grids <- lapply(seq_len(p), function(k) list(mu = gamma_grid,
    omega = vector("list", length(gamma_grid)), solver_records = vector("list", length(gamma_grid)),
    residual = rep(NA_real_, length(gamma_grid)), l1_norm = rep(NA_real_, length(gamma_grid)),
    status = rep("not_required_after_global_infeasibility", length(gamma_grid)),
    feasible = rep(NA, length(gamma_grid)), solver = "compiled_HiGHS_complete_primal_dual", coordinate = k))
  records <- list()
  cutoff <- NA_real_
  models <- 0L
  path_cache <- solver_opts$path_cache
  if (is.null(path_cache)) path_cache <- new.env(parent = emptyenv())
  for (index in order(gamma_grid, decreasing = TRUE)) {
    prepared <- solver_opts$computed_native
    if (!is.null(prepared) && (length(gamma_grid) != 1L || prepared$gamma != gamma_grid[index] ||
        !identical(prepared$H, H) || isTRUE(prepared$value$pruned)))
      stop("The completed native matrix belongs to a different CLIME problem.")
    if (is.null(prepared) && isTRUE(solver_opts$native_guard) && gamma_grid[index] < .25)
      submission_clime_path_guard(H, gamma_grid[index], path_cache)
    value <- if (is.null(prepared)) submission_clime_highs_cpp(H, gamma_grid[index], isTRUE(solver_opts$presolve), NULL, Inf) else prepared$value
    models <- models + sum(value$programmes)
    reached <- which(value$status %in% c("optimal", "infeasible"))
    for (k in reached) {
      feasible <- value$status[k] == "optimal"
      grids[[k]]$status[index] <- value$status[k]
      grids[[k]]$feasible[index] <- feasible
      grids[[k]]$omega[[index]] <- if (feasible) value$omega_raw[k, ] else rep(NA_real_, p)
      grids[[k]]$residual[index] <- value$residual[k]
      grids[[k]]$l1_norm[index] <- value$objective[k]
      grids[[k]]$solver_records[[index]] <- list(dual = value$dual[k, ],
        dual_gap = value$dual_gap[k], dual_residual = value$dual_residual[k],
        programmes = value$programmes[k], simplex_iterations = value$simplex_iterations[k],
        active_rows = value$rows[k], active_columns = value$columns[k])
    }
    records[[length(records) + 1L]] <- data.frame(gamma = gamma_grid[index], coordinate = reached,
      status = value$status[reached], residual = value$residual[reached],
      feasible = value$status[reached] == "optimal", dual_gap = value$dual_gap[reached],
      programmes = value$programmes[reached], simplex_iterations = value$simplex_iterations[reached])
    if (!value$feasible) {
      cutoff <- gamma_grid[index]
      break
    }
  }
  list(grids = grids, attempts = do.call(rbind, records), global_infeasible_at = cutoff,
    model_constructions = models, basis_policy = "compiled_complete_primal_and_dual_pricing",
    solver_control = list(engine = "HiGHS", version = value$highs_version, threads = 1L,
      native_path_guard = isTRUE(solver_opts$native_guard), native_path_maximum_pivots = "4p",
      primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8),
    feasibility_certificate = list(verified = FALSE, definition = "per_row_primal_dual_certificates"),
    source_precision_definition = "common_gamma_full_matrix_CLIME")
}

submission_clime_bounded_cv <- function(H, gamma, validation_H, loss_limit, solver_opts, coordinates) {
  H <- (H + t(H)) / 2
  cache <- solver_opts$path_cache
  if (!is.environment(cache)) stop("The bounded CLIME calculation requires its fixed-fold cache.")
  if (is.null(cache$hessian)) cache$hessian <- H else
    if (!identical(cache$hessian, H)) stop("The CLIME certificate belongs to a different Hessian.")
  proof <- cache$infeasibility
  if (!is.null(proof) && gamma < proof$gamma_bound - 1e-7)
    return(list(pruned = FALSE, feasible = FALSE, record = list(
      feasible = FALSE, gamma = gamma, status = "infeasible_by_full_matrix_certificate",
      feasibility_certificate = proof)))
  if (isTRUE(solver_opts$native_guard) && gamma < .1)
    submission_clime_path_guard(H, gamma, cache)
  value <- submission_clime_highs_cpp(H, gamma, isTRUE(solver_opts$presolve), validation_H, loss_limit)
  if (isTRUE(value$pruned)) return(list(pruned = TRUE, lower_bound = value$loss_lower_bound,
    record = list(gamma = gamma, status = "certified_nonnegative_column_loss_bound",
      loss_lower_bound = value$loss_lower_bound, loss_limit = loss_limit,
      evaluated_columns = value$evaluated_columns, solved_rows = value$solved_rows,
      column_losses = value$column_losses, column_lower_bounds = value$column_lower_bounds,
      primal_residuals = value$residual,
      dual_residuals = value$dual_residual, dual_gaps = value$dual_gap,
      solver = "compiled_HiGHS_complete_primal_dual")))
  if (!isTRUE(value$feasible)) {
    coordinate <- value$failed_coordinate
    ray <- value$dual[coordinate, ]
    cache$infeasibility <- list(coordinate = coordinate, ray = ray,
      gamma_bound = abs(ray[coordinate]) / sum(abs(ray)),
      null_residual = value$ray_residual, verified_gamma = gamma,
      definition = "complete_Hessian_Farkas_certificate")
  } else if (isTRUE(solver_opts$native_guard) && gamma < .25) {
    submission_clime_path_guard(H, gamma, cache)
  }
  solver_opts$computed_native <- list(H = H, gamma = gamma, value = value)
  precision <- solve_clime_matrix_iid_v2(H, gamma, solver_opts = solver_opts, coordinates = coordinates)
  record <- sqr_cv_precision_record_v2(precision, coordinates)
  list(pruned = FALSE, feasible = isTRUE(precision$feasible), record = record,
    loss = if (isTRUE(precision$feasible)) submission_sqr_sparse_loss_cpp(validation_H, precision$omega_symmetric) else NA_real_)
}

dantzig_iid_solve_matrix_grid_v2 <- local({
  reference <- dantzig_iid_solve_matrix_grid_v2
  function(H, gamma_grid, solver_opts = list()) {
    if (identical(solver_opts$algorithm, "compiled_certified"))
      return(submission_clime_compiled_grid(H, gamma_grid, solver_opts))
    reference(H, gamma_grid, solver_opts)
  }
})

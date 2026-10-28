profile_pair_geometry_grid_20261022 <- function() {
  seeds <- 20261022L + seq_len(20L)
  grid <- expand.grid(seed_index = seq_along(seeds), tau_index = 1:4)
  grid <- grid[order(grid$seed_index, grid$tau_index), ]
  grid$seed <- seeds[grid$seed_index]
  grid$tau <- c(0.1, 0.25, 0.5, 0.75)[grid$tau_index]
  grid$q <- 1L + (floor((grid$seed_index - 1L) / 4L) + grid$tau_index) %% 2L
  grid$h_multiplier <- c(0.75, 1, 1.25)[1L + (grid$seed_index + 2L * grid$tau_index) %% 3L]
  grid$balanced <- (grid$seed_index + grid$tau_index) %% 3L == 0L
  grid$n_clusters <- 8L + grid$seed_index %% 3L
  grid$non_diagonal_lambda <- grid$q == 2L & grid$seed_index %% 2L == 1L
  grid$case_id <- sprintf("seed%02d_tau%02d", grid$seed_index, grid$tau_index)
  rownames(grid) <- NULL
  grid
}

profile_pair_geometry_limits_20261022 <- function() {
  list(gradient = 2e-4, hessian = 3e-3, identity = 1e-8,
       nuisance = 1e-7, dantzig = 1e-6, score_combination = 1e-12,
       sandwich = 1e-10, assembly = 1e-8)
}

profile_pair_minimax_check_20261022 <- function(H, coordinate) {
  variable <- CVXR::Variable(nrow(H))
  e <- replace(numeric(nrow(H)), coordinate, 1)
  problem <- CVXR::Problem(CVXR::Minimize(CVXR::norm_inf(H %*% variable - e)))
  result <- CVXR::solve(problem, solver = "ECOS", verbose = FALSE,
                        feastol = 1e-9, reltol = 1e-9, abstol = 1e-9)
  if (!identical(result$status, "optimal") || !is.finite(result$value)) {
    stop("Independent minimax feasibility LP did not return optimal status.")
  }
  omega <- as.numeric(result$getValue(variable))
  list(minimum_defect = result$value,
       achieved_defect = max(abs(as.numeric(H %*% omega) - e)),
       solver = "ECOS", status = result$status, omega = omega)
}

profile_pair_geometry_check_20261022 <- function(dat, beta, tau, h, Lambda,
                                                coordinate = 1L,
                                                check_assembly = FALSE) {
  nuisance <- list(grad_tol = 1e-7, reltol = 1e-12, maxit = 600L)
  evaluate <- function(value, assembly = NULL) {
    profile_bandwidth_pair_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, value, tau, h,
      lambda_gamma = Lambda, nuisance_control = nuisance,
      identity_check = TRUE, assembly = assembly, need_fold_sums = TRUE
    )
  }
  comp <- evaluate(beta)
  full <- profile_components_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, h, Lambda,
    identity_check = TRUE, nuisance_control = nuisance
  )
  half <- profile_components_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, h / 2, Lambda,
    identity_check = TRUE, nuisance_control = nuisance
  )
  combination_score <- max(abs(comp$score - (2 * half$score - full$score)))
  combination_H <- max(abs(comp$hessian - (2 * half$hessian - full$hessian)))
  steps <- c(1e-5, 2e-5)
  gradient_error <- hessian_error <- numeric(length(steps))
  fd_records <- vector("list", length(steps))
  for (a in seq_along(steps)) {
    delta <- steps[a]
    gradient <- numeric(ncol(dat$X))
    Hessian <- matrix(0, ncol(dat$X), ncol(dat$X))
    for (k in seq_len(ncol(dat$X))) {
      increment <- replace(numeric(ncol(dat$X)), k, delta)
      plus <- evaluate(beta + increment)
      minus <- evaluate(beta - increment)
      gradient[k] <- (plus$objective - minus$objective) / (2 * delta)
      Hessian[, k] <- (plus$score - minus$score) / (2 * delta)
    }
    gradient_error[a] <- max(abs(gradient - comp$score))
    hessian_error[a] <- max(abs(Hessian - comp$hessian))
    fd_records[[a]] <- list(step = delta, gradient = gradient, Hessian = Hessian)
  }
  precision <- profile_pair_solve_row_v2(comp$hessian, coordinate, 0.05)
  statuses <- gsub("_", "", tolower(precision$attempts$status))
  minimax <- NULL
  if (isTRUE(precision$feasible)) {
    dantzig_consistent <- all(statuses %in% c("solved", "almostsolved", "optimal")) &&
      precision$residual <= precision$mu + 1e-6
  } else {
    minimax <- profile_pair_minimax_check_20261022(comp$hessian, coordinate)
    dantzig_consistent <- all(statuses %in% c("primalinfeasible", "almostprimalinfeasible",
                                             "unbounded", "unboundedorinfeasible",
                                             "infeasiblebymonotonicity")) &&
      minimax$minimum_defect > 0.05 + 1e-6
  }
  canonical <- replace(numeric(ncol(dat$X)), coordinate, 1)
  interval <- profile_bandwidth_pair_wald_v2(
    comp, beta, matrix(canonical, ncol = 1L), coordinate
  )
  projected <- as.numeric(comp$cluster_scores %*% canonical)
  matched_variance <- mean((projected - mean(projected))^2) / comp$n_clusters
  dantzig_sandwich_error <- NA_real_
  if (isTRUE(precision$feasible)) {
    dantzig_interval <- profile_bandwidth_pair_wald_v2(
      comp, beta, matrix(precision$omega, ncol = 1L), coordinate
    )
    projection <- as.numeric(comp$cluster_scores %*% precision$omega)
    dantzig_sandwich_error <- abs(dantzig_interval$se^2 -
      mean((projection - mean(projection))^2) / comp$n_clusters)
  }
  moments <- profile_pair_moments_v2(comp)
  parity <- c(full = 0, folds = 0, identity = 0)
  if (check_assembly) {
    sequential <- evaluate(beta, "sequential")
    batch <- evaluate(beta, "batch")
    parity <- c(
      full = max(abs(sequential$hessian - batch$hessian)),
      folds = max(vapply(seq_len(4), function(f) {
        max(abs(sequential$hessian_fold_sums[[f]] - batch$hessian_fold_sums[[f]]))
      }, numeric(1))),
      identity = max(abs(sequential$hessian_identity - batch$hessian_identity))
    )
  }
  metrics <- data.frame(
    gradient_error_small_step = gradient_error[1],
    gradient_error_large_step = gradient_error[2],
    gradient_error = max(gradient_error),
    hessian_error_small_step = hessian_error[1],
    hessian_error_large_step = hessian_error[2],
    hessian_error = max(hessian_error),
    score_combination_error = combination_score,
    hessian_combination_error = combination_H,
    profile_identity_error = max(comp$profile_identity_error,
                                 full$profile_identity_error, half$profile_identity_error),
    max_nuisance_gradient = max(comp$max_nuisance_gradient,
                                full$max_nuisance_gradient, half$max_nuisance_gradient),
    dantzig_mu = 0.05, dantzig_residual = precision$residual,
    dantzig_solver = precision$solver, dantzig_status = precision$solver_status,
    dantzig_attempt_status = paste(precision$attempts$status, collapse = "|"),
    dantzig_feasible = precision$feasible, dantzig_status_consistent = dantzig_consistent,
    independent_minimum_defect = if (is.null(minimax)) NA_real_ else minimax$minimum_defect,
    independent_feasibility_solver = if (is.null(minimax)) "not_required" else minimax$solver,
    numerical_rank = qr(comp$hessian, tol = 1e-12)$rank,
    smallest_singular_value = min(svd(comp$hessian, nu = 0L, nv = 0L)$d),
    sandwich_identity_error = abs(interval$se^2 - matched_variance),
    variance_identity_direction = "prespecified_canonical_coordinate",
    dantzig_sandwich_error = dantzig_sandwich_error,
    covariance_identity_error = moments$covariance_identity_error,
    assembly_hessian_error = parity[["full"]], assembly_fold_error = parity[["folds"]],
    assembly_identity_error = parity[["identity"]],
    min_paired_eigenvalue = comp$hessian_min_eigenvalue,
    step_small = steps[1], step_large = steps[2], assembly_checked = check_assembly
  )
  limits <- profile_pair_geometry_limits_20261022()
  metrics$pass <- with(metrics,
    gradient_error < limits$gradient && hessian_error < limits$hessian &&
      score_combination_error < limits$score_combination &&
      hessian_combination_error < limits$score_combination &&
      profile_identity_error < limits$identity && max_nuisance_gradient < limits$nuisance &&
      dantzig_status_consistent &&
      sandwich_identity_error < limits$sandwich && covariance_identity_error < limits$identity &&
      (is.na(dantzig_sandwich_error) || dantzig_sandwich_error < limits$sandwich) &&
      max(parity) < limits$assembly
  )
  list(metrics = metrics, data = dat, beta = beta, Lambda = Lambda,
       components = comp, finite_differences = fd_records,
       precision = precision, independent_minimax = minimax, interval = interval)
}

profile_pair_audit_configuration_20261022 <- function(mode = "reference") {
  if (!mode %in% c("reference", "production")) stop("Unknown audit configuration.")
  production <- identical(mode, "production")
  list(dantzig_solver = if (production) "dual_warm" else "clarabel_native",
       cluster_solver = if (production) "compiled" else "newton",
       assembly = if (production) "batch" else "sequential",
       cv_workers = 1L, hessian_refresh = if (production) 20L else 1L,
       shared_first_stage = production,
       matrix_backend = if (production) "Eigen_single_thread" else "R_BLAS",
       identity_detail = if (production) "aggregate" else "cluster",
       proximal_scaling = if (production) "cluster_rms" else "none",
       unpenalised_optimizer = if (production) "newton_armijo" else "proximal_gradient")
}

run_profile_pair_geometry_audit_20261022 <- function(root, run_id, geometry_source_run = NULL,
                                                    mode = "reference") {
  root <- normalizePath(root, mustWork = TRUE)
  configuration <- profile_pair_audit_configuration_20261022(mode)
  environment_fields <- c(dantzig_solver = "JPDMEQR_DANTZIG_SOLVER",
                          cluster_solver = "JPDMEQR_CLUSTER_SOLVER",
                          assembly = "JPDMEQR_ASSEMBLY", cv_workers = "JPDMEQR_CV_CORES",
                          hessian_refresh = "JPDMEQR_HESSIAN_REFRESH",
                          shared_first_stage = "JPDMEQR_SHARED_FIRST_STAGE",
                          matrix_backend = "JPDMEQR_MATRIX_BACKEND",
                          identity_detail = "JPDMEQR_IDENTITY_DETAIL",
                          proximal_scaling = "JPDMEQR_PROXIMAL_SCALING",
                          unpenalised_optimizer = "JPDMEQR_UNPENALISED_OPTIMIZER")
  do.call(Sys.setenv, setNames(lapply(configuration, as.character), environment_fields))
  local_library <- file.path(root, "results", "scratch", "profile_pair_package_20261022", "library")
  if (dir.exists(local_library)) .libPaths(c(normalizePath(local_library), .libPaths()))
  if (identical(mode, "production") && !requireNamespace("jpDMEQR", quietly = TRUE)) {
    stop("The compiled production audit requires the installed jpDMEQR package.")
  }
  out <- file.path(root, "results", "scratch", "profile_pair_geometry_development_20261022", run_id)
  if (dir.exists(out)) stop("The geometry run directory already exists: ", out)
  dir.create(file.path(out, "objects"), recursive = TRUE)
  source(file.path(root, "scripts", "00_source_v2.R"), local = TRUE)
  module_environment <- environment(profile_pair_geometry_check_20261022)
  source_v2_module(root, "profile_v2", module_environment)
  source_v2_module(root, "simulation_v2", module_environment)
  grid <- profile_pair_geometry_grid_20261022()
  source_out <- NULL
  if (!is.null(geometry_source_run)) {
    source_out <- file.path(root, "results", "scratch", "profile_pair_geometry_development_20261022",
                            geometry_source_run)
    source_manifest <- jsonlite::read_json(file.path(source_out, "manifest.json"), simplifyVector = TRUE)
    source_grid <- utils::read.csv(file.path(source_out, "coverage_matrix.csv"))
    if (!isTRUE(source_manifest$pass) || source_manifest$case_count != 81L ||
        !configuration_matches_v2(source_manifest$solver_configuration, configuration) ||
        !isTRUE(all.equal(source_grid, grid, check.attributes = FALSE, tolerance = 0))) {
      stop("The geometry source run does not match the current audit design.")
    }
  }
  utils::write.csv(grid, file.path(out, "coverage_matrix.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  stage <- "geometry"
  current <- NULL
  globalCallingHandlers(error = function(error) {
    saveRDS(list(condition = error, stage = stage, context = current),
            file.path(out, "failure.rds"))
    writeLines(c(paste("stage:", stage), conditionMessage(error)), file.path(out, "FAILED.txt"))
  })
  options(error = function() quit(save = "no", status = 1L))
  geometry_records <- list()
  fit_records <- list()
  for (r in seq_len(nrow(grid))) {
    row <- grid[r, ]
    stage <- paste("geometry", row$case_id)
    m <- if (row$balanced) rep(4L, row$n_clusters) else
      rep(c(2L, 3L, 5L, 4L), length.out = row$n_clusters)
    dat <- generate_profile_qr_data_v2(
      n_clusters = row$n_clusters, p = 3L, s = 1L, tau = row$tau,
      q = row$q, m_multiset = m, rho_x = 0.3, signal = 0.75, seed = row$seed
    )
    Lambda <- if (row$q == 1L) matrix(0.8, 1L) else if (row$non_diagonal_lambda) {
      matrix(c(1, 0.25, 0.25, 0.8), 2L)
    } else diag(c(1, 0.8))
    beta <- as.numeric(dat$beta0) + c(0.08, -0.04, 0.03)
    h <- row$h_multiplier * row$n_clusters^(-3 / 10)
    current <- list(row = row, data = dat, beta = beta, h = h, Lambda = Lambda)
    result <- if (!is.null(source_out)) readRDS(
      file.path(source_out, "objects", paste0(row$case_id, ".rds"))
    ) else profile_pair_geometry_check_20261022(
      dat, beta, row$tau, h, Lambda,
      coordinate = 1L + (row$seed_index + row$tau_index) %% 3L,
      check_assembly = row$tau_index == 1L
    )
    saveRDS(result, file.path(out, "objects", paste0(row$case_id, ".rds")))
    geometry_records[[r]] <- cbind(row, h = h, result$metrics)
    utils::write.csv(do.call(rbind, geometry_records), file.path(out, "geometry.csv"), row.names = FALSE)
    if (!result$metrics$pass) stop("Geometry thresholds failed at ", row$case_id, ".")
    if (r %% 8L == 0L) message("Geometry cases completed: ", r, "/", nrow(grid))
  }

  # 此例的配对 Hessian 为负，检查仍使用微分与同 score 方差关系。
  stage <- "indefinite_paired_curvature"
  dat <- list(y = rep(c(0.75, 0.8, 0.85), each = 2L),
              X = matrix(1, 6L, 1L, dimnames = list(NULL, "intercept")),
              Z = matrix(1, 6L, 1L), cluster_id = rep(1:3, each = 2L))
  current <- list(data = dat, beta = 0, tau = 0.5, h = 0.4, Lambda = matrix(1, 1L))
  special <- if (!is.null(source_out)) readRDS(
    file.path(source_out, "objects", "indefinite_paired_curvature.rds")
  ) else profile_pair_geometry_check_20261022(dat, 0, 0.5, 0.4, matrix(1, 1L),
                                                 check_assembly = TRUE)
  saveRDS(special, file.path(out, "objects", "indefinite_paired_curvature.rds"))
  utils::write.csv(special$metrics, file.path(out, "indefinite_curvature.csv"), row.names = FALSE)
  if (!special$metrics$pass || special$metrics$min_paired_eigenvalue >= 0) {
    stop("The indefinite paired-curvature geometry check failed.")
  }
  jsonlite::write_json(list(
    development = TRUE, method_version = "profile_pair_20261022",
    score_definition = "profile_bandwidth_pair_2_half_minus_full",
    seed_count = 20L, case_count = length(geometry_records) + 1L,
    pass = TRUE, thresholds = profile_pair_geometry_limits_20261022(),
    solver_configuration = solver_configuration_v2()
  ), file.path(out, "geometry_completed_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = NA)

  for (i in seq_len(20L)) {
    row <- grid[grid$seed_index == i & grid$tau_index == 1L + (i - 1L) %% 4L, ]
    stage <- paste("fit_inference", row$case_id)
    dat <- generate_profile_qr_data_v2(
      n_clusters = 32L, p = 3L, s = 1L, tau = row$tau, q = row$q,
      m_multiset = rep(c(2L, 3L, 5L, 4L), 8L), rho_x = 0.3,
      signal = 0.75, seed = row$seed + 10000L
    )
    Lambda <- if (row$q == 1L) matrix(0.8, 1L) else matrix(c(1, 0.25, 0.25, 0.8), 2L)
    penalised <- floor((i - 1L) / 4L) %% 2L == 0L
    current <- list(row = row, data = dat, Lambda = Lambda, penalised = penalised)
    fit <- fit_profile_pair_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, row$tau,
      h_inf = row$h_multiplier * 32^(-3 / 10), lambda_gamma = Lambda,
      penalty_factor = rep(if (penalised) 1 else 0, 3L),
      target_coordinates = c("x00001", "x00003"),
      control = list(nuisance_control = list(grad_tol = 1e-8, reltol = 1e-12, maxit = 600L))
    )
    saveRDS(fit, file.path(out, "objects", paste0("fit_", row$case_id, ".rds")))
    current$fit <- fit
    mu_multipliers <- c(0.02, 0.05, 0.10, 0.25, 0.50, 1, 2, 4)
    mu_grid <- mu_multipliers * (sqrt(log(3) / (32 * fit$h_inf)) + fit$h_inf^2)
    interval <- infer_profile_pair_v2(fit, c("x00001", "x00003"), mu_grid)
    saveRDS(interval, file.path(out, "objects", paste0("inference_", row$case_id, ".rds")))
    checks <- vapply(seq_len(nrow(interval$table)), function(k) {
      omega <- interval$directions[[k]]$omega
      projection <- as.numeric(fit$components$cluster_scores %*% omega)
      abs(interval$table$se[k]^2 - mean((projection - mean(projection))^2) / 32L)
    }, numeric(1))
    inference_pass <- all(interval$table$feasible) &&
      all(interval$table$dantzig_residual <= interval$table$mu + 1e-6) &&
      max(checks) < 1e-10 && isTRUE(fit$converged) &&
      fit$post_refit_gradient_max <= 1e-7 &&
      fit$components$max_nuisance_gradient < 1e-7 &&
      fit$components$profile_identity_error < 1e-8
    fit_records[[i]] <- data.frame(
      case_id = row$case_id, seed = row$seed + 10000L, tau = row$tau, q = row$q,
      h_multiplier = row$h_multiplier, penalised = penalised,
      mu_candidates = length(mu_grid),
      kkt = fit$post_refit_gradient_max,
      max_nuisance_gradient = fit$components$max_nuisance_gradient,
      identity_error = fit$components$profile_identity_error,
      sandwich_error = max(checks),
      dantzig_max_excess = max(interval$table$dantzig_residual - interval$table$mu),
      pass = inference_pass
    )
    utils::write.csv(do.call(rbind, fit_records), file.path(out, "fit_inference.csv"), row.names = FALSE)
    if (!inference_pass) stop("Fitted paired inference acceptance failed at ", row$case_id, ".")
  }
  geometry <- do.call(rbind, geometry_records)
  fitting <- do.call(rbind, fit_records)
  complete_metrics <- rbind(geometry[, names(special$metrics)], special$metrics)
  commit <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
  manifest <- list(
    run_id = run_id, development = TRUE, final_authorisation = FALSE,
    method_version = "profile_pair_20261022",
    score_definition = "profile_bandwidth_pair_2_half_minus_full",
    implementation_commit = commit, commit_sha = commit,
    strict = length(unique(grid$seed)) >= 20L && all(complete_metrics$pass),
    dantzig_checked = nrow(complete_metrics) == 81L && nrow(fitting) == 20L,
    dantzig_pass = all(complete_metrics$dantzig_status_consistent) &&
      all(fitting$dantzig_max_excess <= 1e-6),
    assembly_modes_checked = sum(complete_metrics$assembly_checked) >= 21L,
    assembly_pass = all(complete_metrics$assembly_hessian_error < 1e-8) &&
      all(complete_metrics$assembly_fold_error < 1e-8) &&
      all(complete_metrics$assembly_identity_error < 1e-8),
    audit_configuration = mode,
    seed_count = length(unique(grid$seed)),
    case_count = nrow(geometry) + 1L, fit_inference_count = nrow(fitting),
    tau_values = sort(unique(grid$tau)), q_values = sort(unique(grid$q)),
    h_multipliers = sort(unique(grid$h_multiplier)),
    non_diagonal_lambda_cases = sum(grid$non_diagonal_lambda),
    unbalanced_cases = sum(!grid$balanced), balanced_cases = sum(grid$balanced),
    assembly_comparisons = sum(geometry$assembly_checked) + 1L,
    dantzig_feasible_cases = sum(complete_metrics$dantzig_feasible),
    dantzig_infeasible_status_verified_cases = sum(!complete_metrics$dantzig_feasible),
    dantzig_solver_used = sort(unique(complete_metrics$dantzig_solver)),
    geometry_source_run = geometry_source_run,
    fit_mu_multipliers = mu_multipliers,
    geometry_nuisance_control = list(grad_tol = 1e-7, reltol = 1e-12, maxit = 600L),
    fit_nuisance_control = list(grad_tol = 1e-8, reltol = 1e-12, maxit = 600L),
    indefinite_case_min_eigenvalue = special$metrics$min_paired_eigenvalue,
    thresholds = profile_pair_geometry_limits_20261022(),
    maxima = lapply(complete_metrics[c("gradient_error", "hessian_error", "profile_identity_error",
                                      "max_nuisance_gradient", "sandwich_identity_error",
                                      "covariance_identity_error", "assembly_fold_error")], max),
    solver_configuration = solver_configuration_v2(),
    pass = all(complete_metrics$pass) && all(fitting$pass),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  jsonlite::write_json(manifest, file.path(out, "manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA)
  writeLines(c(
    paste0("运行标识为 `", run_id, "`，`development=TRUE`。"),
    sprintf("%d 个 deterministic seeds 完成 %d 项 geometry 检查，另有 %d 项真实 fit/infer 检查。",
            manifest$seed_count, manifest$case_count, manifest$fit_inference_count),
    sprintf("gradient 最大误差为 %.6g，Hessian 最大误差为 %.6g。",
            max(complete_metrics$gradient_error), max(complete_metrics$hessian_error)),
    sprintf("profile identity 最大误差为 %.6g，nuisance gradient 最大值为 %.6g。",
            max(complete_metrics$profile_identity_error), max(complete_metrics$max_nuisance_gradient)),
    sprintf("负曲率用例的 paired Hessian 特征值为 %.6g，微分与 sandwich 检查通过。",
            special$metrics$min_paired_eigenvalue),
    "全部指标及逐例对象保存在本目录。本次运行未启动最终规模模拟。"
  ), file.path(out, "REPORT.md"))
  writeLines("COMPLETE", file.path(out, "COMPLETE.txt"))
  print(manifest)
  invisible(manifest)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  root <- if (length(args)) args[[1]] else "."
  run_id <- if (length(args) > 1L) args[[2]] else format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
  geometry_source_run <- if (length(args) > 2L) args[[3]] else NULL
  if (identical(geometry_source_run, "none")) geometry_source_run <- NULL
  mode <- if (length(args) > 3L) args[[4]] else "reference"
  run_profile_pair_geometry_audit_20261022(root, run_id, geometry_source_run, mode)
}

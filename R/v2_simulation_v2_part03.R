# Population target and population effective-direction utilities.

population_target_columns_v2 <- function(config, extra_null = 3L) {
  active <- seq_len(config$s)
  null <- seq.int(config$s + 1L, min(config$p, config$s + as.integer(extra_null)))
  sort(unique(c(active, null)))
}

fit_working_nuisance_design_v2 <- function(dat, fit_random_effects) {
  out <- dat
  if (identical(fit_random_effects, "intercept") && ncol(out$Z) > 1L) {
    out$Z <- out$Z[, 1L, drop = FALSE]
  } else if (!fit_random_effects %in% c("intercept", "intercept+slope")) {
    stop("Unknown fit_random_effects: ", fit_random_effects)
  }
  out
}

approximate_profile_target_v2 <- function(config,
                                          target_columns = population_target_columns_v2(config),
                                          n_population = 100000L,
                                          seed = 87001L,
                                          fit_control = list(
                                            max_iter = 500L,
                                            beta_tol = 1e-8,
                                            kkt_tol = 1e-7
                                          ),
                                          retain_fit = FALSE) {
  if (!exists("fit_profile_lasso_v2", mode = "function")) {
    stop("Source the profile-v2 module before approximating the profile target.")
  }
  target_columns <- sort(unique(as.integer(target_columns)))
  if (!length(target_columns) || min(target_columns) < 1L || max(target_columns) > config$p) {
    stop("Invalid target_columns.")
  }
  if (!all(seq_len(config$s) %in% target_columns)) {
    stop("Population target submodel must contain every active coordinate.")
  }

  cfg <- config
  cfg$n_clusters <- as.integer(n_population)
  cfg$p <- max(target_columns)
  dat <- simulate_from_config_v2(cfg, seed = seed)
  dat <- fit_working_nuisance_design_v2(dat, cfg$fit_random_effects)
  X_target <- dat$X[, target_columns, drop = FALSE]
  h_target <- min(
    config$h_multiplier * config$n_clusters^(-1 / 3),
    as.integer(n_population)^(-1 / 3)
  )

  fit <- fit_profile_lasso_v2(
    y = dat$y,
    X = X_target,
    Z = dat$Z,
    cluster_id = dat$cluster_id,
    tau = cfg$tau,
    h = h_target,
    lambda_beta = 0,
    lambda_gamma = cfg$lambda_gamma,
    penalty_factor = rep(0, length(target_columns)),
    control = fit_control
  )
  beta <- fit$beta
  names(beta) <- sprintf("x%05d", target_columns)
  list(
    beta_star_mc = beta,
    converged = fit$converged,
    kkt_residual = fit$kkt_residual,
    n_population = as.integer(n_population),
    seed = as.integer(seed),
    target_columns = target_columns,
    h_target = h_target,
    fit_random_effects = cfg$fit_random_effects,
    profile_identity_error = fit$components$profile_identity_error,
    max_nuisance_gradient = fit$components$max_nuisance_gradient,
    nuisance_converged = fit$components$nuisance_converged,
    fit = if (isTRUE(retain_fit)) fit else NULL
  )
}

summarise_profile_target_repeats_v2 <- function(fits, config, target_columns,
                                                implementation_commit = "unknown") {
  if (length(fits) < 2L) stop("At least two target repeats are required.")
  target_names <- sprintf("x%05d", target_columns)
  estimates <- do.call(rbind, lapply(fits, function(x) {
    out <- setNames(rep(NA_real_, length(target_names)), target_names)
    out[names(x$beta_star_mc)] <- x$beta_star_mc
    out
  }))
  if (any(!is.finite(estimates))) stop("Non-finite population target estimate.")
  mean_target <- colMeans(estimates)
  sd_by_coordinate <- apply(estimates, 2, stats::sd)
  se_by_coordinate <- sd_by_coordinate / sqrt(nrow(estimates))
  range_by_coordinate <- apply(estimates, 2, function(x) diff(range(x)))
  list(
    experiment_id = config$experiment_id,
    beta_star_mc = mean_target,
    target_mc_sd_by_coordinate = sd_by_coordinate,
    target_mc_se_by_coordinate = se_by_coordinate,
    target_mc_se = max(se_by_coordinate),
    estimates_by_repeat = estimates,
    target_columns = target_columns,
    n_population = unique(vapply(fits, `[[`, numeric(1), "n_population")),
    repeats = length(fits),
    between_run_max_difference = max(range_by_coordinate),
    between_run_range_by_coordinate = range_by_coordinate,
    all_converged = all(vapply(fits, `[[`, logical(1), "converged")),
    max_kkt_residual = max(vapply(fits, `[[`, numeric(1), "kkt_residual")),
    config = config,
    fits = fits,
    implementation_commit = implementation_commit,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}

stable_inverse_symmetric_v2 <- function(H, eigen_floor = 1e-8) {
  H <- (as.matrix(H) + t(as.matrix(H))) / 2
  eig <- eigen(H, symmetric = TRUE)
  floor_value <- max(eigen_floor, max(abs(eig$values)) * eigen_floor)
  values <- pmax(eig$values, floor_value)
  inv <- eig$vectors %*% diag(1 / values, nrow = length(values)) %*% t(eig$vectors)
  list(
    inverse = (inv + t(inv)) / 2,
    eigenvalues = eig$values,
    eigen_floor = floor_value,
    condition_number = max(values) / min(values)
  )
}

approximate_population_direction_v2 <- function(config,
                                                beta_target,
                                                target_columns = population_target_columns_v2(config),
                                                requested_coordinates = NULL,
                                                n_population = 100000L,
                                                seed = 97001L,
                                                nuisance_control = list(
                                                  reltol = 1e-11,
                                                  maxit = 500L,
                                                  grad_tol = 1e-8
                                                )) {
  target_columns <- sort(unique(as.integer(target_columns)))
  target_names <- sprintf("x%05d", target_columns)
  beta_target <- as.numeric(beta_target[target_names])
  names(beta_target) <- target_names
  if (any(!is.finite(beta_target))) stop("beta_target does not cover target_columns.")

  cfg <- config
  cfg$n_clusters <- as.integer(n_population)
  cfg$p <- max(target_columns)
  dat <- simulate_from_config_v2(cfg, seed = seed)
  dat <- fit_working_nuisance_design_v2(dat, cfg$fit_random_effects)
  X_target <- dat$X[, target_columns, drop = FALSE]
  h_analysis <- config$h_multiplier * config$n_clusters^(-1 / 3)
  comp <- profile_components_v2(
    y = dat$y,
    X = X_target,
    Z = dat$Z,
    cluster_id = dat$cluster_id,
    beta = beta_target,
    tau = config$tau,
    h = h_analysis,
    lambda_gamma = config$lambda_gamma,
    need_hessian = TRUE,
    nuisance_control = nuisance_control
  )
  H <- comp$hessian
  colnames(H) <- rownames(H) <- target_names
  inv <- stable_inverse_symmetric_v2(H)
  Omega <- inv$inverse
  colnames(Omega) <- rownames(Omega) <- target_names

  if (is.null(requested_coordinates)) {
    requested_coordinates <- target_names
  }
  requested_coordinates <- intersect(as.character(requested_coordinates), target_names)
  if (!length(requested_coordinates)) stop("No requested coordinate is in the population submodel.")
  directions <- lapply(requested_coordinates, function(nm) {
    k <- match(nm, target_names)
    omega <- Omega[, k]
    e <- numeric(length(target_names)); e[k] <- 1
    list(
      coordinate = nm,
      omega = setNames(as.numeric(omega), target_names),
      residual = max(abs(as.numeric(H %*% omega - e))),
      l1_norm = sum(abs(omega)),
      l2_norm = sqrt(sum(omega^2))
    )
  })
  names(directions) <- requested_coordinates

  list(
    experiment_id = config$experiment_id,
    H_population = H,
    Omega_population = Omega,
    directions = directions,
    target_columns = target_columns,
    target_names = target_names,
    beta_target = setNames(beta_target, target_names),
    h_analysis = h_analysis,
    n_population = as.integer(n_population),
    seed = as.integer(seed),
    condition_number = inv$condition_number,
    raw_eigenvalues = inv$eigenvalues,
    eigen_floor = inv$eigen_floor,
    profile_identity_error = comp$profile_identity_error,
    max_nuisance_gradient = comp$max_nuisance_gradient,
    nuisance_converged = comp$nuisance_converged,
    config = config,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}

target_asset_path_v2 <- function(root, experiment_id) {
  file.path(root, "results", "intermediate", "targets", paste0(experiment_id, ".rds"))
}

population_direction_asset_path_v2 <- function(root, experiment_id) {
  file.path(root, "results", "intermediate", "population_directions",
            paste0(experiment_id, ".rds"))
}

validate_profile_target_asset_v2 <- function(obj,
                                             min_population = 100000L,
                                             min_repeats = 4L,
                                             expected_commit = NULL,
                                             final = TRUE) {
  problems <- character()
  if (is.null(obj$beta_star_mc) || any(!is.finite(obj$beta_star_mc))) {
    problems <- c(problems, "target estimate is missing or non-finite")
  }
  if (isTRUE(final) && (is.null(obj$n_population) || min(obj$n_population) < min_population)) {
    problems <- c(problems, "population sample is below the final minimum")
  }
  if (isTRUE(final) && (is.null(obj$repeats) || obj$repeats < min_repeats)) {
    problems <- c(problems, "independent repeat count is below the final minimum")
  }
  if (!isTRUE(obj$all_converged)) problems <- c(problems, "one or more target fits did not converge")
  if (!is.null(expected_commit) && !identical(obj$implementation_commit, expected_commit)) {
    problems <- c(problems, "target asset implementation commit is stale")
  }
  list(valid = length(problems) == 0L, problems = problems)
}

submission_population_cache <- new.env(parent = emptyenv())

submission_registry <- function(phase = "main") {
  specification <- yaml::read_yaml("config/joint_profile_submission_20261027.yml")
  path <- switch(phase, main = specification$generating_registry, pilot = specification$pilot_registry,
                 preflight = specification$preflight_registry, stop("Unknown simulation phase."))
  registry <- read_simulation_registry_v2(path)
  mapping <- unlist(specification$method_mapping)
  registry$methods <- vapply(strsplit(registry$methods, "|", fixed = TRUE), function(methods) {
    matched <- match(methods, names(mapping))
    methods[!is.na(matched)] <- mapping[matched[!is.na(matched)]]
    paste(methods, collapse = "|")
  }, character(1))
  for (id in names(specification$target_mode_overrides)) registry$target_mode[registry$experiment_id == id] <- specification$target_mode_overrides[[id]]
  registry$budget_multiplier <- 1
  for (id in names(specification$budget_multipliers)) registry$budget_multiplier[registry$experiment_id == id] <- specification$budget_multipliers[[id]]
  registry$method_version <- "joint_profile_design_20261026"
  registry$seed_family <- registry$experiment_id
  registry$seed_family[registry$experiment_id %in% c("G06", "G07")] <- specification$parallel_seed_family
  registry
}

submission_target_path <- function(config) {
  specification <- yaml::read_yaml("config/joint_profile_submission_20261027.yml")
  override <- specification$population_targets[[config$experiment_id]]
  if (!is.null(override)) return(override)
  file.path(specification$population_targets$default, paste0(config$experiment_id, ".rds"))
}

submission_target <- function(config, dat) {
  if (config$target_mode == "structural") return(list(beta = dat$beta0,
    mcse = setNames(numeric(length(dat$beta0)), names(dat$beta0)), type = "full_structural_slope", path = NA_character_))
  path <- submission_target_path(config)
  if (!identical(submission_population_cache$target_path, path)) {
    submission_population_cache$target <- readRDS(path)
    submission_population_cache$target_path <- path
  }
  asset <- submission_population_cache$target
  if (asset$development || asset$n_population < 100000L || asset$repeats < 4L || !asset$all_converged ||
      !identical(asset$target_names, names(dat$beta0)) || any(!is.finite(asset$beta_star_mc))) stop("Invalid complete profile target.")
  list(beta = asset$beta_star_mc, mcse = asset$target_mc_se_by_coordinate,
       type = "full_unsmoothed_profile", path = path)
}

submission_direction <- function(config) {
  specification <- yaml::read_yaml("config/joint_profile_submission_20261027.yml")
  path <- file.path(specification$population_directions, paste0(config$experiment_id, ".rds"))
  if (!identical(submission_population_cache$direction_path, path)) {
    submission_population_cache$direction <- readRDS(path)
    submission_population_cache$direction_path <- path
  }
  asset <- submission_population_cache$direction
  if (asset$development || asset$n_population < 100000L || asset$repeats < 4L || !asset$all_converged ||
      abs(asset$h_analysis - config$h_multiplier * config$n_clusters^(-.3)) > 1e-12) stop("Invalid population direction.")
  asset
}

submission_core <- function(dat, config, coordinates) {
  fitted <- ds_fit_profile(dat$y, dat$X, dat$Z, dat$cluster_id, config$tau, coordinates,
    Lambda = diag(config$lambda_gamma, ncol(dat$Z)), h = config$h_multiplier * dat$n_clusters^(-.3),
    keep_prepared = TRUE, budget_multiplier = config$budget_multiplier)
  beta <- Reduce(`+`, lapply(fitted$training, function(fold) Reduce(`+`, fold$initial$coefficients) / 3)) / 2
  names(beta) <- colnames(dat$X)
  list(beta = beta, table = fitted$table, fit = fitted, converged = TRUE,
       global_estimator = "cross_fitted_outcome_learner", method_version = fitted$method_version)
}

submission_support_oracle <- function(dat, config, coordinates) {
  columns <- sort(unique(c(dat$active, match(coordinates, colnames(dat$X)))))
  fitted <- jp_fit_fixed(dat$y, dat$X[, columns, drop = FALSE], dat$Z, dat$cluster_id, config$tau,
    coordinates = coordinates, Lambda = diag(config$lambda_gamma, ncol(dat$Z)),
    h = config$h_multiplier * dat$n_clusters^(-.3))
  table <- data.frame(coordinate = fitted$table$coordinate, beta_tilde = fitted$table$beta_hat,
    se = fitted$table$se, lower = fitted$table$lower, upper = fitted$table$upper,
    max_kkt = pmax(fitted$diagnostics$full_kkt, fitted$diagnostics$half_kkt))
  beta <- map_beta_to_full_v2(fitted$beta, colnames(dat$X))
  list(beta = beta, table = table, fit = fitted, converged = TRUE,
       global_estimator = "oracle_support_paired_coefficients", method_version = "joint_profile_true_support_20261027")
}

submission_fold_oracle <- function(core, dat, config, coordinates, variant, population = NULL) {
  if (!variant %in% c("known_nuisance", "population_direction")) stop("Unknown oracle construction.")
  prepared <- core$fit$prepared
  folds_used <- 1:2
  sizes <- vapply(prepared$folds[folds_used], function(fold) length(fold$validation_subjects), integer(1))
  total <- sum(sizes)
  details <- vector("list", length(coordinates))
  for (j in seq_along(coordinates)) {
    estimates <- numeric(length(folds_used))
    influences <- vector("list", length(folds_used))
    diagnostics <- vector("list", length(folds_used))
    for (g in seq_along(folds_used)) {
      f <- folds_used[g]
      model <- prepared$folds[[f]]$models[[j]]
      fitted <- core$fit$details[[j]]$fits[[f]]
      if (variant == "known_nuisance") {
        small <- model$data
        effects <- dat$random_effects[as.character(small$cluster_id), , drop = FALSE]
        if (ncol(effects) != ncol(small$Z)) stop("The known nuisance dimensions differ.")
        small$y <- small$y - rowSums(small$Z * effects)
        small$Z <- matrix(0, length(small$y), 1L)
        fitted <- ds_pair_small(small, model$initial, config$tau, core$fit$h, diag(1), gram = model$gram)
        estimates[g] <- fitted$estimate
        influences[[g]] <- fitted$influence
        diagnostics[[g]] <- fitted
      } else if (variant == "population_direction") {
        rows <- unlist(cluster_index_v2(dat$cluster_id)$rows[prepared$folds[[f]]$validation_subjects], use.names = FALSE)
        full_data <- sb_data(dat$y[rows], dat$X[rows, , drop = FALSE], dat$Z[rows, , drop = FALSE], dat$cluster_id[rows])
        if (!identical(full_data$cluster_id, model$data$cluster_id) ||
            max(abs(full_data$X[, model$support, drop = FALSE] - model$data$X)) > 1e-12) stop("Oracle estimation rows differ.")
        omega <- population$directions[[coordinates[j]]]$omega[colnames(dat$X)]
        if (anyNA(omega)) stop("The population direction omits a fitted coordinate.")
        components <- lapply(c("full", "half"), function(band) {
          beta <- numeric(ncol(dat$X))
          beta[model$support] <- fitted[[band]]$beta
          bandwidth <- if (band == "full") core$fit$h else core$fit$h / 2
          profile <- ds_profile(full_data, beta, config$tau, bandwidth, diag(config$lambda_gamma, ncol(dat$Z)))
          contribution <- -drop(rowsum(drop(full_data$X %*% omega) * as.numeric(profile$psi) * full_data$weight,
                                       full_data$group, reorder = FALSE))
          augmented <- cbind(model$data$X, projected = drop(full_data$X %*% omega))
          curvature <- sb_geometry_cpp(augmented, full_data$Z, profile$phi, full_data$order,
            full_data$pointers, diag(config$lambda_gamma, ncol(dat$Z)), FALSE)
          dimension <- ncol(model$data$X)
          canonical <- c(1, rep(0, dimension - 1L))
          adjustment <- drop(solve(fitted[[band]]$hessian, canonical - curvature[seq_len(dimension), dimension + 1L]))
          adjustment_score <- -drop(rowsum(drop(model$data$X %*% adjustment) * as.numeric(profile$psi) * full_data$weight,
                                           full_data$group, reorder = FALSE))
          list(estimate = fitted[[band]]$beta[1L] - mean(contribution), contribution = contribution,
               linearised_score = contribution + adjustment_score, adjustment = adjustment)
        })
        estimates[g] <- 2 * components[[2L]]$estimate - components[[1L]]$estimate
        influences[[g]] <- -2 * components[[2L]]$linearised_score + components[[1L]]$linearised_score
        diagnostics[[g]] <- components
      } else {
        estimates[g] <- fitted$estimate
        influences[[g]] <- fitted$influence
        diagnostics[[g]] <- fitted
      }
    }
    estimate <- sum(sizes / total * estimates)
    se <- sqrt(sum(vapply(influences, function(value) sum((value - mean(value))^2), numeric(1)))) / total
    if (!is.finite(se) || se <= 0) stop("Invalid oracle cluster variance.")
    details[[j]] <- list(table = data.frame(coordinate = coordinates[j], beta_tilde = estimate,
      se = se, lower = estimate - qnorm(.975) * se, upper = estimate + qnorm(.975) * se),
      influences = influences, diagnostics = diagnostics)
  }
  list(beta = core$beta, table = do.call(rbind, lapply(details, `[[`, "table")),
       fit = list(details = details, variant = variant, folds_used = folds_used, n_clusters = total),
       converged = TRUE, global_estimator = "shared_outcome_learner", method_version = paste0("joint_profile_", variant, "_20261027"))
}

submission_benchmark <- function(method_id, train, tau, target_coords, tuning, seed, context, control) {
  started <- proc.time()[[3L]]
  set.seed(seed)
  config <- context$config
  if (method_id == "JOINT-PROFILE-TRUE-SUPPORT") {
    full <- subset_clusters_v2(context$raw_training, unique(train$cluster_id))
    full <- fit_working_nuisance_design_v2(full, config$fit_random_effects)
    result <- submission_support_oracle(full, config, target_coords)
  } else {
    cache <- context$shared_first_stage
    if (!is.environment(cache)) stop("A replication cache is required for the joint procedure.")
    if (!exists("joint_profile", envir = cache, inherits = FALSE)) {
      assign("joint_profile", submission_core(train, config, target_coords), envir = cache)
    }
    core <- get("joint_profile", envir = cache, inherits = FALSE)
    result <- switch(method_id,
      `JOINT-PROFILE` = core,
      `JOINT-PROFILE-SPLIT` = core,
      `JOINT-PROFILE-TRUE-NUISANCE` = submission_fold_oracle(core, train, config, target_coords, "known_nuisance"),
      `JOINT-PROFILE-POP-H` = submission_fold_oracle(core, train, config, target_coords, "population_direction", context$population_direction_asset),
      stop("Unknown joint procedure."))
  }
  table <- result$table
  table$index <- match(table$coordinate, colnames(train$X))
  table$feasible <- TRUE
  table$beta_hat <- result$beta[table$coordinate]
  table$ci_lower <- table$lower
  table$ci_upper <- table$upper
  for (name in c("mu", "dantzig_residual", "omega_l1", "omega_l2")) table[[name]] <- NA_real_
  kkt <- if ("max_kkt" %in% names(table)) max(table$max_kkt) else NA_real_
  retained <- result$fit
  retained$prepared <- NULL
  list(method_id = method_id, status = "ok", failure_stage = "", failure_message = "",
    beta_hat = result$beta, beta_tilde = setNames(table$beta_tilde, table$coordinate),
    se = setNames(table$se, table$coordinate), ci_lower = table$ci_lower, ci_upper = table$ci_upper,
    selected = if (method_id == "JOINT-PROFILE-TRUE-SUPPORT") train$active else which(abs(result$beta) > 1e-8),
    fit_object = list(beta = result$beta), inference_object = list(table = table, ci_level = .95),
    runtime_sec = proc.time()[[3L]] - started, converged = TRUE, kkt_residual = kkt,
    method_version = result$method_version, implementation_version = result$method_version,
    score_definition = "two_positive_profile_fits_paired_coefficients",
    reference_identifier = "docs/METHOD_SPECIFICATION.md", adapter_fidelity_status = "verified_joint_profile_definition",
    target_scope = if (method_id == "JOINT-PROFILE-TRUE-NUISANCE") "known_nuisance_structural_slope" else "complete_unsmoothed_profile",
    lambda_base = NA_real_, warning_messages = character(),
    submission_diagnostics = list(global_estimator = result$global_estimator, fit = retained))
}

submission_interval_target <- function(method_id, config, target) {
  if (method_id == "SQR-DEBIASED-IID" && config$experiment_id %in% c("B16", "B19")) {
    path <- file.path("results/intermediate/joint_profile_submission_20261027/marginal_targets",
                       paste0(config$experiment_id, ".rds"))
    asset <- readRDS(path)
    if (!asset$deterministic_target || asset$max_kkt_residual > 1e-8 ||
        !identical(names(asset$beta_marginal), names(target$beta))) stop("Invalid marginal quantile target.")
    target$beta <- asset$beta_marginal
    target$target_mc_se_by_coordinate <- asset$target_mc_se_by_coordinate
    target$target_mc_se <- 0
    target$target_type <- "own_marginal_quantile_projection"
    target$asset_path <- path
  }
  target
}

submission_pinball <- function(y, prediction, cluster_id, tau) {
  residual <- y - prediction
  mean(tapply(residual * (tau - (residual < 0)), cluster_id, mean))
}

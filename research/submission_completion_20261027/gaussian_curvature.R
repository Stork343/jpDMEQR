submission_gaussian_curvature <- function(tau, m_values, lambda, sigma_b) {
  if (tau <= 0 || tau >= 1 || lambda <= 0 || sigma_b <= 0 || any(m_values < 1))
    stop("Invalid Gaussian profile curvature parameters.")
  integrate_score <- function(order) {
    gaussian <- statmod::gauss.quad.prob(order, dist = "normal", sigma = sigma_b)
    legendre <- statmod::gauss.quad(order / 2L, kind = "legendre")
    vapply(m_values, function(m) {
      total <- numeric(length(gaussian$nodes))
      for (j in 0:m) {
        z <- (tau - j / m) / lambda - gaussian$nodes + qnorm(tau)
        probability <- pnorm(z)
        total <- total + dnorm(z) * ((1 - tau) * dbinom(j - 1L, m - 1L, probability) +
          tau * dbinom(j, m - 1L, probability))
      }
      for (j in seq_len(m)) {
        lower <- (tau - j / m) / lambda
        upper <- (tau - (j - 1L) / m) / lambda
        residual <- (lower + upper) / 2 + (upper - lower) / 2 * legendre$nodes
        z <- outer(gaussian$nodes, residual, function(effect, value) value - effect + qnorm(tau))
        probability <- pnorm(z)
        density <- dnorm(z)
        order_density <- m * dbinom(j - 1L, m - 1L, probability) * density
        dual <- lambda * residual - ((m - 1L) * tau - j + 1L) / m
        integrand <- sweep(order_density * z, 2L, dual, "*")
        if (m > 1L) integrand <- integrand + (m - 1L) * density^2 *
          ((1 - tau) * dbinom(j - 2L, m - 2L, probability) + tau * dbinom(j - 1L, m - 2L, probability))
        total <- total + drop(integrand %*% legendre$weights) * (upper - lower) / 2
      }
      sum(gaussian$weights * total)
    }, numeric(1))
  }
  integrated <- integrate_score(256L)
  checked <- integrate_score(512L)
  if (max(abs(integrated - checked)) > 1e-9 || any(integrated <= 0))
    stop("Independent Gaussian curvature integrals disagree.")
  list(curvature = mean(checked), curvature_by_cluster_size = checked,
    cluster_sizes = m_values, independent_difference = max(abs(integrated - checked)),
    tau = tau, lambda = lambda, sigma_b = sigma_b,
    definition = "Gaussian_Stein_expected_innovation_times_profile_dual")
}

submission_ar1_inverse_row <- function(p, rho, coordinate) {
  if (p < 1 || abs(rho) >= 1 || !coordinate %in% seq_len(p)) stop("Invalid AR1 inverse row.")
  result <- numeric(p)
  if (p == 1L) return(1)
  result[coordinate] <- if (coordinate %in% c(1L, p)) 1 else 1 + rho^2
  neighbours <- c(coordinate - 1L, coordinate + 1L)
  neighbours <- neighbours[neighbours >= 1L & neighbours <= p]
  result[neighbours] <- -rho
  result / (1 - rho^2)
}

submission_gaussian_diagnostics <- function(result, curvature) {
  config <- result$config
  if (config$module != "A_theory_scaling" || config$q != 1L || config$screen_method != "none" ||
      config$error_dist != "gaussian" || config$random_effect_dist != "normal" || config$heteroskedastic ||
      config$design_type != "ar1" || config$x_b_corr != 0 || config$informative_size)
    stop("The exact curvature calculation requires the declared Gaussian scaling design.")
  metric <- result$metrics[result$metrics$method_id == "JOINT-PROFILE", ]
  if (nrow(metric) != 1L) stop("The Gaussian diagnostic lacks its joint procedure.")
  if (!isTRUE(metric$converged)) return(list(available = FALSE, failure_stage = metric$failure_stage,
    failure_message = metric$failure_message))
  dat <- simulate_from_config_v2(config, result$seed)
  data <- sb_data(dat$y, dat$X, dat$Z, dat$cluster_id)
  fit <- result$method_diagnostics[["JOINT-PROFILE"]]$submission$fit
  if (abs(fit$h - config$h_multiplier * data$n^(-.3)) > 1e-12) stop("Saved diagnostic bandwidth differs.")
  profiles <- lapply(c(fit$h, fit$h / 2), function(h) ds_profile(data, dat$beta0, config$tau, h, diag(config$lambda_gamma)))
  paired_score <- (2 * profiles[[2L]]$psi - profiles[[1L]]$psi) * data$weight
  fold_records <- list()
  coordinate_records <- list()
  for (j in seq_len(nrow(fit$table))) {
    coordinate <- match(fit$table$coordinate[j], colnames(dat$X))
    omega <- submission_ar1_inverse_row(ncol(dat$X), config$rho_x, coordinate) / curvature$curvature
    projected <- drop(dat$X %*% omega)
    influence <- drop(rowsum(projected * paired_score, data$group, reorder = FALSE))
    estimate <- fit$table$beta_tilde[j]
    reconstruction <- 0
    variance <- 0
    shifts <- numeric(2L)
    for (f in 1:2) {
      training <- fit$training[[f]]
      graph <- training$graphs[[j]]
      support <- graph$support
      if (support[1L] != coordinate) stop("The selected target coordinate is not first.")
      rows <- which(as.character(dat$cluster_id) %in% training$validation_subjects)
      small <- sb_data(dat$y[rows], dat$X[rows, support, drop = FALSE], dat$Z[rows, , drop = FALSE], dat$cluster_id[rows])
      covariance <- config$rho_x^abs(outer(support, support, "-"))
      population <- curvature$curvature * covariance
      projected_target <- drop(solve(covariance,
        config$rho_x^abs(outer(support, dat$active, "-")) %*% dat$beta0[dat$active]))
      shifts[f] <- projected_target[1L] - dat$beta0[coordinate]
      saved <- fit$details[[j]]$fits[[f]]
      reconstruction <- reconstruction + small$n / data$n * (2 * saved$half$beta[1L] - saved$full$beta[1L])
      variance <- variance + sum((saved$influence - mean(saved$influence))^2)
      for (b in 1:2) {
        band <- c("full", "half")[b]
        fitted <- saved[[band]]
        H <- sb_geometry_cpp(small$X, small$Z, profiles[[b]]$phi[rows], small$order, small$pointers,
          diag(config$lambda_gamma), FALSE)
        inverse <- numeric(ncol(dat$X))
        inverse[support] <- fitted$inverse
        fold_records[[length(fold_records) + 1L]] <- data.frame(
          coordinate = fit$table$coordinate[j], fold = f, bandwidth = fit$h / (2^(b - 1L)),
          n_estimation = small$n, selected_dimension = length(support),
          curvature_max_error = max(abs(H - population)),
          curvature_operator_error = max(abs(eigen(H - population, symmetric = TRUE, only.values = TRUE)$values)),
          fitted_curvature_max_error = max(abs(fitted$hessian - population)),
          inverse_row_l2_error = sqrt(sum((inverse - omega)^2)),
          coefficient_l2_error = sqrt(sum((fitted$beta - projected_target)^2)),
          target_shift = shifts[f], kkt = fitted$kkt)
      }
    }
    if (abs(reconstruction - estimate) > 1e-12 || abs(sqrt(variance) / data$n - fit$table$se[j]) > 1e-12)
      stop("Saved paired coefficients or their covariance fail reconstruction.")
    coordinate_records[[j]] <- data.frame(coordinate = fit$table$coordinate[j],
      coordinate_type = if (coordinate %in% dat$active) "active" else "null",
      estimate = estimate, target = dat$beta0[coordinate], estimated_se = fit$table$se[j],
      leading_mean = mean(influence), scaled_remainder = sqrt(data$n) *
        (estimate - dat$beta0[coordinate] - mean(influence)),
      remainder_over_se = (estimate - dat$beta0[coordinate] - mean(influence)) / fit$table$se[j],
      maximum_target_shift = max(abs(shifts)))
  }
  list(available = TRUE, coordinates = do.call(rbind, coordinate_records),
    folds = do.call(rbind, fold_records), curvature = curvature,
    seed = result$seed, method_version = fit$method_version)
}

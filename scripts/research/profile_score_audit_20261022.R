args <- base::commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
out <- file.path(root, "results", "scratch", "profile_score_audit_20261022")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (identical(Sys.getenv("JPDMEQR_REUSE_RESEARCH_RESULTS"), "TRUE")) {
  required <- c("repeated_geometry.csv", "population_targets.csv",
                "population_sandwich.csv", "cluster_dgp.csv", "checks.csv",
                "COMPLETE.txt")
  stopifnot(all(file.exists(file.path(out, required))))
  checks <- utils::read.csv(file.path(out, "checks.csv"))
  stopifnot(all(checks$pass))
  source(file.path(root, "scripts", "research", "build_profile_score_summary_20261022.R"))
  quit(save = "no", status = 0L)
}
stage <- "initialisation"
options(error = function() {
  writeLines(c(paste("stage:", stage), geterrmessage()), file.path(out, "FAILED.txt"))
  traceback(12)
  quit(save = "no", status = 1L)
})

write_csv <- function(x, name) {
  utils::write.csv(x, file.path(out, name), row.names = FALSE, na = "")
}
packages <- c("MASS", "quantreg", "CVXR", "Matrix", "testthat", "jsonlite", "ps", "highs")
installed <- utils::installed.packages()
write_csv(data.frame(
  package = packages,
  installed = packages %in% rownames(installed),
  version = vapply(packages, function(x) {
    if (x %in% rownames(installed)) installed[x, "Version"] else NA_character_
  }, character(1))
), "dependencies.csv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2")
source_v2_module(root, "simulation_v2")

# 使用单调方程求解重复观测的平滑 profile score。
repeat_score <- function(w, tau, lambda, h) {
  if (h == 0) return(pmax(tau - 1, pmin(tau, lambda * w)))
  lower <- (tau - 1) / lambda - h
  upper <- tau / lambda + h
  vapply(w, function(x) {
    if (x <= lower) return(tau - 1)
    if (x >= upper) return(tau)
    t <- stats::uniroot(function(t) {
      (tau - 1 + kernel_epan_cdf_v2(t)) / lambda + h * t - x
    }, c(-1, 1), tol = 2e-13)$root
    tau - 1 + kernel_epan_cdf_v2(t)
  }, numeric(1))
}
repeat_hessian <- function(w, tau, lambda, h) {
  s <- repeat_score(w, tau, lambda, h)
  t <- (w - s / lambda) / h
  k <- kernel_epan_pdf_v2(t)
  lambda * k / (k + lambda * h)
}

# 变量变换使积分只涉及固定区间，并明确保留两端的概率。
repeat_moment <- function(beta, tau, lambda, h, power = 1L) {
  lo <- tau - 1
  hi <- tau
  if (h == 0) {
    middle <- stats::integrate(function(w) {
      (lambda * w)^power * stats::dnorm(beta + w)
    }, lo / lambda, hi / lambda, abs.tol = 2e-12, rel.tol = 2e-11)$value
  } else {
    middle <- stats::integrate(function(t) {
      s <- lo + kernel_epan_cdf_v2(t)
      s^power * stats::dnorm(beta + s / lambda + h * t) *
        (kernel_epan_pdf_v2(t) / lambda + h)
    }, -1, 1, abs.tol = 2e-12, rel.tol = 2e-11)$value
  }
  lo^power * stats::pnorm(beta + lo / lambda - h) + middle +
    hi^power * stats::pnorm(beta + hi / lambda + h, lower.tail = FALSE)
}
repeat_population_hessian <- function(beta, tau, lambda, h) {
  if (h == 0) {
    return(lambda * (stats::pnorm(beta + tau / lambda) -
                      stats::pnorm(beta + (tau - 1) / lambda)))
  }
  stats::integrate(function(t) {
    s <- tau - 1 + kernel_epan_cdf_v2(t)
    kernel_epan_pdf_v2(t) * stats::dnorm(beta + s / lambda + h * t)
  }, -1, 1, abs.tol = 2e-12, rel.tol = 2e-11)$value
}
population_root <- function(tau, lambda, h, extrapolated = FALSE) {
  equation <- function(beta) {
    if (extrapolated) {
      2 * repeat_moment(beta, tau, lambda, h / 2) -
        repeat_moment(beta, tau, lambda, h)
    } else repeat_moment(beta, tau, lambda, h)
  }
  stats::uniroot(equation, c(-5, 5), tol = 2e-12)$root
}

stage <- "repeated_geometry"
geometry <- list()
index <- 0L
for (tau in c(0.25, 0.5, 0.75)) {
  for (lambda in c(0.5, 1, 2)) {
    for (h in c(0.2, 0.1, 0.05, 0.025)) {
      w_values <- c((tau - 1) / lambda - 0.3,
                    (tau - 0.85) / lambda,
                    (tau - 0.5) / lambda - 0.03,
                    (tau - 0.5) / lambda + 0.03,
                    (tau - 0.15) / lambda,
                    tau / lambda + 0.3)
      for (w in w_values) {
        dat <- list(y = rep(w, 2), X = matrix(1, 2, 1),
                    Z = matrix(1, 2, 1), cluster_id = rep(1, 2))
        component <- function(beta) profile_components_v2(
          dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, h,
          lambda_gamma = lambda, identity_check = TRUE,
          nuisance_control = list(grad_tol = 1e-11, reltol = 1e-13)
        )
        comp <- component(0)
        eps <- 2e-6
        numeric_hessian <- (component(eps)$score - component(-eps)$score) / (2 * eps)
        exact_fd <- (-repeat_score(w - eps, tau, lambda, h) +
                       repeat_score(w + eps, tau, lambda, h)) / (2 * eps)
        s <- repeat_score(w, tau, lambda, h)
        H <- repeat_hessian(w, tau, lambda, h)
        index <- index + 1L
        geometry[[index]] <- data.frame(
          tau = tau, lambda = lambda, h = h, w = w,
          score = comp$score, analytic_score = -s,
          score_error = abs(comp$score + s),
          hessian = comp$hessian[1, 1], analytic_hessian = H,
          hessian_error = abs(comp$hessian[1, 1] - H),
          hessian_fd_error = abs(numeric_hessian - H),
          hessian_exact_fd_error = abs(exact_fd - H),
          identity_error = comp$identity_error %||% comp$profile_identity_error,
          ordinary_score = comp$unsmoothed_cluster_scores[1, 1]
        )
      }
    }
  }
}
geometry <- do.call(rbind, geometry)
write_csv(geometry, "repeated_geometry.csv")
stopifnot(max(geometry$score_error) < 2e-8,
          max(geometry$hessian_error) < 2e-7,
          max(geometry$hessian_fd_error) < 3e-3,
          max(geometry$hessian_exact_fd_error) < 1e-6)

stage <- "population_integration"
targets <- list()
sandwich <- list()
index <- 0L
h_values <- 0.4 / 2^(0:6)
for (tau in c(0.25, 0.5, 0.75)) {
  for (lambda in c(0.5, 1, 2)) {
    beta0 <- population_root(tau, lambda, 0)
    H0 <- repeat_population_hessian(beta0, tau, lambda, 0)
    V0 <- repeat_moment(beta0, tau, lambda, 0, 2)
    score_linear <- -stats::integrate(function(t) {
      s <- tau - 1 + kernel_epan_cdf_v2(t)
      t * kernel_epan_pdf_v2(t) * stats::dnorm(beta0 + s / lambda)
    }, -1, 1, abs.tol = 2e-12, rel.tol = 2e-11)$value
    probability <- stats::pnorm(beta0 + (tau - 0.5) / lambda)
    Vraw0 <- probability * (1 - probability)
    for (h in h_values) {
      beta_h <- population_root(tau, lambda, h)
      beta_r <- population_root(tau, lambda, h, TRUE)
      Hh <- repeat_population_hessian(beta_h, tau, lambda, h)
      mean_h <- repeat_moment(beta_h, tau, lambda, h)
      Vh <- repeat_moment(beta_h, tau, lambda, h, 2) - mean_h^2
      probability_h <- stats::pnorm(beta_h + (tau - 0.5) / lambda)
      Vraw <- probability_h * (1 - probability_h)
      index <- index + 1L
      targets[[index]] <- data.frame(
        tau = tau, lambda = lambda, h = h, beta_star = beta0,
        beta_h = beta_h, beta_extrapolated = beta_r,
        error_h = beta_h - beta0, error_extrapolated = beta_r - beta0,
        target_linear_coefficient = score_linear / H0,
        target_error_over_h = (beta_h - beta0) / h,
        score_at_beta_star = repeat_moment(beta0, tau, lambda, h),
        extrapolated_score_at_beta_star =
          2 * repeat_moment(beta0, tau, lambda, h / 2) -
          repeat_moment(beta0, tau, lambda, h)
      )
      sandwich[[index]] <- data.frame(
        tau = tau, lambda = lambda, h = h, beta_h = beta_h,
        atom_probability_limit = stats::pnorm(beta0 + tau / lambda) -
          stats::pnorm(beta0 + (tau - 1) / lambda),
        hessian = Hh, hessian_limit = H0,
        meat_profile = Vh, meat_profile_limit = V0,
        meat_ordinary = Vraw, meat_ordinary_limit = Vraw0,
        root_n_variance_profile = Vh / Hh^2,
        root_n_variance_ordinary = Vraw / Hh^2,
        root_n_variance_profile_limit = V0 / H0^2,
        root_n_variance_ordinary_limit = Vraw0 / H0^2
      )
    }
  }
}
targets <- do.call(rbind, targets)
sandwich <- do.call(rbind, sandwich)
targets$order_h <- NA_real_
targets$order_extrapolated <- NA_real_
for (ix in split(seq_len(nrow(targets)), interaction(targets$tau, targets$lambda))) {
  for (j in seq_len(length(ix) - 1L)) {
    a <- ix[j]
    b <- ix[j + 1L]
    targets$order_h[b] <- log(abs(targets$error_h[a] / targets$error_h[b])) / log(2)
    targets$order_extrapolated[b] <-
      log(abs(targets$error_extrapolated[a] / targets$error_extrapolated[b])) / log(2)
  }
}
write_csv(targets, "population_targets.csv")
write_csv(sandwich, "population_sandwich.csv")
stopifnot(all(is.finite(targets$beta_h)), all(is.finite(targets$beta_extrapolated)))

# q=1 截距的非平滑优化同时比较数值最优点与观测折点。
exact_cluster <- function(y, X, beta, tau, lambda) {
  u <- as.numeric(y - X %*% beta)
  objective <- function(gamma) mean(check_loss_v2(u - gamma, tau)) +
    lambda * gamma^2 / 2
  opt <- stats::optimize(objective, c((tau - 1) / lambda - 1, tau / lambda + 1),
                         tol = 1e-12)
  candidates <- c(unique(u), opt$minimum)
  gamma <- candidates[which.min(vapply(candidates, objective, numeric(1)))]
  r <- u - gamma
  zero <- abs(r) < 1e-7
  s <- psi_quantile_v2(r, tau)
  if (any(zero)) {
    s[zero] <- (length(u) * lambda * gamma - sum(s[!zero])) / sum(zero)
  }
  stopifnot(all(s >= tau - 1 - 1e-6), all(s <= tau + 1e-6),
            abs(mean(s) - lambda * gamma) < 2e-6)
  list(score = -as.numeric(crossprod(X, s)) / length(u),
       residual = r, subgradient = s,
       stationarity = abs(mean(s) - lambda * gamma))
}

stage <- "cluster_dgp"
cluster_results <- list()
index <- 0L
for (tau in c(0.25, 0.5, 0.75)) {
  dat <- generate_profile_qr_data_v2(
    n_clusters = 80, p = 2, s = 1, tau = tau, q = 1,
    m_values = 2:4, rho_x = 0.3, sigma_b0 = 0.8, seed = 20261022
  )
  beta <- as.numeric(dat$beta0)
  rows <- cluster_index_v2(dat$cluster_id)$rows
  reference <- lapply(rows, function(ix) {
    exact_cluster(dat$y[ix], dat$X[ix, , drop = FALSE], beta, tau, 1)
  })
  Gref <- do.call(rbind, lapply(reference, `[[`, "score"))
  atom <- mean(unlist(lapply(reference, `[[`, "residual")) == 0)
  atom_cluster <- mean(vapply(reference, function(x) any(abs(x$residual) < 1e-7), logical(1)))
  saveRDS(list(data = dat, reference = reference),
          file.path(out, paste0("cluster_reference_tau", tau, ".rds")))
  for (h in c(0.2, 0.1, 0.05, 0.025)) {
    comp <- profile_components_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, h, lambda_gamma = 1,
      identity_check = TRUE,
      nuisance_control = list(grad_tol = 1e-10, reltol = 1e-12)
    )
    fit <- fit_profile_lasso_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, tau, h,
      lambda_beta = 0, lambda_gamma = 1, penalty_factor = c(0, 0),
      control = list(max_iter = 500, beta_tol = 1e-8,
                     kkt_normalized_tol = 1e-7)
    )
    saveRDS(list(component = comp, fit = fit),
            file.path(out, paste0("cluster_tau", tau, "_h", h, ".rds")))
    for (coordinate in seq_len(2)) {
      smooth <- comp$cluster_scores[, coordinate]
      ordinary <- comp$unsmoothed_cluster_scores[, coordinate]
      actual <- Gref[, coordinate]
      index <- index + 1L
      cluster_results[[index]] <- data.frame(
        tau = tau, h = h, coordinate = coordinate,
        residual_atom_fraction = atom, clusters_with_atom_fraction = atom_cluster,
        nuisance_kkt_max = max(vapply(reference, `[[`, numeric(1), "stationarity")),
        smooth_to_profile_rmse = sqrt(mean((smooth - actual)^2)),
        ordinary_to_profile_rmse = sqrt(mean((ordinary - actual)^2)),
        profile_meat = mean((actual - mean(actual))^2),
        smooth_meat = mean((smooth - mean(smooth))^2),
        ordinary_meat = mean((ordinary - mean(ordinary))^2),
        fit_converged = fit$converged, fit_kkt = fit$kkt_residual,
        fit_beta = fit$beta[coordinate],
        fit_smooth_meat = mean((fit$components$cluster_scores[, coordinate] -
                                 mean(fit$components$cluster_scores[, coordinate]))^2),
        fit_ordinary_meat = mean((fit$components$unsmoothed_cluster_scores[, coordinate] -
                                   mean(fit$components$unsmoothed_cluster_scores[, coordinate]))^2)
      )
    }
    write_csv(do.call(rbind, cluster_results), "cluster_dgp.csv")
  }
}
cluster_results <- do.call(rbind, cluster_results)
write_csv(cluster_results, "cluster_dgp.csv")
write_csv(data.frame(
  check = c("geometry_score", "geometry_hessian", "geometry_hessian_fd",
            "population_integral", "small_dgp_optimisation"),
  pass = c(max(geometry$score_error) < 2e-8,
           max(geometry$hessian_error) < 2e-7,
           max(geometry$hessian_fd_error) < 3e-3,
           all(is.finite(targets$beta_h)), all(cluster_results$fit_converged))
), "checks.csv")
writeLines("COMPLETE", file.path(out, "COMPLETE.txt"))
print(utils::tail(subset(targets, tau == 0.25 & lambda == 1), 3))
print(subset(sandwich, tau == 0.25 & lambda == 1 & h == min(h_values)))
print(utils::tail(cluster_results, 4))
source(file.path(root, "scripts", "research", "build_profile_score_summary_20261022.R"))

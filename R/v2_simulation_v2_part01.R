# Data-generating mechanisms and low-level simulation utilities for the SJS
# reconstruction. Source the profile-v2 module before fitting any estimator.

parse_semicolon_numeric_v2 <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(as.character(x))) return(numeric())
  as.numeric(strsplit(as.character(x), ";", fixed = TRUE)[[1]])
}

parse_pipe_character_v2 <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(as.character(x))) return(character())
  strsplit(as.character(x), "|", fixed = TRUE)[[1]]
}

seed_from_id_v2 <- function(base_seed, experiment_id, replicate, role = "replicate") {
  key <- paste(experiment_id, replicate, role, sep = "::")
  ints <- utf8ToInt(key)
  hash <- 0
  for (v in ints) hash <- (hash * 131 + v) %% 2147483000
  seed <- (as.numeric(base_seed) + hash + 104729 * as.numeric(replicate)) %% 2147483000
  as.integer(max(seed, 1))
}

ar1_covariance_v2 <- function(p, rho) {
  outer(seq_len(p), seq_len(p), function(i, j) rho^abs(i - j))
}

block_covariance_v2 <- function(p, block_size = 25L, rho = 0.7) {
  Sigma <- diag(p)
  blocks <- split(seq_len(p), ceiling(seq_len(p) / block_size))
  for (idx in blocks) Sigma[idx, idx] <- rho
  diag(Sigma) <- 1
  Sigma
}

sparse_precision_covariance_v2 <- function(p, offdiag = -0.25) {
  Omega <- diag(1, p)
  if (p > 1L) {
    Omega[cbind(seq_len(p - 1L), 2:p)] <- offdiag
    Omega[cbind(2:p, seq_len(p - 1L))] <- offdiag
  }
  Sigma <- solve(Omega)
  d <- sqrt(diag(Sigma))
  Sigma / tcrossprod(d)
}

dense_precision_covariance_v2 <- function(p, rho = 0.35) {
  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  Sigma
}

factor_covariance_v2 <- function(p, n_factors = 3L, loading_scale = 0.35) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(918273 + p + n_factors)
  L <- matrix(stats::rnorm(p * n_factors, sd = loading_scale), p, n_factors)
  Sigma <- tcrossprod(L) + diag(0.75, p)
  d <- sqrt(diag(Sigma))
  Sigma / tcrossprod(d)
}

make_design_covariance_v2 <- function(p, design_type = "ar1", rho_x = 0.5) {
  switch(
    as.character(design_type),
    ar1 = ar1_covariance_v2(p, rho_x),
    independent = diag(p),
    block = block_covariance_v2(p, block_size = 25L, rho = 0.7),
    sparse_precision = sparse_precision_covariance_v2(p),
    factor = factor_covariance_v2(p),
    dense_precision = dense_precision_covariance_v2(p, rho = 0.35),
    stop("Unknown design_type: ", design_type)
  )
}


tridiagonal_precision_sampler_v2 <- function(p, offdiag = -0.25) {
  if (p < 1L) stop("p must be positive.")
  if (!is.finite(offdiag) || abs(offdiag) >= 0.5) {
    stop("For unit diagonal, abs(offdiag) must be below 0.5.")
  }
  diag_l <- numeric(p)
  sub_l <- numeric(p)
  diag_l[1] <- 1
  if (p > 1L) {
    for (k in 2:p) {
      sub_l[k] <- offdiag / diag_l[k - 1L]
      value <- 1 - sub_l[k]^2
      if (value <= 0) stop("Precision Cholesky recurrence lost positive definiteness.")
      diag_l[k] <- sqrt(value)
    }
  }

  # Exact diagonal of the inverse by the Usmani tridiagonal recurrence.
  theta <- numeric(p + 1L)
  theta[1] <- 1
  theta[2] <- 1
  if (p > 1L) {
    for (k in 2:p) theta[k + 1L] <- theta[k] - offdiag^2 * theta[k - 1L]
  }
  phi <- numeric(p + 2L)
  phi[p + 2L] <- 1
  phi[p + 1L] <- 1
  if (p > 1L) {
    for (k in seq.int(p - 1L, 1L)) {
      phi[k + 1L] <- phi[k + 2L] - offdiag^2 * phi[k + 3L]
    }
  }
  inverse_diag <- vapply(seq_len(p), function(k) {
    theta[k] * phi[k + 2L] / theta[p + 1L]
  }, numeric(1))
  if (any(!is.finite(inverse_diag)) || any(inverse_diag <= 0)) {
    stop("Failed to obtain a positive inverse-precision diagonal.")
  }

  draw <- function(n) {
    n <- as.integer(n)
    z <- matrix(stats::rnorm(n * p), n, p)
    x <- matrix(0, n, p)
    x[, p] <- z[, p] / diag_l[p]
    if (p > 1L) {
      for (k in seq.int(p - 1L, 1L)) {
        x[, k] <- (z[, k] - sub_l[k + 1L] * x[, k + 1L]) / diag_l[k]
      }
    }
    sweep(x, 2, sqrt(inverse_diag), "/")
  }
  list(
    draw = draw,
    covariance = NULL,
    metadata = list(type = "tridiagonal_precision", offdiag = offdiag,
                    inverse_diagonal = inverse_diag)
  )
}

factor_loadings_v2 <- function(p, n_factors = 3L, loading_scale = 0.35) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(918273 + p + n_factors)
  matrix(stats::rnorm(p * n_factors, sd = loading_scale), p, n_factors)
}

make_design_sampler_v2 <- function(p, design_type = "ar1", rho_x = 0.5) {
  p <- as.integer(p)
  design_type <- as.character(design_type)
  if (p < 1L) stop("p must be positive.")

  if (design_type %in% c("ar1", "independent")) {
    rho <- if (design_type == "independent") 0 else rho_x
    if (!is.finite(rho) || abs(rho) >= 1) stop("AR(1) rho must lie in (-1,1).")
    draw <- function(n) {
      n <- as.integer(n)
      z <- matrix(stats::rnorm(n * p), n, p)
      if (p == 1L || rho == 0) return(z)
      x <- z
      innovation_sd <- sqrt(1 - rho^2)
      for (k in 2:p) x[, k] <- rho * x[, k - 1L] + innovation_sd * z[, k]
      x
    }
    covariance <- if (p <= 2000L) ar1_covariance_v2(p, rho) else NULL
    return(list(draw = draw, covariance = covariance,
                metadata = list(type = design_type, rho = rho)))
  }

  if (design_type == "block") {
    block_size <- 25L
    rho <- 0.7
    blocks <- split(seq_len(p), ceiling(seq_len(p) / block_size))
    draw <- function(n) {
      n <- as.integer(n)
      x <- matrix(0, n, p)
      for (idx in blocks) {
        common <- stats::rnorm(n)
        noise <- matrix(stats::rnorm(n * length(idx)), n, length(idx))
        x[, idx] <- sqrt(rho) * common + sqrt(1 - rho) * noise
      }
      x
    }
    covariance <- if (p <= 2000L) block_covariance_v2(p, block_size, rho) else NULL
    return(list(draw = draw, covariance = covariance,
                metadata = list(type = "block", block_size = block_size, rho = rho)))
  }

  if (design_type == "sparse_precision") {
    return(tridiagonal_precision_sampler_v2(p, offdiag = -0.25))
  }

  if (design_type == "factor") {
    n_factors <- 3L
    idiosyncratic_variance <- 0.75
    loadings <- factor_loadings_v2(p, n_factors = n_factors, loading_scale = 0.35)
    scale <- sqrt(rowSums(loadings^2) + idiosyncratic_variance)
    draw <- function(n) {
      n <- as.integer(n)
      factors <- matrix(stats::rnorm(n * n_factors), n, n_factors)
      noise <- matrix(stats::rnorm(n * p, sd = sqrt(idiosyncratic_variance)), n, p)
      sweep(factors %*% t(loadings) + noise, 2, scale, "/")
    }
    covariance <- if (p <= 1500L) {
      raw <- tcrossprod(loadings) + diag(idiosyncratic_variance, p)
      raw / tcrossprod(scale)
    } else NULL
    return(list(draw = draw, covariance = covariance,
                metadata = list(type = "factor", n_factors = n_factors,
                                loadings = loadings, scale = scale)))
  }

  if (design_type == "dense_precision") {
    # Equicorrelation covariance has a dense inverse and admits an O(np) sampler.
    rho <- 0.35
    draw <- function(n) {
      n <- as.integer(n)
      common <- stats::rnorm(n)
      noise <- matrix(stats::rnorm(n * p), n, p)
      sqrt(rho) * common + sqrt(1 - rho) * noise
    }
    covariance <- if (p <= 2000L) dense_precision_covariance_v2(p, rho) else NULL
    return(list(draw = draw, covariance = covariance,
                metadata = list(type = "dense_precision", rho = rho)))
  }

  stop("Unknown design_type: ", design_type)
}

rmvnorm_chol_v2 <- function(n, Sigma, mean = NULL) {
  Sigma <- as.matrix(Sigma)
  p <- ncol(Sigma)
  mean <- mean %||% rep(0, p)
  if (length(mean) != p) stop("mean dimension mismatch.")
  R <- chol((Sigma + t(Sigma)) / 2)
  Z <- matrix(stats::rnorm(n * p), n, p)
  sweep(Z %*% R, 2, mean, "+")
}

laplace_quantile_v2 <- function(prob, scale = 1 / sqrt(2)) {
  ifelse(
    prob < 0.5,
    scale * log(2 * prob),
    -scale * log(2 * (1 - prob))
  )
}

r_laplace_v2 <- function(n, scale = 1 / sqrt(2)) {
  laplace_quantile_v2(stats::runif(n), scale = scale)
}

contaminated_normal_cdf_v2 <- function(x, contamination = 0.1,
                                       sd_clean = 1, sd_outlier = 5) {
  (1 - contamination) * stats::pnorm(x, sd = sd_clean) +
    contamination * stats::pnorm(x, sd = sd_outlier)
}

contaminated_normal_quantile_v2 <- function(prob, contamination = 0.1,
                                             sd_clean = 1, sd_outlier = 5) {
  solve_one <- function(pp) {
    if (pp <= 0) return(-Inf)
    if (pp >= 1) return(Inf)
    stats::uniroot(
      function(x) contaminated_normal_cdf_v2(
        x, contamination = contamination,
        sd_clean = sd_clean, sd_outlier = sd_outlier
      ) - pp,
      interval = c(-50, 50), tol = 1e-11
    )$root
  }
  vapply(prob, solve_one, numeric(1))
}

asymmetric_laplace_quantile_zero_v2 <- function(prob, tau) {
  raw <- ifelse(
    prob < tau,
    log(prob / tau) / (1 - tau),
    -log((1 - prob) / (1 - tau)) / tau
  )
  sd_raw <- sqrt((1 - 2 * tau + 2 * tau^2) /
                   (tau^2 * (1 - tau)^2))
  raw / sd_raw
}

r_asymmetric_laplace_zero_quantile_v2 <- function(n, tau) {
  asymmetric_laplace_quantile_zero_v2(stats::runif(n), tau)
}

error_quantile_function_v2 <- function(prob, error_dist, tau) {
  prob <- pmin(pmax(as.numeric(prob), 1e-12), 1 - 1e-12)
  switch(
    as.character(error_dist),
    gaussian = stats::qnorm(prob),
    t3 = stats::qt(prob, df = 3) / sqrt(3),
    laplace = laplace_quantile_v2(prob),
    skew_chisq3 = (stats::qchisq(prob, df = 3) - 3) / sqrt(6),
    contaminated_normal = contaminated_normal_quantile_v2(prob),
    asymmetric_laplace = asymmetric_laplace_quantile_zero_v2(prob, tau),
    random_scale = stats::qnorm(prob),
    stop("Unknown error_dist: ", error_dist)
  )
}

r_error_marginal_v2 <- function(n, error_dist, tau) {
  if (error_dist == "contaminated_normal") {
    is_outlier <- stats::runif(n) < 0.1
    return(stats::rnorm(n, sd = ifelse(is_outlier, 5, 1)))
  }
  error_quantile_function_v2(stats::runif(n), error_dist, tau)
}

r_quantile_centered_error_v2 <- function(n,
                                         tau,
                                         error_dist = "gaussian",
                                         x1 = NULL,
                                         random_intercept = NULL,
                                         heteroskedastic = FALSE) {
  if (tau <= 0 || tau >= 1) stop("tau must be in (0,1).")
  if (is.null(x1)) x1 <- rep(0, n)
  if (is.null(random_intercept)) random_intercept <- rep(0, n)
  raw <- r_error_marginal_v2(n, error_dist, tau)
  qtau <- if (error_dist == "asymmetric_laplace") {
    0
  } else {
    error_quantile_function_v2(tau, error_dist, tau)
  }
  centered <- raw - qtau
  scale <- rep(1, n)
  if (isTRUE(heteroskedastic)) {
    scale <- scale * exp(0.35 * pmax(pmin(x1, 3), -3))
  }
  if (error_dist == "random_scale") {
    scale <- scale * exp(0.35 * pmax(pmin(random_intercept, 3), -3))
  }
  scale * centered
}

gaussian_copula_uniform_v2 <- function(m, rho = 0.4) {
  if (m < 1L) return(numeric())
  if (!is.finite(rho) || abs(rho) >= 1) stop("copula rho must lie in (-1,1).")
  R <- ar1_covariance_v2(m, rho)
  latent <- as.numeric(rmvnorm_chol_v2(1L, R))
  pmin(pmax(stats::pnorm(latent), 1e-12), 1 - 1e-12)
}

r_cluster_quantile_error_v2 <- function(m,
                                        tau,
                                        error_dist = "gaussian",
                                        x1 = NULL,
                                        random_intercept = NULL,
                                        heteroskedastic = FALSE,
                                        dependence = "independent",
                                        copula_rho = 0.4) {
  if (is.null(x1)) x1 <- rep(0, m)
  if (is.null(random_intercept)) random_intercept <- rep(0, m)
  if (is.na(dependence) || dependence %in% c("", "independent")) {
    return(r_quantile_centered_error_v2(
      m, tau, error_dist = error_dist, x1 = x1,
      random_intercept = random_intercept,
      heteroskedastic = heteroskedastic
    ))
  }
  if (!identical(dependence, "gaussian_copula_ar1")) {
    stop("Unknown error dependence: ", dependence)
  }
  u <- gaussian_copula_uniform_v2(m, rho = copula_rho)
  raw <- error_quantile_function_v2(u, error_dist, tau)
  qtau <- if (error_dist == "asymmetric_laplace") 0 else
    error_quantile_function_v2(tau, error_dist, tau)
  centered <- raw - qtau
  scale <- rep(1, m)
  if (isTRUE(heteroskedastic)) {
    scale <- scale * exp(0.35 * pmax(pmin(x1, 3), -3))
  }
  if (error_dist == "random_scale") {
    scale <- scale * exp(0.35 * pmax(pmin(random_intercept, 3), -3))
  }
  scale * centered
}

make_random_effect_covariance_v2 <- function(q, sigma_b0 = 1,
                                             sigma_b1 = 0.5, rho_b = 0.4) {
  if (q == 1L) return(matrix(sigma_b0^2, 1, 1))
  if (q != 2L) stop("Reference DGP currently supports q=1 or q=2.")
  matrix(c(
    sigma_b0^2, rho_b * sigma_b0 * sigma_b1,
    rho_b * sigma_b0 * sigma_b1, sigma_b1^2
  ), 2, 2)
}

r_random_effects_v2 <- function(n, Psi, distribution = "normal") {
  q <- ncol(Psi)
  if (distribution == "normal") {
    return(rmvnorm_chol_v2(n, Psi))
  }
  if (distribution == "t5") {
    Z <- rmvnorm_chol_v2(n, Psi)
    scale <- sqrt(stats::rchisq(n, df = 5) / 5)
    return(Z / scale * sqrt(3 / 5))
  }
  if (distribution == "mixture") {
    # Symmetric two-component mixture with exactly covariance Psi and
    # independent clusters. The component shift is along the first axis.
    delta <- numeric(q)
    delta[1] <- 0.5 * sqrt(Psi[1, 1])
    base_cov <- Psi - tcrossprod(delta)
    min_eig <- min(eigen((base_cov + t(base_cov)) / 2,
                         symmetric = TRUE, only.values = TRUE)$values)
    if (min_eig <= 1e-10) stop("Mixture base covariance is not positive definite.")
    sign_group <- ifelse(stats::rbinom(n, 1, 0.5) == 1L, 1, -1)
    return(rmvnorm_chol_v2(n, base_cov) + sign_group *
             matrix(delta, nrow = n, ncol = q, byrow = TRUE))
  }
  if (distribution == "skew_lognormal") {
    if (q != 1L) stop("skew_lognormal is defined only for q=1.")
    log_sd <- 0.65
    raw <- exp(stats::rnorm(n, sd = log_sd))
    mean_raw <- exp(log_sd^2 / 2)
    sd_raw <- sqrt((exp(log_sd^2) - 1) * exp(log_sd^2))
    x <- (raw - mean_raw) / sd_raw * sqrt(Psi[1, 1])
    return(matrix(x, n, 1))
  }
  stop("Unknown random_effect_dist: ", distribution)
}

make_beta_v2 <- function(p, s, signal = 0.75) {
  if (s > p) stop("s cannot exceed p.")
  beta <- numeric(p)
  beta[seq_len(s)] <- signal * (-1)^(seq_len(s) + 1)
  names(beta) <- sprintf("x%05d", seq_len(p))
  beta
}

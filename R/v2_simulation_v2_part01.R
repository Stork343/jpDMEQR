# Data-generating mechanisms and simulation utilities for the SJS reconstruction.
# Source R/profile_v2.R before using functions that fit the proposed estimator.

parse_semicolon_numeric_v2 <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) return(numeric())
  as.numeric(strsplit(as.character(x), ";", fixed = TRUE)[[1]])
}

parse_pipe_character_v2 <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) return(character())
  strsplit(as.character(x), "|", fixed = TRUE)[[1]]
}

seed_from_id_v2 <- function(base_seed, experiment_id, replicate) {
  ints <- utf8ToInt(paste0(experiment_id, "::", replicate))
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
  if (p > 1) {
    Omega[cbind(seq_len(p - 1), 2:p)] <- offdiag
    Omega[cbind(2:p, seq_len(p - 1))] <- offdiag
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
  set.seed(918273 + p + n_factors)
  L <- matrix(stats::rnorm(p * n_factors, sd = loading_scale), p, n_factors)
  Sigma <- tcrossprod(L) + diag(0.75, p)
  d <- sqrt(diag(Sigma))
  Sigma / tcrossprod(d)
}

make_design_covariance_v2 <- function(p, design_type = "ar1", rho_x = 0.5) {
  switch(
    design_type,
    ar1 = ar1_covariance_v2(p, rho_x),
    block = block_covariance_v2(p, block_size = 25L, rho = 0.7),
    sparse_precision = sparse_precision_covariance_v2(p),
    factor = factor_covariance_v2(p),
    dense_precision = dense_precision_covariance_v2(p, rho = 0.35),
    stop("Unknown design_type: ", design_type)
  )
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
  ifelse(prob < 0.5,
         scale * log(2 * prob),
         -scale * log(2 * (1 - prob)))
}

r_laplace_v2 <- function(n, scale = 1 / sqrt(2)) {
  laplace_quantile_v2(stats::runif(n), scale = scale)
}

contaminated_normal_quantile_v2 <- function(prob, contamination = 0.1,
                                             sd_clean = 1, sd_outlier = 5) {
  cdf <- function(x) {
    (1 - contamination) * stats::pnorm(x, sd = sd_clean) +
      contamination * stats::pnorm(x, sd = sd_outlier)
  }
  stats::uniroot(function(x) cdf(x) - prob, interval = c(-50, 50),
                 tol = 1e-12)$root
}

r_asymmetric_laplace_zero_quantile_v2 <- function(n, tau) {
  u <- stats::runif(n)
  raw <- ifelse(
    u < tau,
    log(u / tau) / (1 - tau),
    -log((1 - u) / (1 - tau)) / tau
  )
  # The variance for scale one is (1 - 2*tau + 2*tau^2)/(tau^2(1-tau)^2).
  sd_raw <- sqrt((1 - 2 * tau + 2 * tau^2) / (tau^2 * (1 - tau)^2))
  raw / sd_raw
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

  if (error_dist == "gaussian") {
    raw <- stats::rnorm(n)
    qtau <- stats::qnorm(tau)
  } else if (error_dist == "t3") {
    raw <- stats::rt(n, df = 3) / sqrt(3)
    qtau <- stats::qt(tau, df = 3) / sqrt(3)
  } else if (error_dist == "laplace") {
    raw <- r_laplace_v2(n)
    qtau <- laplace_quantile_v2(tau)
  } else if (error_dist == "skew_chisq3") {
    raw <- (stats::rchisq(n, df = 3) - 3) / sqrt(6)
    qtau <- (stats::qchisq(tau, df = 3) - 3) / sqrt(6)
  } else if (error_dist == "contaminated_normal") {
    is_outlier <- stats::runif(n) < 0.1
    raw <- stats::rnorm(n, sd = ifelse(is_outlier, 5, 1))
    qtau <- contaminated_normal_quantile_v2(tau)
  } else if (error_dist == "asymmetric_laplace") {
    raw <- r_asymmetric_laplace_zero_quantile_v2(n, tau)
    qtau <- 0
  } else if (error_dist == "random_scale") {
    raw <- stats::rnorm(n)
    qtau <- stats::qnorm(tau)
  } else {
    stop("Unknown error_dist: ", error_dist)
  }

  centered <- raw - qtau
  scale <- rep(1, n)
  if (isTRUE(heteroskedastic)) scale <- scale * exp(0.35 * pmax(pmin(x1, 3), -3))
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
    # t5 covariance inflation is 5/(5-2); correct it.
    return(Z / scale * sqrt(3 / 5))
  }
  if (distribution == "mixture") {
    group <- stats::rbinom(n, 1, 0.5)
    means <- matrix(0, n, q)
    means[, 1] <- ifelse(group == 1, 0.65, -0.65)
    base <- rmvnorm_chol_v2(n, 0.6 * Psi)
    out <- base + means
    out <- sweep(out, 2, colMeans(out), "-")
    # Rescale approximately to target marginal SDs.
    current_sd <- apply(out, 2, stats::sd)
    target_sd <- sqrt(diag(Psi))
    sweep(out, 2, current_sd / target_sd, "/")
  } else if (distribution == "skew_lognormal") {
    if (q != 1L) stop("skew_lognormal is defined only for q=1.")
    x <- exp(stats::rnorm(n, sd = 0.65))
    x <- (x - mean(x)) / stats::sd(x) * sqrt(Psi[1, 1])
    matrix(x, n, 1)
  } else {
    stop("Unknown random_effect_dist: ", distribution)
  }
}

make_beta_v2 <- function(p, s, signal = 0.75) {
  if (s > p) stop("s cannot exceed p.")
  beta <- numeric(p)
  beta[seq_len(s)] <- signal * (-1)^(seq_len(s) + 1)
  names(beta) <- sprintf("x%05d", seq_len(p))
  beta
}

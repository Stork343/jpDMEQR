check_loss <- function(u, tau) {
  u * (tau - as.numeric(u < 0))
}

.safe_mean <- function(x) {
  if (length(x) == 0) {
    return(NA_real_)
  }
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

.safe_sd <- function(x) {
  if (length(x) <= 1) {
    return(NA_real_)
  }
  if (all(is.na(x))) {
    return(NA_real_)
  }
  sd(x, na.rm = TRUE)
}

.make_signal_beta <- function(p, s, b0) {
  s <- min(s, p)
  beta0 <- rep(0, p)
  if (s > 0) {
    signs <- ifelse(seq_len(s) %% 2 == 0, -1, 1)
    beta0[seq_len(s)] <- b0 * signs
  }
  beta0
}

.scenario_q <- function(scenario) {
  if (scenario %in% c("S2", "S6")) 2L else 1L
}

.qmixnorm <- function(tau, w = 0.9, sd1 = 1, sd2 = 5) {
  if (tau <= 0 || tau >= 1) {
    stop("tau must be in (0,1).")
  }
  f <- function(x) {
    w * pnorm(x, mean = 0, sd = sd1) + (1 - w) * pnorm(x, mean = 0, sd = sd2) - tau
  }
  uniroot(f, interval = c(-20 * sd2, 20 * sd2), tol = 1e-9)$root
}

#' Generate Clustered Mixed-Effects Quantile-Regression Data
#'
#' Simulates clustered data under one of six error scenarios used in the
#' manuscript simulation study.
#'
#' @param n Number of clusters.
#' @param p Number of fixed-effect covariates.
#' @param s Number of non-zero fixed effects.
#' @param b0 Signal magnitude for active fixed effects.
#' @param rho_x Correlation parameter for Toeplitz covariance in `X`.
#' @param tau Target quantile.
#' @param scenario Scenario label in `c("S1","S2","S3","S4","S5","S6")`.
#' @param sigma_gamma2 Random-effect variance.
#' @param delta Heteroscedasticity strength used in scenario `S3`.
#' @param sigma_eps Base error scale.
#' @param m_choices Candidate cluster sizes when `cluster_sizes` is not supplied.
#' @param cluster_sizes Optional explicit cluster sizes.
#' @param beta0 Optional true fixed-effect vector.
#' @param gamma0 Optional true random-effect matrix.
#' @param seed Optional random seed.
#'
#' @return A list containing simulated response, design matrices, and truth.
#' @export
generate_mixed_qr_data <- function(n = 200,
                                   p = 500,
                                   s = 10,
                                   b0 = 0.75,
                                   rho_x = 0.5,
                                   tau = 0.5,
                                   scenario = c("S1", "S2", "S3", "S4", "S5", "S6"),
                                   sigma_gamma2 = 1,
                                   delta = 0.3,
                                   sigma_eps = 1,
                                   m_choices = c(4, 5, 6),
                                   cluster_sizes = NULL,
                                   beta0 = NULL,
                                   gamma0 = NULL,
                                   seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  scenario <- match.arg(scenario)
  if (is.null(cluster_sizes)) {
    cluster_sizes <- sample(m_choices, size = n, replace = TRUE)
  } else {
    cluster_sizes <- as.integer(cluster_sizes)
    n <- length(cluster_sizes)
  }

  N <- sum(cluster_sizes)
  cluster_id <- rep(seq_len(n), times = cluster_sizes)

  Sigma_x <- rho_x^abs(outer(seq_len(p), seq_len(p), "-"))
  X <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = Sigma_x)

  if (is.null(beta0)) {
    beta0 <- .make_signal_beta(p = p, s = s, b0 = b0)
  } else {
    beta0 <- as.numeric(beta0)
    if (length(beta0) != p) {
      stop("beta0 length must equal p.")
    }
  }

  q <- .scenario_q(scenario)
  if (q == 1L) {
    Z <- matrix(1, nrow = N, ncol = 1)
    colnames(Z) <- "intercept"
  } else {
    Z <- cbind(1, X[, 1])
    colnames(Z) <- c("intercept", "x1")
  }

  if (is.null(gamma0)) {
    Psi0 <- diag(rep(sigma_gamma2, q), nrow = q)
    gamma0 <- MASS::mvrnorm(n = n, mu = rep(0, q), Sigma = Psi0)
  } else {
    gamma0 <- as.matrix(gamma0)
    if (nrow(gamma0) != n || ncol(gamma0) != q) {
      stop("gamma0 has incompatible shape.")
    }
  }
  rownames(gamma0) <- as.character(seq_len(n))

  gamma_obs <- gamma0[cluster_id, , drop = FALSE]
  z_gamma <- rowSums(Z * gamma_obs)

  if (scenario == "S1") {
    eps_raw <- rnorm(N, mean = 0, sd = sigma_eps)
    eps <- eps_raw - qnorm(tau, mean = 0, sd = sigma_eps)
  } else if (scenario == "S2") {
    eps_raw <- rt(N, df = 3) / sqrt(3)
    eps <- eps_raw - qt(tau, df = 3) / sqrt(3)
  } else if (scenario == "S3") {
    sigma_i <- 1 + delta * gamma0[, 1]
    sigma_i <- pmax(0.1, sigma_i)
    sigma_obs <- sigma_i[cluster_id]
    eps <- sigma_obs * (rnorm(N) - qnorm(tau))
  } else if (scenario == "S4") {
    # Laplace errors, scale chosen so Var(eps)=sigma_eps^2.
    b_lap <- sigma_eps / sqrt(2)
    eps_raw <- rexp(N, rate = 1 / b_lap) - rexp(N, rate = 1 / b_lap)
    q_tau <- if (tau < 0.5) {
      b_lap * log(2 * tau)
    } else {
      -b_lap * log(2 * (1 - tau))
    }
    eps <- eps_raw - q_tau
  } else if (scenario == "S5") {
    # Skewed chi-square errors, standardized then tau-centered.
    eps_raw <- (rchisq(N, df = 3) - 3) / sqrt(6)
    q_tau <- (qchisq(tau, df = 3) - 3) / sqrt(6)
    eps <- eps_raw - q_tau
  } else {
    # Contaminated normal mixture (outlier-prone).
    is_outlier <- rbinom(N, size = 1, prob = 0.1)
    eps_raw <- (1 - is_outlier) * rnorm(N, mean = 0, sd = 1) +
      is_outlier * rnorm(N, mean = 0, sd = 5)
    q_tau <- .qmixnorm(tau = tau, w = 0.9, sd1 = 1, sd2 = 5)
    eps <- eps_raw - q_tau
  }

  y <- as.numeric(X %*% beta0 + z_gamma + eps)

  list(
    y = y,
    X = X,
    Z = Z,
    cluster_id = cluster_id,
    beta0 = beta0,
    support = which(beta0 != 0),
    gamma0 = gamma0,
    cluster_sizes = cluster_sizes,
    n = n,
    N = N,
    p = p,
    q = q,
    scenario = scenario,
    tau = tau,
    rho_x = rho_x,
    sigma_gamma2 = sigma_gamma2,
    delta = delta,
    b0 = b0,
    s = length(which(beta0 != 0))
  )
}

make_test_data <- function(train_data, n_test = 2000, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n_test_clusters <- max(50L, ceiling(n_test / mean(train_data$cluster_sizes)))
  generate_mixed_qr_data(
    n = n_test_clusters,
    p = train_data$p,
    s = train_data$s,
    b0 = train_data$b0,
    rho_x = train_data$rho_x,
    tau = train_data$tau,
    scenario = train_data$scenario,
    sigma_gamma2 = train_data$sigma_gamma2,
    delta = train_data$delta,
    m_choices = unique(train_data$cluster_sizes),
    beta0 = train_data$beta0,
    seed = seed
  )
}

.empty_method_result <- function(method, p, status, reason) {
  list(
    method = method,
    status = status,
    reason = reason,
    beta_hat = rep(NA_real_, p),
    beta_infer = rep(NA_real_, p),
    se = rep(NA_real_, p),
    screened_set = integer(0),
    selected_set = integer(0),
    runtime_sec = NA_real_,
    diagnostics = list(
      bcd_iterations = NA_real_,
      average_bfgs_steps = NA_real_,
      hess_cond_eff = NA_real_,
      hess_cond_gamma = NA_real_
    )
  )
}

run_method_jp_dme <- function(train_data,
                              tau,
                              c_screen = 3,
                              T_iter = 2,
                              h = NULL,
                              lambda1 = NULL,
                              lambda2 = 1,
                              lambda_clime = NULL,
                              h_mult = 1,
                              lambda1_mult = 0.3,
                              lambda_clime_mult = 1,
                              inference_mode = c("auto", "CLIME", "RIDGE"),
                              verbose = FALSE) {
  inference_mode <- match.arg(inference_mode)
  p <- train_data$p
  n_clusters <- length(unique(train_data$cluster_id))
  N <- train_data$N

  # Use n/log n screening scale to retain enough candidates before penalization.
  d_target <- min(
    p,
    max(1L, floor(c_screen * n_clusters / log(max(3, n_clusters))))
  )
  screened_idx <- CA_IQR_SIS(
    y = train_data$y,
    X = train_data$X,
    cluster_id = train_data$cluster_id,
    tau_grid = sort(unique(c(0.25, 0.5, 0.75, tau))),
    tau_target = tau,
    d_target = d_target,
    T_iter = T_iter,
    verbose = verbose
  )

  if (length(screened_idx) == 0) {
    return(.empty_method_result("JP-DME-QR", p, "failed", "screening returned empty set"))
  }

  X_red <- train_data$X[, screened_idx, drop = FALSE]
  p_red <- ncol(X_red)
  h_base <- if (is.null(h)) N^(-1 / 3) else h
  h <- h_mult * h_base
  lambda1_base <- if (is.null(lambda1)) sqrt(log(max(2, p_red)) / N) else lambda1
  lambda1 <- lambda1_mult * lambda1_base
  lambda_clime_base <- if (is.null(lambda_clime)) default_lambda_clime(p = p_red, N = N, h = h) else lambda_clime
  lambda_clime <- lambda_clime_mult * lambda_clime_base

  fit <- tryCatch(
    JP_DME_QR(
      y = train_data$y,
      X = X_red,
      Z = train_data$Z,
      cluster_id = train_data$cluster_id,
      tau = tau,
      h = h,
      lambda1 = lambda1,
      lambda2 = lambda2,
      verbose = verbose
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(.empty_method_result("JP-DME-QR", p, "failed", fit$message))
  }
  if (any(!is.finite(fit$beta)) || max(abs(fit$beta), na.rm = TRUE) > 100) {
    return(.empty_method_result("JP-DME-QR", p, "failed", "penalised_estimate_unstable"))
  }

  inf_clime <- NULL
  inf_ridge <- NULL
  if (inference_mode %in% c("auto", "CLIME")) {
    inf_clime <- tryCatch(
      Debias_Inference(
        y = train_data$y,
        X = X_red,
        Z = train_data$Z,
        cluster_id = train_data$cluster_id,
        beta_hat = fit$beta,
        gamma_hat = fit$gamma,
        tau = tau,
        h = h,
        lambda2 = lambda2,
        lambda_clime = lambda_clime,
        method = "CLIME",
        verbose = FALSE
      ),
      error = function(e) e
    )
  }
  if (inference_mode %in% c("auto", "RIDGE")) {
    inf_ridge <- tryCatch(
      Debias_Inference(
        y = train_data$y,
        X = X_red,
        Z = train_data$Z,
        cluster_id = train_data$cluster_id,
        beta_hat = fit$beta,
        gamma_hat = fit$gamma,
        tau = tau,
        h = h,
        lambda2 = lambda2,
        method = "RIDGE"
      ),
      error = function(e) e
    )
  }

  inf <- NULL
  inf_reason <- ""
  if (inference_mode == "CLIME") {
    if (inherits(inf_clime, "error")) {
      return(.empty_method_result("JP-DME-QR", p, "failed", inf_clime$message))
    }
    inf <- inf_clime
    inf_reason <- "clime"
  } else if (inference_mode == "RIDGE") {
    if (inherits(inf_ridge, "error")) {
      return(.empty_method_result("JP-DME-QR", p, "failed", inf_ridge$message))
    }
    inf <- inf_ridge
    inf_reason <- "ridge"
  } else {
    clime_ok <- !inherits(inf_clime, "error")
    ridge_ok <- !inherits(inf_ridge, "error")
    if (!clime_ok && !ridge_ok) {
      return(.empty_method_result("JP-DME-QR", p, "failed", "both CLIME and RIDGE inference failed"))
    }
    if (clime_ok && ridge_ok) {
      se_c <- as.numeric(inf_clime$se)
      se_r <- as.numeric(inf_ridge$se)
      se_c <- se_c[is.finite(se_c) & se_c > 0]
      se_r <- se_r[is.finite(se_r) & se_r > 0]
      ratio_med <- if (length(se_c) > 0 && length(se_r) > 0) median(se_c) / median(se_r) else Inf
      if (is.finite(ratio_med) && ratio_med < 0.5) {
        inf <- inf_ridge
        inf_reason <- "auto_ridge_fallback_for_underdispersion"
      } else {
        inf <- inf_clime
        inf_reason <- "auto_clime"
      }
    } else if (clime_ok) {
      inf <- inf_clime
      inf_reason <- "auto_clime_only"
    } else {
      inf <- inf_ridge
      inf_reason <- "auto_ridge_only"
    }
  }

  beta_hat_full <- rep(0, p)
  beta_hat_full[screened_idx] <- fit$beta

  beta_tilde_use <- as.numeric(inf$beta_tilde)
  debias_unstable <- any(!is.finite(beta_tilde_use))
  if (!debias_unstable) {
    ref_scale <- max(1, max(abs(fit$beta), na.rm = TRUE))
    debias_unstable <- max(abs(beta_tilde_use), na.rm = TRUE) > 25 * ref_scale
  }
  if (debias_unstable) {
    beta_tilde_use <- as.numeric(fit$beta)
    se_use <- rep(NA_real_, length(beta_tilde_use))
  } else {
    se_use <- as.numeric(inf$se)
  }

  beta_infer_full <- rep(0, p)
  beta_infer_full[screened_idx] <- beta_tilde_use

  se_full <- rep(NA_real_, p)
  se_full[screened_idx] <- se_use

  selected_set <- which(abs(beta_hat_full) > 1e-8)
  gamma_cond <- inf$H_gg_condition
  gamma_cond <- gamma_cond[is.finite(gamma_cond)]

  list(
    method = "JP-DME-QR",
    status = "ok",
    reason = paste(
      c(
        if (debias_unstable) "debias_unstable_fallback_to_penalised" else NULL,
        inf_reason
      ),
      collapse = ";"
    ),
    beta_hat = beta_hat_full,
    beta_infer = beta_infer_full,
    se = se_full,
    screened_set = screened_idx,
    selected_set = selected_set,
    runtime_sec = NA_real_,
    diagnostics = list(
      bcd_iterations = fit$iterations,
      average_bfgs_steps = fit$average_bfgs_steps,
      hess_cond_eff = inf$H_eff_condition,
      hess_cond_gamma = if (length(gamma_cond) == 0) NA_real_ else mean(gamma_cond)
    )
  )
}

run_method_naive_qr_lasso <- function(train_data, tau, lambda = NULL) {
  p <- train_data$p
  X_aug <- cbind(1, train_data$X)
  colnames(X_aug) <- c("(Intercept)", paste0("x", seq_len(p)))
  if (is.null(lambda)) {
    lambda <- tryCatch(
      quantreg::LassoLambdaHat(X_aug, tau = tau),
      error = function(e) sqrt(log(max(2, p)) / train_data$N)
    )
  }
  lambda <- as.numeric(lambda)[1]
  lambda_vec <- c(0, rep(lambda, p))

  fit <- tryCatch(
    quantreg::rq.fit.lasso(x = X_aug, y = train_data$y, tau = tau, lambda = lambda_vec),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(.empty_method_result("Naive QR-Lasso", p, "failed", fit$message))
  }

  coef_full <- as.numeric(fit$coefficients)
  beta_hat <- coef_full[-1]
  selected_set <- which(abs(beta_hat) > 1e-8)

  list(
    method = "Naive QR-Lasso",
    status = "ok",
    reason = "",
    beta_hat = beta_hat,
    beta_infer = beta_hat,
    se = rep(NA_real_, p),
    screened_set = seq_len(p),
    selected_set = selected_set,
    runtime_sec = NA_real_,
    diagnostics = list(
      bcd_iterations = NA_real_,
      average_bfgs_steps = NA_real_,
      hess_cond_eff = NA_real_,
      hess_cond_gamma = NA_real_
    )
  )
}

run_method_oracle_qr <- function(train_data,
                                 tau,
                                 h = NULL,
                                 lambda2 = 1,
                                 lambda_clime = NULL,
                                 verbose = FALSE) {
  p <- train_data$p
  support <- train_data$support
  if (length(support) == 0) {
    return(.empty_method_result("Oracle QR", p, "failed", "true support is empty"))
  }

  X_oracle <- train_data$X[, support, drop = FALSE]
  oracle_df <- data.frame(y = train_data$y, X_oracle)
  fit <- tryCatch(
    quantreg::rq(y ~ ., data = oracle_df, tau = tau, method = "fn"),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(.empty_method_result("Oracle QR", p, "failed", fit$message))
  }
  coef_fit <- stats::coef(fit)
  beta_oracle <- as.numeric(coef_fit[-1])
  if (length(beta_oracle) != length(support)) {
    return(.empty_method_result("Oracle QR", p, "failed", "oracle coefficient length mismatch"))
  }

  beta_hat_full <- rep(0, p)
  beta_hat_full[support] <- beta_oracle
  beta_infer_full <- rep(0, p)
  beta_infer_full[support] <- beta_oracle
  se_full <- rep(NA_real_, p)

  list(
    method = "Oracle QR",
    status = "ok",
    reason = "",
    beta_hat = beta_hat_full,
    beta_infer = beta_infer_full,
    se = se_full,
    screened_set = support,
    selected_set = support,
    runtime_sec = NA_real_,
    diagnostics = list(
      bcd_iterations = NA_real_,
      average_bfgs_steps = NA_real_,
      hess_cond_eff = NA_real_,
      hess_cond_gamma = NA_real_
    )
  )
}

run_method_two_step_lqmm <- function(train_data, tau) {
  p <- train_data$p
  if (!requireNamespace("lqmm", quietly = TRUE)) {
    return(.empty_method_result("Two-step LQMM", p, "not_available", "package 'lqmm' is not installed"))
  }

  n_clusters <- length(unique(train_data$cluster_id))
  d_target <- min(40L, p, max(5L, floor(2 * sqrt(n_clusters))))
  screened_idx <- CA_IQR_SIS(
    y = train_data$y,
    X = train_data$X,
    cluster_id = train_data$cluster_id,
    tau_grid = sort(unique(c(0.25, 0.5, 0.75, tau))),
    tau_target = tau,
    d_target = d_target,
    T_iter = 2,
    verbose = FALSE
  )
  if (length(screened_idx) == 0) {
    return(.empty_method_result("Two-step LQMM", p, "failed", "screening returned empty set"))
  }

  X_sub <- train_data$X[, screened_idx, drop = FALSE]
  colnames(X_sub) <- paste0("x", seq_len(ncol(X_sub)))
  dat <- data.frame(
    y = train_data$y,
    id = factor(train_data$cluster_id),
    X_sub,
    stringsAsFactors = FALSE
  )
  fixed_formula <- as.formula(paste("y ~", paste(colnames(X_sub), collapse = " + ")))
  lqmm_fit <- getFromNamespace("lqmm", "lqmm")

  fit <- tryCatch(
    lqmm_fit(
      fixed = fixed_formula,
      random = ~ 1,
      group = dat[["id"]],
      tau = tau,
      data = dat
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(.empty_method_result("Two-step LQMM", p, "failed", fit$message))
  }

  coef_fix <- tryCatch(
    {
      cf <- stats::coef(fit)
      if (is.list(cf) && !is.null(cf$fixed)) {
        as.numeric(cf$fixed)
      } else {
        as.numeric(cf)
      }
    },
    error = function(e) NULL
  )
  if (is.null(coef_fix) && !is.null(fit$theta)) {
    coef_fix <- as.numeric(fit$theta)
  }
  if (is.null(coef_fix) || length(coef_fix) < (ncol(X_sub) + 1L)) {
    return(.empty_method_result("Two-step LQMM", p, "failed", "cannot extract fixed-effect coefficients"))
  }

  beta_sub <- coef_fix[seq_len(ncol(X_sub) + 1L)][-1]
  beta_hat <- rep(0, p)
  beta_hat[screened_idx] <- beta_sub

  list(
    method = "Two-step LQMM",
    status = "ok",
    reason = "",
    beta_hat = beta_hat,
    beta_infer = beta_hat,
    se = rep(NA_real_, p),
    screened_set = screened_idx,
    selected_set = which(abs(beta_hat) > 1e-8),
    runtime_sec = NA_real_,
    diagnostics = list(
      bcd_iterations = NA_real_,
      average_bfgs_steps = NA_real_,
      hess_cond_eff = NA_real_,
      hess_cond_gamma = NA_real_
    )
  )
}

run_method_glmm_lasso <- function(train_data, tau) {
  p <- train_data$p
  if (!requireNamespace("glmmLasso", quietly = TRUE)) {
    return(.empty_method_result("GLMM-Lasso", p, "not_available", "package 'glmmLasso' is not installed"))
  }

  n_clusters <- length(unique(train_data$cluster_id))
  d_target <- min(40L, p, max(5L, floor(2 * sqrt(n_clusters))))
  screened_idx <- CA_IQR_SIS(
    y = train_data$y,
    X = train_data$X,
    cluster_id = train_data$cluster_id,
    tau_grid = c(0.25, 0.5, 0.75),
    tau_target = tau,
    d_target = d_target,
    T_iter = 2,
    verbose = FALSE
  )
  if (length(screened_idx) == 0) {
    return(.empty_method_result("GLMM-Lasso", p, "failed", "screening returned empty set"))
  }

  X_sub <- train_data$X[, screened_idx, drop = FALSE]
  colnames(X_sub) <- paste0("x", seq_len(ncol(X_sub)))
  dat <- data.frame(
    y = train_data$y,
    id = factor(train_data$cluster_id),
    X_sub,
    stringsAsFactors = FALSE
  )
  fixed_formula <- as.formula(paste("y ~", paste(colnames(X_sub), collapse = " + ")))
  lambda_glmm <- sqrt(log(max(2, ncol(X_sub))) / nrow(dat))
  glmm_fit <- getFromNamespace("glmmLasso", "glmmLasso")

  fit <- tryCatch(
    glmm_fit(
      fix = fixed_formula,
      rnd = list(id = ~ 1),
      data = dat,
      lambda = lambda_glmm,
      family = gaussian(link = "identity")
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(.empty_method_result("GLMM-Lasso", p, "failed", fit$message))
  }

  coef_full <- tryCatch(as.numeric(fit$coefficients), error = function(e) NULL)
  if (is.null(coef_full) || length(coef_full) < (ncol(X_sub) + 1L)) {
    return(.empty_method_result("GLMM-Lasso", p, "failed", "cannot extract fixed-effect coefficients"))
  }
  beta_sub <- coef_full[seq_len(ncol(X_sub) + 1L)][-1]

  beta_hat <- rep(0, p)
  beta_hat[screened_idx] <- beta_sub
  list(
    method = "GLMM-Lasso",
    status = "ok",
    reason = "",
    beta_hat = beta_hat,
    beta_infer = beta_hat,
    se = rep(NA_real_, p),
    screened_set = screened_idx,
    selected_set = which(abs(beta_hat) > 1e-8),
    runtime_sec = NA_real_,
    diagnostics = list(
      bcd_iterations = NA_real_,
      average_bfgs_steps = NA_real_,
      hess_cond_eff = NA_real_,
      hess_cond_gamma = NA_real_
    )
  )
}

run_method_bayes_mixed_lasso <- function(train_data,
                                         tau,
                                         nIter = 600,
                                         burnIn = 200,
                                         thin = 3) {
  p <- train_data$p
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    return(.empty_method_result("Bayesian Mixed-Lasso", p, "not_available", "package 'BGLR' is not installed"))
  }

  n_clusters <- length(unique(train_data$cluster_id))
  d_target <- min(40L, p, max(5L, floor(2 * sqrt(n_clusters))))
  screened_idx <- CA_IQR_SIS(
    y = train_data$y,
    X = train_data$X,
    cluster_id = train_data$cluster_id,
    tau_grid = c(0.25, 0.5, 0.75),
    tau_target = tau,
    d_target = d_target,
    T_iter = 2,
    verbose = FALSE
  )
  if (length(screened_idx) == 0) {
    return(.empty_method_result("Bayesian Mixed-Lasso", p, "failed", "screening returned empty set"))
  }

  X_sub <- train_data$X[, screened_idx, drop = FALSE]
  # Bayesian lasso part for fixed effects, random intercept part for clusters.
  id_fac <- factor(train_data$cluster_id)
  Z_id <- stats::model.matrix(~ id_fac - 1)
  bglr_fit <- getFromNamespace("BGLR", "BGLR")
  fit <- tryCatch(
    bglr_fit(
      y = train_data$y,
      ETA = list(
        list(X = X_sub, model = "BL"),
        list(X = Z_id, model = "BRR")
      ),
      nIter = as.integer(nIter),
      burnIn = as.integer(burnIn),
      thin = as.integer(thin),
      verbose = FALSE,
      saveAt = tempfile(pattern = "bglr_")
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(.empty_method_result("Bayesian Mixed-Lasso", p, "failed", fit$message))
  }

  beta_sub <- as.numeric(fit$ETA[[1]]$b)
  if (length(beta_sub) != ncol(X_sub)) {
    return(.empty_method_result("Bayesian Mixed-Lasso", p, "failed", "cannot extract BL fixed effects"))
  }

  beta_hat <- rep(0, p)
  beta_hat[screened_idx] <- beta_sub

  sd_sub <- fit$ETA[[1]]$SD.b
  if (!is.null(sd_sub) && length(sd_sub) == length(beta_sub)) {
    selected_local <- which(abs(beta_sub) > 1.96 * as.numeric(sd_sub))
  } else {
    selected_local <- which(abs(beta_sub) > 1e-4)
  }
  selected_set <- screened_idx[selected_local]

  list(
    method = "Bayesian Mixed-Lasso",
    status = "ok",
    reason = "",
    beta_hat = beta_hat,
    beta_infer = beta_hat,
    se = rep(NA_real_, p),
    screened_set = screened_idx,
    selected_set = selected_set,
    runtime_sec = NA_real_,
    diagnostics = list(
      bcd_iterations = NA_real_,
      average_bfgs_steps = NA_real_,
      hess_cond_eff = NA_real_,
      hess_cond_gamma = NA_real_
    )
  )
}

evaluate_method <- function(method_out,
                            train_data,
                            test_data,
                            tau,
                            alpha_levels = c(0.1, 0.05),
                            target_beta = NULL,
                            target_label = "beta0") {
  p <- train_data$p
  if (is.null(target_beta)) {
    target_beta <- train_data$beta0
    target_label <- "beta0"
  } else {
    target_beta <- as.numeric(target_beta)
    if (length(target_beta) != p) {
      stop("target_beta length mismatch.")
    }
  }
  support <- which(abs(target_beta) > 1e-12)
  inactive <- setdiff(seq_len(p), support)
  null_eval <- inactive[seq_len(min(length(inactive), max(1L, length(support))))]

  base_row <- data.frame(
    method = method_out$method,
    target = target_label,
    status = method_out$status,
    reason = method_out$reason,
    l2_error = NA_real_,
    l1_error = NA_real_,
    prediction_risk = NA_real_,
    tpr = NA_real_,
    fpr = NA_real_,
    sure_screen = NA_real_,
    screen_size = NA_real_,
    model_size = NA_real_,
    runtime_sec = method_out$runtime_sec,
    coverage_signal_90 = NA_real_,
    coverage_signal_95 = NA_real_,
    coverage_signal_90_cond = NA_real_,
    coverage_signal_95_cond = NA_real_,
    coverage_signal_90_full = NA_real_,
    coverage_signal_95_full = NA_real_,
    coverage_null_90 = NA_real_,
    coverage_null_95 = NA_real_,
    avg_ci_length_90 = NA_real_,
    avg_ci_length_95 = NA_real_,
    bcd_iterations = method_out$diagnostics$bcd_iterations,
    average_bfgs_steps = method_out$diagnostics$average_bfgs_steps,
    hess_cond_eff = method_out$diagnostics$hess_cond_eff,
    hess_cond_gamma = method_out$diagnostics$hess_cond_gamma,
    stringsAsFactors = FALSE
  )

  if (method_out$status != "ok") {
    return(base_row)
  }

  beta_hat_eval <- method_out$beta_hat
  beta_infer <- method_out$beta_infer
  if (length(beta_hat_eval) != p || length(beta_infer) != p) {
    base_row$status <- "failed"
    base_row$reason <- "beta dimension mismatch"
    return(base_row)
  }

  # Estimation metrics align with the penalized estimator target.
  base_row$l2_error <- sqrt(sum((beta_hat_eval - target_beta)^2))
  base_row$l1_error <- sum(abs(beta_hat_eval - target_beta))

  screened <- unique(as.integer(method_out$screened_set))
  screened <- screened[screened >= 1L & screened <= p]
  if (length(screened) > 0) {
    base_row$screen_size <- length(screened)
    base_row$sure_screen <- as.numeric(all(support %in% screened))
  }

  pred_test <- as.numeric(test_data$X %*% beta_hat_eval)
  base_row$prediction_risk <- mean(check_loss(test_data$y - pred_test, tau))

  selected <- method_out$selected_set
  base_row$model_size <- length(selected)
  base_row$tpr <- if (length(support) == 0) NA_real_ else length(intersect(selected, support)) / length(support)
  base_row$fpr <- if (length(inactive) == 0) NA_real_ else length(setdiff(selected, support)) / length(inactive)

  se <- method_out$se
  if (length(se) == p && any(is.finite(se))) {
    for (alpha in alpha_levels) {
      zcrit <- qnorm(1 - alpha / 2)
      ci_low <- beta_infer - zcrit * se
      ci_up <- beta_infer + zcrit * se

      sig_mask <- is.finite(ci_low[support]) & is.finite(ci_up[support])
      null_mask <- is.finite(ci_low[null_eval]) & is.finite(ci_up[null_eval])

      cov_sig_cond <- if (any(sig_mask)) {
        mean(ci_low[support][sig_mask] <= target_beta[support][sig_mask] &
               ci_up[support][sig_mask] >= target_beta[support][sig_mask])
      } else {
        NA_real_
      }
      # Full-signal coverage: missing/non-finite CI counts as non-coverage.
      cov_sig_full <- {
        hit <- rep(FALSE, length(support))
        hit[sig_mask] <- ci_low[support][sig_mask] <= target_beta[support][sig_mask] &
          ci_up[support][sig_mask] >= target_beta[support][sig_mask]
        mean(hit)
      }
      cov_null <- if (any(null_mask)) {
        mean(ci_low[null_eval][null_mask] <= 0 & ci_up[null_eval][null_mask] >= 0)
      } else {
        NA_real_
      }

      ci_len <- 2 * zcrit * se
      ci_len_mask <- is.finite(ci_len[support])
      avg_len <- if (any(ci_len_mask)) mean(ci_len[support][ci_len_mask]) else NA_real_

      if (abs(alpha - 0.1) < 1e-12) {
        base_row$coverage_signal_90 <- cov_sig_full
        base_row$coverage_signal_90_cond <- cov_sig_cond
        base_row$coverage_signal_90_full <- cov_sig_full
        base_row$coverage_null_90 <- cov_null
        base_row$avg_ci_length_90 <- avg_len
      } else if (abs(alpha - 0.05) < 1e-12) {
        base_row$coverage_signal_95 <- cov_sig_full
        base_row$coverage_signal_95_cond <- cov_sig_cond
        base_row$coverage_signal_95_full <- cov_sig_full
        base_row$coverage_null_95 <- cov_null
        base_row$avg_ci_length_95 <- avg_len
      }
    }
  }

  base_row
}

build_simulation_grid <- function(p_grid = c(200, 500, 1000, 2000),
                                  tau_grid = c(0.25, 0.5, 0.75),
                                  scenarios = c("S1", "S2", "S3", "S4", "S5", "S6"),
                                  s_grid = c(10, 20, 40),
                                  b0_grid = c(0.5, 0.75),
                                  rho_grid = c(0.3, 0.7),
                                  sigma_gamma2_grid = c(0.25, 1, 4),
                                  delta_grid = c(0.3, 0.5),
                                  full_factorial = FALSE) {
  scenarios <- unique(scenarios)
  out <- list()
  idx <- 1L

  for (sc in scenarios) {
    if (full_factorial) {
      local_grid <- expand.grid(
        p = p_grid,
        tau = tau_grid,
        s = s_grid,
        b0 = b0_grid,
        rho_x = rho_grid,
        sigma_gamma2 = sigma_gamma2_grid,
        delta = if (sc == "S3") delta_grid else 0,
        stringsAsFactors = FALSE
      )
    } else {
      local_grid <- expand.grid(
        p = p_grid,
        tau = tau_grid,
        stringsAsFactors = FALSE
      )
      local_grid$s <- s_grid[1]
      local_grid$b0 <- b0_grid[1]
      local_grid$rho_x <- rho_grid[1]
      sigma_idx <- min(2L, length(sigma_gamma2_grid))
      local_grid$sigma_gamma2 <- sigma_gamma2_grid[sigma_idx]
      local_grid$delta <- if (sc == "S3") delta_grid[1] else 0
    }
    local_grid$scenario <- sc
    out[[idx]] <- local_grid
    idx <- idx + 1L
  }

  grid <- do.call(rbind, out)
  rownames(grid) <- NULL
  grid[, c("scenario", "p", "tau", "s", "b0", "rho_x", "sigma_gamma2", "delta")]
}

run_single_replication <- function(config_row,
                                   replication_id,
                                   n = 200,
                                   m_choices = c(4, 5, 6),
                                   n_test = 2000,
                                   c_screen = 3,
                                   T_iter = 2,
                                   h_mult = 1,
                                   lambda1_mult = 1,
                                   lambda_clime_mult = 1,
                                   inference_mode = c("auto", "CLIME", "RIDGE"),
                                   jp_target_beta = NULL,
                                   jp_target_label = "beta0",
                                   methods = c("JP-DME-QR", "Naive QR-Lasso", "Two-step LQMM", "GLMM-Lasso", "Bayesian Mixed-Lasso", "Oracle QR"),
                                   seed = 2026,
                                   verbose = FALSE) {
  inference_mode <- match.arg(inference_mode)
  base_seed <- seed + 100000L * as.integer(replication_id)
  train_data <- generate_mixed_qr_data(
    n = n,
    p = as.integer(config_row$p),
    s = as.integer(config_row$s),
    b0 = as.numeric(config_row$b0),
    rho_x = as.numeric(config_row$rho_x),
    tau = as.numeric(config_row$tau),
    scenario = as.character(config_row$scenario),
    sigma_gamma2 = as.numeric(config_row$sigma_gamma2),
    delta = as.numeric(config_row$delta),
    m_choices = m_choices,
    seed = base_seed
  )
  test_data <- make_test_data(train_data = train_data, n_test = n_test, seed = base_seed + 1L)
  tau <- as.numeric(config_row$tau)

  timed_eval <- function(expr) {
    t0 <- proc.time()[["elapsed"]]
    out <- expr
    out$runtime_sec <- as.numeric(proc.time()[["elapsed"]] - t0)
    out
  }

  method_results <- list()
  if ("JP-DME-QR" %in% methods) {
    method_results[[length(method_results) + 1L]] <- timed_eval(run_method_jp_dme(
      train_data = train_data,
      tau = tau,
      c_screen = c_screen,
      T_iter = T_iter,
      h_mult = h_mult,
      lambda1_mult = lambda1_mult,
      lambda_clime_mult = lambda_clime_mult,
      inference_mode = inference_mode,
      verbose = verbose
    ))
  }
  if ("Naive QR-Lasso" %in% methods) {
    method_results[[length(method_results) + 1L]] <- timed_eval(run_method_naive_qr_lasso(
      train_data = train_data,
      tau = tau
    ))
  }
  if ("Two-step LQMM" %in% methods) {
    method_results[[length(method_results) + 1L]] <- timed_eval(run_method_two_step_lqmm(
      train_data = train_data,
      tau = tau
    ))
  }
  if ("GLMM-Lasso" %in% methods) {
    method_results[[length(method_results) + 1L]] <- timed_eval(run_method_glmm_lasso(
      train_data = train_data,
      tau = tau
    ))
  }
  if ("Bayesian Mixed-Lasso" %in% methods) {
    method_results[[length(method_results) + 1L]] <- timed_eval(run_method_bayes_mixed_lasso(
      train_data = train_data,
      tau = tau
    ))
  }
  if ("Oracle QR" %in% methods) {
    method_results[[length(method_results) + 1L]] <- timed_eval(run_method_oracle_qr(
      train_data = train_data,
      tau = tau,
      verbose = verbose
    ))
  }

  eval_rows <- lapply(
    method_results,
    function(x) {
      if (x$method == "JP-DME-QR" && !is.null(jp_target_beta)) {
        evaluate_method(
          method_out = x,
          train_data = train_data,
          test_data = test_data,
          tau = tau,
          target_beta = jp_target_beta,
          target_label = jp_target_label
        )
      } else {
        evaluate_method(
          method_out = x,
          train_data = train_data,
          test_data = test_data,
          tau = tau
        )
      }
    }
  )
  eval_df <- do.call(rbind, eval_rows)
  eval_df$replication <- replication_id
  eval_df$scenario <- as.character(config_row$scenario)
  eval_df$p <- as.integer(config_row$p)
  eval_df$tau <- as.numeric(config_row$tau)
  eval_df$s <- as.integer(config_row$s)
  eval_df$b0 <- as.numeric(config_row$b0)
  eval_df$rho_x <- as.numeric(config_row$rho_x)
  eval_df$sigma_gamma2 <- as.numeric(config_row$sigma_gamma2)
  eval_df$delta <- as.numeric(config_row$delta)
  eval_df
}

estimate_beta_star_reference <- function(config_row,
                                         n_ref = 400,
                                         b_ref = 3,
                                         m_choices = c(4, 5, 6),
                                         lambda2 = 1,
                                         seed = 90210) {
  p <- as.integer(config_row$p)
  s <- as.integer(config_row$s)
  s <- min(s, p)
  support <- seq_len(s)
  beta_store <- matrix(NA_real_, nrow = b_ref, ncol = s)

  for (b in seq_len(b_ref)) {
    d_ref <- generate_mixed_qr_data(
      n = n_ref,
      p = p,
      s = s,
      b0 = as.numeric(config_row$b0),
      rho_x = as.numeric(config_row$rho_x),
      tau = as.numeric(config_row$tau),
      scenario = as.character(config_row$scenario),
      sigma_gamma2 = as.numeric(config_row$sigma_gamma2),
      delta = as.numeric(config_row$delta),
      m_choices = m_choices,
      seed = seed + 1000L * b
    )

    X_ref <- d_ref$X[, support, drop = FALSE]
    h_ref <- d_ref$N^(-1 / 3)
    fit_ref <- tryCatch(
      JP_DME_QR(
        y = d_ref$y,
        X = X_ref,
        Z = d_ref$Z,
        cluster_id = d_ref$cluster_id,
        tau = d_ref$tau,
        h = h_ref,
        lambda1 = 0,
        lambda2 = lambda2,
        verbose = FALSE
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_ref) && all(is.finite(fit_ref$beta))) {
      beta_store[b, ] <- fit_ref$beta
    }
  }

  beta_target <- rep(0, p)
  if (all(apply(beta_store, 2, function(x) all(!is.finite(x))))) {
    beta_target <- .make_signal_beta(p = p, s = s, b0 = as.numeric(config_row$b0))
    return(beta_target)
  }
  beta_target[support] <- colMeans(beta_store, na.rm = TRUE)
  beta_target
}

#' Summarise Simulation Replications
#'
#' Aggregates replication-level results into mean and standard deviation
#' by scenario and method.
#'
#' @param raw_results Data frame returned by [Simulation_Study()].
#'
#' @return A summary data frame with `_mean` and `_sd` columns.
#' @export
summarise_simulation_results <- function(raw_results) {
  id_cols_pref <- c("scenario", "p", "tau", "s", "b0", "rho_x", "sigma_gamma2", "delta", "method", "target")
  id_cols <- id_cols_pref[id_cols_pref %in% colnames(raw_results)]

  key <- interaction(
    raw_results[, id_cols, drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )
  split_list <- split(raw_results, key)

  metric_cols <- c(
    "l2_error", "l1_error", "prediction_risk", "tpr", "fpr", "sure_screen",
    "screen_size", "model_size", "runtime_sec",
    "coverage_signal_90", "coverage_signal_95",
    "coverage_signal_90_cond", "coverage_signal_95_cond",
    "coverage_signal_90_full", "coverage_signal_95_full",
    "coverage_null_90", "coverage_null_95",
    "avg_ci_length_90", "avg_ci_length_95",
    "bcd_iterations", "average_bfgs_steps", "hess_cond_eff", "hess_cond_gamma"
  )

  summary_rows <- lapply(split_list, function(df) {
    out <- df[1, id_cols, drop = FALSE]
    out$n_runs <- nrow(df)
    out$n_ok <- sum(df$status == "ok")
    ok_df <- df[df$status == "ok", , drop = FALSE]
    for (m in metric_cols) {
      out[[paste0(m, "_mean")]] <- .safe_mean(ok_df[[m]])
      out[[paste0(m, "_sd")]] <- .safe_sd(ok_df[[m]])
    }
    out
  })
  out <- do.call(rbind, summary_rows)
  rownames(out) <- NULL
  out
}

#' Run the Full Monte Carlo Simulation Study
#'
#' Executes the manuscript simulation design across scenario/`p`/`tau` grids
#' and returns both replication-level and aggregated outputs.
#'
#' @param B Number of replications per configuration.
#' @param n Number of clusters in each replication.
#' @param m_choices Candidate cluster sizes.
#' @param n_test Test-set size for prediction-risk evaluation.
#' @param p_grid Grid of fixed-effect dimensions.
#' @param tau_grid Grid of quantile levels.
#' @param scenarios Error scenarios.
#' @param s_grid,b0_grid,rho_grid,sigma_gamma2_grid,delta_grid Scenario parameters.
#' @param full_factorial Logical; whether to use full factorial design.
#' @param methods Methods to benchmark.
#' @param c_screen Screening scaling constant.
#' @param T_iter Number of screening iterations.
#' @param h_mult,lambda1_mult,lambda_clime_mult Multipliers for tuning levels.
#' @param inference_mode Inference backend for `JP-DME-QR`.
#' @param jp_target Target used for JP-DME-QR evaluation.
#' @param n_ref_pseudotrue,b_ref_pseudotrue Settings for pseudo-true reference.
#' @param seed Base random seed.
#' @param verbose Logical; print progress if `TRUE`.
#'
#' @return A list with `config_grid`, `raw_results`, and `summary`.
#' @export
Simulation_Study <- function(B = 500,
                             n = 200,
                             m_choices = c(4, 5, 6),
                             n_test = 2000,
                             p_grid = c(200, 500, 1000, 2000),
                             tau_grid = c(0.25, 0.5, 0.75),
                             scenarios = c("S1", "S2", "S3", "S4", "S5", "S6"),
                             s_grid = c(10, 20, 40),
                             b0_grid = c(0.5, 0.75),
                             rho_grid = c(0.3, 0.7),
                             sigma_gamma2_grid = c(0.25, 1, 4),
                             delta_grid = c(0.3, 0.5),
                             full_factorial = FALSE,
                             methods = c("JP-DME-QR", "Naive QR-Lasso", "Two-step LQMM", "GLMM-Lasso", "Bayesian Mixed-Lasso", "Oracle QR"),
                             c_screen = 3,
                             T_iter = 2,
                             h_mult = 1,
                             lambda1_mult = 0.3,
                             lambda_clime_mult = 1,
                             inference_mode = c("auto", "CLIME", "RIDGE"),
                             jp_target = c("beta_star_ref", "beta0"),
                             n_ref_pseudotrue = 400,
                             b_ref_pseudotrue = 3,
                             seed = 2026,
                             verbose = TRUE) {
  inference_mode <- match.arg(inference_mode)
  jp_target <- match.arg(jp_target)
  config_grid <- build_simulation_grid(
    p_grid = p_grid,
    tau_grid = tau_grid,
    scenarios = scenarios,
    s_grid = s_grid,
    b0_grid = b0_grid,
    rho_grid = rho_grid,
    sigma_gamma2_grid = sigma_gamma2_grid,
    delta_grid = delta_grid,
    full_factorial = full_factorial
  )

  total_jobs <- nrow(config_grid) * B
  job_id <- 1L
  raw_list <- vector("list", total_jobs)

  for (cfg_idx in seq_len(nrow(config_grid))) {
    cfg <- config_grid[cfg_idx, , drop = FALSE]
    if (verbose) {
      cat(sprintf(
        "\n[Config %d/%d] scenario=%s, p=%d, tau=%.2f, s=%d, b0=%.2f, rho=%.2f, sigma_gamma2=%.2f, delta=%.2f\n",
        cfg_idx, nrow(config_grid),
        cfg$scenario, cfg$p, cfg$tau, cfg$s, cfg$b0, cfg$rho_x, cfg$sigma_gamma2, cfg$delta
      ))
    }

    jp_target_beta_cfg <- NULL
    jp_target_label_cfg <- "beta0"
    if ("JP-DME-QR" %in% methods && jp_target == "beta_star_ref") {
      if (verbose) cat("  computing beta_star reference...\n")
      jp_target_beta_cfg <- estimate_beta_star_reference(
        config_row = cfg,
        n_ref = n_ref_pseudotrue,
        b_ref = b_ref_pseudotrue,
        m_choices = m_choices,
        seed = seed + cfg_idx * 10000L
      )
      jp_target_label_cfg <- "beta_star_ref"
    }

    for (b in seq_len(B)) {
      if (verbose && (b %% max(1L, floor(B / 10))) == 0L) {
        cat(sprintf("  replication %d/%d\n", b, B))
      }
      raw_list[[job_id]] <- run_single_replication(
        config_row = cfg,
        replication_id = b,
        n = n,
        m_choices = m_choices,
        n_test = n_test,
        c_screen = c_screen,
        T_iter = T_iter,
        h_mult = h_mult,
        lambda1_mult = lambda1_mult,
        lambda_clime_mult = lambda_clime_mult,
        inference_mode = inference_mode,
        jp_target_beta = jp_target_beta_cfg,
        jp_target_label = jp_target_label_cfg,
        methods = methods,
        seed = seed + cfg_idx * 1000000L,
        verbose = FALSE
      )
      job_id <- job_id + 1L
    }
  }

  raw_results <- do.call(rbind, raw_list)
  summary <- summarise_simulation_results(raw_results)

  list(
    config_grid = config_grid,
    raw_results = raw_results,
    summary = summary
  )
}

#' Export Simulation Outputs to CSV
#'
#' @param sim_out Output list returned by [Simulation_Study()].
#' @param out_dir Output directory.
#' @param prefix File prefix for generated CSV files.
#'
#' @return Invisibly returns paths for raw and summary files.
#' @export
export_simulation_outputs <- function(sim_out, out_dir = "simulation_outputs", prefix = "jp_dme_qr") {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  raw_path <- file.path(out_dir, paste0(prefix, "_raw.csv"))
  summary_path <- file.path(out_dir, paste0(prefix, "_summary.csv"))
  utils::write.csv(sim_out$raw_results, raw_path, row.names = FALSE)
  utils::write.csv(sim_out$summary, summary_path, row.names = FALSE)
  invisible(list(raw = raw_path, summary = summary_path))
}

# Optional quick plotting helpers (base graphics)
#' Plot Coverage Against Dimension
#'
#' @param summary_df Summary data frame from [summarise_simulation_results()].
#' @param method Method name.
#' @param tau Quantile level.
#' @param coverage_col Coverage column name in `summary_df`.
#'
#' @return Invisibly returns `NULL`.
#' @export
plot_coverage_vs_p <- function(summary_df,
                               method = "JP-DME-QR",
                               tau = 0.5,
                               coverage_col = "coverage_signal_95_mean") {
  sub <- summary_df[summary_df$method == method & abs(summary_df$tau - tau) < 1e-12, , drop = FALSE]
  if (nrow(sub) == 0) {
    warning("No rows matched method/tau.")
    return(invisible(NULL))
  }
  sub <- sub[order(sub$p), , drop = FALSE]
  plot(
    x = sub$p,
    y = sub[[coverage_col]],
    type = "b",
    pch = 19,
    xlab = "p",
    ylab = coverage_col,
    main = sprintf("Coverage vs p (%s, tau=%.2f)", method, tau)
  )
  graphics::abline(h = if (grepl("95", coverage_col)) 0.95 else 0.90, lty = 2, col = "gray40")
}

#' Plot `L2` Error Against Dimension
#'
#' @param summary_df Summary data frame from [summarise_simulation_results()].
#' @param method Method name.
#' @param tau Quantile level.
#'
#' @return Invisibly returns `NULL`.
#' @export
plot_l2_error_vs_p <- function(summary_df, method = "JP-DME-QR", tau = 0.5) {
  sub <- summary_df[summary_df$method == method & abs(summary_df$tau - tau) < 1e-12, , drop = FALSE]
  if (nrow(sub) == 0) {
    warning("No rows matched method/tau.")
    return(invisible(NULL))
  }
  sub <- sub[order(sub$p), , drop = FALSE]
  plot(
    x = sub$p,
    y = sub$l2_error_mean,
    type = "b",
    pch = 19,
    xlab = "p",
    ylab = "L2 error (mean)",
    main = sprintf("L2 error vs p (%s, tau=%.2f)", method, tau)
  )
}

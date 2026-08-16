generate_profile_qr_data_v2 <- function(n_clusters,
                                        p,
                                        s,
                                        tau,
                                        q = 1L,
                                        m_values = 3:8,
                                        design_type = "ar1",
                                        rho_x = 0.5,
                                        signal = 0.75,
                                        error_dist = "gaussian",
                                        heteroskedastic = FALSE,
                                        random_effect_dist = "normal",
                                        sigma_b0 = 1,
                                        sigma_b1 = 0.5,
                                        rho_b = 0.4,
                                        seed = 1L) {
  set.seed(seed)
  n_clusters <- as.integer(n_clusters)
  p <- as.integer(p)
  s <- as.integer(s)
  q <- as.integer(q)
  m_values <- as.integer(m_values)
  if (n_clusters < 2 || p < 1 || s < 1 || q < 1) stop("Invalid dimensions.")

  Sigma_x <- make_design_covariance_v2(p, design_type, rho_x)
  Psi_b <- make_random_effect_covariance_v2(q, sigma_b0, sigma_b1, rho_b)
  beta0 <- make_beta_v2(p, s, signal)
  b_mat <- r_random_effects_v2(n_clusters, Psi_b, random_effect_dist)
  colnames(b_mat) <- if (q == 1) "b_intercept" else c("b_intercept", "b_time")

  m_i <- sample(m_values, n_clusters, replace = TRUE)
  total_rows <- sum(m_i)
  X <- matrix(NA_real_, total_rows, p)
  Z <- matrix(NA_real_, total_rows, q)
  y <- numeric(total_rows)
  epsilon <- numeric(total_rows)
  time <- numeric(total_rows)
  cluster_id <- character(total_rows)
  row_cursor <- 1L

  for (i in seq_len(n_clusters)) {
    idx <- row_cursor:(row_cursor + m_i[i] - 1L)
    X_i <- rmvnorm_chol_v2(m_i[i], Sigma_x)
    t_i <- sort(stats::runif(m_i[i], min = -1, max = 1))
    if (q == 1L) {
      Z_i <- matrix(1, m_i[i], 1)
    } else if (q == 2L) {
      Z_i <- cbind(1, t_i)
    } else {
      stop("q>2 is not implemented in the reference DGP.")
    }
    b_i <- as.numeric(b_mat[i, ])
    eps_i <- r_quantile_centered_error_v2(
      n = m_i[i], tau = tau, error_dist = error_dist,
      x1 = X_i[, 1], random_intercept = rep(b_i[1], m_i[i]),
      heteroskedastic = heteroskedastic
    )
    y_i <- as.numeric(X_i %*% beta0 + Z_i %*% b_i + eps_i)

    X[idx, ] <- X_i
    Z[idx, ] <- Z_i
    y[idx] <- y_i
    epsilon[idx] <- eps_i
    time[idx] <- t_i
    cluster_id[idx] <- sprintf("c%05d", i)
    row_cursor <- max(idx) + 1L
  }

  colnames(X) <- names(beta0)
  colnames(Z) <- if (q == 1L) "intercept" else c("intercept", "time")

  list(
    y = y,
    X = X,
    Z = Z,
    cluster_id = cluster_id,
    time = time,
    epsilon = epsilon,
    beta0 = beta0,
    active = seq_len(s),
    null_targets = seq.int(s + 1L, min(p, s + 3L)),
    random_effects = b_mat,
    m_i = m_i,
    Sigma_x = Sigma_x,
    Psi_b = Psi_b,
    n_clusters = n_clusters,
    total_rows = total_rows,
    p = p,
    s = s,
    q = q,
    tau = tau,
    seed = seed,
    dgp = list(
      design_type = design_type,
      rho_x = rho_x,
      signal = signal,
      error_dist = error_dist,
      heteroskedastic = heteroskedastic,
      random_effect_dist = random_effect_dist,
      sigma_b0 = sigma_b0,
      sigma_b1 = sigma_b1,
      rho_b = rho_b,
      m_values = m_values
    )
  )
}

subset_clusters_v2 <- function(dat, clusters) {
  keep <- dat$cluster_id %in% clusters
  X <- dat$X[keep, , drop = FALSE]
  Z <- dat$Z[keep, , drop = FALSE]
  out <- dat
  out$y <- dat$y[keep]
  out$X <- X
  out$Z <- Z
  out$cluster_id <- dat$cluster_id[keep]
  out$time <- dat$time[keep]
  out$epsilon <- dat$epsilon[keep]
  out$n_clusters <- length(unique(out$cluster_id))
  out$total_rows <- sum(keep)
  out
}

make_subject_folds_v2 <- function(cluster_id, k = 5L, seed = 1L) {
  clusters <- unique(as.character(cluster_id))
  set.seed(seed)
  shuffled <- sample(clusters)
  fold <- rep(seq_len(k), length.out = length(shuffled))
  data.frame(cluster_id = shuffled, fold = fold, stringsAsFactors = FALSE)
}

split_clusters_v2 <- function(cluster_id, fraction = 0.5, seed = 1L) {
  clusters <- unique(as.character(cluster_id))
  set.seed(seed)
  n_a <- max(1L, min(length(clusters) - 1L, floor(length(clusters) * fraction)))
  a <- sample(clusters, n_a)
  list(a = a, b = setdiff(clusters, a))
}

quantile_score_screen_v2 <- function(y, X, cluster_id, tau, d,
                                     center_within_cluster = TRUE) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  d <- min(as.integer(d), ncol(X))
  if (center_within_cluster) {
    y_center <- y - ave(y, cluster_id, FUN = stats::median)
    X_center <- X
    for (j in seq_len(ncol(X))) {
      X_center[, j] <- X[, j] - ave(X[, j], cluster_id, FUN = mean)
    }
  } else {
    y_center <- y
    X_center <- X
  }
  q0 <- as.numeric(stats::quantile(y_center, probs = tau, type = 8, names = FALSE))
  psi0 <- tau - as.numeric(y_center <= q0)
  denom <- sqrt(colMeans(X_center^2))
  denom[denom <= 1e-12] <- Inf
  scores <- abs(colMeans(X_center * psi0)) / denom
  selected <- order(scores, decreasing = TRUE)[seq_len(d)]
  list(selected = selected, scores = scores)
}

read_simulation_registry_v2 <- function(path) {
  cfg <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("experiment_id", "n_clusters", "p", "s", "tau", "q",
                "m_values", "error_dist", "replications", "methods")
  missing <- setdiff(required, names(cfg))
  if (length(missing)) stop("Config is missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(cfg$experiment_id)) stop("experiment_id values must be unique.")
  cfg
}

config_row_to_list_v2 <- function(row) {
  row <- as.list(row)
  numeric_names <- c("n_clusters", "p", "s", "tau", "q", "rho_x", "signal",
                     "sigma_b0", "sigma_b1", "rho_b", "lambda_gamma",
                     "h_multiplier", "screen_dim", "replications")
  for (nm in intersect(numeric_names, names(row))) {
    if (!is.na(row[[nm]]) && nzchar(as.character(row[[nm]]))) row[[nm]] <- as.numeric(row[[nm]])
  }
  row$m_values <- parse_semicolon_numeric_v2(row$m_values)
  row$lambda_beta_multipliers <- parse_semicolon_numeric_v2(row$lambda_beta_multipliers)
  row$dantzig_multipliers <- parse_semicolon_numeric_v2(row$dantzig_multipliers)
  row$methods <- parse_pipe_character_v2(row$methods)
  row$heteroskedastic <- tolower(as.character(row$heteroskedastic)) %in% c("true", "1", "yes")
  row
}

simulate_from_config_v2 <- function(config, seed) {
  generate_profile_qr_data_v2(
    n_clusters = config$n_clusters,
    p = config$p,
    s = config$s,
    tau = config$tau,
    q = config$q,
    m_values = config$m_values,
    design_type = config$design_type,
    rho_x = config$rho_x,
    signal = config$signal,
    error_dist = config$error_dist,
    heteroskedastic = config$heteroskedastic,
    random_effect_dist = config$random_effect_dist,
    sigma_b0 = config$sigma_b0,
    sigma_b1 = config$sigma_b1,
    rho_b = config$rho_b,
    seed = seed
  )
}

reference_tuning_values_v2 <- function(config) {
  n <- config$n_clusters
  p <- config$p
  h <- config$h_multiplier * n^(-1 / 3)
  lambda_beta <- config$lambda_beta_multipliers * sqrt(log(max(p, 2)) / n)
  mu <- config$dantzig_multipliers * (sqrt(log(max(p, 2)) / (n * h)) + h^2)
  list(h = h, lambda_beta = lambda_beta, mu = mu,
       lambda_gamma = config$lambda_gamma)
}

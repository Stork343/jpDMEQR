# Higher-level DGP construction, screening utilities, and registry parsing.

draw_cluster_sizes_v2 <- function(n_clusters,
                                  m_values = 3:8,
                                  m_rule = "uniform",
                                  random_intercept = NULL,
                                  latent_size = NULL,
                                  geometric_prob = 0.35) {
  m_values <- sort(unique(as.integer(m_values)))
  if (!length(m_values) || any(m_values < 1L)) stop("m_values must contain positive integers.")
  if (m_rule %in% c("", "uniform")) {
    return(sample(m_values, n_clusters, replace = TRUE))
  }
  if (m_rule == "geometric_imbalance") {
    weights <- geometric_prob * (1 - geometric_prob)^(seq_along(m_values) - 1L)
    weights <- weights / sum(weights)
    return(sample(m_values, n_clusters, replace = TRUE, prob = weights))
  }
  if (m_rule == "informative") {
    if (is.null(random_intercept) || is.null(latent_size)) {
      stop("informative cluster size requires standardised random_intercept and latent_size.")
    }
    b_std <- as.numeric(random_intercept)
    u_std <- as.numeric(latent_size)
    if (length(b_std) != n_clusters || length(u_std) != n_clusters ||
        any(!is.finite(b_std)) || any(!is.finite(u_std))) {
      stop("Informative-size latent variables must be finite and cluster aligned.")
    }
    prob <- stats::plogis(-0.4 + 0.8 * (b_std + u_std) / sqrt(2))
    lower <- min(m_values)
    upper <- max(m_values)
    return(lower + stats::rbinom(n_clusters, size = upper - lower, prob = prob))
  }
  stop("Unknown m_rule: ", m_rule)
}

cluster_size_summary_v2 <- function(m_i) {
  q <- stats::quantile(m_i, probs = c(0.25, 0.5, 0.75), names = FALSE, type = 1)
  list(
    m_min = min(m_i),
    m_q25 = q[1],
    m_median = q[2],
    m_mean = mean(m_i),
    m_q75 = q[3],
    m_max = max(m_i),
    m_sd = if (length(m_i) > 1L) stats::sd(m_i) else 0,
    m_frequency = table(factor(m_i, levels = seq.int(min(m_i), max(m_i))))
  )
}

generate_profile_qr_data_v2 <- function(n_clusters,
                                        p,
                                        s,
                                        tau,
                                        q = 1L,
                                        m_values = 3:8,
                                        m_rule = "uniform",
                                        design_type = "ar1",
                                        rho_x = 0.5,
                                        signal = 0.75,
                                        error_dist = "gaussian",
                                        heteroskedastic = FALSE,
                                        error_dependence = "independent",
                                        copula_rho = 0.4,
                                        random_effect_dist = "normal",
                                        sigma_b0 = 1,
                                        sigma_b1 = 0.5,
                                        rho_b = 0.4,
                                        x_b_corr = 0,
                                        informative_size = FALSE,
                                        nonlinear_re_strength = 0,
                                        seed = 1L) {
  set.seed(seed)
  n_clusters <- as.integer(n_clusters)
  p <- as.integer(p)
  s <- as.integer(s)
  q <- as.integer(q)
  m_values <- as.integer(m_values)
  if (n_clusters < 2L || p < 1L || s < 1L || q < 1L) stop("Invalid dimensions.")
  if (s > p) stop("s cannot exceed p.")
  if (!is.finite(x_b_corr) || abs(x_b_corr) >= 1) stop("x_b_corr must lie in (-1,1).")

  design_sampler <- make_design_sampler_v2(p, design_type, rho_x)
  Psi_b <- make_random_effect_covariance_v2(q, sigma_b0, sigma_b1, rho_b)
  beta0 <- make_beta_v2(p, s, signal)
  b_mat <- r_random_effects_v2(n_clusters, Psi_b, random_effect_dist)
  cluster_names <- sprintf("c%05d", seq_len(n_clusters))
  rownames(b_mat) <- cluster_names
  colnames(b_mat) <- if (q == 1L) "b_intercept" else c("b_intercept", "b_time")
  latent_size <- stats::rnorm(n_clusters)

  effective_m_rule <- if (isTRUE(informative_size)) "informative" else m_rule
  m_i <- draw_cluster_sizes_v2(
    n_clusters = n_clusters,
    m_values = m_values,
    m_rule = effective_m_rule,
    random_intercept = b_mat[, 1] / max(sigma_b0, 1e-8),
    latent_size = latent_size
  )
  total_rows <- sum(m_i)
  cluster_id <- rep(cluster_names, times = m_i)
  cluster_index <- rep(seq_len(n_clusters), times = m_i)
  X <- design_sampler$draw(total_rows)
  b0_by_row <- b_mat[cluster_index, 1]
  latent_by_row <- latent_size[cluster_index]
  if (abs(x_b_corr) > 0) {
    X[, 1] <- sqrt(1 - x_b_corr^2) * X[, 1] +
      x_b_corr * b0_by_row / max(sigma_b0, 1e-8)
  }
  if (isTRUE(informative_size)) {
    X[, 1] <- (X[, 1] + 0.5 * latent_by_row) / sqrt(1.25)
  }

  time_list <- lapply(m_i, function(m) sort(stats::runif(m, min = -1, max = 1)))
  time <- unlist(time_list, use.names = FALSE)
  if (q == 1L) {
    Z <- matrix(1, total_rows, 1L)
  } else if (q == 2L) {
    Z <- cbind(1, time)
  } else {
    stop("q>2 is not implemented in the reference DGP.")
  }

  epsilon <- numeric(total_rows)
  row_start <- cumsum(c(1L, head(m_i, -1L)))
  for (i in seq_len(n_clusters)) {
    idx <- row_start[i]:(row_start[i] + m_i[i] - 1L)
    epsilon[idx] <- r_cluster_quantile_error_v2(
      m = m_i[i], tau = tau, error_dist = error_dist,
      x1 = X[idx, 1], random_intercept = rep(b_mat[i, 1], m_i[i]),
      heteroskedastic = heteroskedastic,
      dependence = error_dependence,
      copula_rho = copula_rho
    )
  }
  nonlinear_component <- nonlinear_re_strength *
    ((b0_by_row / max(sigma_b0, 1e-8))^2 - 1) * time
  nuisance_rows <- rowSums(Z * b_mat[cluster_index, seq_len(q), drop = FALSE])
  y <- as.numeric(X %*% beta0 + nuisance_rows + nonlinear_component + epsilon)

  colnames(X) <- names(beta0)
  colnames(Z) <- if (q == 1L) "intercept" else c("intercept", "time")
  m_summary <- cluster_size_summary_v2(m_i)

  list(
    y = y,
    X = X,
    Z = Z,
    cluster_id = cluster_id,
    time = time,
    epsilon = epsilon,
    nonlinear_component = nonlinear_component,
    beta0 = beta0,
    active = seq_len(s),
    null_targets = if (s < p) seq.int(s + 1L, min(p, s + 6L)) else integer(),
    random_effects = b_mat,
    latent_size = setNames(latent_size, cluster_names),
    m_i = setNames(m_i, cluster_names),
    m_summary = m_summary,
    Sigma_x = design_sampler$covariance,
    design_metadata = design_sampler$metadata,
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
      error_dependence = error_dependence,
      copula_rho = copula_rho,
      random_effect_dist = random_effect_dist,
      sigma_b0 = sigma_b0,
      sigma_b1 = sigma_b1,
      rho_b = rho_b,
      x_b_corr = x_b_corr,
      informative_size = informative_size,
      nonlinear_re_strength = nonlinear_re_strength,
      m_values = m_values,
      m_rule = effective_m_rule
    )
  )
}

subset_clusters_v2 <- function(dat, clusters) {
  clusters <- unique(as.character(clusters))
  keep <- dat$cluster_id %in% clusters
  out <- dat
  out$y <- dat$y[keep]
  out$X <- dat$X[keep, , drop = FALSE]
  out$Z <- dat$Z[keep, , drop = FALSE]
  out$cluster_id <- dat$cluster_id[keep]
  out$time <- dat$time[keep]
  out$epsilon <- dat$epsilon[keep]
  if (!is.null(dat$nonlinear_component)) {
    out$nonlinear_component <- dat$nonlinear_component[keep]
  }
  kept_levels <- unique(out$cluster_id)
  if (!is.null(dat$random_effects)) {
    out$random_effects <- dat$random_effects[kept_levels, , drop = FALSE]
  }
  if (!is.null(dat$latent_size)) out$latent_size <- dat$latent_size[kept_levels]
  out$m_i <- table(factor(out$cluster_id, levels = kept_levels))
  out$m_i <- setNames(as.integer(out$m_i), kept_levels)
  out$m_summary <- cluster_size_summary_v2(out$m_i)
  out$n_clusters <- length(kept_levels)
  out$total_rows <- sum(keep)
  out
}

subset_features_v2 <- function(dat, selected) {
  selected <- sort(unique(as.integer(selected)))
  selected <- selected[selected >= 1L & selected <= ncol(dat$X)]
  if (!length(selected)) stop("Feature subset is empty.")
  original_active <- dat$active
  original_null <- dat$null_targets
  out <- dat
  out$X <- dat$X[, selected, drop = FALSE]
  out$beta0 <- dat$beta0[selected]
  out$active <- which(selected %in% original_active)
  out$null_targets <- which(selected %in% original_null)
  out$p <- ncol(out$X)
  out$feature_map_global <- selected
  out
}

make_subject_folds_v2 <- function(cluster_id, k = 5L, seed = 1L) {
  clusters <- unique(as.character(cluster_id))
  if (length(clusters) < k) stop("Number of clusters is smaller than k.")
  set.seed(seed)
  shuffled <- sample(clusters)
  fold <- rep(seq_len(k), length.out = length(shuffled))
  data.frame(cluster_id = shuffled, fold = fold, stringsAsFactors = FALSE)
}

split_clusters_v2 <- function(cluster_id, fraction = 0.5, seed = 1L) {
  clusters <- unique(as.character(cluster_id))
  if (length(clusters) < 2L) stop("At least two clusters are required for splitting.")
  if (!is.finite(fraction) || fraction <= 0 || fraction >= 1) {
    stop("fraction must lie in (0,1).")
  }
  set.seed(seed)
  n_a <- max(1L, min(length(clusters) - 1L, floor(length(clusters) * fraction)))
  a <- sample(clusters, n_a)
  list(a = a, b = setdiff(clusters, a))
}

quantile_score_screen_v2 <- function(y, X, cluster_id, tau, d,
                                     center_within_cluster = TRUE) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  d <- min(max(as.integer(d), 1L), ncol(X))
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
  list(selected = selected, scores = scores, method = "quantile_score")
}

cluster_weighted_variance_screen_v2 <- function(X, cluster_id, d) {
  X <- as.matrix(X)
  d <- min(max(as.integer(d), 1L), ncol(X))
  clusters <- unique(as.character(cluster_id))
  second_moment <- matrix(0, nrow = length(clusters), ncol = ncol(X))
  for (ii in seq_along(clusters)) {
    Xi <- X[cluster_id == clusters[ii], , drop = FALSE]
    Xi <- sweep(Xi, 2, colMeans(Xi), "-")
    second_moment[ii, ] <- colMeans(Xi^2)
  }
  scores <- colMeans(second_moment)
  selected <- order(scores, decreasing = TRUE)[seq_len(d)]
  list(selected = selected, scores = scores, method = "variance_filter")
}

ca_iqr_sis_screen_v2 <- function(y, X, cluster_id, tau, d,
                                 quantiles = c(0.25, 0.5, 0.75),
                                 rounds = 2L) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  p <- ncol(X)
  d <- min(max(as.integer(d), 1L), p)
  rounds <- max(1L, as.integer(rounds))
  yc <- y - ave(y, cluster_id, FUN = stats::median)
  Xc <- X
  for (j in seq_len(p)) Xc[, j] <- X[, j] - ave(X[, j], cluster_id, FUN = mean)
  residual <- yc
  selected <- integer()
  score_history <- vector("list", rounds)
  engine <- character(rounds)

  for (rr in seq_len(rounds)) {
    available <- setdiff(seq_len(p), selected)
    if (!length(available)) break
    denom <- sqrt(colMeans(Xc[, available, drop = FALSE]^2))
    denom[denom <= 1e-12] <- Inf
    score_mat <- sapply(quantiles, function(qq) {
      q0 <- as.numeric(stats::quantile(residual, probs = qq, type = 8, names = FALSE))
      psi0 <- qq - as.numeric(residual <= q0)
      abs(colMeans(Xc[, available, drop = FALSE] * psi0)) / denom
    })
    if (is.null(dim(score_mat))) score_mat <- matrix(score_mat, ncol = 1L)
    scores <- apply(score_mat, 1, max)
    add_n <- min(length(available), ceiling((d - length(selected)) / (rounds - rr + 1L)))
    if (add_n <= 0L) break
    add <- available[order(scores, decreasing = TRUE)[seq_len(add_n)]]
    selected <- sort(unique(c(selected, add)))
    score_history[[rr]] <- setNames(scores, colnames(X)[available])

    if (requireNamespace("quantreg", quietly = TRUE) && length(selected) < length(y)) {
      fit <- tryCatch(
        quantreg::rq.fit(x = cbind(1, Xc[, selected, drop = FALSE]), y = yc, tau = tau),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        residual <- yc - as.numeric(cbind(1, Xc[, selected, drop = FALSE]) %*% fit$coefficients)
        engine[rr] <- "quantreg::rq.fit"
      } else {
        residual <- stats::lm.fit(cbind(1, Xc[, selected, drop = FALSE]), yc)$residuals
        engine[rr] <- "lm.fit_fallback"
      }
    } else {
      residual <- stats::lm.fit(cbind(1, Xc[, selected, drop = FALSE]), yc)$residuals
      engine[rr] <- "lm.fit"
    }
  }

  if (length(selected) < d) {
    remaining <- setdiff(seq_len(p), selected)
    extra <- head(remaining, d - length(selected))
    selected <- c(selected, extra)
  }
  list(
    selected = selected[seq_len(min(d, length(selected)))],
    score_history = score_history,
    fit_engine = engine,
    method = "ca_iqr_sis"
  )
}

read_simulation_registry_v2 <- function(path) {
  cfg <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "experiment_id", "module", "n_clusters", "p", "s", "tau", "q",
    "m_values", "m_rule", "error_dist", "error_dependence",
    "fit_random_effects", "lambda_gamma", "h_multiplier",
    "lambda_beta_multipliers", "dantzig_multipliers", "screen_method",
    "target_mode", "replications", "methods"
  )
  missing <- setdiff(required, names(cfg))
  if (length(missing)) stop("Config is missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(cfg$experiment_id)) stop("experiment_id values must be unique.")
  cfg
}

config_row_to_list_v2 <- function(row) {
  row <- as.list(row)
  numeric_names <- c(
    "n_clusters", "p", "s", "tau", "q", "rho_x", "signal",
    "copula_rho", "sigma_b0", "sigma_b1", "rho_b", "x_b_corr",
    "nonlinear_re_strength", "lambda_gamma", "h_multiplier", "screen_dim",
    "screen_fraction", "target_coordinate_count", "workers", "replications"
  )
  for (nm in intersect(numeric_names, names(row))) {
    value <- as.character(row[[nm]])
    row[[nm]] <- if (is.na(value) || !nzchar(value)) NA_real_ else as.numeric(value)
  }
  row$m_values <- parse_semicolon_numeric_v2(row$m_values)
  row$lambda_beta_multipliers <- parse_semicolon_numeric_v2(row$lambda_beta_multipliers)
  row$dantzig_multipliers <- parse_semicolon_numeric_v2(row$dantzig_multipliers)
  row$methods <- parse_pipe_character_v2(row$methods)
  bool_names <- c("heteroskedastic", "informative_size", "force_target_inclusion")
  for (nm in intersect(bool_names, names(row))) {
    row[[nm]] <- tolower(as.character(row[[nm]])) %in% c("true", "1", "yes", "y")
  }
  defaults <- list(
    m_rule = "uniform",
    design_type = "ar1",
    rho_x = 0.5,
    signal = 0.75,
    heteroskedastic = FALSE,
    error_dependence = "independent",
    copula_rho = 0,
    random_effect_dist = "normal",
    sigma_b0 = 1,
    sigma_b1 = 0.5,
    rho_b = 0.4,
    x_b_corr = 0,
    informative_size = FALSE,
    nonlinear_re_strength = 0,
    fit_random_effects = "intercept",
    lambda_gamma = 1,
    h_multiplier = 1,
    screen_method = "none",
    screen_dim = NA_real_,
    screen_fraction = 0.5,
    force_target_inclusion = TRUE,
    target_coordinate_count = 4,
    workers = 1,
    target_mode = "structural"
  )
  for (nm in names(defaults)) {
    if (is.null(row[[nm]]) || length(row[[nm]]) == 0L ||
        (length(row[[nm]]) == 1L && (is.na(row[[nm]]) || identical(row[[nm]], "")))) {
      row[[nm]] <- defaults[[nm]]
    }
  }
  row$n_clusters <- as.integer(row$n_clusters)
  row$p <- as.integer(row$p)
  row$s <- as.integer(row$s)
  row$q <- as.integer(row$q)
  row$screen_dim <- if (is.na(row$screen_dim)) NA_integer_ else as.integer(row$screen_dim)
  row$target_coordinate_count <- as.integer(row$target_coordinate_count)
  row$workers <- as.integer(row$workers)
  row$replications <- as.integer(row$replications)
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
    m_rule = config$m_rule,
    design_type = config$design_type,
    rho_x = config$rho_x,
    signal = config$signal,
    error_dist = config$error_dist,
    heteroskedastic = config$heteroskedastic,
    error_dependence = config$error_dependence,
    copula_rho = config$copula_rho,
    random_effect_dist = config$random_effect_dist,
    sigma_b0 = config$sigma_b0,
    sigma_b1 = config$sigma_b1,
    rho_b = config$rho_b,
    x_b_corr = config$x_b_corr,
    informative_size = config$informative_size,
    nonlinear_re_strength = config$nonlinear_re_strength,
    seed = seed
  )
}

reference_tuning_values_v2 <- function(config, p = config$p, n = config$n_clusters) {
  h <- config$h_multiplier * n^(-3 / 10)
  lambda_beta <- config$lambda_beta_multipliers * sqrt(log(max(p, 2)) / n)
  mu <- config$dantzig_multipliers *
    (sqrt(log(max(p, 2)) / (n * h)) + h^2)
  list(
    h = h,
    lambda_beta = lambda_beta,
    mu = mu,
    lambda_gamma = config$lambda_gamma
  )
}

resolve_target_coordinates_v2 <- function(config, dat) {
  count <- min(max(as.integer(config$target_coordinate_count %||% 4L), 1L), ncol(dat$X))
  active_names <- names(dat$beta0)[dat$active]
  null_names <- names(dat$beta0)[dat$null_targets]
  n_active <- min(length(active_names), ceiling(count / 2))
  chosen <- head(active_names, n_active)
  chosen <- unique(c(chosen, head(null_names, count - length(chosen))))
  if (length(chosen) < count) chosen <- unique(c(chosen, head(names(dat$beta0), count)))
  head(chosen, count)
}

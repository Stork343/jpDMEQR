# Metrics, Monte Carlo summaries and output validation.

pinball_loss_v2 <- function(y, qhat, tau, na_rm = FALSE) {
  y <- as.numeric(y)
  qhat <- as.numeric(qhat)
  if (length(y) != length(qhat)) stop("y and qhat length mismatch.")
  u <- y - qhat
  loss <- u * (tau - as.numeric(u < 0))
  if (na_rm) mean(loss, na.rm = TRUE) else mean(loss)
}

# 区间、覆盖指示、簇级损失分位数与内存采样：记录层共用的辅助函数。
normal_interval_multiplier_v2 <- function(level = 0.95) {
  stats::qnorm(1 - (1 - level) / 2)
}

interval_bounds_v2 <- function(estimate, se, level = 0.95) {
  estimate <- as.numeric(estimate)
  se <- as.numeric(se)
  if (length(estimate) != 1L || length(se) != 1L ||
      !is.finite(estimate) || !is.finite(se) || se < 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  half <- normal_interval_multiplier_v2(level) * se
  c(lower = estimate - half, upper = estimate + half)
}

coverage_indicator_v2 <- function(lower, upper, value) {
  if (!is.finite(lower) || !is.finite(upper) || !is.finite(value)) return(NA)
  isTRUE(lower <= value && upper >= value)
}

cluster_loss_quantiles_v2 <- function(y, qhat, cluster_id, tau,
                                      probs = c(0.1, 0.5, 0.9)) {
  y <- as.numeric(y)
  qhat <- as.numeric(qhat)
  if (length(y) != length(qhat)) stop("y and qhat length mismatch.")
  if (length(cluster_id) != length(y)) stop("cluster_id and y length mismatch.")
  labels <- paste0("q", formatC(probs * 100, format = "d"))
  if (any(!is.finite(qhat))) {
    out <- rep(NA_real_, length(probs))
    names(out) <- labels
    return(out)
  }
  u <- y - qhat
  loss <- u * (tau - as.numeric(u < 0))
  per_cluster <- tapply(loss, cluster_id, mean)
  per_cluster <- per_cluster[is.finite(per_cluster)]
  if (!length(per_cluster)) {
    out <- rep(NA_real_, length(probs))
    names(out) <- labels
    return(out)
  }
  out <- stats::quantile(per_cluster, probs = probs, names = FALSE, type = 8)
  names(out) <- labels
  out
}

memory_snapshot_v2 <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) {
    stop("The ps package is required for peak-memory recording.")
  }
  info <- ps::ps_memory_full_info()
  c(rss_mb = as.numeric(info[["rss"]]) / 1048576,
    peak_wset_mb = if ("peak_wset" %in% names(info))
      as.numeric(info[["peak_wset"]]) / 1048576 else NA_real_,
    peak_rss_mb = if ("maxrss" %in% names(info))
      as.numeric(info[["maxrss"]]) / 1048576 else NA_real_)
}

# 推断摘要辅助：极值比例、零假设 p 值分布、studentised 分位数、对数斜率与
# 选择稳定性。汇总只读取记录层已有的列。
extremes_summary_v2 <- function(p_value, alpha) {
  p_value <- p_value[is.finite(p_value)]
  total <- length(p_value)
  if (!total) {
    return(c(scheduled = 0L, evaluable = 0L, rate = NA_real_, mcse = NA_real_,
             binom_lower = NA_real_, binom_upper = NA_real_))
  }
  success <- sum(p_value <= alpha)
  rate <- success / total
  ci <- binomial_interval_v2(success, total)
  c(scheduled = total, evaluable = total, rate = rate,
    mcse = sqrt(rate * (1 - rate) / total),
    binom_lower = ci[1], binom_upper = ci[2])
}

replicate_coverage_summary_v2 <- function(covered, replicate) {
  keep <- is.finite(covered)
  if (!any(keep)) {
    return(c(replications = 0L, rate = NA_real_, mcse = NA_real_))
  }
  per_replicate <- tapply(covered[keep], replicate[keep], mean)
  per_replicate <- per_replicate[is.finite(per_replicate)]
  if (!length(per_replicate)) {
    return(c(replications = 0L, rate = NA_real_, mcse = NA_real_))
  }
  c(replications = length(per_replicate),
    rate = mean(per_replicate),
    mcse = if (length(per_replicate) > 1L) {
      stats::sd(per_replicate) / sqrt(length(per_replicate))
    } else NA_real_)
}

replicate_rate_summary_v2 <- function(p_value, replicate, alpha) {
  keep <- is.finite(p_value)
  if (!any(keep)) {
    return(c(replications = 0L, rate = NA_real_, mcse = NA_real_))
  }
  per_replicate <- tapply(p_value[keep] <= alpha, replicate[keep], mean)
  per_replicate <- per_replicate[is.finite(per_replicate)]
  if (!length(per_replicate)) {
    return(c(replications = 0L, rate = NA_real_, mcse = NA_real_))
  }
  c(replications = length(per_replicate),
    rate = mean(per_replicate),
    mcse = if (length(per_replicate) > 1L) {
      stats::sd(per_replicate) / sqrt(length(per_replicate))
    } else NA_real_)
}

null_pvalue_bins_v2 <- function(p_value, breaks = seq(0, 1, by = 0.1)) {
  p_value <- p_value[is.finite(p_value)]
  if (!length(p_value)) {
    return(data.frame(bin_lower = head(breaks, -1), bin_upper = tail(breaks, -1),
                      count = rep(0L, length(breaks) - 1L),
                      fraction = rep(0, length(breaks) - 1L),
                      expected_fraction = diff(breaks), stringsAsFactors = FALSE))
  }
  counts <- graphics::hist(p_value, breaks = breaks, plot = FALSE)$counts
  data.frame(
    bin_lower = head(breaks, -1),
    bin_upper = tail(breaks, -1),
    count = counts,
    fraction = counts / length(p_value),
    expected_fraction = diff(breaks),
    stringsAsFactors = FALSE
  )
}

ks_uniform_v2 <- function(p_value) {
  p_value <- p_value[is.finite(p_value)]
  if (length(p_value) < 2L) {
    return(c(statistic = NA_real_, mcse_bound = NA_real_, n = length(p_value)))
  }
  statistic <- unname(stats::ks.test(p_value, "punif")$statistic)
  c(statistic = statistic, mcse_bound = 0.5 / sqrt(length(p_value)),
    n = length(p_value))
}

studentised_quantiles_v2 <- function(z, probs = seq(0.01, 0.99, by = 0.01)) {
  z <- z[is.finite(z)]
  if (!length(z)) {
    return(data.frame(probability = probs, empirical = NA_real_,
                      normal = stats::qnorm(probs), stringsAsFactors = FALSE))
  }
  data.frame(
    probability = probs,
    empirical = unname(stats::quantile(z, probs, type = 8)),
    normal = stats::qnorm(probs),
    stringsAsFactors = FALSE
  )
}

loglog_slope_v2 <- function(n, value, min_points = 3L) {
  keep <- is.finite(n) & is.finite(value) & n > 1 & value > 0
  points <- sum(keep)
  if (points < min_points) {
    return(c(slope = NA_real_, slope_mcse = NA_real_, points = points))
  }
  fit <- stats::lm(log(value[keep]) ~ log(n[keep]))
  coefficients <- summary(fit)$coefficients
  c(slope = unname(coefficients[2, 1]),
    slope_mcse = unname(coefficients[2, 2]),
    points = points)
}

selection_stability_v2 <- function(selected_names, active_names,
                                   null_names = character()) {
  selection_sets <- lapply(
    strsplit(selected_names, "|", fixed = TRUE),
    function(x) x[nzchar(x)]
  )
  replications <- length(selection_sets)
  frequency_table <- if (replications) {
    counts <- table(unlist(selection_sets, use.names = FALSE))
    counts / replications
  } else numeric()
  frequencies <- function(names) {
    if (!length(names)) return(numeric())
    # 选择集合全部为空时 table() 没有条目，此时每个候选坐标的频率为 0。
    if (!length(frequency_table)) return(rep(0, length(names)))
    out <- as.numeric(frequency_table[names])
    out[is.na(out)] <- 0
    out
  }
  compact <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(mean = NA_real_, q90 = NA_real_, max = NA_real_))
    c(mean = mean(x),
      q90 = unname(stats::quantile(x, 0.9, names = FALSE)),
      max = max(x))
  }
  active_frequency <- frequencies(active_names)
  null_frequency <- frequencies(null_names)
  active <- compact(active_frequency)
  null <- compact(null_frequency)
  c(active_selection_frequency_mean = active[["mean"]],
    active_selection_frequency_q90 = active[["q90"]],
    active_selection_frequency_max = active[["max"]],
    stable_active_fraction = if (length(active_frequency)) {
      mean(active_frequency >= 0.5)
    } else NA_real_,
    null_selection_frequency_mean = null[["mean"]],
    null_selection_frequency_q90 = null[["q90"]],
    null_selection_frequency_max = null[["max"]],
    unstable_null_count = if (length(null_frequency)) {
      sum(null_frequency >= 0.5)
    } else NA_real_)
}

selection_metrics_v2 <- function(selected, active, p) {
  selected <- sort(unique(as.integer(selected)))
  active <- sort(unique(as.integer(active)))
  selected <- selected[selected >= 1 & selected <= p]
  tp <- length(intersect(selected, active))
  fp <- length(setdiff(selected, active))
  fn <- length(setdiff(active, selected))
  tn <- p - tp - fp - fn
  data.frame(
    selected_size = length(selected),
    true_positive = tp,
    false_positive = fp,
    false_negative = fn,
    true_negative = tn,
    tpr = if (length(active)) tp / length(active) else NA_real_,
    fdp = if (length(selected)) fp / length(selected) else 0,
    fpr = if ((p - length(active)) > 0) fp / (p - length(active)) else NA_real_,
    exact_support = identical(selected, active),
    stringsAsFactors = FALSE
  )
}

estimation_metrics_v2 <- function(beta_hat, beta_target, active = NULL) {
  beta_hat <- as.numeric(beta_hat)
  beta_target <- as.numeric(beta_target)
  if (length(beta_hat) != length(beta_target)) stop("Coefficient length mismatch.")
  err <- beta_hat - beta_target
  active <- active %||% which(beta_target != 0)
  data.frame(
    l1_error = sum(abs(err)),
    l2_error = sqrt(sum(err^2)),
    max_error = max(abs(err)),
    max_active_error = if (length(active)) max(abs(err[active])) else NA_real_,
    active_bias_mean = if (length(active)) mean(err[active]) else NA_real_,
    null_bias_mean = if (length(setdiff(seq_along(err), active))) {
      mean(err[setdiff(seq_along(err), active)])
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}

mean_or_na_v2 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

median_or_na_v2 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

max_or_na_v2 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

coordinate_inference_metrics_v2 <- function(inference_table,
                                            beta_target,
                                            ci_level = 0.95,
                                            target_type = "structural",
                                            target_mc_se = 0) {
  tab <- inference_table
  if (!all(c("coordinate", "beta_hat", "beta_tilde", "se", "ci_lower", "ci_upper") %in% names(tab))) {
    stop("inference_table has missing columns.")
  }
  if (is.null(names(beta_target))) stop("beta_target must be named.")
  target <- as.numeric(beta_target[tab$coordinate])
  if (anyNA(target)) stop("Target missing for one or more coordinates.")
  data.frame(
    coordinate = tab$coordinate,
    target_type = target_type,
    target_value = target,
    target_mc_se = target_mc_se,
    beta_hat = tab$beta_hat,
    beta_tilde = tab$beta_tilde,
    estimated_se = tab$se,
    ci_level = ci_level,
    ci_lower = tab$ci_lower,
    ci_upper = tab$ci_upper,
    covered = tab$ci_lower <= target & tab$ci_upper >= target,
    bias = tab$beta_tilde - target,
    wald_z = (tab$beta_tilde - target) / tab$se,
    interval_length = tab$ci_upper - tab$ci_lower,
    feasible = tab$feasible,
    dantzig_mu = tab$mu,
    dantzig_residual = tab$dantzig_residual,
    omega_l1 = tab$omega_l1,
    omega_l2 = tab$omega_l2,
    stringsAsFactors = FALSE
  )
}

mcse_mean_v2 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

mcse_proportion_v2 <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  p <- mean(x)
  sqrt(p * (1 - p) / length(x))
}

binomial_interval_v2 <- function(success, total, conf_level = 0.95) {
  if (total <= 0) return(c(lower = NA_real_, upper = NA_real_))
  unname(stats::binom.test(success, total, conf.level = conf_level)$conf.int)
}

summarise_numeric_v2 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(data.frame(n = 0L, mean = NA_real_, sd = NA_real_, median = NA_real_,
                      q25 = NA_real_, q75 = NA_real_, mcse = NA_real_))
  }
  data.frame(
    n = length(x),
    mean = mean(x),
    sd = if (length(x) > 1) stats::sd(x) else NA_real_,
    median = stats::median(x),
    q25 = unname(stats::quantile(x, 0.25, type = 8)),
    q75 = unname(stats::quantile(x, 0.75, type = 8)),
    mcse = mcse_mean_v2(x),
    stringsAsFactors = FALSE
  )
}

summarise_coverage_v2 <- function(covered, scheduled = NULL, conf_level = 0.95) {
  valid <- !is.na(covered)
  success <- sum(covered[valid])
  total <- sum(valid)
  scheduled <- scheduled %||% length(covered)
  ci <- binomial_interval_v2(success, total, conf_level)
  data.frame(
    scheduled = scheduled,
    evaluable = total,
    covered = success,
    coverage = if (total) success / total else NA_real_,
    mcse = if (total) sqrt((success / total) * (1 - success / total) / total) else NA_real_,
    binom_lower = ci[1],
    binom_upper = ci[2],
    evaluable_fraction = if (scheduled) total / scheduled else NA_real_,
    stringsAsFactors = FALSE
  )
}

paired_difference_summary_v2 <- function(data, value, method, reference,
                                         id = c("experiment_id", "replicate")) {
  if (!all(c(value, method, id) %in% names(data))) stop("Missing pairing columns.")
  methods <- unique(data[[method]])
  methods <- setdiff(methods, reference)
  out <- vector("list", length(methods))
  for (ii in seq_along(methods)) {
    m <- methods[ii]
    left <- data[data[[method]] == m, c(id, value), drop = FALSE]
    right <- data[data[[method]] == reference, c(id, value), drop = FALSE]
    names(left)[ncol(left)] <- "value_method"
    names(right)[ncol(right)] <- "value_reference"
    paired <- merge(left, right, by = id, all = FALSE)
    d <- paired$value_method - paired$value_reference
    out[[ii]] <- data.frame(
      method = m,
      reference = reference,
      n_pairs = sum(is.finite(d)),
      mean_difference = mean(d, na.rm = TRUE),
      sd_difference = stats::sd(d, na.rm = TRUE),
      mcse_difference = mcse_mean_v2(d),
      median_difference = stats::median(d, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

bahadur_remainder_v2 <- function(beta_tilde, beta_target, influence_sum, n_clusters) {
  sqrt(n_clusters) * (beta_tilde - beta_target) - influence_sum / sqrt(n_clusters)
}

validate_replication_output_v2 <- function(x) {
  required <- c(
    "experiment_id", "replicate", "seed", "method_id", "status",
    "n_clusters", "p", "tau", "runtime_sec", "implementation_commit"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Replication output missing: ", paste(missing, collapse = ", "))
  allowed_status <- c("ok", "warning", "failed", "not_implemented")
  if (!all(x$status %in% allowed_status)) stop("Unknown status value.")
  invisible(TRUE)
}

failure_row_v2 <- function(experiment_id, replicate, seed, method_id,
                           stage, message, metadata = list()) {
  base <- data.frame(
    experiment_id = experiment_id,
    replicate = replicate,
    seed = seed,
    method_id = method_id,
    status = "failed",
    failure_stage = stage,
    failure_message = as.character(message),
    stringsAsFactors = FALSE
  )
  for (nm in names(metadata)) base[[nm]] <- metadata[[nm]]
  base
}

write_atomic_csv_v2 <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.csv(x, tmp, row.names = row.names, na = "")
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
  invisible(path)
}

write_session_info_v2 <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  capture.output(utils::sessionInfo(), file = path)
  invisible(path)
}

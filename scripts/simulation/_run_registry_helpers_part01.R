# Shared registry runner. Source from pilot/main scripts after all R/*_v2.R files.

parse_cli_args_v2 <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- pair[1]
    value <- if (length(pair) > 1) paste(pair[-1], collapse = "=") else "true"
    out[[key]] <- value
  }
  out
}

as_bool_cli_v2 <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}

current_commit_v2 <- function(root = ".") {
  ans <- tryCatch(system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
                  error = function(e) character())
  if (length(ans)) trimws(ans[1]) else "unknown"
}

prepare_analysis_data_v2 <- function(dat, config, seed) {
  if (identical(config$fit_random_effects, "intercept") && ncol(dat$Z) > 1L) {
    dat$Z <- dat$Z[, 1, drop = FALSE]
  }

  selected_global <- seq_len(ncol(dat$X))
  fit_clusters <- unique(dat$cluster_id)
  screening_info <- list(method = config$screen_method, selected_global = selected_global,
                         screen_clusters = character(), fit_clusters = fit_clusters)

  if (config$screen_method == "split_quantile_score") {
    split <- split_clusters_v2(dat$cluster_id, fraction = 0.5, seed = seed + 17L)
    screen_dat <- subset_clusters_v2(dat, split$a)
    d <- as.integer(config$screen_dim)
    scr <- quantile_score_screen_v2(screen_dat$y, screen_dat$X,
                                    screen_dat$cluster_id, config$tau, d)
    selected_global <- scr$selected
    dat <- subset_clusters_v2(dat, split$b)
    dat$X <- dat$X[, selected_global, drop = FALSE]
    dat$beta0 <- dat$beta0[selected_global]
    dat$active <- which(selected_global %in% seq_len(config$s))
    dat$null_targets <- which(selected_global %in% seq.int(config$s + 1L,
                                                            min(config$p, config$s + 3L)))
    screening_info <- list(method = config$screen_method,
                           selected_global = selected_global,
                           scores = scr$scores,
                           screen_clusters = split$a,
                           fit_clusters = split$b)
  } else if (config$screen_method == "ca_iqr_sis") {
    stop("CA-IQR-SIS v2 adapter is not implemented; do not substitute another screen.")
  } else if (!config$screen_method %in% c("none", "")) {
    stop("Unknown screen_method: ", config$screen_method)
  }

  list(data = dat, screening = screening_info)
}

load_profile_target_v2 <- function(root, config, dat) {
  if (identical(config$target_mode, "structural")) {
    return(list(beta = dat$beta0, target_type = "structural", target_mc_se = 0))
  }
  if (identical(config$target_mode, "profile_mc")) {
    path <- file.path(root, "results", "intermediate", "targets",
                      paste0(config$experiment_id, ".rds"))
    if (!file.exists(path)) {
      stop("Profile-MC target is required but missing: ", path)
    }
    obj <- readRDS(path)
    beta <- obj$beta_star_mc
    if (!all(names(dat$beta0) %in% names(beta))) {
      stop("Profile-MC target does not cover the analysis design.")
    }
    return(list(beta = beta[names(dat$beta0)], target_type = "profile_mc",
                target_mc_se = obj$target_mc_se %||% NA_real_))
  }
  stop("Unknown target_mode: ", config$target_mode)
}

common_lqmm_subset_v2 <- function(dat, target_names, tau, max_features = 20L) {
  scr <- quantile_score_screen_v2(dat$y, dat$X, dat$cluster_id, tau,
                                  d = min(max_features, ncol(dat$X)))
  target_idx <- match(target_names, colnames(dat$X))
  target_idx <- target_idx[!is.na(target_idx)]
  keep <- sort(unique(c(scr$selected, target_idx)))
  if (length(keep) > max_features) {
    non_target <- setdiff(keep, target_idx)
    keep <- unique(c(target_idx, head(non_target, max_features - length(target_idx))))
  }
  sub <- dat
  sub$X <- dat$X[, keep, drop = FALSE]
  list(data = sub, keep = keep)
}

map_beta_to_full_v2 <- function(beta, full_names) {
  out <- setNames(rep(0, length(full_names)), full_names)
  if (!is.null(beta)) {
    if (is.null(names(beta))) stop("Benchmark beta_hat must be named.")
    common <- intersect(names(beta), full_names)
    out[common] <- beta[common]
  }
  out
}

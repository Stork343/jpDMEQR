.find_existing_file <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[1]
}

.numeric_like <- function(x, min_prop = 0.5) {
  v <- suppressWarnings(as.numeric(as.character(x)))
  mean(is.finite(v)) >= min_prop
}

.first_numeric_column <- function(df, preferred = character(0)) {
  nm <- names(df)
  for (x in preferred) {
    if (x %in% nm && .numeric_like(df[[x]], min_prop = 0.5)) {
      return(x)
    }
  }
  num_nm <- nm[vapply(df, .numeric_like, logical(1))]
  if (length(num_nm) == 0) {
    return(NA_character_)
  }
  num_nm[1]
}

.first_cluster_column <- function(df, preferred = character(0)) {
  nm <- names(df)
  for (x in preferred) {
    if (x %in% nm) {
      return(x)
    }
  }
  char_nm <- nm[vapply(df, function(z) is.character(z) || is.factor(z), logical(1))]
  if (length(char_nm) == 0) {
    return(NA_character_)
  }
  char_nm[1]
}

.detect_time_slope <- function(pheno, valid_mask) {
  cand <- grep("time|visit|month|week|day|year", names(pheno), ignore.case = TRUE, value = TRUE)
  if (length(cand) == 0) {
    return(NULL)
  }
  for (nm in cand) {
    x <- suppressWarnings(as.numeric(pheno[[nm]]))
    x <- x[valid_mask]
    ok <- is.finite(x)
    if (sum(ok) > 5 && stats::sd(x[ok]) > 1e-8) {
      out <- rep(0, length(valid_mask))
      out[valid_mask] <- as.numeric(scale(x))
      return(list(name = nm, value = out[valid_mask]))
    }
  }
  NULL
}

#' Run the Real-Data Workflow
#'
#' Executes screening, penalised fitting, and debiased inference on
#' GEO series matrix data.
#'
#' @param gse_name Default GEO series matrix filename.
#' @param input_file Optional explicit file path to GEO matrix.
#' @param tau_grid Quantile levels to fit.
#' @param random_effect Random-effect structure: intercept only or intercept+slope.
#' @param d_prefilter Number of probes kept after variance prefiltering.
#' @param c_screen Screening scaling constant.
#' @param T_iter Number of screening iterations.
#' @param seed Random seed.
#' @param max_iter Maximum BCD iterations for model fitting.
#' @param tol Outer-loop convergence tolerance.
#' @param verbose Logical; print progress if `TRUE`.
#'
#' @return A list with metadata, per-quantile results, and diagnostics.
#' @export
run_real_data_analysis <- function(
    gse_name = "GSE121239_series_matrix.txt.gz",
    input_file = NULL,
    tau_grid = c(0.25, 0.5, 0.75),
    random_effect = c("intercept", "intercept_slope"),
    d_prefilter = 5000L,
    c_screen = 3,
    T_iter = 2L,
    seed = 2026,
    max_iter = 200,
    tol = 1e-4,
    verbose = TRUE) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Package 'GEOquery' is required for real-data analysis.")
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package 'Biobase' is required for real-data analysis.")
  }

  random_effect <- match.arg(random_effect)
  tau_grid <- sort(unique(as.numeric(tau_grid)))
  if (any(tau_grid <= 0 | tau_grid >= 1)) {
    stop("tau_grid must be in (0,1).")
  }

  if (is.null(input_file)) {
    file_candidates <- c(
      gse_name,
      file.path("..", "..", "data", gse_name),
      file.path("data", gse_name)
    )
    input_file <- .find_existing_file(file_candidates)
    if (is.na(input_file)) {
      stop(sprintf("Cannot find '%s'. Place it in current dir, ../../data, or pass input_file.", gse_name))
    }
  }

  if (verbose) {
    cat(sprintf("Reading GEO matrix from: %s\n", input_file))
  }
  geo_get <- getFromNamespace("getGEO", "GEOquery")
  biobase_pdata <- getFromNamespace("pData", "Biobase")
  biobase_exprs <- getFromNamespace("exprs", "Biobase")
  eset_obj <- geo_get(filename = input_file, GSEMatrix = TRUE, getGPL = FALSE)
  if (is.list(eset_obj)) {
    eset_obj <- eset_obj[[1]]
  }

  pheno <- biobase_pdata(eset_obj)
  expr_mat <- biobase_exprs(eset_obj)

  y_col <- .first_numeric_column(pheno, preferred = c("sledai:ch1", "SLEDAI:ch1"))
  if (is.na(y_col)) {
    stop("Cannot find a numeric response column in phenotype data.")
  }

  id_col <- .first_cluster_column(pheno, preferred = c("patient id:ch1", "patient_id", "subject_id", "id"))
  if (is.na(id_col)) {
    stop("Cannot find a cluster/subject id column in phenotype data.")
  }

  y_raw <- suppressWarnings(as.numeric(as.character(pheno[[y_col]])))
  cluster_raw <- trimws(as.character(pheno[[id_col]]))
  bad_cluster <- is.na(cluster_raw) |
    cluster_raw == "" |
    tolower(cluster_raw) %in% c("na", "nan", "null")

  d_prefilter <- min(as.integer(d_prefilter), nrow(expr_mat))
  probe_var <- apply(expr_mat, 1, var)
  keep_idx <- order(probe_var, decreasing = TRUE)[seq_len(d_prefilter)]
  expr_filt <- expr_mat[keep_idx, , drop = FALSE]
  X_raw <- t(expr_filt)
  X_raw <- scale(X_raw)

  valid <- is.finite(y_raw) & !bad_cluster
  valid <- valid & rowSums(is.na(X_raw)) == 0

  y <- as.numeric(y_raw[valid])
  cluster_id <- as.character(cluster_raw[valid])
  X <- as.matrix(X_raw[valid, , drop = FALSE])

  if (random_effect == "intercept_slope") {
    slope_info <- .detect_time_slope(pheno, valid)
  } else {
    slope_info <- NULL
  }

  if (!is.null(slope_info)) {
    Z <- cbind(1, as.numeric(slope_info$value))
    colnames(Z) <- c("intercept", "time_std")
    if (verbose) {
      cat(sprintf("Random-effects design: intercept + slope(%s)\n", slope_info$name))
    }
  } else {
    Z <- matrix(1, nrow = length(y), ncol = 1)
    colnames(Z) <- "intercept"
    if (verbose && random_effect == "intercept_slope") {
      cat("Random slope requested but no valid time-like variable was detected; using intercept-only.\n")
    }
  }

  N <- length(y)
  n_clusters <- length(unique(cluster_id))
  stopifnot(N == nrow(X), N == nrow(Z), N == length(cluster_id))

  tau_target_screen <- 0.5
  d_target <- min(300L, floor(c_screen * sqrt(n_clusters)), ncol(X))
  set.seed(seed)
  screened_idx <- CA_IQR_SIS(
    y = y,
    X = X,
    cluster_id = cluster_id,
    tau_grid = tau_grid,
    tau_target = tau_target_screen,
    d_target = d_target,
    T_iter = as.integer(T_iter),
    verbose = verbose
  )

  if (length(screened_idx) == 0) {
    stop("CA-IQR-SIS returned an empty set on real data.")
  }

  X_reduced <- X[, screened_idx, drop = FALSE]
  p_red <- ncol(X_reduced)
  h <- N^(-1 / 3)
  lambda1 <- sqrt(log(max(2, p_red)) / N)
  lambda2 <- 1
  lambda_clime <- default_lambda_clime(p = p_red, N = N, h = h)

  probe_ids_prefilt <- rownames(expr_filt)
  probe_ids_sel <- probe_ids_prefilt[screened_idx]

  fit_list <- vector("list", length(tau_grid))
  inf_list <- vector("list", length(tau_grid))
  result_tables <- vector("list", length(tau_grid))
  names(fit_list) <- names(inf_list) <- names(result_tables) <- sprintf("tau_%.2f", tau_grid)

  for (k in seq_along(tau_grid)) {
    tau <- tau_grid[k]
    if (verbose) {
      cat(sprintf("\n[Real-data] Fitting tau = %.2f ...\n", tau))
    }

    fit_k <- JP_DME_QR(
      y = y,
      X = X_reduced,
      Z = Z,
      cluster_id = cluster_id,
      tau = tau,
      h = h,
      lambda1 = lambda1,
      lambda2 = lambda2,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose
    )
    inf_k <- Debias_Inference(
      y = y,
      X = X_reduced,
      Z = Z,
      cluster_id = cluster_id,
      beta_hat = fit_k$beta,
      gamma_hat = fit_k$gamma,
      tau = tau,
      h = h,
      lambda2 = lambda2,
      lambda_clime = lambda_clime,
      method = "CLIME",
      verbose = FALSE
    )

    beta_tilde <- inf_k$beta_tilde
    se_beta <- inf_k$se
    z_score <- beta_tilde / se_beta
    pval <- 2 * (1 - pnorm(abs(z_score)))

    result_table <- data.frame(
      tau = tau,
      probe = probe_ids_sel,
      beta = beta_tilde,
      se = se_beta,
      z = z_score,
      pval = pval,
      stringsAsFactors = FALSE
    )
    result_table$adj_pval <- p.adjust(result_table$pval, method = "BH")
    z90 <- qnorm(0.95)
    z95 <- qnorm(0.975)
    result_table$ci90_low <- result_table$beta - z90 * result_table$se
    result_table$ci90_high <- result_table$beta + z90 * result_table$se
    result_table$ci95_low <- result_table$beta - z95 * result_table$se
    result_table$ci95_high <- result_table$beta + z95 * result_table$se
    result_table <- result_table[order(result_table$pval), ]

    fit_list[[k]] <- fit_k
    inf_list[[k]] <- inf_k
    result_tables[[k]] <- result_table

    if (verbose) {
      cat(sprintf("Top findings at tau=%.2f (first 10 rows):\n", tau))
      print(utils::head(result_table, 10))
      cat(sprintf("#Significant probes (BH-FDR < 0.05): %d\n", sum(result_table$adj_pval < 0.05, na.rm = TRUE)))
    }
  }

  combined_table <- do.call(rbind, result_tables)
  rownames(combined_table) <- NULL

  diagnostics <- data.frame(
    tau = tau_grid,
    bcd_iterations = vapply(fit_list, function(x) x$iterations, numeric(1)),
    avg_bfgs_steps = vapply(fit_list, function(x) x$average_bfgs_steps, numeric(1)),
    hess_cond_eff = vapply(inf_list, function(x) x$H_eff_condition, numeric(1)),
    hess_cond_gamma = vapply(
      inf_list,
      function(x) {
        gg <- x$H_gg_condition
        gg <- gg[is.finite(gg)]
        if (length(gg) == 0) NA_real_ else mean(gg)
      },
      numeric(1)
    )
  )

  list(
    metadata = list(
      input_file = input_file,
      response_column = y_col,
      cluster_column = id_col,
      random_effect = if (ncol(Z) == 1) "intercept" else "intercept_slope",
      screened_size = p_red,
      N = N,
      n_clusters = n_clusters,
      h = h,
      lambda1 = lambda1,
      lambda2 = lambda2,
      lambda_clime = lambda_clime
    ),
    screened_idx = screened_idx,
    probe_ids_selected = probe_ids_sel,
    fit = fit_list,
    inference = inf_list,
    result_tables = result_tables,
    combined_table = combined_table,
    diagnostics = diagnostics
  )
}

# Example:
# setwd("code/code_v2/real_data")
# res_real <- run_real_data_analysis()

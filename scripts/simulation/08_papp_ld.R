#!/usr/bin/env Rscript
# PAPP-LD runner v2 (optimised): full delete-one-cluster jackknife computed ONCE
# per replication (all coordinates from the same leave-out fits), warm-started
# from the full-data fit, parallelised over leave-out clusters on Linux.
# Usage: Rscript scripts/simulation/08_papp_ld.R [shard_index] [n_shards]
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/08_papp_ld.R"
root <- Sys.getenv("JPDMEQR_ROOT", unset = NA)
if (is.na(root) || !nzchar(root)) {
  root <- if (file.exists("R/profile_v2.R")) "." else
    normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
}
setwd(root)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source_v2_module(root, "benchmark_adapters_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

argv <- commandArgs(trailingOnly = TRUE)
shard <- as.integer(argv[1] %||% "1")
n_shards <- as.integer(argv[2] %||% "1")
B <- 200L
taus <- c(0.25, 0.50, 0.75)
multiset <- read.csv(file.path(root, "config", "gse65391_visit_multiset.csv"))$cluster_size
out_dir <- file.path(root, "results", "raw", "simulation", "papp_ld_B200")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("papp_ld_s%d.csv", shard))
jack_cores <- as.integer(Sys.getenv("JPDMEQR_LD_CORES", unset = "2"))

rows <- list()
k <- 0L
for (tau in taus) {
  reps <- seq.int(from = shard, to = B, by = n_shards)
  for (r in reps) {
    seed <- seed_from_id_v2(20260817L, paste0("PAPP-LD", tau * 100), r, "training")
    dat <- generate_profile_qr_data_v2(
      n_clusters = 129L, p = 20L, s = 5L, tau = tau, q = 1L,
      m_values = 2:22, m_multiset = multiset,
      design_type = "ar1", rho_x = 0.5, signal = 0.75,
      error_dist = "gaussian", seed = seed
    )
    y <- dat$y; X <- dat$X; Z <- dat$Z; cid <- dat$cluster_id
    beta0 <- dat$beta0
    h_inf <- 129^(-3 / 10)
    nc <- list(reltol = 1e-11, maxit = 400L, grad_tol = 1e-8)
    fit <- tryCatch(
      fit_profile_lasso_v2(y, X, Z, cid, tau, h_inf, lambda_beta = 0,
                           lambda_gamma = 1, penalty_factor = rep(0, 20L),
                           control = list(max_iter = 1000L, max_backtrack = 50L,
                                          beta_tol = 1e-8, kkt_normalized_tol = 1e-6,
                                          nuisance_control = nc)),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      for (coord in colnames(X)[1:4]) {
        k <- k + 1L
        rows[[k]] <- data.frame(tau = tau, replicate = r, seed = seed, coordinate = coord,
                                beta_hat = NA_real_, sandwich_se = NA_real_, jackknife_se = NA_real_,
                                sand_cov = NA, jack_cov = NA, status = paste("error:", conditionMessage(fit)),
                                stringsAsFactors = FALSE)
      }
      cat(sprintf("shard %d: tau=%.2f rep %d ERROR %s\n", shard, tau, r,
                  conditionMessage(fit)), file = stderr())
      next
    }
    # ONE jackknife pass over leave-out clusters, warm-started from the full fit.
    rows_by_cluster <- split(seq_along(cid), cid)
    jack_fun <- function(ix) {
      f <- tryCatch(fit_profile_lasso_v2(
        y[-ix], X[-ix, , drop = FALSE], Z[-ix, , drop = FALSE], cid[-ix],
        tau, h_inf, lambda_beta = 0, lambda_gamma = 1, penalty_factor = rep(0, 20L),
        beta_init = fit$beta,
        control = list(max_iter = 300L, beta_tol = 1e-7, kkt_normalized_tol = 1e-5)
      ), error = function(e) NULL)
      if (is.null(f)) rep(NA_real_, 20L) else f$beta
    }
    jack_beta <- if (jack_cores > 1L && .Platform$OS.type != "windows") {
      do.call(rbind, parallel::mclapply(rows_by_cluster, jack_fun, mc.cores = jack_cores))
    } else {
      do.call(rbind, lapply(rows_by_cluster, jack_fun))
    }
    G0 <- fit$components$unsmoothed_cluster_scores
    Hf <- fit$components$hessian
    for (coord in colnames(X)[1:4]) {
      kk <- match(coord, colnames(X))
      e <- numeric(20L); e[kk] <- 1
      omega <- tryCatch(as.numeric(solve((Hf + t(Hf)) / 2, e)), error = function(e) NULL)
      jvals <- jack_beta[, kk]
      jok <- sum(is.finite(jvals)) > 1L
      jack_se <- if (jok) {
        j <- jvals[is.finite(jvals)]
        sqrt((length(j) - 1) / length(j) * sum((j - mean(j))^2))
      } else NA_real_
      sand_se <- if (!is.null(omega)) {
        proj <- as.numeric(G0 %*% omega)
        sqrt(mean((proj - mean(proj))^2) / 129L)
      } else NA_real_
      beta_hat <- fit$beta[kk]
      tgt <- as.numeric(beta0[kk])
      sand_cov <- if (is.finite(sand_se)) abs(beta_hat - tgt) <= 1.96 * sand_se else NA
      jack_cov <- if (is.finite(jack_se)) abs(beta_hat - tgt) <= 1.96 * jack_se else NA
      k <- k + 1L
      rows[[k]] <- data.frame(tau = tau, replicate = r, seed = seed, coordinate = coord,
                              beta_hat = beta_hat, sandwich_se = sand_se, jackknife_se = jack_se,
                              sand_cov = sand_cov, jack_cov = jack_cov, status = "ok",
                              stringsAsFactors = FALSE)
    }
    cat(sprintf("shard %d: tau=%.2f rep %d done\n", shard, tau, r), file = stderr())
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, out_file, row.names = FALSE)
cat("PAPP-LD shard", shard, "written:", nrow(tab), "rows\n", file = stderr())
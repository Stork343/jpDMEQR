#!/usr/bin/env Rscript
# 每次replication完成全部leave-one-subject拟合后计算jackknife。
# 使用方式 Rscript scripts/simulation/08_papp_ld.R [shard_index] [n_shards] [run_id]
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
run_id <- if (length(argv) >= 3L) argv[3] else
  paste0("papp_ld_full_jackknife_", format(Sys.time(), "%Y%m%dT%H%M%S"))
out_dir <- file.path(root, "results", "raw", "simulation", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("papp_ld_s%d.csv", shard))
jack_cores <- as.integer(Sys.getenv("JPDMEQR_LD_CORES", unset = "2"))

blocks <- list()
for (tau in taus) {
  rows <- list()
  k <- 0L
  # 每个tau保存独立检查点，恢复时验证所有预定replication。
  ck_file <- sub("\\.csv$", sprintf("_tau%.0f.csv", tau * 100), out_file)
  if (file.exists(ck_file)) {
    cat(sprintf("shard %d: tau=%.2f checkpoint exists, skipping\n", shard, tau),
        file = stderr())
    checkpoint <- read.csv(ck_file, stringsAsFactors = FALSE)
    expected_reps <- seq.int(from = shard, to = B, by = n_shards)
    if (!all(checkpoint$tau == tau) || nrow(checkpoint) != 4L * length(expected_reps) ||
        !setequal(checkpoint$replicate, expected_reps) ||
        anyDuplicated(checkpoint[c("tau", "replicate", "coordinate")]) ||
        !"jackknife_completed" %in% names(checkpoint) ||
        !all(checkpoint$jackknife_completed == 129L |
               (grepl("^error:", checkpoint$status) & checkpoint$jackknife_completed == 0L))) {
      stop("PAPP_CHECKPOINT_INVALID: checkpoint does not contain the complete frozen tau block.")
    }
    blocks[[length(blocks) + 1L]] <- checkpoint
    next
  }
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
                                sand_cov = NA, jack_cov = NA, jackknife_completed = 0L,
                                score_source = "smoothed", h = h_inf,
                                target_type = "structural",
                                score_definition = "gradient_of_complete_smoothed_ridge_profile_criterion",
                                status = paste("error:", conditionMessage(fit)),
                                stringsAsFactors = FALSE)
      }
      cat(sprintf("shard %d: tau=%.2f rep %d ERROR %s\n", shard, tau, r,
                  conditionMessage(fit)), file = stderr())
      next
    }
    if (!isTRUE(fit$converged) || any(!is.finite(fit$beta))) {
      stop("PAPP_FIT_FAILED: full-data unpenalised fit did not converge.")
    }
    # 每次删除一个完整subject，所有leave-one拟合都必须收敛。
    rows_by_cluster <- split(seq_along(cid), cid)
    jack_fun <- function(ix) {
      f <- fit_profile_lasso_v2(
        y[-ix], X[-ix, , drop = FALSE], Z[-ix, , drop = FALSE], cid[-ix],
        tau, h_inf, lambda_beta = 0, lambda_gamma = 1, penalty_factor = rep(0, 20L),
        beta_init = fit$beta,
        control = list(max_iter = 1000L, max_backtrack = 50L,
                       beta_tol = 1e-8, kkt_normalized_tol = 1e-6,
                       nuisance_control = nc)
      )
      if (!isTRUE(f$converged) || any(!is.finite(f$beta)) ||
          f$components$max_nuisance_gradient > 1e-7) {
        stop("JACKKNIFE_FIT_FAILED: leave-one-subject profile fit did not converge.")
      }
      f$beta
    }
    jack_beta <- if (jack_cores > 1L && .Platform$OS.type != "windows") {
      do.call(rbind, parallel::mclapply(rows_by_cluster, jack_fun, mc.cores = jack_cores))
    } else {
      do.call(rbind, lapply(rows_by_cluster, jack_fun))
    }
    if (nrow(jack_beta) != 129L || ncol(jack_beta) != 20L ||
        any(!is.finite(jack_beta))) {
      stop("JACKKNIFE_INCOMPLETE: all 129 leave-one-subject fits are required.")
    }
    G0 <- fit$components$cluster_scores
    Hf <- fit$components$hessian
    for (coord in colnames(X)[1:4]) {
      kk <- match(coord, colnames(X))
      e <- numeric(20L); e[kk] <- 1
      omega <- as.numeric(solve((Hf + t(Hf)) / 2, e))
      jvals <- jack_beta[, kk]
      jack_se <- sqrt(128 / 129 * sum((jvals - mean(jvals))^2))
      proj <- as.numeric(G0 %*% omega)
      sand_se <- sqrt(mean((proj - mean(proj))^2) / 129L)
      beta_hat <- fit$beta[kk]
      tgt <- as.numeric(beta0[kk])
      sand_cov <- if (is.finite(sand_se)) abs(beta_hat - tgt) <= 1.96 * sand_se else NA
      jack_cov <- if (is.finite(jack_se)) abs(beta_hat - tgt) <= 1.96 * jack_se else NA
      k <- k + 1L
      rows[[k]] <- data.frame(tau = tau, replicate = r, seed = seed, coordinate = coord,
                              beta_hat = beta_hat, sandwich_se = sand_se, jackknife_se = jack_se,
                              sand_cov = sand_cov, jack_cov = jack_cov,
                              jackknife_completed = nrow(jack_beta),
                              score_source = "smoothed", h = h_inf,
                              target_type = "structural",
                              score_definition = "gradient_of_complete_smoothed_ridge_profile_criterion",
                              status = "ok",
                              stringsAsFactors = FALSE)
    }
    cat(sprintf("shard %d: tau=%.2f rep %d done\n", shard, tau, r), file = stderr())
  }
  # 当前tau完成后保存检查点，既有结果文件保持原状。
  block <- do.call(rbind, rows)
  write.csv(block, ck_file, row.names = FALSE)
  blocks[[length(blocks) + 1L]] <- block
  cat(sprintf("shard %d: tau=%.2f checkpoint written (%d rows)\n", shard, tau,
              length(rows)), file = stderr())
}
tab <- do.call(rbind, blocks)
if (file.exists(out_file)) stop("PAPP_OUTPUT_EXISTS: completed shard cannot be overwritten.")
write.csv(tab, out_file, row.names = FALSE)
cat("PAPP-LD shard", shard, "written:", nrow(tab), "rows\n", file = stderr())

# -----------------------------------------------------------------------------
# Unsmoothed quantile score (primary Wald meat)
# -----------------------------------------------------------------------------

# psi_tau(u) = tau - I(u < 0), with the deterministic convention I(u < 0) so
# that psi(0) = tau. Used for the primary cluster meat of the unsmoothed
# regularised profile target (theory decision: PILOT_GATE_THEORY_DECISION.md).
psi_quantile_v2 <- function(u, tau) {
  tau - as.numeric(u < 0)
}

profile_components_v2 <- function(y,
                                  X,
                                  Z,
                                  cluster_id,
                                  beta,
                                  tau,
                                  h,
                                  lambda_gamma = 1,
                                  gamma_init = NULL,
                                  need_hessian = TRUE,
                                  identity_check = FALSE,
                                  nuisance_control = list(),
                                  assembly = NULL,
                                  need_fold_sums = TRUE) {
  validate_tau_h_v2(tau, h)
  y <- as.numeric(y)
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  beta <- as.numeric(beta)
  N <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  if (nrow(X) != N || nrow(Z) != N || length(cluster_id) != N) {
    stop("y, X, Z and cluster_id dimensions do not match.")
  }
  if (length(beta) != p) stop("beta has wrong length.")

  Lambda <- make_lambda_matrix_v2(lambda_gamma, q)
  ci <- cluster_index_v2(cluster_id)
  n <- ci$n

  # 两式 Hessian 的比较值只服务于几何审计：默认不构造，中间迭代因此省下一次
  # p x p 的 crossprod 和一份 p x p 的驻留矩阵。
  want_identity <- isTRUE(need_hessian) && isTRUE(identity_check)

  # 装配模式 sequential 逐簇累加 p x p 矩阵，求和次序与冻结结果相同；
  # 装配模式 batch 逐簇只保留 w、X'WZ 与 A 的行，随后用整块 crossprod 装配，
  # 大 p 下 BLAS 的效率显著提高，代价是求和次序改变，Hessian 出现 1e-15 量级的
  # 相对差异。默认 sequential，由 jpDMEQR.assembly 选项或 JPDMEQR_ASSEMBLY
  # 环境变量选择。
  assembly_mode <- assembly %||% Sys.getenv(
    "JPDMEQR_ASSEMBLY", unset = getOption("jpDMEQR.assembly", "sequential"))
  if (!assembly_mode %in% c("sequential", "batch")) {
    stop("Unknown assembly mode: ", assembly_mode,
         ". Allowed values: sequential, batch.")
  }
  batch_assembly <- isTRUE(need_hessian) && identical(assembly_mode, "batch")
  identity_detail <- getOption("jpDMEQR.identity_detail",
                               Sys.getenv("JPDMEQR_IDENTITY_DETAIL", unset = "cluster"))
  if (!identity_detail %in% c("aggregate", "cluster")) stop("Unknown identity detail.")
  want_cluster_identity <- want_identity &&
    (!batch_assembly || identical(identity_detail, "cluster"))

  # 折 Hessian 只被 mu 交叉验证读取。拟合迭代中的 Hessian 刷新随后即被丢弃，
  # 折内交叉乘积因此由调用者声明是否需要，缺省保留退回原行为。
  want_fold_sums <- isTRUE(need_hessian) && isTRUE(need_fold_sums)

  cluster_scores <- matrix(NA_real_, nrow = n, ncol = p,
                           dimnames = list(ci$levels, colnames(X)))
  unsmoothed_cluster_scores <- cluster_scores
  H_sum <- if (need_hessian) matrix(0, p, p) else NULL
  H_identity_sum <- if (want_identity) matrix(0, p, p) else NULL
  hessian_fold_sums <- if (want_fold_sums) rep(list(matrix(0, p, p)), 4L) else NULL
  hessian_fold_counts <- if (want_fold_sums) integer(4L) else NULL
  if (batch_assembly) {
    # 逐簇返回的量：w（长度 N）、X'WZ 与 A 的行（各 q 份 n x p）、折编号，
    # 以及两式比较所需的 X_tilde 与 L A 的行。
    batch_w <- numeric(N)
    batch_U <- lapply(seq_len(q), function(cc) matrix(0, n, p))
    batch_V <- lapply(seq_len(q), function(cc) matrix(0, n, p))
    batch_fold <- if (want_fold_sums) integer(n) else NULL
    batch_Xt <- if (want_identity) matrix(0, N, p) else NULL
    batch_AL <- if (want_identity) matrix(0, n * q, p) else NULL
    chol_Lambda <- if (want_identity) chol(Lambda) else NULL
  }
  gamma_list <- setNames(vector("list", n), ci$levels)
  residual_list <- setNames(vector("list", n), ci$levels)
  diag_m <- integer(n)
  diag_converged <- logical(n)
  diag_gradient <- numeric(n)
  diag_optim_convergence <- integer(n)
  diag_function_evaluations <- integer(n)
  diag_gradient_evaluations <- integer(n)
  diag_newton_iterations <- integer(n)
  objective_values <- numeric(n)
  identity_errors <- rep(NA_real_, n)
  nuisance_engine <- Sys.getenv(
    "JPDMEQR_CLUSTER_SOLVER", unset = getOption("jpDMEQR.cluster_solver", "bfgs")
  )
  if (!nuisance_engine %in% c("bfgs", "newton", "compiled")) {
    stop("Unknown nuisance solver: ", nuisance_engine)
  }
  compiled_profiles <- if (identical(nuisance_engine, "compiled")) {
    compiled_cluster_profiles_v2(
      y - as.numeric(X %*% beta), Z, ci$rows, tau, h, Lambda,
      gamma_init = gamma_init, control = nuisance_control
    )
  } else NULL

  # Per-cluster worker: pure function of (y_i, X_i, Z_i, beta, ...) with no
  # shared state, so parallel execution is numerically identical to the
  # sequential loop. Each cluster's nuisance solve is deterministic (BFGS +
  # damped Newton; no global RNG dependence beyond the supplied seed).
  #
  # accumulate=TRUE 时 worker 把 p x p 的 Hessian 直接累加进调用环境，只回传
  # 标量与向量，于是同时驻留的 p x p 矩阵从 n 份降到 0 份。累加顺序与父循环
  # 逐项累加相同，和值与参考实现按位相同。
  cluster_worker <- function(ii, accumulate = FALSE) {
    idx <- ci$rows[[ii]]
    cname <- ci$levels[ii]
    y_i <- y[idx]
    X_i <- X[idx, , drop = FALSE]
    Z_i <- Z[idx, , drop = FALSE]
    m_i <- length(idx)
    g0 <- resolve_gamma_init_v2(gamma_init, cname, ii, q)

    ng <- if (is.null(compiled_profiles)) {
      cluster_gamma_dispatch_v2(
        y_i = y_i, X_i = X_i, Z_i = Z_i, beta = beta,
        tau = tau, h = h, Lambda = Lambda,
        gamma_init = g0, control = nuisance_control
      )
    } else list(
      gamma = as.numeric(compiled_profiles$gamma[ii, ]),
      converged = TRUE, gradient_max = compiled_profiles$gradient_max[ii],
      optim_convergence = NA_integer_,
      function_evaluations = compiled_profiles$evaluations[ii],
      gradient_evaluations = compiled_profiles$evaluations[ii],
      newton_iterations = compiled_profiles$iterations[ii]
    )
    gamma_i <- ng$gamma
    r_i <- if (is.null(compiled_profiles)) {
      y_i - as.numeric(X_i %*% beta) - as.numeric(Z_i %*% gamma_i)
    } else compiled_profiles$residual[idx]
    psi_i <- if (is.null(compiled_profiles)) psi_smooth_v2(r_i, tau, h) else
      compiled_profiles$psi[idx]

    objective_i <- if (is.null(compiled_profiles)) {
      mean(rho_smooth_v2(r_i, tau, h)) +
        0.5 * drop(crossprod(gamma_i, Lambda %*% gamma_i))
    } else compiled_profiles$objective[ii]
    score_i <- as.numeric(-crossprod(X_i, psi_i) / m_i)
    score0_i <- as.numeric(-crossprod(X_i, psi_quantile_v2(r_i, tau)) / m_i)
    fold_i <- ((ii - 1L) %% 4L) + 1L

    H_i <- NULL
    H_identity_i <- NULL
    batch_row <- NULL
    identity_err_i <- NA_real_
    if (need_hessian) {
      w <- if (is.null(compiled_profiles)) phi_smooth_v2(r_i, tau, h) / m_i else
        compiled_profiles$phi[idx] / m_i
      WX <- X_i * w
      WZ <- Z_i * w
      B_i <- crossprod(Z_i, WZ) + Lambda
      A_i <- solve_linear_pd_v2(B_i, crossprod(Z_i, WX),
                                name = paste0("B_i for cluster ", cname))
      H_xz <- crossprod(X_i, WZ)
      if (batch_assembly) {
        X_tilde <- if (want_identity) X_i - Z_i %*% A_i else NULL
        AL_i <- if (want_identity) chol_Lambda %*% A_i else NULL
        batch_row <- list(w = w, U = H_xz, V = A_i, Xt = X_tilde, AL = AL_i)
        if (want_cluster_identity) {
          # 两式比较的逐簇诊断按逐簇表达式计算，与装配模式无关
          H_i <- crossprod(X_i, WX) - H_xz %*% A_i
          H_identity_i <- crossprod(X_tilde, X_tilde * w) +
            t(A_i) %*% Lambda %*% A_i
          identity_err_i <- max(abs(H_i - H_identity_i))
        }
      } else {
        H_i <- crossprod(X_i, WX) - H_xz %*% A_i

        if (want_identity) {
          X_tilde <- X_i - Z_i %*% A_i
          H_identity_i <- crossprod(X_tilde, X_tilde * w) +
            t(A_i) %*% Lambda %*% A_i
          identity_err_i <- max(abs(H_i - H_identity_i))
        }
      }
    }
    accumulated <- FALSE
    if (accumulate && !batch_assembly) {
      H_sum <<- H_sum + H_i
      if (want_fold_sums) {
        hessian_fold_sums[[fold_i]] <<- hessian_fold_sums[[fold_i]] + H_i
        hessian_fold_counts[fold_i] <<- hessian_fold_counts[fold_i] + 1L
      }
      if (want_identity) H_identity_sum <<- H_identity_sum + H_identity_i
      accumulated <- TRUE
    }
    list(
      ii = ii,
      gamma = gamma_i,
      residual = r_i,
      m_i = m_i,
      converged = ng$converged,
      gradient_max = ng$gradient_max,
      optim_convergence = ng$optim_convergence,
      function_evaluations = ng$function_evaluations,
      gradient_evaluations = ng$gradient_evaluations,
      newton_iterations = ng$newton_iterations %||% NA_integer_,
      objective = objective_i,
      score = score_i,
      score0 = score0_i,
      fold = fold_i,
      H = if (accumulated) NULL else H_i,
      H_identity = if (accumulated) NULL else H_identity_i,
      identity_error = identity_err_i,
      batch = batch_row
    )
  }

  consume <- function(item) {
    ii <- item$ii
    gamma_list[[ii]] <<- item$gamma
    residual_list[[ii]] <<- item$residual
    diag_m[ii] <<- item$m_i
    diag_converged[ii] <<- item$converged
    diag_gradient[ii] <<- item$gradient_max
    diag_optim_convergence[ii] <<- item$optim_convergence
    diag_function_evaluations[ii] <<- item$function_evaluations
    diag_gradient_evaluations[ii] <<- item$gradient_evaluations
    diag_newton_iterations[ii] <<- item$newton_iterations
    objective_values[ii] <<- item$objective
    cluster_scores[ii, ] <<- item$score
    unsmoothed_cluster_scores[ii, ] <<- item$score0
    if (need_hessian && !batch_assembly) {
      if (!is.null(item$H)) {
        H_sum <<- H_sum + item$H
        if (want_fold_sums) {
          f <- item$fold
          hessian_fold_sums[[f]] <<- hessian_fold_sums[[f]] + item$H
          hessian_fold_counts[f] <<- hessian_fold_counts[f] + 1L
        }
      }
      if (want_identity && !is.null(item$H_identity)) {
        H_identity_sum <<- H_identity_sum + item$H_identity
      }
    }
    if (need_hessian) {
      if (batch_assembly && !is.null(item$batch)) {
        b <- item$batch
        idx_ii <- as.integer(unlist(ci$rows[[ii]]))
        batch_w[idx_ii] <<- b$w
        for (cc in seq_len(q)) {
          batch_U[[cc]][ii, ] <<- b$U[, cc]
          batch_V[[cc]][ii, ] <<- b$V[cc, ]
        }
        if (want_fold_sums) batch_fold[ii] <<- item$fold
        if (want_identity) {
          batch_Xt[idx_ii, ] <<- b$Xt
          batch_AL[((ii - 1L) * q + seq_len(q)), ] <<- b$AL
        }
      }
      identity_errors[ii] <<- item$identity_error
    }
    invisible(NULL)
  }

  cluster_cores <- getOption("jpDMEQR.cluster_cores", 1L)
  use_parallel <- cluster_cores > 1L && .Platform$OS.type != "windows" && n > 1L
  if (use_parallel) {
    # 分叉子进程无法写回父进程的累加器，并行分支沿用逐项累加，并在用完后立即
    # 释放每个 cluster 的矩阵。
    out <- parallel::mclapply(seq_len(n), cluster_worker, accumulate = FALSE,
                              mc.cores = min(cluster_cores, n),
                              mc.preschedule = TRUE)
    for (k in seq_along(out)) {
      consume(out[[k]])
      out[[k]]$H <- NULL
      out[[k]]$H_identity <- NULL
    }
    rm(out)
  } else {
    # batch 模式的逐簇数据由 consume 统一写入，两种进程路径的写入位置一致
    for (ii in seq_len(n)) {
      consume(cluster_worker(ii, accumulate = need_hessian && !batch_assembly))
    }
  }

  if (batch_assembly) {
    # 整块装配：全样本与每个折各用一次大 crossprod，Schur 修正项按 q 列分别用
    # U'V 一次装配。逐簇只收集 X'WZ 与 A 的行，p x p 矩阵不再逐簇构造。
    correction <- matrix(0, p, p)
    for (cc in seq_len(q)) {
      correction <- correction + profile_dense_crossprod_v2(batch_U[[cc]], batch_V[[cc]])
    }
    H_sum <- profile_dense_crossprod_v2(X * batch_w, X) - correction
    if (want_fold_sums) {
      for (f in seq_len(4L)) {
        sel <- which(batch_fold == f)
        if (!length(sel)) next
        idx_f <- as.integer(unlist(ci$rows[sel]))
        Xf <- X[idx_f, , drop = FALSE]
        fold_correction <- matrix(0, p, p)
        for (cc in seq_len(q)) {
          fold_correction <- fold_correction +
            profile_dense_crossprod_v2(batch_U[[cc]][sel, , drop = FALSE],
                                      batch_V[[cc]][sel, , drop = FALSE])
        }
        hessian_fold_sums[[f]] <- profile_dense_crossprod_v2(Xf * batch_w[idx_f], Xf) -
          fold_correction
        hessian_fold_counts[f] <- length(sel)
      }
    }
    if (want_identity) {
      H_identity_sum <- profile_dense_crossprod_v2(batch_Xt * batch_w, batch_Xt) +
        profile_dense_crossprod_v2(batch_AL)
    }
  }

  nuisance_df <- data.frame(
    cluster = ci$levels,
    m_i = diag_m,
    converged = diag_converged,
    gradient_max = diag_gradient,
    optim_convergence = diag_optim_convergence,
    function_evaluations = diag_function_evaluations,
    gradient_evaluations = diag_gradient_evaluations,
    newton_iterations = diag_newton_iterations,
    stringsAsFactors = FALSE
  )
  gamma_mat <- do.call(rbind, gamma_list)
  rownames(gamma_mat) <- ci$levels
  colnames(gamma_mat) <- colnames(Z) %||% paste0("z", seq_len(q))

  score <- colMeans(cluster_scores)
  score0 <- colMeans(unsmoothed_cluster_scores)
  hessian <- if (need_hessian) (H_sum / n + t(H_sum / n)) / 2 else NULL
  hessian_identity <- if (want_identity) (H_identity_sum / n + t(H_identity_sum / n)) / 2 else NULL

  list(
    objective = mean(objective_values),
    score = as.numeric(score),
    score0 = as.numeric(score0),
    hessian = hessian,
    hessian_identity = hessian_identity,
    profile_identity_error = if (want_identity) max(abs(hessian - hessian_identity)) else NA_real_,
    identity_checked = want_identity,
    cluster_identity_checked = want_cluster_identity,
    matrix_backend = getOption("jpDMEQR.matrix_backend",
                                Sys.getenv("JPDMEQR_MATRIX_BACKEND", unset = "R_BLAS")),
    cluster_identity_errors = identity_errors,
    cluster_scores = cluster_scores,
    unsmoothed_cluster_scores = unsmoothed_cluster_scores,
    hessian_fold_sums = hessian_fold_sums,
    hessian_fold_counts = hessian_fold_counts,
    gamma = gamma_mat,
    gamma_list = gamma_list,
    residuals = residual_list,
    nuisance_diagnostics = nuisance_df,
    nuisance_converged = all(nuisance_df$converged),
    max_nuisance_gradient = max(nuisance_df$gradient_max),
    nuisance_iterations = sum(nuisance_df$newton_iterations, na.rm = TRUE),
    Lambda = Lambda,
    n_clusters = n,
    total_rows = N,
    p = p,
    q = q,
    tau = tau,
    h = h
  )
}

profile_objective_v2 <- function(y, X, Z, cluster_id, beta, tau, h,
                                 lambda_gamma = 1, gamma_init = NULL,
                                 nuisance_control = list()) {
  profile_components_v2(
    y = y, X = X, Z = Z, cluster_id = cluster_id,
    beta = beta, tau = tau, h = h, lambda_gamma = lambda_gamma,
    gamma_init = gamma_init, need_hessian = FALSE,
    nuisance_control = nuisance_control
  )$objective
}

# -----------------------------------------------------------------------------
# Penalised profile estimator
# -----------------------------------------------------------------------------

soft_threshold_weighted_v2 <- function(x, threshold) {
  sign(x) * pmax(abs(x) - threshold, 0)
}

largest_eigenvalue_power_v2 <- function(A, maxit = 200L, tol = 1e-8) {
  A <- (A + t(A)) / 2
  p <- nrow(A)
  if (p == 1L) return(max(A[1, 1], 0))
  v <- rep(1 / sqrt(p), p)
  value_old <- 0
  for (iter in seq_len(maxit)) {
    Av <- as.numeric(A %*% v)
    nv <- sqrt(sum(Av^2))
    if (!is.finite(nv) || nv == 0) return(0)
    v <- Av / nv
    value <- drop(crossprod(v, A %*% v))
    if (abs(value - value_old) <= tol * (1 + abs(value))) break
    value_old <- value
  }
  max(as.numeric(value), 0)
}

kkt_residual_profile_v2 <- function(beta, score, lambda_beta, penalty_factor,
                                    active_tol = 1e-10) {
  beta <- as.numeric(beta)
  score <- as.numeric(score)
  pf <- as.numeric(penalty_factor)
  out <- numeric(length(beta))
  unpen <- pf == 0
  out[unpen] <- abs(score[unpen])

  pen <- !unpen
  active <- pen & abs(beta) > active_tol
  inactive <- pen & !active
  out[active] <- abs(score[active] + lambda_beta * pf[active] * sign(beta[active]))
  out[inactive] <- pmax(0, abs(score[inactive]) - lambda_beta * pf[inactive])
  max(out)
}

# -----------------------------------------------------------------------------
# Round-3 cluster self-normalised first-stage penalty (METHOD_SPECIFICATION_
# ROUND3_AMENDMENT.md sections 1-3).
# -----------------------------------------------------------------------------

# ell_j(b) = [ n^{-1} sum_i { s_ij(b) - sbar_j(b) }^2 ]^{1/2}: the standard
# deviation across INDEPENDENT CLUSTERS of the centred smoothed cluster score
# contributions (n denominator, exactly as specified -- not the n-1 sample sd).
score_loadings_v2 <- function(cluster_scores) {
  S <- as.matrix(cluster_scores)
  if (!nrow(S) || ncol(S) < 1L) stop("cluster_scores must be a non-empty matrix.")
  if (any(!is.finite(S))) {
    stop("Non-finite cluster score contributions in score_loadings_v2.")
  }
  if (nrow(S) == 1L) return(rep(0, ncol(S)))
  S <- S - matrix(colMeans(S), nrow(S), ncol(S), byrow = TRUE)
  sqrt(colMeans(S^2))
}

# Freeze: alpha_lambda,n = 0.10 / log(max(n,3));
#        q_lambda,n      = Phi^{-1}(1 - alpha/(2 p_P));
#        lambda_0,n      = 1.10 q_lambda,n / sqrt(n).
lambda_0_n_v2 <- function(n, p_penalized, safety_constant = 1.10) {
  n <- as.integer(n)
  p_penalized <- as.integer(p_penalized)
  if (n < 2L || p_penalized < 1L) stop("lambda_0_n_v2 requires n>=2 and p_P>=1.")
  alpha <- 0.10 / log(max(n, 3))
  q <- stats::qnorm(1 - alpha / (2 * p_penalized))
  lambda_0 <- safety_constant * q / sqrt(n)
  list(
    alpha = alpha,
    normal_quantile = q,
    safety_constant = safety_constant,
    lambda_0 = lambda_0
  )
}

# Scale-normalised weighted KKT residual R_KKT (round-3 Q1(b)):
#   r_j = |g_j + p_j sign(beta_j)|          j in P, |beta_j| > 1e-10
#       = max(0, |g_j| - p_j)              j in P, |beta_j| <= 1e-10
#       = |g_j|                            j in U (p_j == 0)
#   p_ref = median positive p_j;
#   R_KKT = max( max_{j in P} r_j/p_j, max_{j in U} r_j/p_ref ).
# `penalty` is the EFFECTIVE coordinate penalty vector p_j (zero on U).
kkt_residual_normalized_v2 <- function(beta, score, penalty,
                                       active_tol = 1e-10) {
  beta <- as.numeric(beta)
  score <- as.numeric(score)
  p <- as.numeric(penalty)
  if (length(beta) != length(score) || length(score) != length(p)) {
    stop("beta, score and penalty must have equal length.")
  }
  unpen <- p == 0
  pen <- !unpen
  r <- numeric(length(p))
  active <- pen & abs(beta) > active_tol
  inactive <- pen & !active
  r[active] <- abs(score[active] + p[active] * sign(beta[active]))
  r[inactive] <- pmax(0, abs(score[inactive]) - p[inactive])
  r[unpen] <- abs(score[unpen])
  p_pos <- p[pen & p > 0]
  p_ref <- if (length(p_pos)) stats::median(p_pos) else NA_real_
  if (!is.finite(p_ref) || p_ref <= 0) {
    # No positive penalty anywhere: pure unpenalised KKT is the max |score|.
    return(max(r))
  }
  max(c(
    if (any(pen)) max(r[pen] / p[pen]) else 0,
    if (any(unpen)) max(r[unpen]) / p_ref else 0
  ))
}

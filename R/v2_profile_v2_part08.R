# ---------------------------------------------------------------------------
# 对偶形式下的行规划求解：同一个多面体，只有目标随 mu 与坐标变化
#
# 行规划的原始形式是
#   min ||omega||_1  s.t.  ||H omega - e_k||_inf <= mu
# 它的对偶是
#   min e_k' w + mu ||w||_1  s.t.  ||H' w||_inf <= 1
# 拆成 w = wp - wm 之后是一个标准线性规划
#   min (e_k + mu)' wp + (mu - e_k)' wm
#   s.t. H'wp - H'wm <= 1,  -H'wp + H'wm <= 1,  wp, wm >= 0
# 约束矩阵与右端项只由 H 决定，八个候选共用同一个多面体，变化的只有目标系数。
# 原形式的 omega 由约束的对偶价格给出，omega = alpha - beta。
#
# 因此模型只建一次，基在候选之间带下去。实测在注册单元 A02 的第 1 折上，
# mu=0.43/0.86/1.73 的重解需要 0 到 18 次单纯形迭代，mu=0.108 需要约 1000 次。
# 可行性随 mu 单调：某个候选报无界，比它小的候选全部不可行，直接判定而不再
# 求解，判定记录在 attempts 表里。
#
# 折与折相互独立，每个折的模型各自成链，因此折可以分给多个进程；同一折内的
# 候选共享一个模型，跨坐标共享模型则更慢（小 mu 的不可行求解留下的基会让下
# 一个坐标多付出迭代），所以并行单位是折，不是坐标。
#
# 状态码取自 HighsModelStatus：7 最优，8 不可行，9 无界或不可行，10 无界。
# 只有 9 与 10 作为不可行的判定依据；残差超过 mu + 1e-6 的候选记为失败但不
# 据此推断更小的候选，因为那属于数值容差而不是几何结论。
# ---------------------------------------------------------------------------

# 对偶形式的约束矩阵（2p x 2p，对称）
dantzig_dual_matrix_build_v2 <- function(H) {
  p <- nrow(H)
  Ht <- t(H)
  rbind(cbind(Ht, -Ht), cbind(-Ht, Ht))
}

# 对偶形式在给定 mu 与坐标下的目标系数
dantzig_dual_cost_v2 <- function(p, coordinate, mu) {
  e <- numeric(p)
  e[coordinate] <- 1
  c(e + mu, mu - e)
}

# 从对偶解恢复 omega：omega = alpha - beta
dantzig_dual_omega_v2 <- function(row_dual, p) {
  as.numeric(row_dual[seq_len(p)]) - as.numeric(row_dual[(p + 1):(2 * p)])
}

# 一次调用解完整个 mu 网格，基在候选之间复用。返回网格上每个候选的解，
# 不可行的候选返回 NA。降序扫描配单调判定，避免在无界状态之后继续求解。
dantzig_dual_solve_grid_v2 <- function(H, coordinate, mu_grid, solver_opts = list()) {
  if (!requireNamespace("highs", quietly = TRUE)) {
    stop("Package 'highs' is required for the dual warm-start path.")
  }

  H <- assert_numeric_matrix_v2(H, "H")
  H <- (H + t(H)) / 2
  p <- nrow(H)
  if (ncol(H) != p) stop("H must be square.")

  if (is.character(coordinate)) {
    if (is.null(colnames(H)) || !coordinate %in% colnames(H)) {
      stop("Unknown coordinate name: ", coordinate)
    }
    k <- match(coordinate, colnames(H))
  } else {
    k <- as.integer(coordinate)
    if (length(k) != 1L || k < 1L || k > p) stop("Invalid coordinate index.")
  }

  mu_grid <- sort(unique(as.numeric(mu_grid)))
  if (length(mu_grid) == 0L || any(!is.finite(mu_grid)) || any(mu_grid <= 0)) {
    stop("mu_grid must contain positive finite values.")
  }

  e <- numeric(p)
  e[k] <- 1
  n_var <- 2L * p
  # solver_opts 里可以给 solver / presolve / simplex_strategy 等 HiGHS 选项，
  # 也可以直接给一个现成的 control 对象。实测 presolve 与 simplex_strategy
  # 在注册单元的折矩阵上都不改变结果，默认组合最快。
  control <- solver_opts$control
  if (is.null(control)) {
    control <- do.call(highs::highs_control,
                       c(list(threads = 1L, solver = "simplex", presolve = "off",
                              log_to_console = FALSE), solver_opts))
  }
  model <- highs::highs_model(L = dantzig_dual_cost_v2(p, k, max(mu_grid)),
                              lower = rep(0, n_var), upper = rep(Inf, n_var),
                              A = dantzig_dual_matrix_build_v2(H),
                              lhs = rep(-Inf, 2L * p), rhs = rep(1, 2L * p),
                              maximum = FALSE)
  solver <- highs::highs_solver(model, control = control)

  n_mu <- length(mu_grid)
  out <- list(mu = mu_grid,
              omega = vector("list", n_mu),
              residual = rep(NA_real_, n_mu),
              l1_norm = rep(NA_real_, n_mu),
              status = rep(NA_character_, n_mu),
              feasible = rep(FALSE, n_mu),
              unbounded_at = NA_real_)

  for (ii in order(mu_grid, decreasing = TRUE)) {
    mu <- mu_grid[ii]
    if (!is.na(out$unbounded_at)) {
      out$status[ii] <- "infeasible_by_monotonicity"
      next
    }
    solver$L(seq_len(n_var), dantzig_dual_cost_v2(p, k, mu))
    solver$solve(log_to_console = FALSE)
    status <- highs_model_status_label_v2(solver$status())
    out$status[ii] <- status
    if (status %in% c("unbounded", "unbounded_or_infeasible", "infeasible")) {
      out$unbounded_at <- mu
      next
    }
    if (!identical(status, "optimal")) next
    sol <- solver$solution()
    if (is.null(sol$row_dual)) next
    omega <- dantzig_dual_omega_v2(sol$row_dual, p)
    if (length(omega) != p || any(!is.finite(omega))) next
    residual <- max(abs(as.numeric(H %*% omega) - e))
    out$omega[[ii]] <- omega
    out$residual[ii] <- residual
    out$l1_norm[ii] <- sum(abs(omega))
    if (residual <= mu + 1e-6) out$feasible[ii] <- TRUE
  }
  out$solver <- "HiGHS-dual-simplex"
  out$coordinate <- k
  out
}

# 冻结语义下的行规划结果：网格上最小的可行 mu，与 solve_dantzig_row_v2 的结构一致
solve_dantzig_row_dual_warm_v2 <- function(H, coordinate, mu_grid, solver_opts = list()) {
  grid <- dantzig_dual_solve_grid_v2(H, coordinate, mu_grid, solver_opts)
  H <- assert_numeric_matrix_v2(H, "H")
  p <- nrow(H)
  attempts <- data.frame(
    mu = grid$mu,
    status = grid$status,
    residual = grid$residual,
    message = ifelse(is.na(grid$unbounded_at), "",
                     sprintf("unbounded at mu=%.10g", grid$unbounded_at)),
    stringsAsFactors = FALSE
  )
  hit <- which(grid$feasible)[1]
  if (is.na(hit)) {
    return(list(
      omega = rep(NA_real_, p),
      coordinate = grid$coordinate,
      coordinate_name = colnames(H)[grid$coordinate] %||% as.character(grid$coordinate),
      mu = NA_real_,
      residual = NA_real_,
      solver = grid$solver,
      solver_status = "infeasible_on_grid",
      l1_norm = NA_real_,
      l2_norm = NA_real_,
      attempts = attempts,
      feasible = FALSE
    ))
  }
  omega <- grid$omega[[hit]]
  list(
    omega = omega,
    coordinate = grid$coordinate,
    coordinate_name = colnames(H)[grid$coordinate] %||% as.character(grid$coordinate),
    mu = grid$mu[hit],
    residual = grid$residual[hit],
    solver = grid$solver,
    solver_status = grid$status[hit],
    l1_norm = sum(abs(omega)),
    l2_norm = sqrt(sum(omega^2)),
    attempts = attempts,
    feasible = TRUE
  )
}

# 折与折相互独立：每个折有自己的训练 Hessian，因此各建一个模型、各走一遍网格，
# 结果与逐个折串行调用逐位相同。cv_cores 大于 1 时按折分给多个进程。
dantzig_dual_solve_grids_v2 <- function(folds_train, coordinate, mu_grid,
                                        solver_opts = list(), cv_cores = 1L) {
  one <- function(Hf) {
    if (is.null(Hf)) return(NULL)
    dantzig_dual_solve_grid_v2(Hf, coordinate, mu_grid, solver_opts)
  }
  n_used <- sum(!vapply(folds_train, is.null, logical(1)))
  cores <- max(1L, as.integer(cv_cores))
  if (cores <= 1L || n_used <= 1L) {
    return(lapply(folds_train, one))
  }
  if (.Platform$OS.type != "windows") {
    return(parallel::mclapply(folds_train, one, mc.cores = min(cores, n_used)))
  }
  cl <- parallel::makeCluster(min(cores, n_used), type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  # 求解链按函数对象导出到工作进程的 globalenv，工作函数定义在 .GlobalEnv，
  # 于是每个进程只序列化一次模型构建所需的函数与少量参数，折矩阵由
  # parLapply 按块分发。所有包调用都写了命名空间，工作进程不需要 attach。
  fn_env <- new.env(parent = .GlobalEnv)
  for (nm in c("dantzig_dual_solve_grid_v2", "dantzig_dual_matrix_build_v2",
               "dantzig_dual_cost_v2", "dantzig_dual_omega_v2",
               "highs_model_status_label_v2", "assert_numeric_matrix_v2",
               "%||%")) {
    assign(nm, get(nm, envir = environment(), inherits = TRUE), envir = fn_env)
  }
  fn_env$dantzig_dual_coordinate_v2 <- coordinate
  fn_env$dantzig_dual_mu_grid_v2 <- mu_grid
  fn_env$dantzig_dual_solver_opts_v2 <- solver_opts
  parallel::clusterExport(cl, ls(fn_env), envir = fn_env)
  eval(quote(dantzig_dual_worker_v2 <- function(Hf) {
    if (is.null(Hf)) return(NULL)
    dantzig_dual_solve_grid_v2(
      Hf,
      get0("dantzig_dual_coordinate_v2", envir = .GlobalEnv, inherits = FALSE),
      get0("dantzig_dual_mu_grid_v2", envir = .GlobalEnv, inherits = FALSE),
      get0("dantzig_dual_solver_opts_v2", envir = .GlobalEnv, inherits = FALSE)
    )
  }), envir = .GlobalEnv)
  on.exit(rm("dantzig_dual_worker_v2", envir = .GlobalEnv), add = TRUE)
  parallel::parLapply(cl, folds_train, dantzig_dual_worker_v2)
}

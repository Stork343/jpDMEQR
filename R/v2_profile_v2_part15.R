unpenalised_optimizer_v2 <- function() {
  value <- Sys.getenv(
    "JPDMEQR_UNPENALISED_OPTIMIZER",
    unset = getOption("jpDMEQR.unpenalised_optimizer", "proximal_gradient")
  )
  if (length(value) != 1L || is.na(value) ||
      !value %in% c("proximal_gradient", "newton_armijo")) {
    stop("JPDMEQR_UNPENALISED_OPTIMIZER must be proximal_gradient or newton_armijo.")
  }
  value
}

proximal_scaling_v2 <- function() {
  value <- Sys.getenv(
    "JPDMEQR_PROXIMAL_SCALING",
    unset = getOption("jpDMEQR.proximal_scaling", "none")
  )
  if (length(value) != 1L || is.na(value) || !value %in% c("none", "cluster_rms")) {
    stop("JPDMEQR_PROXIMAL_SCALING must be none or cluster_rms.")
  }
  value
}

fit_profile_unpenalised_newton_v2 <- function(
    y, X, Z, cluster_id, tau, h, lambda_beta = 0, lambda_gamma = 1,
    penalty_factor = rep(0, ncol(X)), beta_init = rep(0, ncol(X)),
    control = list()) {
  validate_tau_h_v2(tau, h)
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  y <- as.numeric(y)
  p <- ncol(X)
  clusters <- cluster_index_v2(cluster_id)
  beta <- as.numeric(beta_init)
  if (!length(y) || p < 1L || ncol(Z) < 1L ||
      nrow(X) != length(y) || nrow(Z) != length(y) ||
      length(cluster_id) != length(y) || any(!is.finite(y)) ||
      length(beta) != p || any(!is.finite(beta)) ||
      length(penalty_factor) != p || any(!is.finite(penalty_factor)) ||
      any(penalty_factor < 0) || length(lambda_beta) != 1L ||
      !is.finite(lambda_beta) || lambda_beta < 0 ||
      any(lambda_beta * penalty_factor != 0)) {
    stop("The Newton profile solver requires finite unpenalised data and coefficients.")
  }
  max_iter <- control$max_iter %||% 2000L
  max_backtrack <- control$max_backtrack %||% 50L
  beta_tol <- control$beta_tol %||% 1e-7
  kkt_tol <- control$kkt_normalized_tol %||% 1e-3
  backtrack_factor <- control$backtrack_factor %||% 0.5
  min_step <- control$min_step %||% 1e-12
  nuisance_control <- control$nuisance_control %||% list()
  requested_nuisance_grad_tol <- nuisance_control$grad_tol %||% 1e-8
  nuisance_control$grad_tol <- min(requested_nuisance_grad_tol, 1e-11)
  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1L ||
      max_iter != as.integer(max_iter) || length(max_backtrack) != 1L ||
      !is.finite(max_backtrack) || max_backtrack < 0L ||
      max_backtrack != as.integer(max_backtrack) ||
      length(beta_tol) != 1L || !is.finite(beta_tol) || beta_tol <= 0 ||
      length(kkt_tol) != 1L || !is.finite(kkt_tol) || kkt_tol <= 0 ||
      length(backtrack_factor) != 1L || !is.finite(backtrack_factor) ||
      backtrack_factor <= 0 || backtrack_factor >= 1 ||
      length(min_step) != 1L || !is.finite(min_step) || min_step <= 0) {
    stop("Invalid Newton profile solver control.")
  }
  # 固定坐标尺度只用于线性方程；目标函数和接受检验保持原始单位。
  rows <- split(seq_along(y), as.character(cluster_id))
  second_moment <- Reduce(`+`, lapply(rows, function(index) {
    colMeans(X[index, , drop = FALSE]^2)
  })) / clusters$n
  coordinate_scale <- sqrt(second_moment)
  if (any(!is.finite(coordinate_scale)) || any(coordinate_scale <= 0)) {
    stop("Every unpenalised coordinate must have positive finite measurement scale.")
  }
  evaluations <- 0L
  evaluate <- function(value, gamma, hessian = TRUE, identity = FALSE, folds = FALSE) {
    evaluations <<- evaluations + 1L
    result <- profile_components_v2(
      y, X, Z, cluster_id, value, tau, h, lambda_gamma,
      gamma_init = gamma, need_hessian = hessian, identity_check = identity,
      need_fold_sums = folds, nuisance_control = nuisance_control
    )
    if (!isTRUE(result$nuisance_converged) ||
        !is.finite(result$max_nuisance_gradient) ||
        result$max_nuisance_gradient > (nuisance_control$grad_tol %||% 1e-8) ||
        any(!is.finite(result$score)) || !is.finite(result$objective)) {
      stop("Newton profile nuisance or objective acceptance failed.")
    }
    result
  }
  current <- evaluate(beta, NULL)
  history <- vector("list", max_iter)
  converged <- FALSE
  failure_stage <- ""
  for (iter in seq_len(max_iter)) {
    H <- (current$hessian + t(current$hessian)) / 2
    scaled_H <- H / outer(coordinate_scale, coordinate_scale)
    eig <- eigen(scaled_H, symmetric = TRUE)
    curvature_scale <- max(1, max(abs(eig$values)))
    if (min(eig$values) < -1e-10 * curvature_scale) {
      stop("The convex profile Hessian has negative numerical curvature.")
    }
    # 有限带宽的零曲率区间采用明确的正定 Newton 方程，Armijo 检查完整目标。
    curvature_floor <- 1e-8 * curvature_scale
    denominator <- pmax(eig$values, curvature_floor)
    scaled_score <- current$score / coordinate_scale
    direction <- -as.numeric(eig$vectors %*%
                               (crossprod(eig$vectors, scaled_score) / denominator)) /
      coordinate_scale
    # 限制一次 Newton 更新的预测量变化，有限支撑区间里的零曲率仍可开始迭代。
    prediction_change <- max(abs(as.numeric(X %*% direction)))
    prediction_radius <- max(1, max(abs(y - as.numeric(X %*% beta))))
    direction_multiplier <- min(1, prediction_radius / max(prediction_change, .Machine$double.eps))
    direction <- direction_multiplier * direction
    derivative <- sum(current$score * direction)
    kkt_current <- kkt_residual_normalized_v2(beta, current$score, rep(0, p))
    if (kkt_current <= kkt_tol &&
        max(abs(direction)) <= beta_tol * max(1, max(abs(beta)))) {
      beta_new <- beta
      candidate <- current
      step <- 0
      bt <- 0L
      accepted <- TRUE
    } else {
      if (!is.finite(derivative) || derivative >= 0 || any(!is.finite(direction))) {
        stop("The Newton profile direction must be a finite descent direction.")
      }
      step <- 1
      accepted <- FALSE
      allowance <- 64 * .Machine$double.eps * max(1, abs(current$objective))
      for (bt in 0:max_backtrack) {
        beta_new <- beta + step * direction
        candidate <- evaluate(beta_new, current$gamma, hessian = FALSE)
        if (candidate$objective <= current$objective +
            1e-4 * step * derivative + allowance) {
          accepted <- TRUE
          break
        }
        step <- step * backtrack_factor
        if (step < min_step) break
      }
    }
    if (!accepted) {
      failure_stage <- "profile_backtracking"
      break
    }
    beta_change <- max(abs(beta_new - beta))
    kkt <- kkt_residual_profile_v2(beta_new, candidate$score, lambda_beta, penalty_factor)
    kkt_norm <- kkt_residual_normalized_v2(beta_new, candidate$score, rep(0, p))
    history[[iter]] <- data.frame(
      iteration = iter, profile_objective = candidate$objective,
      composite_objective = candidate$objective, beta_change = beta_change,
      kkt_residual = kkt, kkt_normalized = kkt_norm, step = step,
      backtracks = bt, prox_refresh = TRUE,
      nuisance_converged = candidate$nuisance_converged,
      max_nuisance_gradient = candidate$max_nuisance_gradient,
      profile_identity_error = candidate$profile_identity_error,
      scaled_hessian_min_eigenvalue = min(eig$values),
      newton_curvature_floor = curvature_floor,
      newton_direction_max = max(abs(direction)),
      newton_direction_multiplier = direction_multiplier,
      newton_prediction_radius = prediction_radius,
      stringsAsFactors = FALSE
    )
    beta <- beta_new
    current <- candidate
    if (beta_change <= beta_tol * max(1, max(abs(beta))) && kkt_norm <= kkt_tol) {
      converged <- TRUE
      break
    }
    current <- evaluate(beta, current$gamma)
  }
  current <- evaluate(beta, current$gamma, identity = TRUE, folds = TRUE)
  history <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  if (is.null(history)) history <- data.frame()
  final_kkt <- kkt_residual_profile_v2(beta, current$score, lambda_beta, penalty_factor)
  final_kkt_norm <- kkt_residual_normalized_v2(beta, current$score, rep(0, p))
  converged <- converged && final_kkt_norm <= kkt_tol
  list(
    beta = setNames(beta, colnames(X) %||% paste0("x", seq_len(p))),
    gamma = current$gamma, components = current,
    profile_objective = current$objective, composite_objective = current$objective,
    lambda_beta = lambda_beta, lambda_gamma = lambda_gamma,
    penalty_factor = penalty_factor, converged = converged,
    failure_stage = failure_stage, kkt_residual = final_kkt,
    kkt_normalized = final_kkt_norm, kkt_normalized_tol = kkt_tol,
    beta_tol = beta_tol,
    first_stage_beta_change = if (nrow(history)) tail(history$beta_change, 1) else NA_real_,
    first_stage_nonzero_count = sum(abs(beta) > 1e-10),
    iterations = nrow(history), history = history,
    optimizer = "newton_armijo", coordinate_scale = setNames(coordinate_scale, colnames(X)),
    component_evaluations = evaluations,
    requested_nuisance_grad_tol = requested_nuisance_grad_tol,
    numerical_nuisance_grad_tol = nuisance_control$grad_tol, call = match.call()
  )
}

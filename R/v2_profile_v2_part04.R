cluster_sandwich_coordinate_v2 <- function(cluster_scores, omega, df_correction = FALSE) {
  G <- assert_numeric_matrix_v2(cluster_scores, "cluster_scores")
  omega <- as.numeric(omega)
  if (ncol(G) != length(omega)) stop("omega dimension mismatch.")
  n <- nrow(G)
  projected <- as.numeric(G %*% omega)
  projected_centered <- projected - mean(projected)
  sigma2 <- mean(projected_centered^2)
  if (isTRUE(df_correction) && n > 1) sigma2 <- sigma2 * n / (n - 1)
  list(
    sigma2 = sigma2,
    se = sqrt(max(sigma2, 0) / n),
    projected_scores = projected,
    projected_centered = projected_centered,
    n_clusters = n
  )
}

# -----------------------------------------------------------------------------
# Non-truth-based inverse-Hessian CV selection of the Dantzig tolerance
# (ROUND-2 adjudication: docs/METHOD_SPECIFICATION_ROUND2_AMENDMENT.md).
#
# The round-1 rule (one-SE band, largest-in-band mu) is REVOKED: for Dantzig
# inverse rows a larger mu admits a coarser row, so once the held-out loss is
# flat, largest-in-band systematically biases toward the worst direction.
#
# New rule:
#   - n < 200: deterministic 2-fold cluster CV;
#   - n >= 200: deterministic 4-fold cluster CV;
#   - held-out loss  L_defect_f(mu) = || H_val omega_train(-f)(mu) - e_k ||_inf;
#   - choose argmin over mu of the fold-mean defect;
#   - numerical ties choose the SMALLER mu.
# The inverse-quadratic loss L_quad = 0.5 omega'H_val omega - e'omega is kept
# as a diagnostic only. No simulation truth, coverage, or bias is used.
select_mu_inverse_hessian_cv_v2 <- function(H,
                                            coordinate,
                                            mu_grid,
                                            fold_sums,
                                            fold_counts,
                                            solver_preference = c("CLARABEL", "ECOS", "SCS"),
                                            solver_opts = list(),
                                            n_folds = NULL) {
  mode <- Sys.getenv("JPDMEQR_DANTZIG_SOLVER",
                     unset = getOption("jpDMEQR.dantzig_solver", "cvxr"))
  p <- nrow(H)
  H <- (H + t(H)) / 2
  n <- sum(fold_counts)
  if (is.null(n_folds)) n_folds <- if (n < 200L) 2L else 4L
  n_folds <- max(2L, as.integer(n_folds))
  e <- numeric(p); e[coordinate] <- 1
  mu_grid <- sort(unique(as.numeric(mu_grid)))
  mu_grid <- mu_grid[mu_grid > 0]

  # Build fold indicators deterministically over the cluster indices stored in
  # fold_counts has length = n_folds; if the counts came from a 4-fold build we
  # re-bucket to n_folds deterministically.
  # The fold sums are indexed 1..L (length(fold_sums)); if n_folds differs from
  # length(fold_sums), merge adjacent fold buckets deterministically.
  L <- length(fold_sums)
  if (L < n_folds) stop("fold_sums length is smaller than the requested n_folds")
  # profile_components stores 4 fold buckets (clusters by (ii-1)%%4).
  # For n_folds=2, merging buckets (1,3)->1 and (2,4)->2 reproduces the
  # deterministic ((ii-1)%%2) cluster partition.
  bucket_to_fold <- rep(seq_len(n_folds), length.out = L)
  H_folds <- lapply(seq_len(n_folds), function(f) {
    idx <- which(bucket_to_fold == f)
    s <- Reduce(`+`, fold_sums[idx])
    cnt <- sum(fold_counts[idx])
    if (cnt > 0) (s + t(s)) / 2 / cnt else NULL
  })
  fold_n <- vapply(seq_len(n_folds), function(f) sum(fold_counts[bucket_to_fold == f]), numeric(1))

  results <- if (identical(mode, "osqp")) {
    # Lever 1: one OSQP model per fold, warm-start chained across the mu grid
    # (the constraint matrix is mu-invariant; each subsequent solve is O((n+m)^2)).
    chain_output <- vector("list", length(mu_grid))
    for (f in seq_len(n_folds)) {
      n_f <- fold_n[f]
      if (n_f <= 0L || n_f >= n) next
      H_val_f <- H_folds[[f]]
      if (is.null(H_val_f)) next
      H_train_f <- (n * H - n_f * H_val_f) / (n - n_f)
      H_train_f <- (H_train_f + t(H_train_f)) / 2
      parts <- dantzig_build_osqp_model_v2(H_train_f, p)
      sols <- dantzig_osqp_solve_chain_v2(parts, e, mu_grid)
      for (ii in seq_along(sols)) {
        if (is.null(chain_output[[ii]])) {
          chain_output[[ii]] <- list(mu = mu_grid[ii], defects = numeric(n_folds),
                                     quads = numeric(n_folds), infeasible_folds = 0L)
        }
        s <- sols[[ii]]
        if (isTRUE(s$feasible) && !is.null(s$omega)) {
          omega <- s$omega
          r <- as.numeric(H_val_f %*% omega) - e
          chain_output[[ii]]$defects[f] <- max(abs(r))
          chain_output[[ii]]$quads[f] <-
            0.5 * drop(crossprod(omega, H_val_f %*% omega)) - drop(crossprod(e, omega))
        } else {
          chain_output[[ii]]$infeasible_folds <- chain_output[[ii]]$infeasible_folds + 1L
        }
      }
    }
    lapply(chain_output, function(co) {
      nf <- n_folds
      if (co$infeasible_folds > 0L && nf - co$infeasible_folds < ceiling(nf / 2)) return(NULL)
      list(
        mu = co$mu,
        defect_mean = mean(co$defects[is.finite(co$defects)]),
        defect_se = stats::sd(co$defects[is.finite(co$defects)]) /
          sqrt(sum(is.finite(co$defects))),
        quad_mean = mean(co$quads[is.finite(co$quads)]),
        infeasible_folds = co$infeasible_folds
      )
    })
  } else {
    # CLARABEL reference path (frozen solver). Each (fold, mu) solve is
    # independent; when jpDMEQR.cv_cores > 1 the fold x mu grid is solved in
    # parallel (Linux: fork-based mclapply - keep `jobs` at 1 and use shards to
    # avoid nested-fork memory blowup; Windows: deterministic PSOCK) and
    # aggregated in the original fold order, so the selected mu and defect
    # statistics are bit-identical to the serial reference.
    # Bit-identity guarantee: every (fold, mu) CLARABEL solve is a pure
    # function of its inputs; the defects/quads arrays are refilled fold by
    # fold in the same order, so mean/sd and the argmin-tie rule see exactly
    # the same numbers as the serial path.
    cv_cores <- max(1L, as.integer(getOption("jpDMEQR.cv_cores", 1L)))
    if (cv_cores <= 1L || n_folds * length(mu_grid) <= 1L) {
      results <- lapply(mu_grid, function(mu) {
        defects <- numeric(n_folds)
        quads <- numeric(n_folds)
        infeasible_folds <- 0L
        for (f in seq_len(n_folds)) {
          n_f <- fold_n[f]
          if (n_f <= 0L || n_f >= n) next
          H_val_f <- H_folds[[f]]
          if (is.null(H_val_f)) next
          # H_train = (n H - n_f H_val)/(n - n_f)
          H_train_f <- (n * H - n_f * H_val_f) / (n - n_f)
          H_train_f <- (H_train_f + t(H_train_f)) / 2
          sol <- tryCatch(
            solve_dantzig_row_v2(
              H = H_train_f, coordinate = coordinate, mu_grid = c(mu),
              solver_preference = solver_preference, solver_opts = solver_opts
            ),
            error = function(e) list(feasible = FALSE, omega = NULL)
          )
          if (!isTRUE(sol$feasible) || is.null(sol$omega)) { infeasible_folds <- infeasible_folds + 1L; next }
          omega <- as.numeric(sol$omega)
          r <- as.numeric(H_val_f %*% omega) - e
          defects[f] <- max(abs(r))
          quads[f] <- 0.5 * drop(crossprod(omega, H_val_f %*% omega)) - drop(crossprod(e, omega))
        }
        if (infeasible_folds > 0L && n_folds - infeasible_folds < ceiling(n_folds / 2)) {
          return(NULL)
        }
        list(
          mu = mu,
          defect_mean = mean(defects[is.finite(defects)]),
          defect_se = stats::sd(defects[is.finite(defects)]) / sqrt(sum(is.finite(defects))),
          quad_mean = mean(quads[is.finite(quads)]),
          infeasible_folds = infeasible_folds
        )
      })
    } else {
      results <- dantzig_cv_clarabel_parallel_v2(
        H = H, n = n, n_folds = n_folds, fold_n = fold_n, H_folds = H_folds,
        mu_grid = mu_grid, coordinate = coordinate, e = e,
        solver_preference = solver_preference, solver_opts = solver_opts,
        cv_cores = cv_cores
      )
    }
  }
  ok <- results[!vapply(results, is.null, logical(1))]
  if (!length(ok)) {
    return(list(mu = NA_real_, defect_mean = NA_real_, defect_se = NA_real_,
                quad_mean = NA_real_, selected = FALSE, candidates = data.frame()))
  }
  defect_means <- vapply(ok, `[[`, numeric(1), "defect_mean")
  quad_means <- vapply(ok, `[[`, numeric(1), "quad_mean")
  mus <- vapply(ok, `[[`, numeric(1), "mu")
  # argmin of the held-out defect; ties choose the SMALLER mu (round-2).
  loss_min <- min(defect_means, na.rm = TRUE)
  cand <- which(abs(defect_means - loss_min) <= 1e-12 * pmax(1, abs(loss_min)))
  selected_idx <- cand[which.min(mus[cand])]
  list(
    mu = mus[selected_idx],
    defect_mean = defect_means[selected_idx],
    quad_mean = quad_means[selected_idx],
    min_defect = loss_min,
    selected = TRUE,
    candidates = data.frame(
      mu = mus, defect_mean = defect_means, quad_mean = quad_means,
      stringsAsFactors = FALSE
    )
  )
}

# ---------------------------------------------------------------------------
# Parallel CLARABEL mu-CV (the ~75-80% cost term).
#
# Bit-identity contract: every (fold, mu) CLARABEL solve is a pure function of
# its inputs, so executing the fold x mu grid in any order and then refilling
# the defects/quads vectors fold by fold reproduces the serial reference path
# bit-for-bit (same mean/sd inputs, same argmin-tie rule, same
# infeasible-fold drop rule). cv_cores == 1 never enters this path.
#
# Worker environment: the solver chain (solve_dantzig_row_v2 and its two
# helpers %||% and assert_numeric_matrix_v2) is exported to PSOCK workers as
# function objects - this is the smallest self-contained unit that actually
# runs the frozen CLARABEL program; no module re-source is needed on workers.
# Linux uses fork-based mclapply (no serialisation, inherits everything), but
# keep `jobs` at 1 and use shards for task-level parallelism to avoid
# nested-fork memory blowup (see launch_sharded_pilot.R).
# ---------------------------------------------------------------------------

# Top-level worker for one (mu_idx, fold) cell: a pure function of its
# inputs, so the exact same call runs in serial, forked, and PSOCK contexts
# and the result is independent of execution order.
dantzig_cv_solve_task_v2 <- function(j, ctx, tasks) {
  mu_idx <- tasks$mu_idx[j]
  fold <- tasks$fold[j]
  mu <- ctx$mu_grid[mu_idx]
  n <- ctx$n
  n_f <- ctx$fold_n[fold]
  if (n_f <= 0L || n_f >= n) return(list(mu_idx = mu_idx, fold = fold, skipped = TRUE))
  H_val_f <- ctx$H_folds[[fold]]
  if (is.null(H_val_f)) return(list(mu_idx = mu_idx, fold = fold, skipped = TRUE))
  # H_train = (n H - n_f H_val)/(n - n_f)
  H_train_f <- (n * ctx$H - n_f * H_val_f) / (n - n_f)
  H_train_f <- (H_train_f + t(H_train_f)) / 2
  sol <- tryCatch(
    solve_dantzig_row_v2(
      H = H_train_f, coordinate = ctx$coordinate, mu_grid = c(mu),
      solver_preference = ctx$solver_preference, solver_opts = ctx$solver_opts
    ),
    error = function(e) list(feasible = FALSE, omega = NULL)
  )
  if (!isTRUE(sol$feasible) || is.null(sol$omega)) {
    return(list(mu_idx = mu_idx, fold = fold, infeasible = TRUE))
  }
  omega <- as.numeric(sol$omega)
  r <- as.numeric(H_val_f %*% omega) - ctx$e
  list(
    mu_idx = mu_idx, fold = fold, feasible = TRUE,
    defect = max(abs(r)),
    quad = 0.5 * drop(crossprod(omega, H_val_f %*% omega)) - drop(crossprod(ctx$e, omega))
  )
}

dantzig_cv_clarabel_parallel_v2 <- function(H, n, n_folds, fold_n, H_folds,
                                            mu_grid, coordinate, e,
                                            solver_preference, solver_opts,
                                            cv_cores) {
  tasks <- expand.grid(
    mu_idx = seq_along(mu_grid),
    fold = seq_len(n_folds)
  )
  tasks <- tasks[order(tasks$mu_idx, tasks$fold), , drop = FALSE]
  n_tasks <- nrow(tasks)
  ctx <- list(
    H = H, n = n, fold_n = fold_n, H_folds = H_folds, mu_grid = mu_grid,
    coordinate = coordinate, e = e,
    solver_preference = solver_preference, solver_opts = solver_opts
  )
  out <- if (.Platform$OS.type != "windows") {
    # Forked children inherit the full process image: plain closure over the
    # local ctx/tasks (no serialisation); mclapply returns results in input
    # order regardless of execution order.
    parallel::mclapply(seq_len(n_tasks),
                       function(j) dantzig_cv_solve_task_v2(j, ctx, tasks),
                       mc.cores = min(cv_cores, n_tasks),
                       mc.preschedule = FALSE)
  } else {
    cl <- parallel::makeCluster(min(cv_cores, n_tasks), type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Shared data + the solver chain are exported by value into each worker's
    # globalenv. The worker closure is DEFINED in .GlobalEnv (eval - not
    # assign - so the function's environment is globalenv, which PSOCK
    # serialises as an id reference rather than expanding the enclosing
    # environment with the large H matrices); it resolves the ctx/tasks/solver
    # symbols from the worker globalenv where they were exported. One
    # serialisation per worker, not per task.
    fn_env <- new.env(parent = .GlobalEnv)
    for (nm in c("dantzig_cv_solve_task_v2", "solve_dantzig_row_v2",
                 "assert_numeric_matrix_v2", "%||%")) {
      if (exists(nm, envir = environment(), inherits = TRUE)) {
        assign(nm, get(nm, envir = environment(), inherits = TRUE), envir = fn_env)
      } else if (exists(nm, envir = .GlobalEnv, inherits = TRUE)) {
        assign(nm, get(nm, envir = .GlobalEnv, inherits = TRUE), envir = fn_env)
      }
    }
    fn_env$dantzig_cv_ctx_v2 <- ctx
    fn_env$dantzig_cv_tasks_v2 <- tasks
    parallel::clusterExport(cl, ls(fn_env), envir = fn_env)
    eval(quote(dantzig_cv_worker_v2 <- function(j) {
      dantzig_cv_solve_task_v2(j,
                               get0("dantzig_cv_ctx_v2", envir = .GlobalEnv,
                                    inherits = FALSE),
                               get0("dantzig_cv_tasks_v2", envir = .GlobalEnv,
                                    inherits = FALSE))
    }), envir = .GlobalEnv)
    on.exit(rm("dantzig_cv_worker_v2", envir = .GlobalEnv), add = TRUE)
    # parLapplyLB: dynamic load balancing. CLARABEL wall time varies widely
    # across (fold, mu) cells (observed 2.9-7.4s at p=500 for fixed inputs);
    # static round-robin (parLapply) assigns each fold's whole mu chain to one
    # worker and waits for the slowest worker. LB hands out tasks on demand;
    # results are still returned in input order.
    parallel::parLapplyLB(cl, seq_len(n_tasks), dantzig_cv_worker_v2)
  }
  # Aggregate in the original (mu outer, fold inner) order - bit-identical to
  # the serial reference path.
  lapply(seq_along(mu_grid), function(mi) {
    defects <- numeric(n_folds)
    quads <- numeric(n_folds)
    infeasible_folds <- 0L
    for (f in seq_len(n_folds)) {
      r <- out[[(mi - 1L) * n_folds + f]]
      if (is.null(r) || isTRUE(r$skipped)) next
      if (isTRUE(r$infeasible)) { infeasible_folds <- infeasible_folds + 1L; next }
      defects[f] <- r$defect
      quads[f] <- r$quad
    }
    if (infeasible_folds > 0L && n_folds - infeasible_folds < ceiling(n_folds / 2)) {
      return(NULL)
    }
    list(
      mu = mu_grid[mi],
      defect_mean = mean(defects[is.finite(defects)]),
      defect_se = stats::sd(defects[is.finite(defects)]) / sqrt(sum(is.finite(defects))),
      quad_mean = mean(quads[is.finite(quads)]),
      infeasible_folds = infeasible_folds
    )
  })
}

debias_profile_coordinates_v2 <- function(fit,
                                           coordinates,
                                           mu_grid,
                                           ci_level = 0.95,
                                           solver_preference = c("CLARABEL", "ECOS", "SCS"),
                                           solver_opts = list(),
                                           df_correction = FALSE) {
  if (is.null(fit$components$hessian) || is.null(fit$components$score)) {
    stop("fit does not contain profile score/Hessian components.")
  }
  H <- fit$components$hessian
  colnames(H) <- rownames(H) <- names(fit$beta)
  coords <- coordinates
  if (is.character(coords)) {
    idx <- match(coords, names(fit$beta))
    if (anyNA(idx)) stop("Unknown coordinates: ", paste(coords[is.na(idx)], collapse = ", "))
  } else {
    idx <- as.integer(coords)
  }

  alpha <- 1 - ci_level
  zcrit <- stats::qnorm(1 - alpha / 2)
  rows <- vector("list", length(idx))
  directions <- vector("list", length(idx))

  fold_sums <- fit$components$hessian_fold_sums
  fold_counts <- fit$components$hessian_fold_counts
  use_cv <- !is.null(fold_sums) && length(fold_sums) == 4L &&
    !is.null(fold_counts) && all(fold_counts > 0L)

  for (ii in seq_along(idx)) {
    k <- idx[ii]
    cv <- if (use_cv) {
      select_mu_inverse_hessian_cv_v2(
        H = H, coordinate = k, mu_grid = mu_grid,
        fold_sums = fold_sums, fold_counts = fold_counts,
        solver_preference = solver_preference, solver_opts = solver_opts
      )
    } else {
      list(mu = NA_real_, selected = FALSE,
           candidates = data.frame(mu = NA_real_, loss_mean = NA_real_,
                                   loss_se = NA_real_, in_one_se = NA,
                                   stringsAsFactors = FALSE))
    }

    direction <- if (isTRUE(cv$selected)) {
      solve_dantzig_row_v2(
        H = H, coordinate = k, mu_grid = c(cv$mu),
        solver_preference = solver_preference, solver_opts = solver_opts
      )
    } else {
      # backward-compatible fallback: smallest feasible value on the grid
      solve_dantzig_row_v2(
        H = H, coordinate = k, mu_grid = mu_grid,
        solver_preference = solver_preference, solver_opts = solver_opts
      )
    }
    directions[[ii]] <- direction
    if (!direction$feasible) {
      rows[[ii]] <- data.frame(
        coordinate = names(fit$beta)[k], index = k, feasible = FALSE,
        beta_hat = fit$beta[k], beta_tilde = NA_real_, se = NA_real_,
        se_smoothed = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        mu = NA_real_, mu_cv_min_defect = cv$min_defect,
        mu_cv_defect = cv$defect_mean, mu_cv_quad = cv$quad_mean,
        dantzig_residual = NA_real_,
        omega_l1 = NA_real_, omega_l2 = NA_real_, adjacent_stability = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    omega <- direction$omega
    beta_tilde <- as.numeric(fit$beta[k] - sum(omega * fit$components$score))
    # Primary Wald meat: unsmoothed quantile cluster scores (theory decision
    # section 1.1). Bread remains the smoothed effective Hessian via omega.
    score0 <- fit$components$unsmoothed_cluster_scores
    sand0 <- if (!is.null(score0)) {
      cluster_sandwich_coordinate_v2(score0, omega, df_correction = df_correction)
    } else {
      list(se = NA_real_)
    }
    # Diagnostic/sensitivity meat: smoothed cluster scores (section 1.2).
    sandh <- cluster_sandwich_coordinate_v2(
      fit$components$cluster_scores, omega, df_correction = df_correction
    )
    se0 <- sand0$se
    seh <- sandh$se
    se <- if (is.finite(se0)) se0 else seh
    # adjacent-grid direction stability: relative change between consecutive
    # feasible candidate rows (diagnostic only).
    adj <- NA_real_
    cand <- cv$candidates
    if (isTRUE(cv$selected) && nrow(cand) >= 2L) {
      feas_mus <- cand$mu[is.finite(cand$loss_mean)]
      if (length(feas_mus) >= 2L) {
        sols <- lapply(feas_mus, function(mu) {
          s <- tryCatch(
            solve_dantzig_row_v2(H, k, c(mu), solver_preference, solver_opts),
            error = function(e) list(feasible = FALSE)
          )
          if (isTRUE(s$feasible)) as.numeric(s$omega) else NULL
        })
        sols <- sols[!vapply(sols, is.null, logical(1))]
        if (length(sols) >= 2L) {
          adj <- max(vapply(seq_len(length(sols) - 1L), function(j) {
            d <- sols[[j + 1L]] - sols[[j]]
            sqrt(sum(d^2)) / max(sqrt(sum(sols[[j]]^2)), 1e-12)
          }, numeric(1)))
        }
      }
    }
    rows[[ii]] <- data.frame(
      coordinate = names(fit$beta)[k],
      index = k,
      feasible = TRUE,
      beta_hat = as.numeric(fit$beta[k]),
      beta_tilde = beta_tilde,
      se = se,
      se_smoothed = seh,
      ci_lower = beta_tilde - zcrit * se,
      ci_upper = beta_tilde + zcrit * se,
      mu = direction$mu,
      mu_cv_min_defect = cv$min_defect,
      mu_cv_defect = cv$defect_mean,
      mu_cv_quad = cv$quad_mean,
      dantzig_residual = direction$residual,
      omega_l1 = direction$l1_norm,
      omega_l2 = direction$l2_norm,
      adjacent_stability = adj,
      stringsAsFactors = FALSE
    )
  }

  list(
    table = do.call(rbind, rows),
    directions = directions,
    ci_level = ci_level,
    fit = fit
  )
}

# -----------------------------------------------------------------------------
# Numerical geometry checks
# -----------------------------------------------------------------------------

finite_difference_gradient_v2 <- function(fn, x, eps = 1e-5) {
  x <- as.numeric(x)
  p <- length(x)
  out <- numeric(p)
  for (k in seq_len(p)) {
    step <- numeric(p); step[k] <- eps
    out[k] <- (fn(x + step) - fn(x - step)) / (2 * eps)
  }
  out
}

finite_difference_jacobian_v2 <- function(fn, x, eps = 1e-5) {
  x <- as.numeric(x)
  f0 <- as.numeric(fn(x))
  out <- matrix(NA_real_, nrow = length(f0), ncol = length(x))
  for (k in seq_along(x)) {
    step <- numeric(length(x)); step[k] <- eps
    out[, k] <- (as.numeric(fn(x + step)) - as.numeric(fn(x - step))) / (2 * eps)
  }
  out
}

validate_profile_geometry_v2 <- function(y,
                                         X,
                                         Z,
                                         cluster_id,
                                         beta,
                                         tau,
                                         h,
                                         lambda_gamma = 1,
                                         eps = 1e-5,
                                         nuisance_control = list(
                                           reltol = 1e-13,
                                           maxit = 1000L,
                                           grad_tol = 1e-8
                                         )) {
  analytic <- profile_components_v2(
    y, X, Z, cluster_id, beta, tau, h, lambda_gamma,
    gamma_init = NULL, need_hessian = TRUE,
    nuisance_control = nuisance_control
  )
  obj_fn <- function(b) profile_objective_v2(
    y, X, Z, cluster_id, b, tau, h, lambda_gamma,
    nuisance_control = nuisance_control
  )
  score_fn <- function(b) profile_components_v2(
    y, X, Z, cluster_id, b, tau, h, lambda_gamma,
    need_hessian = FALSE, nuisance_control = nuisance_control
  )$score

  grad_fd <- finite_difference_gradient_v2(obj_fn, beta, eps = eps)
  hess_fd <- finite_difference_jacobian_v2(score_fn, beta, eps = eps)
  hess_fd <- (hess_fd + t(hess_fd)) / 2

  data.frame(
    tau = tau,
    h = h,
    p = ncol(X),
    q = ncol(Z),
    n_clusters = length(unique(cluster_id)),
    gradient_max_error = max(abs(analytic$score - grad_fd)),
    hessian_max_error = max(abs(analytic$hessian - hess_fd)),
    profile_identity_error = analytic$profile_identity_error,
    max_nuisance_gradient = analytic$max_nuisance_gradient,
    nuisance_converged = analytic$nuisance_converged,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# OSQP warm-start chained path for the Dantzig row and its CV (Lever 1).
# Implementation layer only: the frozen Dantzig program, mu grid, smallest-
# feasible rule and CV loss are unchanged. Per fold-Hessian the constraint
# matrix A and objective P are identical across (mu, coordinate); one sparse
# KKT factorisation per fold is reused and every subsequent solve is a warm-
# started O((n+m)^2) step. Selected by env JPDMEQR_DANTZIG_SOLVER=osqp.
# -----------------------------------------------------------------------------

dantzig_build_osqp_model_v2 <- function(H, p) {
  if (!requireNamespace("osqp", quietly = TRUE)) stop("osqp required.")
  osqpSetup <- getFromNamespace("osqpSetup", "osqp")
  n <- 2L * p
  P <- as(Matrix::Matrix(0, n, n, sparse = TRUE), "dgCMatrix")
  q <- c(rep(0, p), rep(1, p))
  Hz <- cbind(H, matrix(0, p, p))
  Iz1 <- cbind(Matrix::Diagonal(p), -Matrix::Diagonal(p))
  Iz2 <- cbind(-Matrix::Diagonal(p), -Matrix::Diagonal(p))
  Iz3 <- cbind(matrix(0, p, p), Matrix::Diagonal(p))
  A <- methods::as(rbind(Hz, Iz1, Iz2, Iz3), "dgCMatrix")
  osqpSettings0 <- getFromNamespace("osqpSettings", "osqp")
  list(P = P, q = q, A = A)
}

dantzig_osqp_solve_chain_v2 <- function(model_parts, e, mu_values,
                                        feas_tol = 1e-6, max_iter = 20000L) {
  # model_parts: list(P, q, A) from dantzig_build_osqp_model_v2; mu_values in
  # ascending order. Returns a list of per-mu solutions (omega, residual,
  # feasible, status) with warm-start chaining.
  osqpSetup <- getFromNamespace("osqpSetup", "osqp")
  osqpSolve <- getFromNamespace("osqpSolve", "osqp")
  osqpUpdate <- getFromNamespace("osqpUpdate", "osqp")
  osqpWarmStart <- getFromNamespace("osqpWarmStart", "osqp")
  osqpSettings0 <- getFromNamespace("osqpSettings", "osqp")
  p <- length(e)
  pars <- osqpSettings0(verbose = FALSE, polishing = FALSE,
                        eps_abs = 1e-7, eps_rel = 1e-7, max_iter = max_iter)
  mdl <- osqpSetup(P = model_parts$P, q = model_parts$q, A = model_parts$A,
                   l = c(e - mu_values[1], rep(-Inf, 3L * p)),
                   u = c(e + mu_values[1], rep(0, p), rep(0, p), rep(Inf, p)),
                   pars = pars)
  out <- vector("list", length(mu_values))
  x_prev <- NULL; y_prev <- NULL
  for (ii in seq_along(mu_values)) {
    mu <- mu_values[ii]
    if (ii > 1L) {
      osqpUpdate(mdl, q_new = model_parts$q,
                 l_new = c(e - mu, rep(-Inf, 3L * p)),
                 u_new = c(e + mu, rep(0, p), rep(0, p), rep(Inf, p)),
                 Px = NULL, Px_idx = NULL, Ax = NULL, Ax_idx = NULL)
      osqpWarmStart(mdl, x = x_prev, y = y_prev)
    }
    sol <- tryCatch(osqpSolve(mdl), error = function(err) err)
    if (inherits(sol, "error")) {
      out[[ii]] <- list(mu = mu, omega = NULL, residual = NA_real_,
                        feasible = FALSE, status = conditionMessage(sol))
      next
    }
    x_prev <- sol$x; y_prev <- sol$y
    omega <- as.numeric(sol$x[seq_len(p)])
    residual <- max(abs(as.numeric(model_parts$A[seq_len(p), seq_len(p)] %*% omega) - e))
    out[[ii]] <- list(mu = mu, omega = omega, residual = residual,
                      feasible = residual <= mu + feas_tol,
                      status = sol$info$status %||% "")
  }
  out
}

# Contract-identical wrapper for the final-row/fallback Dantzig solve via the
# chained OSQP path (same return structure as solve_dantzig_row_v2).
solve_dantzig_row_osqp_contract_v2 <- function(H,
                                               coordinate,
                                               mu_grid,
                                               feas_tol = 1e-6) {
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
  e <- numeric(p); e[k] <- 1
  parts <- dantzig_build_osqp_model_v2(H, p)
  sols <- dantzig_osqp_solve_chain_v2(parts, e, mu_grid, feas_tol = feas_tol)
  attempts <- vector("list", length(mu_grid))
  for (ii in seq_along(sols)) {
    s <- sols[[ii]]
    attempts[[ii]] <- data.frame(mu = mu_grid[ii],
                                 status = if (isTRUE(s$feasible)) "optimal" else "infeasible",
                                 residual = s$residual, message = s$status %||% "")
    if (isTRUE(s$feasible) && !is.null(s$omega)) {
      return(list(
        omega = s$omega,
        coordinate = k,
        coordinate_name = colnames(H)[k] %||% as.character(k),
        mu = mu_grid[ii],
        residual = s$residual,
        solver = "OSQP-chained",
        solver_status = s$status %||% "optimal",
        l1_norm = sum(abs(s$omega)),
        l2_norm = sqrt(sum(s$omega^2)),
        attempts = do.call(rbind, attempts[seq_len(ii)]),
        feasible = TRUE
      ))
    }
  }
  list(
    omega = rep(NA_real_, p),
    coordinate = k,
    coordinate_name = colnames(H)[k] %||% as.character(k),
    mu = NA_real_,
    residual = NA_real_,
    solver = "OSQP-chained",
    solver_status = "infeasible_on_grid",
    l1_norm = NA_real_,
    l2_norm = NA_real_,
    attempts = do.call(rbind, attempts),
    feasible = FALSE
  )
}

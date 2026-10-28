# ---------------------------------------------------------------------------
# Layer-1 computational optimisation: compact native CLARABEL Dantzig path.
#
# WHY: the frozen Dantzig row  min ||omega||_1  s.t.  ||H omega - e_k||_inf <= mu
# is solved 4 coords x 8 mu x 4 folds = 128 times per replication, and that term
# is ~75-80% of task time (see the computational-optimisation plan). Measured on
# this machine (p=200 / p=500, warmup + median of 7/5 reps):
#
#   form                              p=200    p=500   max|domega| vs CVXR
#   A: CVXR reference (frozen)        0.420s   5.590s  0 (baseline)
#   B: native, x=(omega,t), 5p cons   0.460s   4.920s   3.5e-11
#   C: native, x=(omega,t), 4p cons   0.260s   4.360s   2.2e-16  <-- this file
#   D: native L1-split (wp,wm)        1.010s  15.610s   2.2e-10
#
# Form C drops the redundant `-t <= 0` block: with `omega - t <= 0` and
# `-omega - t <= 0` both present, adding them gives `-2t <= 0`, i.e. `t >= 0`.
# The block is therefore implied and its removal changes neither the feasible
# set nor the objective. It also reproduces the CVXR reference to machine
# precision (2.2e-16), which indicates the reference canonicalisation is
# already this compact form.
#
# Parity detail (do not overstate): agreement with the frozen CVXR reference is
# at machine precision, but the compact form is NOT bit-identical to the
# redundant 5-block native form - removing the block changes the
# interior-point iteration path and moves omega by ~3e-10. Both are far below
# any decision threshold (mu + 1e-6), and selected mu / feasibility are asserted
# equal in the parity test.
#
# SCOPE: statistical-object-neutral. The LP is identical (same objective, same
# constraint set, same mu grid, same smallest-feasible rule); only the cone
# formulation handed to the same solver (CLARABEL) differs. Enabled explicitly
# via JPDMEQR_DANTZIG_SOLVER=clarabel_native; the default remains the frozen
# CVXR path. Parity is enforced by tests/testthat/test-dantzig-native-parity.R.
# ---------------------------------------------------------------------------

# CLARABEL solver-status code -> label. The R package returns the status as a
# 1-based integer whose order matches clarabel::solver_status_descriptions()
# (1 = Unsolved, 2 = Solved, 3 = PrimalInfeasible, ...). Code 2 was observed on
# every successful solve in probing, consistent with that order. The label is
# built from the package's own description table when available, so it cannot
# drift from the installed clarabel version. NOTE: feasibility decisions in this
# file deliberately do NOT depend on this label - they use the primal residual
# against the same `mu + 1e-6` tolerance as the frozen CVXR path - so a wrong or
# missing label can never change a selected mu.
clarabel_status_label_v2 <- function(status) {
  s <- suppressWarnings(as.integer(status))
  if (length(s) != 1L || is.na(s)) return(paste0("status_", as.character(status)))
  labs <- tryCatch(names(clarabel::solver_status_descriptions()),
                   error = function(e) NULL)
  if (!is.null(labs) && s >= 1L && s <= length(labs)) {
    lab <- labs[s]
    if (!is.na(lab) && nzchar(lab)) return(tolower(lab))
  }
  switch(as.character(s),
    "1" = "unsolved",
    "2" = "solved",
    "3" = "primal_infeasible",
    "4" = "dual_infeasible",
    "5" = "almost_solved",
    "6" = "almost_primal_infeasible",
    "7" = "almost_dual_infeasible",
    "8" = "max_iterations",
    "9" = "max_time",
    "10" = "numerical_error",
    "11" = "insufficient_progress",
    paste0("status_", s)
  )
}

# Compact cone form for the frozen Dantzig LP (form C above).
#   x = (omega, t) in R^{2p};  min 1't
#   blocks: omega - t <= 0 ; -omega - t <= 0 ; H omega <= e + mu ; -H omega <= -e + mu
# All four blocks live in the nonnegative cone, so A x <= b with b depending on
# (e, mu) only. A itself depends on H alone, so it is mu- and coordinate-invariant.
dantzig_lp_compact_build_v2 <- function(H, p = nrow(H)) {
  Ip <- Matrix::Diagonal(p)
  Zp <- Matrix::Matrix(0, p, p, sparse = TRUE)
  A <- rbind(
    cbind(Ip, -Ip),
    cbind(-Ip, -Ip),
    cbind(H, Zp),
    cbind(-H, Zp)
  )
  methods::as(methods::as(A, "CsparseMatrix"), "dgCMatrix")
}

# Compact-form right-hand side for a given (e, mu); length 4p.
dantzig_lp_compact_rhs_v2 <- function(p, e, mu) {
  c(numeric(p), numeric(p), e + mu, -e + mu)
}

# Native CLARABEL Dantzig row. Returns exactly the structure of
# solve_dantzig_row_v2 so it is a drop-in replacement on the same mu-grid
# semantics (ascending mu, first feasible mu wins).
solve_dantzig_row_native_v2 <- function(H,
                                        coordinate,
                                        mu_grid,
                                        solver_opts = list()) {
  if (!requireNamespace("clarabel", quietly = TRUE)) {
    stop("Package 'clarabel' is required for the native Dantzig path.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for the native Dantzig path.")
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

  A <- dantzig_lp_compact_build_v2(H, p)
  q <- c(numeric(p), rep(1, p))
  cones <- list(l = nrow(A))

  attempts <- vector("list", length(mu_grid))
  for (ii in seq_along(mu_grid)) {
    mu <- mu_grid[ii]
    b <- dantzig_lp_compact_rhs_v2(p, e, mu)
    res <- tryCatch(
      do.call(
        clarabel::clarabel,
        c(list(A = A, b = b, q = q, P = NULL, cones = cones), solver_opts)
      ),
      error = function(err) err
    )
    if (inherits(res, "error")) {
      attempts[[ii]] <- data.frame(mu = mu, status = "solver_error",
                                   residual = NA_real_,
                                   message = conditionMessage(res))
      next
    }
    status <- clarabel_status_label_v2(res$status)
    omega <- if (is.null(res$x)) NULL else as.numeric(res$x[seq_len(p)])
    if (is.null(omega) || length(omega) != p || any(!is.finite(omega))) {
      attempts[[ii]] <- data.frame(mu = mu, status = status, residual = NA_real_,
                                   message = "No finite solution returned")
      next
    }
    # Same acceptance criterion as the frozen CVXR path (residual against
    # mu + 1e-6). Deliberately independent of the status label.
    residual <- max(abs(as.numeric(H %*% omega) - e))
    attempts[[ii]] <- data.frame(mu = mu, status = status, residual = residual,
                                 message = "")
    if (residual <= mu + 1e-6) {
      return(list(
        omega = omega,
        coordinate = k,
        coordinate_name = colnames(H)[k] %||% as.character(k),
        mu = mu,
        residual = residual,
        solver = "CLARABEL-native",
        solver_status = status,
        l1_norm = sum(abs(omega)),
        l2_norm = sqrt(sum(omega^2)),
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
    solver = "CLARABEL-native",
    solver_status = "infeasible_on_grid",
    l1_norm = NA_real_,
    l2_norm = NA_real_,
    attempts = do.call(rbind, attempts),
    feasible = FALSE
  )
}
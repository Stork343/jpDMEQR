args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
destination <- if (length(args) >= 2L) args[2] else
  file.path(root, "results", "preflight", "baseline_reference_20261022")
if (dir.exists(destination)) stop("The reference-fidelity output directory already exists.")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, environment())
}
objects <- list()
records <- list()
record <- function(method, calculation, discrepancy, tolerance, source) {
  records[[length(records) + 1L]] <<- data.frame(
    method_id = method, calculation = calculation, discrepancy = discrepancy,
    tolerance = tolerance,
    passed = is.finite(discrepancy) && discrepancy <= tolerance,
    source = source, stringsAsFactors = FALSE
  )
  utils::write.csv(do.call(rbind, records), file.path(destination, "checks.csv"), row.names = FALSE)
  saveRDS(objects, file.path(destination, "reference_objects.rds"))
}
max_difference <- function(x, y) {
  if (length(x) != length(y) || any(!is.finite(x)) || any(!is.finite(y))) return(Inf)
  max(abs(as.numeric(x) - as.numeric(y)))
}
cvxr_quantile <- function(X, y, tau, penalty) {
  beta <- CVXR::Variable(ncol(X))
  residual <- y - X %*% beta
  objective <- CVXR::sum_entries(tau * CVXR::pos(residual) +
                                  (1 - tau) * CVXR::pos(-residual)) +
    CVXR::p_norm(diag(rep(penalty, length.out = ncol(X)), nrow = ncol(X)) %*% beta, 1)
  problem <- CVXR::Problem(CVXR::Minimize(objective))
  solved <- CVXR::psolve(problem, solver = "CLARABEL", feastol = 1e-9,
                          abstol = 1e-9, reltol = 1e-9)
  if (!identical(solved$status, "optimal")) stop("The independent quantile LP did not solve optimally.")
  list(beta = as.numeric(solved$getValue(beta)), objective = solved$value)
}

dat <- generate_profile_qr_data_v2(n_clusters = 40L, p = 3L, s = 2L,
                                  tau = 0.25, q = 1L, m_values = 3:5,
                                  sigma_b0 = 0.8, seed = 20261120L)
coordinates <- colnames(dat$X)
tuning <- list(lambda_beta = 0.3, h = 0.4, lambda_gamma = 1)
pooled <- fit_benchmark_pooled_qr_lasso_v2(dat, dat$tau, coordinates, tuning, 20261121L)
pooled_reference <- cvxr_quantile(dat$X, dat$y, dat$tau, tuning$lambda_beta)
objects$pooled <- list(data = dat, adapter = pooled, reference = pooled_reference)
record("POOLED-QR-LASSO", "independent_LP_coefficients",
       max_difference(pooled$beta_hat, pooled_reference$beta), 1e-5, "CVXR quantile LP")

frame <- data.frame(y = dat$y, cluster = factor(dat$cluster_id), dat$X,
                     check.names = FALSE)
formula <- stats::reformulate(coordinates, response = "y", intercept = FALSE)
lqmm_adapter <- fit_benchmark_lqmm_v2(dat, dat$tau, coordinates, tuning, 20261122L)
set.seed(20261122L)
lqmm_reference <- lqmm::lqmm(
  fixed = formula, random = ~1, group = "cluster", data = frame,
  tau = dat$tau, covariance = "pdSymm", nK = 7L,
  control = lqmm::lqmmControl(LP_tol_ll = 1e-6)
)
objects$lqmm <- list(adapter = lqmm_adapter, reference = lqmm_reference)
record("LQMM", "maintained_package_coefficients",
       max_difference(lqmm_adapter$beta_hat, stats::coef(lqmm_reference)),
       1e-10, "lqmm 1.5.8 direct call")

official_path <- file.path(root, "data", "raw", "public_methods", "battagliola", "official_code", "ran_int_bias_adj.R")
if (!file.exists(official_path)) stop("The official bias-adjustment reference function is required.")
official <- new.env(parent = environment())
sys.source(official_path, envir = official)
set.seed(20261123L)
bias_reference <- official$ran_int_bias_adj(
  fixed = formula, ID = "cluster", intercept = FALSE, tau = dat$tau,
  nK = 15L, type = "normal", data = frame, B = 100L,
  seMethod = "nid", level = 0.95, plotCI = FALSE
)
bias_adapter <- fit_benchmark_bias_adj_lqmm_v2(
  dat, dat$tau, coordinates, tuning, 20261123L,
  control = list(B = 100L, nK = 15L)
)
objects$bias <- list(adapter = bias_adapter, reference = bias_reference,
                     official_path = normalizePath(official_path))
record("BIAS-ADJ-LQMM", "official_two_step_coefficients",
       max_difference(bias_adapter$beta_hat, bias_reference[[1]]), 1e-8,
       "official ran_int_bias_adj.R")
record("BIAS-ADJ-LQMM", "official_adjusted_coefficients",
       max_difference(bias_adapter$beta_tilde, bias_reference[[2]]), 1e-8,
       "official ran_int_bias_adj.R B100")
record("BIAS-ADJ-LQMM", "official_interval",
       max_difference(c(bias_adapter$ci_lower, bias_adapter$ci_upper),
                        c(as.numeric(bias_reference[[3]][1, ]),
                          as.numeric(bias_reference[[3]][2, ]))),
       1e-8, "official ran_int_bias_adj.R B100")

double_data <- generate_profile_qr_data_v2(n_clusters = 12L, p = 2L, s = 1L,
                                          tau = 0.37, q = 1L, m_values = 3:5,
                                          seed = 20261124L)
ci <- cluster_index_v2(double_data$cluster_id)
lb <- 0.7
la <- 0.5
double_adapter <- fit_double_pen_qlmm_v2(
  double_data$y, double_data$X, double_data$Z, double_data$cluster_id,
  double_data$tau, lb, la, control = list(max_iter = 100L)
)
alpha <- matrix(0, ci$n, ncol(double_data$Z))
beta <- cvxr_quantile(double_data$X, double_data$y, double_data$tau, lb)$beta
reference_converged <- FALSE
for (iteration in seq_len(100L)) {
  previous <- beta
  for (i in seq_len(ci$n)) {
    rows <- ci$rows[[i]]
    response <- double_data$y[rows] - as.numeric(double_data$X[rows, , drop = FALSE] %*% beta)
    alpha[i, ] <- cvxr_quantile(double_data$Z[rows, , drop = FALSE], response,
                                double_data$tau, la)$beta
  }
  response <- double_data$y - rowSums(double_data$Z * alpha[as.integer(ci$factor), , drop = FALSE])
  beta <- cvxr_quantile(double_data$X, response, double_data$tau, lb)$beta
  if (mean(abs(beta - previous)) < 1e-4 && max(abs(beta - previous)) < 1e-6) {
    reference_converged <- TRUE
    break
  }
}
objects$double <- list(data = double_data, adapter = double_adapter,
                       reference = list(beta = beta, alpha = alpha,
                                        converged = reference_converged, iterations = iteration))
record("DOUBLE-PEN-QLMM", "independent_block_algorithm",
       if (double_adapter$converged && reference_converged)
         max_difference(double_adapter$beta, beta) else Inf,
       1e-4, "Li Liu Luo Section 3.1 with independent CVXR block LPs")

if (as.character(utils::packageVersion("geeVerse")) != "0.3.1") {
  stop("The registered geeVerse reference version is 0.3.1.")
}
qgee_adapter <- fit_benchmark_qgee_scad_v2(dat, dat$tau, coordinates, tuning,
                                         20261125L, control = list(lambda = NULL, standardize = FALSE))
reference_X <- dat$X
attr(reference_X, "assign") <- seq_len(ncol(reference_X))
set.seed(20261125L)
qgee_reference <- geeVerse::qpgee(
  reference_X, y = dat$y,
  nobs = as.integer(table(factor(dat$cluster_id, levels = unique(dat$cluster_id)))),
  tau = dat$tau, corstr = "exchangeable", lambda = NULL, method = "HBIC",
  intercept = FALSE, nfold = 5L,
  control = geeVerse::qpgeeControl(epsilon = 1e-4, maxit = 100L, shrinkCutoff = 1e-4, standardize = FALSE)
)
objects$qgee <- list(adapter = qgee_adapter, reference = qgee_reference)
record("QGEE-SCAD", "maintained_package_complete_HBIC_grid",
       max_difference(qgee_adapter$beta_hat,
                        qgee_reference$coefficients %||% qgee_reference$beta),
       1e-10, "geeVerse 0.3.1 qpgee")

set.seed(20261126L)
X <- scale(matrix(stats::rnorm(400L), 100L, 4L))
Y <- 0.5 + as.numeric(X %*% c(0.8, -0.6, 0, 0)) + stats::rnorm(100L)
augmented_X <- cbind(intercept = 1, X)
colnames(augmented_X) <- c("intercept", paste0("x", 1:4))
sqr <- fit_sqr_iid_lasso_v2(Y, augmented_X, 0.35, 0.45, 0.05,
                            penalty_factor = c(0, rep(1, 4L)),
                            control = list(max_iter = 5000L, beta_tol = 1e-10,
                                           kkt_tol = 1e-9))
sqr_reference <- conquer::conquer.reg(X, Y, lambda = 0.05, tau = 0.35,
                                      kernel = "Gaussian", h = 0.45, penalty = "lasso",
                                      epsilon = 1e-10, iteMax = 5000L)
objects$sqr <- list(data = list(X = X, Y = Y), adapter = sqr, reference = sqr_reference)
record("SQR-DEBIASED-IID", "independent_convolution_L1_optimizer",
       if (sqr$converged) max_difference(sqr$beta, sqr_reference$coeff) else Inf,
       1e-5, "conquer 1.3.3 Gaussian lasso")
reference_beta <- as.numeric(sqr_reference$coeff)
reference_residual <- Y - as.numeric(augmented_X %*% reference_beta)
reference_score <- as.numeric(crossprod(augmented_X,
                                        stats::pnorm(-reference_residual / 0.45) - 0.35)) / length(Y)
reference_H <- crossprod(augmented_X,
                         augmented_X * (stats::dnorm(reference_residual / 0.45) / 0.45)) / length(Y)
precision <- solve_clime_matrix_iid_v2(sqr$hessian, 0.15, coordinates = 2L)
direction <- precision$directions[[2L]]
reference_rows <- lapply(seq_len(ncol(augmented_X)), function(k) {
  solve_dantzig_row_native_v2(reference_H, k, 0.15,
                              list(control = list(verbose = FALSE,
                                                   static_regularization_constant = 1e-6)))
})
objects$sqr$reference_rows <- reference_rows
if (!all(vapply(reference_rows, function(row) isTRUE(row$feasible), logical(1)))) {
  saveRDS(objects, file.path(destination, "reference_objects.rds"))
  stop("An independent CLIME reference row failed its native solver checks.")
}
reference_raw <- do.call(rbind, lapply(reference_rows, `[[`, "omega"))
reference_symmetric <- reference_raw
for (i in seq_len(nrow(reference_raw) - 1L)) {
  for (j in seq.int(i + 1L, ncol(reference_raw))) {
    chosen <- if (abs(reference_raw[i, j]) <= abs(reference_raw[j, i]))
      reference_raw[i, j] else reference_raw[j, i]
    reference_symmetric[i, j] <- reference_symmetric[j, i] <- chosen
  }
}
reference_direction <- reference_symmetric[2L, ]
Sigma <- crossprod(augmented_X) / length(Y)
actual <- c(sqr$beta[2] - sum(direction$omega * sqr$score),
             sqrt(0.35 * 0.65 * drop(crossprod(direction$omega, Sigma %*% direction$omega)) / length(Y)))
expected <- c(reference_beta[2] - sum(reference_direction * reference_score),
               sqrt(0.35 * 0.65 * drop(crossprod(reference_direction,
                                                Sigma %*% reference_direction)) / length(Y)))
objects$sqr$direction <- direction
objects$sqr$reference_direction <- reference_direction
objects$sqr$reference_raw <- reference_raw
objects$sqr$reference_symmetric <- reference_symmetric
record("SQR-DEBIASED-IID", "independent_complete_one_step_and_variance",
       max_difference(actual, expected), 1e-5,
       "conquer optimizer plus all native CLARABEL rows Eq12 symmetrisation and iid variance")

qif_data <- generate_profile_qr_data_v2(n_clusters = 80L, p = 3L, s = 3L,
                                        q = 1L, tau = .5, sigma_b0 = 1, seed = 700773L)
qif_beta <- as.numeric(quantreg::rq.fit(qif_data$X, qif_data$y, tau = .5)$coefficients)
qif_density <- rep(stats::dnorm(0) / sqrt(2), length(qif_data$y))
qif_delta <- pmin(1, 1e-5 / abs(qif_beta)^2)
qif_h <- 80^(-1 / 2)
qif_component <- qif_see_extended_source_v2(qif_data$y, qif_data$X, qif_data$cluster_id,
                                            qif_beta, .5, qif_h, qif_density, qif_delta)
qif_rows <- split(seq_along(qif_data$y), factor(qif_data$cluster_id, levels = unique(qif_data$cluster_id)))
qif_reference_score <- function(beta) {
  vapply(qif_rows, function(rows) {
    Xi <- qif_data$X[rows, , drop = FALSE]
    m <- length(rows)
    I <- diag(m)
    J <- matrix(1, m, m) - I
    adjacent <- matrix(0, m, m)
    if (m > 1L) {
      adjacent[cbind(seq_len(m - 1L), seq.int(2L, m))] <- 1
      adjacent[cbind(seq.int(2L, m), seq_len(m - 1L))] <- 1
    }
    psi <- stats::pnorm((qif_data$y[rows] - as.numeric(Xi %*% beta)) / qif_h) - .5
    vapply(seq_len(ncol(Xi)), function(k) {
      vapply(list(I, J, adjacent), function(M) {
        (1 - qif_delta[k]) * sum(Xi[, k] * qif_density[rows] * as.numeric(M %*% psi)) +
          qif_delta[k] * beta[k]
      }, numeric(1))
    }, numeric(3)) |> as.numeric()
  }, numeric(9))
}
qif_reference_G <- qif_reference_score(qif_beta)
qif_reference_C <- tcrossprod(qif_reference_G) / 80L
qif_reference_D <- finite_difference_jacobian_v2(function(beta) rowMeans(qif_reference_score(beta)),
                                                  qif_beta, eps = 1e-6)
qif_reference_mean <- rowMeans(qif_reference_G)
qif_reference_gradient <- 2 * as.numeric(crossprod(qif_reference_D, solve(qif_reference_C, qif_reference_mean)))
qif_reference_H <- 2 * crossprod(qif_reference_D, solve(qif_reference_C, qif_reference_D))
qif_reference_step <- -as.numeric(solve(qif_reference_H, qif_reference_gradient))
qif_actual <- qif_see_source_step_v2(qif_component, qif_beta)
objects$qif <- list(data = qif_data, density = qif_density, delta = qif_delta,
                    actual = qif_actual, component = qif_component,
                    reference_G = qif_reference_G, reference_C = qif_reference_C,
                    reference_D = qif_reference_D, reference_step = qif_reference_step)
record("QIF-SEE", "independent_Algorithm1_modified_score_and_Newton_step",
       max(max_difference(qif_component$G, qif_reference_G),
           max_difference(qif_component$derivative, qif_reference_D),
           max_difference(qif_actual$step, qif_reference_step)), 1e-5,
       "Bhattacharya et al full text Eq3 Eq5 Eq6 Algorithm1 with explicitly supplied Normal density")
qif_highdim <- generate_profile_qr_data_v2(n_clusters = 30L, p = 40L, s = 3L,
                                           q = 1L, tau = .5, seed = 700774L)
qif_domain <- fit_benchmark_qif_see_v2(qif_highdim, .5, colnames(qif_highdim$X)[1:3], list(), 1L)
objects$qif$domain <- qif_domain
record("QIF-SEE", "source_ordinary_inverse_domain",
       if (qif_domain$status == "failed" && qif_domain$fit_object$matrix_dimension == 120L &&
           qif_domain$fit_object$rank_upper_bound == 30L &&
           qif_domain$fit_object$inverse_definition == "ordinary_inverse") 0 else Inf,
       0, "Eq6 covariance rank at most n for the Bp by Bp ordinary inverse")

commit <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
jsonlite::write_json(list(commit_sha = commit,
                           created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                           checks = do.call(rbind, records),
                           QIF_SEE_source_fidelity = "Algorithm1_explicit_density_and_ordinary_inverse_domain"),
                      file.path(destination, "reference_fidelity.json"),
                      auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(capture.output(sessionInfo()), file.path(destination, "sessionInfo.txt"))
if (!all(vapply(records, function(x) x$passed, logical(1)))) stop("A reference fidelity calculation failed.")
print(do.call(rbind, records))

source("research/score_validity_20261026/load_joint.R")
set.seed(2026102809L)
records <- list()
for (q in 1:4) for (m in c(2L, 5L, 8L)) for (pattern in c("separated", "repeated", "close")) {
  times <- seq(-1, 1, length.out = m)
  if (pattern == "repeated") times[] <- 0
  if (pattern == "close" && m >= 2L) times[2L] <- times[1L] + 1e-10
  Z <- vapply(0:(q - 1L), function(power) times^power, numeric(m))
  Lambda <- diag(q) + matrix(.1, q, q)
  X <- matrix(rnorm(m * 5L), m, 5L, dimnames = list(NULL, paste0("x", 1:5)))
  beta <- seq(.1, .5, length.out = 5L)
  for (tau in c(.25, .5, .75)) for (h in c(.4, .15, .05)) for (mode in c("interior", "mixed")) {
    residual <- h * runif(m, -.8, .8)
    if (mode == "mixed") residual[seq.int(1L, m, by = 2L)] <-
      h * sample(c(-2, 2), length(seq.int(1L, m, by = 2L)), replace = TRUE)
    score <- tau - kernel_epan_cdf_v2(-residual / h)
    gamma <- solve(Lambda, crossprod(Z, score) / m)
    y <- drop(X %*% beta + Z %*% gamma + residual)
    data <- sb_data(y, X, Z, rep("cluster", m))
    fitted <- compiled_cluster_profiles_v2(y - drop(X %*% beta), Z, data$rows, tau, h, Lambda,
      control = list(grad_tol = 1e-11, newton_maxit = 150L))
    H <- sb_geometry(data, fitted$phi, Lambda)
    active <- which(fitted$phi > 0)
    chart <- matrix(0, 5L, 5L)
    moment_ratio <- 0
    if (length(active)) {
      ZA <- Z[active, , drop = FALSE]
      Q <- ZA %*% solve(Lambda, t(ZA)) / m + diag(1 / fitted$phi[active], nrow = length(active))
      chart <- crossprod(X[active, , drop = FALSE], solve(Q, X[active, , drop = FALSE])) / m
      eigenvalues <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
      moment_ratio <- max(vapply(2:6, function(order)
        min(eigenvalues)^(-order) / ((.75 / h)^(order - 1L) * sum(1 / eigenvalues)), numeric(1)))
    }
    discrepancy <- max(abs(H - chart))
    nuisance_difference <- max(abs(as.numeric(fitted$gamma) - as.numeric(gamma)))
    if (discrepancy > 1e-7 || nuisance_difference > 1e-8 || moment_ratio > 1 + 1e-9)
      stop("The profile chart geometry or cofactor inequality failed.")
    records[[length(records) + 1L]] <- data.frame(q = q, m = m, pattern = pattern,
      tau = tau, h = h, residual_pattern = mode, active_coordinates = length(active),
      hessian_discrepancy = discrepancy, nuisance_discrepancy = nuisance_difference,
      maximum_moment_ratio = moment_ratio)
  }
}
output <- "results/scratch/submission_completion_20261027/profile_chart_checks.csv"
write.csv(do.call(rbind, records), output, row.names = FALSE)
cat("Verified", length(records), "profile charts; maximum Hessian discrepancy",
    max(vapply(records, `[[`, numeric(1), "hessian_discrepancy")), "\n")

n <- 90L
m <- 4L
p <- 30L
X <- matrix(rnorm(n * m * p), n * m, p, dimnames = list(NULL, paste0("x", seq_len(p))))
time <- rep(c(-1, 0, 0, 1), n)
Z <- cbind(1, time, time^2)
random <- matrix(rnorm(n * 3L), n, 3L)
beta <- c(1, -.8, .5, rep(0, p - 3L))
cluster <- sprintf("c%03d", rep(seq_len(n), each = m))
y <- drop(X %*% beta) + rowSums(Z * random[rep(seq_len(n), each = m), ]) + rnorm(n * m)
fit <- ds_fit_profile(y, X, Z, cluster, .5, c("x1", "x4"), Lambda = diag(3L))
if (any(!is.finite(fit$table$se)) || any(fit$table$max_kkt > 1e-8)) stop("The three-dimensional nuisance fit failed.")
saveRDS(fit, "results/scratch/submission_completion_20261027/three_nuisance_fit.rds")
print(fit$table)

source("research/submission_completion_20261027/load.R")
source("research/submission_completion_20261027/gaussian_curvature.R")
output <- "results/scratch/submission_completion_20261027/gaussian_curvature_checks"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
curvature <- submission_gaussian_curvature(.5, 3:8, 1, 1)
scalar <- submission_gaussian_curvature(.5, 1L, 1, 1)
scalar_exact <- pnorm(.5 / sqrt(2)) - pnorm(-.5 / sqrt(2))
if (abs(scalar$curvature - scalar_exact) > 1e-10) stop("The scalar Gaussian curvature differs from its closed form.")
set.seed(202610281L)
records <- list()
for (m in 3:8) {
  for (iteration in 1:50) {
    y <- rnorm(m) + rnorm(1L)
    X <- matrix(rnorm(m * 3L), m)
    fitted <- pc_exact_profile(y, matrix(1, m, 1L), list(seq_len(m)), diag(1), .5)
    active <- which(abs(fitted$dual - (-.5 / m)) > 1e-8 & abs(fitted$dual - (.5 / m)) > 1e-8)
    if (length(active) > 1L) stop("An ordinary Gaussian cluster has more than one interior score.")
    expected <- if (length(active)) crossprod(X[active, , drop = FALSE]) else matrix(0, 3L, 3L)
    derivative <- vapply(1:3, function(j) {
      plus <- pc_exact_profile(y - 1e-6 * X[, j], matrix(1, m, 1L), list(seq_len(m)), diag(1), .5)
      minus <- pc_exact_profile(y + 1e-6 * X[, j], matrix(1, m, 1L), list(seq_len(m)), diag(1), .5)
      drop(-crossprod(X, plus$dual - minus$dual)) / 2e-6
    }, numeric(3))
    error <- max(abs(expected - derivative))
    if (error > 1e-7) stop("The exact unsmoothed cluster Hessian differs from the score derivative.")
    records[[length(records) + 1L]] <- data.frame(cluster_size = m, iteration = iteration,
      active = length(active), hessian_error = error)
  }
}
sizes <- sample(3:8, 100000L, replace = TRUE)
group <- rep(seq_along(sizes), sizes)
innovations <- rnorm(length(group))
response <- rep(rnorm(length(sizes)), sizes) + innovations
rows <- split(seq_along(response), group)
exact <- pc_exact_profile(response, matrix(1, length(response), 1L), rows, diag(1), .5)
stein <- drop(rowsum(innovations * exact$dual, group, reorder = FALSE))
mcse <- sd(stein) / sqrt(length(stein))
if (abs(mean(stein) - curvature$curvature) > 5 * mcse) stop("The independent Gaussian population check differs from quadrature.")
write.csv(data.frame(quadrature = curvature$curvature, population_mean = mean(stein),
  population_mcse = mcse, clusters = length(sizes)), file.path(output, "population_check.csv"), row.names = FALSE)
for (id in c("A01", "A02", "A03")) {
  result <- readRDS(file.path("results/raw/simulation/submission_20261027_main_joint/replications", id, "replicate_00001/result.rds"))
  diagnostic <- submission_gaussian_diagnostics(result, curvature)
  if (!diagnostic$available || any(!is.finite(diagnostic$coordinates$scaled_remainder)))
    stop("The Gaussian scaling diagnostic is unavailable.")
  saveRDS(diagnostic, file.path(output, paste0(id, ".rds")))
}
write.csv(do.call(rbind, records), file.path(output, "checks.csv"), row.names = FALSE)
saveRDS(curvature, file.path(output, "population_curvature.rds"))
print(curvature)
cat("Verified", length(records), "cluster derivatives and three saved inference calculations\n")

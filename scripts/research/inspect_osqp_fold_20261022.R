.libPaths(c(normalizePath("results/scratch/profile_pair_package_20261022/library"), .libPaths()))
library(jpDMEQR)
source("scripts/00_source_v2.R")
source_v2_module(".", "profile_v2", environment())
input_path <- "results/scratch/repository_package_checks_20261022/actual_failure_diagnostics_02/fold_actual.rds"
object <- readRDS(input_path)
folds <- object$input
output <- "results/scratch/repository_package_checks_20261022/osqp_fold_diagnostics_01"
if (dir.exists(output)) stop("The diagnostic output already exists.")
dir.create(output, recursive = TRUE)
n <- sum(folds$fold_counts)
p <- nrow(folds$H)
e <- c(1, rep(0, p - 1L))
mus <- c(0.05, 0.25)
setup <- getFromNamespace("osqpSetup", "osqp")
solve <- getFromNamespace("osqpSolve", "osqp")
update <- getFromNamespace("osqpUpdate", "osqp")
warm_start <- getFromNamespace("osqpWarmStart", "osqp")
settings <- osqp::osqpSettings(verbose = FALSE, polish = FALSE,
                                eps_abs = 1e-7, eps_rel = 1e-7, max_iter = 20000L)
records <- list()
objects <- list()
row <- 0L
for (fold in which(folds$fold_counts > 0 & folds$fold_counts < n)) {
  count <- folds$fold_counts[[fold]]
  H_validation <- folds$fold_sums[[fold]] / count
  H_training <- (n * folds$H - count * H_validation) / (n - count)
  H_training <- (H_training + t(H_training)) / 2
  parts <- dantzig_build_osqp_model_v2(H_training, p)
  make_model <- function(mu) setup(
    P = parts$P, q = parts$q, A = parts$A,
    l = c(e - mu, rep(-Inf, 3L * p)),
    u = c(e + mu, rep(0, 2L * p), rep(Inf, p)), pars = settings
  )
  warm_model <- make_model(mus[[1]])
  previous <- NULL
  for (j in seq_along(mus)) {
    mu <- mus[[j]]
    if (j > 1L) {
      update(warm_model, q_new = parts$q,
             l_new = c(e - mu, rep(-Inf, 3L * p)),
             u_new = c(e + mu, rep(0, 2L * p), rep(Inf, p)),
             Px = NULL, Px_idx = NULL, Ax = NULL, Ax_idx = NULL)
      warm_start(warm_model, x = previous$x, y = previous$y)
    }
    warm <- solve(warm_model)
    cold <- solve(make_model(mu))
    previous <- warm
    for (mode in c("cold", "warm")) {
      result <- if (mode == "cold") cold else warm
      row <- row + 1L
      omega <- as.numeric(result$x[seq_len(p)])
      residual <- max(abs(as.numeric(H_training %*% omega) - e))
      records[[row]] <- data.frame(
        fold = fold, mode = mode, mu = mu, status = result$info$status,
        iterations = result$info$iter, residual = residual,
        accepted_residual = residual <= mu + 1e-6,
        min_eigenvalue = min(eigen(H_training, symmetric = TRUE, only.values = TRUE)$values),
        stringsAsFactors = FALSE
      )
      objects[[row]] <- list(training_H = H_training, parts = parts, solution = result)
    }
  }
}
table <- do.call(rbind, records)
utils::write.csv(table, file.path(output, "cold_warm.csv"), row.names = FALSE)
saveRDS(objects, file.path(output, "raw_solutions.rds"))
mass_error <- max(abs(n * folds$H - Reduce(`+`, folds$fold_sums)))
empty_fold_mass <- max(abs(folds$fold_sums[[3L]]))
jsonlite::write_json(list(n = n, original_n = sum(folds$real_counts),
                          total_mass_error = mass_error, empty_fold_nonzero_mass = empty_fold_mass,
                          osqp_version = as.character(packageVersion("osqp"))),
                     file.path(output, "input_diagnostics.json"), auto_unbox = TRUE, pretty = TRUE)
print(table)
print(c(total_mass_error = mass_error, empty_fold_nonzero_mass = empty_fold_mass))

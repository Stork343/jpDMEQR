.libPaths(c(normalizePath("results/scratch/profile_pair_package_20261022/library"), .libPaths()))
library(jpDMEQR)
source("scripts/00_source_v2.R")
source_v2_module(".", "profile_v2", environment())
source_v2_module(".", "simulation_v2", environment())
previous <- readRDS("results/scratch/repository_package_checks_20261022/actual_failure_diagnostics_02/fold_actual.rds")$input
output <- "results/scratch/repository_package_checks_20261022/consistent_fold_diagnostics_01"
if (dir.exists(output)) stop("The diagnostic output already exists.")
dir.create(output, recursive = TRUE)
levels <- unique(as.character(previous$dat$cluster_id))
bucket <- rep(seq_len(4L), length.out = length(levels))
kept <- levels[previous$fold_counts[bucket] > 0]
dat <- subset_clusters_v2(previous$dat, kept)
sums <- previous$fold_sums
for (j in which(previous$fold_counts == 0L)) sums[[j]][] <- 0
H <- Reduce(`+`, sums) / sum(previous$fold_counts)
direct <- profile_components_v2(dat$y, dat$X, dat$Z, dat$cluster_id,
                                 previous$beta, 0.5, 0.4, need_hessian = TRUE)
identity_error <- max(abs(H - direct$hessian))
if (identity_error > 1e-12) stop("The subset profile does not reproduce the retained fold sum.")
results <- lapply(c("osqp", "dual_warm", "cvxr"), function(mode) {
  withr::with_envvar(c(JPDMEQR_DANTZIG_SOLVER = mode), select_mu_inverse_hessian_cv_v2(
    H = H, coordinate = 1L, mu_grid = c(0.05, 0.25), fold_sums = sums,
    fold_counts = previous$fold_counts, n_folds = 4L
  ))
})
names(results) <- c("osqp", "dual_warm", "cvxr")
saveRDS(list(dat = dat, H = H, fold_sums = sums, fold_counts = previous$fold_counts,
             results = results, profile_identity_error = identity_error),
        file.path(output, "actual_consistent_folds.rds"))
print(identity_error)
for (mode in names(results)) {
  cat(mode, "\n")
  print(results[[mode]][c("selected", "mu", "folds_used", "candidates")])
}

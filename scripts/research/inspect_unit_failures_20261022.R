root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(normalizePath("results/scratch/profile_pair_package_20261022/library"), .libPaths()))
library(jpDMEQR)
source("scripts/00_source_v2.R")
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "simulation_v2", environment())
output <- "results/scratch/repository_package_checks_20261022/actual_failure_diagnostics_02"
if (dir.exists(output)) stop("The diagnostic directory already exists.")
dir.create(output, recursive = TRUE)
previous <- readRDS("results/scratch/repository_package_checks_20261022/actual_failure_diagnostics_01/dantzig_actual.rds")
H <- previous$H
mu_grid <- previous$mu_grid
reference <- withr::with_envvar(c(JPDMEQR_DANTZIG_SOLVER = "cvxr"),
                               solve_dantzig_row_v2(H, 1L, mu_grid))
native <- withr::with_envvar(c(JPDMEQR_DANTZIG_SOLVER = "clarabel_native"),
                            solve_dantzig_row_v2(H, 1L, mu_grid))
dual <- withr::with_envvar(c(JPDMEQR_DANTZIG_SOLVER = "dual_warm"),
                          solve_dantzig_row_v2(H, 1L, mu_grid))
saveRDS(list(H = H, mu_grid = mu_grid, reference = reference, native = native, dual = dual),
        file.path(output, "dantzig_actual.rds"))
print(reference)
print(native[c("feasible", "mu", "residual", "solver")])
print(dual[c("feasible", "mu", "residual", "solver")])
print(packageVersion("CVXR"))
print(packageVersion("osqp"))
print(formals(osqp::osqpSettings))
declarations <- parse("tests/testthat/test-cv-fold-failure-semantics.R")
for (expression in declarations) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      as.character(expression[[2L]]) %in% c("cv_fold_fixture", "held_out_defects")) {
    eval(expression, envir = environment())
  }
}
initial <- cv_fold_fixture()
counts <- initial$real_counts
counts[3L] <- 0L
folds <- cv_fold_fixture(counts = counts)
osqp_cv <- withr::with_envvar(c(JPDMEQR_DANTZIG_SOLVER = "osqp"),
  select_mu_inverse_hessian_cv_v2(
    H = folds$H, coordinate = 1L, mu_grid = c(0.05, 0.25),
    fold_sums = folds$fold_sums, fold_counts = folds$fold_counts, n_folds = 4L
  ))
dual_cv <- withr::with_envvar(c(JPDMEQR_DANTZIG_SOLVER = "dual_warm"),
  select_mu_inverse_hessian_cv_v2(
    H = folds$H, coordinate = 1L, mu_grid = c(0.05, 0.25),
    fold_sums = folds$fold_sums, fold_counts = folds$fold_counts, n_folds = 4L
  ))
saveRDS(list(input = folds, osqp = osqp_cv, dual = dual_cv), file.path(output, "fold_actual.rds"))
print(osqp_cv)
print(dual_cv[c("selected", "reason", "mu", "folds_used", "candidates")])

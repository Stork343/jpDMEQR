args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("需要提供待运行的 R 脚本及其参数。")
root <- normalizePath(".", mustWork = TRUE)
script <- normalizePath(args[[1L]], mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "research", "audit_profile_pair_geometry_20261022.R"))
configuration <- profile_pair_audit_configuration_20261022("production")
fields <- c(dantzig_solver = "JPDMEQR_DANTZIG_SOLVER", cluster_solver = "JPDMEQR_CLUSTER_SOLVER",
            assembly = "JPDMEQR_ASSEMBLY", cv_workers = "JPDMEQR_CV_CORES",
            hessian_refresh = "JPDMEQR_HESSIAN_REFRESH", shared_first_stage = "JPDMEQR_SHARED_FIRST_STAGE",
            matrix_backend = "JPDMEQR_MATRIX_BACKEND", identity_detail = "JPDMEQR_IDENTITY_DETAIL",
            proximal_scaling = "JPDMEQR_PROXIMAL_SCALING",
            unpenalised_optimizer = "JPDMEQR_UNPENALISED_OPTIMIZER")
if (!setequal(names(configuration), names(fields))) stop("计算配置字段不完整。")
do.call(Sys.setenv, setNames(lapply(configuration[names(fields)], as.character), fields))
Sys.setenv(JPDMEQR_RESULT_STORAGE = "streaming", JPDMEQR_INITIAL_ACCELERATION = "clinical_kkt",
           R_LIBS_USER = .libPaths()[1L])
result <- processx::run(file.path(R.home("bin"), "Rscript"), c(script, args[-1L]),
                        wd = root, env = Sys.getenv(), echo = TRUE,
                        error_on_status = FALSE, timeout = Inf)
quit(status = result$status)

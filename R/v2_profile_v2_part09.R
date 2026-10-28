# 运行配置记录：求解器开关、折级并行的进程数与第一阶段共享开关。冻结身份文件
# 与单次运行目录用它保留一次运行实际采用的求解路径，使同一提交下的两次运行
# 可以区分开。
solver_configuration_v2 <- function() {
  list(
    dantzig_solver = Sys.getenv(
      "JPDMEQR_DANTZIG_SOLVER",
      unset = getOption("jpDMEQR.dantzig_solver", "cvxr")
    ),
    cluster_solver = Sys.getenv(
      "JPDMEQR_CLUSTER_SOLVER",
      unset = getOption("jpDMEQR.cluster_solver", "bfgs")
    ),
    assembly = Sys.getenv(
      "JPDMEQR_ASSEMBLY",
      unset = getOption("jpDMEQR.assembly", "sequential")
    ),
    cv_workers = as.integer(Sys.getenv(
      "JPDMEQR_CV_CORES",
      unset = as.character(getOption("jpDMEQR.cv_cores", 1L))
    )),
    hessian_refresh = hessian_refresh_v2(),
    unpenalised_optimizer = unpenalised_optimizer_v2(),
    proximal_scaling = Sys.getenv("JPDMEQR_PROXIMAL_SCALING", unset = "none"),
    shared_first_stage = shared_first_stage_v2(),
    matrix_backend = getOption("jpDMEQR.matrix_backend",
                                Sys.getenv("JPDMEQR_MATRIX_BACKEND", unset = "R_BLAS")),
    identity_detail = getOption("jpDMEQR.identity_detail",
                                 Sys.getenv("JPDMEQR_IDENTITY_DETAIL", unset = "cluster"))
  )
}

# 第一阶段 Hessian 的刷新间隔：环境变量、选项、冻结值 5 三级退让，取值 5 时一条
# 拟合路径与冻结实现相同。fit_profile_lasso_calibrated_v2 把它注入前两阶段的控制
# 字段；round-5 重拟合、刀切、总体目标资产与几何审计都保持 5。
hessian_refresh_v2 <- function() {
  value <- Sys.getenv(
    "JPDMEQR_HESSIAN_REFRESH",
    unset = as.character(getOption("jpDMEQR.hessian_refresh", 5L))
  )
  number <- suppressWarnings(as.integer(value))
  if (length(number) != 1L || is.na(number) || number < 1L) {
    stop("JPDMEQR_HESSIAN_REFRESH must be a positive integer, got: ", value)
  }
  number
}

# 第一阶段共享开关（Lever 2）：打开时，同一次重复内校准第一阶段请求逐项相同的
# 方法（PROFILE-DQR、POP-H、POSTREFIT-EXACT-H）复用同一份拟合对象。拟合是确定性
# 计算，共享结果与独立拟合逐位相同；默认关闭保持冻结参考路径。
shared_first_stage_v2 <- function() {
  value <- Sys.getenv(
    "JPDMEQR_SHARED_FIRST_STAGE",
    unset = as.character(getOption("jpDMEQR.shared_first_stage", FALSE))
  )
  # as.logical 对 "1"/"0" 这类字符串返回 NA，取值集合显式列出。
  flag <- if (value %in% c("1", "T", "TRUE", "true", "True")) {
    TRUE
  } else if (value %in% c("0", "F", "FALSE", "false", "False")) {
    FALSE
  } else {
    NA
  }
  if (is.na(flag)) {
    stop("JPDMEQR_SHARED_FIRST_STAGE must be TRUE or FALSE, got: ", value)
  }
  flag
}

# 身份文件里记录的求解路径是否与当前配置相同。一次在默认配置下取得的通过结论
# 不能为另一次采用加速开关的运行担保。shared_first_stage 字段出现之前的记录
# 按默认关闭对待：默认配置下旧身份文件仍然有效，打开共享后旧记录不再匹配，
# 需要在相同配置下重跑几何门。
configuration_matches_v2 <- function(recorded, current = solver_configuration_v2()) {
  if (is.null(recorded)) return(FALSE)
  identical(as.character(recorded$dantzig_solver), as.character(current$dantzig_solver)) &&
    identical(as.character(recorded$cluster_solver), as.character(current$cluster_solver)) &&
    identical(as.character(recorded$assembly), as.character(current$assembly)) &&
    identical(as.integer(recorded$cv_workers %||% NA_integer_), as.integer(current$cv_workers)) &&
    identical(as.integer(recorded$hessian_refresh %||% NA_integer_),
              as.integer(current$hessian_refresh)) &&
    identical(recorded$unpenalised_optimizer %||% "proximal_gradient",
              current$unpenalised_optimizer %||% "proximal_gradient") &&
    identical(recorded$proximal_scaling %||% "none", current$proximal_scaling %||% "none") &&
    identical(as.logical(recorded$shared_first_stage %||% FALSE),
              as.logical(current$shared_first_stage)) &&
    identical(recorded$matrix_backend %||% "R_BLAS", current$matrix_backend) &&
    identical(recorded$identity_detail %||% "cluster", current$identity_detail)
}

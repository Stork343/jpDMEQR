# Lever 2 共享第一阶段拟合的奇偶性测试
# （COMPUTATIONAL_OPTIMISATION_PLAN.md 第 2 节）：
#   开关默认关闭，记录进运行配置；
#   缓存对同一请求返回逐位相同的对象，对不同请求分别保存；
#   PROFILE-DQR / POP-H / POSTREFIT-EXACT-H 在共享缓存与独立拟合下
#   除 runtime_sec 外逐位相同；
#   调度器只在开关打开时注入缓存。

root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
source_v2_module(root, "benchmark_adapters_v2", envir = .GlobalEnv)

shared_fixture <- function(seed = 404L) {
  generate_profile_qr_data_v2(
    n_clusters = 40, p = 12, s = 3, tau = 0.5, q = 1L,
    m_values = 3:5, m_rule = "uniform", sigma_b0 = 1,
    error_dist = "gaussian", seed = seed
  )
}

shared_tuning <- list(h = 0.5, lambda_beta = 0.15, lambda_gamma = 1,
                      mu_grid = c(0.3, 0.6, 1.2))

strip_runtime <- function(ans) {
  timing_fields <- c("runtime_sec", "first_stage_sec", "refit_sec", "inference_sec")
  ans[timing_fields] <- NULL
  if (!is.null(ans$fit_object)) ans$fit_object[timing_fields] <- NULL
  if (!is.null(ans$inference_object$fit)) {
    ans$inference_object$fit[timing_fields] <- NULL
  }
  ans
}

make_shared_pop_asset <- function(coords) {
  p_sub <- length(coords)
  H_pop <- diag(2, p_sub)
  directions <- lapply(coords, function(nm) {
    k <- match(nm, coords)
    omega <- rep(0, p_sub)
    omega[k] <- 1 / 2
    e <- numeric(p_sub); e[k] <- 1
    list(omega = setNames(omega, coords),
         residual = max(abs(H_pop %*% omega - e)),
         l1_norm = sum(abs(omega)),
         l2_norm = sqrt(sum(omega^2)))
  })
  names(directions) <- coords
  list(experiment_id = "TEST",
       H_population = H_pop,
       directions = directions,
       n_population = 100000L,
       repeats = 4L,
       implementation_commit = "test",
       config_sha256 = "test",
       nuisance_converged = TRUE)
}

unset_shared_env <- function(code) {
  withr::with_envvar(c(JPDMEQR_SHARED_FIRST_STAGE = NA), code)
}

unset_shared_option <- function(code) {
  withr::with_options(list(jpDMEQR.shared_first_stage = NULL), code)
}

# --- 1. 开关解析：默认关闭，环境变量与选项两级退让 ---------------------------
testthat::test_that("the sharing switch defaults to off and follows env and option", {
  testthat::expect_false(unset_shared_env(unset_shared_option(shared_first_stage_v2())))
  testthat::expect_identical(
    withr::with_envvar(c(JPDMEQR_SHARED_FIRST_STAGE = "1"), shared_first_stage_v2()),
    TRUE
  )
  testthat::expect_identical(
    withr::with_envvar(c(JPDMEQR_SHARED_FIRST_STAGE = "TRUE"), shared_first_stage_v2()),
    TRUE
  )
  testthat::expect_identical(
    withr::with_envvar(c(JPDMEQR_SHARED_FIRST_STAGE = "0"), shared_first_stage_v2()),
    FALSE
  )
  testthat::expect_identical(
    unset_shared_env(withr::with_options(
      list(jpDMEQR.shared_first_stage = TRUE), shared_first_stage_v2()
    )),
    TRUE
  )
  testthat::expect_identical(
    withr::with_options(
      list(jpDMEQR.shared_first_stage = TRUE),
      withr::with_envvar(c(JPDMEQR_SHARED_FIRST_STAGE = "0"), shared_first_stage_v2())
    ),
    FALSE
  )
  testthat::expect_error(
    withr::with_envvar(c(JPDMEQR_SHARED_FIRST_STAGE = "maybe"), shared_first_stage_v2()),
    "JPDMEQR_SHARED_FIRST_STAGE"
  )
})

# --- 2. 运行配置记录与旧身份文件的兼容 ---------------------------------------
testthat::test_that("the configuration record carries the switch and legacy records match the default", {
  rec_default <- unset_shared_env(unset_shared_option(solver_configuration_v2()))
  testthat::expect_identical(rec_default$shared_first_stage, FALSE)
  rec_on <- unset_shared_option(withr::with_envvar(
    c(JPDMEQR_SHARED_FIRST_STAGE = "1"), solver_configuration_v2()
  ))
  testthat::expect_identical(rec_on$shared_first_stage, TRUE)

  legacy <- rec_default
  legacy$shared_first_stage <- NULL
  testthat::expect_true(configuration_matches_v2(legacy, rec_default))
  testthat::expect_false(configuration_matches_v2(legacy, rec_on))
  testthat::expect_true(configuration_matches_v2(rec_on, rec_on))
  testthat::expect_false(configuration_matches_v2(rec_default, rec_on))
})

# --- 3. 缓存机制：同一请求复用，不同请求分条 ---------------------------------
testthat::test_that("the cache reuses one fit per request and separates distinct requests", {
  train <- shared_fixture()
  coords <- colnames(train$X)[1:3]
  cache <- new.env(parent = emptyenv())
  control <- list()
  v1 <- obtain_calibrated_first_stage_v2(cache, train, 0.5, coords,
                                         shared_tuning, control)
  v2 <- obtain_calibrated_first_stage_v2(cache, train, 0.5, coords,
                                         shared_tuning, control)
  testthat::expect_identical(v2, v1)
  testthat::expect_identical(length(cache$entries), 1L)
  # 共享的前提是确定性：独立重算的对象与缓存逐位相同（时间列除外）。
  fresh <- obtain_calibrated_first_stage_v2(NULL, train, 0.5, coords,
                                            shared_tuning, control)
  testthat::expect_identical(strip_runtime(fresh), strip_runtime(v1))
  # 罚因子不同的请求是新条目，不改写已有缓存。
  control0 <- list(penalty_factor = rep(0, ncol(train$X)))
  v3 <- obtain_calibrated_first_stage_v2(cache, train, 0.5, coords,
                                         shared_tuning, control0)
  testthat::expect_identical(length(cache$entries), 2L)
  testthat::expect_false(identical(v3, v1))
  # 缓存参数类型错误就地停止。
  testthat::expect_error(
    obtain_calibrated_first_stage_v2(list(), train, 0.5, coords,
                                     shared_tuning, control),
    "environment"
  )
})

# --- 4. PROFILE-DQR：共享与独立拟合逐位相同 ----------------------------------
testthat::test_that("PROFILE-DQR is bit-identical with and without the shared cache", {
  train <- shared_fixture()
  coords <- colnames(train$X)[1:3]
  cache <- new.env(parent = emptyenv())
  ans_shared <- fit_benchmark_profile_dqr_v2(
    train, 0.5, coords, shared_tuning, seed = 9L,
    control = list(shared_first_stage = cache)
  )
  ans_plain <- fit_benchmark_profile_dqr_v2(
    train, 0.5, coords, shared_tuning, seed = 9L
  )
  testthat::expect_identical(strip_runtime(ans_shared), strip_runtime(ans_plain))
  testthat::expect_identical(length(cache$entries), 1L)
})

# --- 5. POP-H 与 POSTREFIT：复用同一份拟合并保持逐位相同 ---------------------
testthat::test_that("POP-H and POSTREFIT reuse the shared fit bit-identically", {
  train <- shared_fixture()
  coords <- colnames(train$X)[1:3]
  asset <- make_shared_pop_asset(coords)

  pop_plain <- fit_benchmark_profile_pop_h_v2(
    train, 0.5, coords, shared_tuning, seed = 9L,
    population_direction_asset = asset
  )
  post_plain <- fit_benchmark_postrefit_exact_h_v2(
    train, 0.5, coords, shared_tuning, seed = 9L
  )

  cache <- new.env(parent = emptyenv())
  dqr_shared <- fit_benchmark_profile_dqr_v2(
    train, 0.5, coords, shared_tuning, seed = 9L,
    control = list(shared_first_stage = cache)
  )
  testthat::expect_identical(length(cache$entries), 1L)
  pop_shared <- fit_benchmark_profile_pop_h_v2(
    train, 0.5, coords, shared_tuning, seed = 9L,
    population_direction_asset = asset,
    control = list(shared_first_stage = cache)
  )
  post_shared <- fit_benchmark_postrefit_exact_h_v2(
    train, 0.5, coords, shared_tuning, seed = 9L,
    control = list(shared_first_stage = cache)
  )
  # 三个方法命中同一条缓存，条目数保持 1。
  testthat::expect_identical(length(cache$entries), 1L)
  testthat::expect_identical(strip_runtime(pop_shared), strip_runtime(pop_plain))
  testthat::expect_identical(strip_runtime(post_shared), strip_runtime(post_plain))
  # 共享的拟合对象就是各方法携带的 fit_object。
  testthat::expect_identical(pop_shared$fit_object, dqr_shared$fit_object)
  testthat::expect_identical(post_shared$fit_object, dqr_shared$fit_object)
})

# --- 6. 调度器：只在开关打开时注入缓存 ---------------------------------------
testthat::test_that("the dispatcher shares only when the switch is on", {
  train <- shared_fixture()
  coords <- colnames(train$X)[1:3]
  asset <- make_shared_pop_asset(coords)

  cache_off <- new.env(parent = emptyenv())
  ans_off <- unset_shared_env(fit_benchmark_by_id_v2(
    "PROFILE-DQR", train, 0.5, coords, shared_tuning, 9L,
    context = list(shared_first_stage = cache_off)
  ))
  testthat::expect_null(cache_off$entries)

  cache_on <- new.env(parent = emptyenv())
  ctx_on <- list(shared_first_stage = cache_on,
                 population_direction_asset = asset)
  ans_on <- withr::with_envvar(
    c(JPDMEQR_SHARED_FIRST_STAGE = "1"),
    fit_benchmark_by_id_v2("PROFILE-DQR", train, 0.5, coords, shared_tuning, 9L,
                           context = ctx_on)
  )
  pop_on <- withr::with_envvar(
    c(JPDMEQR_SHARED_FIRST_STAGE = "1"),
    fit_benchmark_by_id_v2("PROFILE-DQR-POP-H", train, 0.5, coords,
                           shared_tuning, 9L, context = ctx_on)
  )
  testthat::expect_identical(length(cache_on$entries), 1L)
  testthat::expect_identical(strip_runtime(ans_on), strip_runtime(ans_off))

  pop_off <- unset_shared_env(fit_benchmark_by_id_v2(
    "PROFILE-DQR-POP-H", train, 0.5, coords, shared_tuning, 9L,
    context = ctx_on
  ))
  # 开关关闭时不写入缓存，结果与打开时逐位相同。
  testthat::expect_identical(length(cache_on$entries), 1L)
  testthat::expect_identical(strip_runtime(pop_on), strip_runtime(pop_off))
})

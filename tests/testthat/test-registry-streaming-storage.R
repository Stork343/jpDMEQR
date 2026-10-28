root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, .GlobalEnv)
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"), local = TRUE)

testthat::test_that("streaming storage preserves every actual generated replication table and failure", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "compiled", JPDMEQR_DANTZIG_SOLVER = "dual_warm",
                        JPDMEQR_UNPENALISED_OPTIMIZER = "newton_armijo", JPDMEQR_PROXIMAL_SCALING = "cluster_rms"))
  library(jpDMEQR)
  registry <- file.path(root, "tests", "fixtures", "simulation_streaming_pilot_20261022.yml")
  tag <- paste0("storage_actual_", format(Sys.time(), "%Y%m%dT%H%M%OS6"))
  memory <- run_registry_v2(root, registry, paste0(tag, "_memory"), max_reps = 2L, storage_mode = "memory")
  streaming <- run_registry_v2(root, registry, paste0(tag, "_stream"), jobs = 2L,
                               max_reps = 2L, storage_mode = "streaming")
  raw_memory <- readRDS(file.path(memory$run_dir, "raw_objects.rds"))
  index <- readRDS(file.path(streaming$run_dir, "raw_objects.rds"))
  testthat::expect_identical(index$format, "replication_file_index_v1")
  testthat::expect_equal(nrow(index$files), 2L)
  testthat::expect_lt(as.numeric(object.size(index)), as.numeric(object.size(raw_memory)))
  testthat::expect_equal(nrow(memory$metrics), 6L)
  testthat::expect_equal(nrow(streaming$metrics), 6L)
  testthat::expect_equal(sum(streaming$metrics$method_id == "QIF-SEE" & streaming$metrics$status == "failed"), 2L)
  stable <- function(table) {
    omit <- names(table)[grepl("run_id|runtime|_sec$|memory|_rss|_wset", names(table))]
    table[, setdiff(names(table), omit), drop = FALSE]
  }
  for (name in c("metrics", "coordinates", "theory", "screening")) {
    testthat::expect_equal(stable(memory[[name]]), stable(streaming[[name]]), tolerance = 1e-10)
  }
  testthat::expect_equal(length(streaming$errors), 0L)
  compact_configuration <- registry_storage_configuration_v2("streaming")
  compact_configuration$matrix_inline_entries <- 16L
  for (i in seq_along(raw_memory)) {
    complete <- raw_memory[[i]]
    compact <- compact_registry_replication_v2(complete, compact_configuration, root)
    table <- registry_table_record_v2(compact)
    # 同一次实际计算的全部表格逐项比较，包含实际runtime。
    testthat::expect_identical(table$metrics, complete$metrics)
    testthat::expect_identical(table$coordinates, complete$coordinates)
    testthat::expect_identical(table$theory, complete$theory)
    testthat::expect_identical(table$screening_table, screening_record_v2(complete))
    testthat::expect_identical(compact$method_diagnostics$`PROFILE-DQR`$directions,
                               complete$method_diagnostics$`PROFILE-DQR`$directions)
    testthat::expect_identical(compact$method_diagnostics$`PROFILE-DQR`$cv,
                               complete$method_diagnostics$`PROFILE-DQR`$cv)
    recovered <- read_registry_replication_v2(file.path(streaming$run_dir, "raw_objects.rds"), index$files$key[i])
    testthat::expect_equal(recovered$config$replicate, i)
    testthat::expect_identical(recovered$storage$format, "compact_diagnostics_v1")
  }
})

testthat::test_that("compact diagnostics retain actual failed fits and exact coefficient score and history values", {
  withr::local_envvar(c(JPDMEQR_CLUSTER_SOLVER = "compiled", JPDMEQR_PROXIMAL_SCALING = "cluster_rms"))
  library(jpDMEQR)
  dat <- generate_profile_qr_data_v2(n_clusters = 24L, p = 12L, s = 2L, q = 1L,
                                      tau = .5, seed = 22624L)
  limited <- fit_profile_lasso_v2(dat$y, dat$X, dat$Z, dat$cluster_id, .5,
                                   (log(12) / 24)^(1 / 4), lambda_beta = .001,
                                   control = list(max_iter = 1L))
  testthat::expect_false(limited$converged)
  failure <- tryCatch(profile_pair_check_fit_v2(limited, "penalised_fit", 1e-3, 1e-7),
                       error = function(condition) condition)
  testthat::expect_s3_class(failure, "profile_pair_error")
  compact <- compact_registry_diagnostics_v2(failure, 16L)
  testthat::expect_identical(compact$message, conditionMessage(failure))
  testthat::expect_identical(compact$fields$failure_stage, failure$failure_stage)
  testthat::expect_identical(compact$fields$fit$beta, failure$fit$beta)
  testthat::expect_length(failure$fit$components$score, 12L)
  testthat::expect_equal(nrow(failure$fit$history), 1L)
  testthat::expect_identical(compact$fields$fit$components$score, failure$fit$components$score)
  testthat::expect_identical(compact$fields$fit$history, failure$fit$history)
  component <- profile_components_v2(dat$y, dat$X, dat$Z, dat$cluster_id,
                                       failure$fit$beta, .5, failure$fit$components$h, need_hessian = TRUE)
  stored_component <- compact_registry_diagnostics_v2(component, 16L)
  H <- component$hessian
  testthat::expect_identical(stored_component$hessian$dimensions, dim(H))
  testthat::expect_equal(stored_component$hessian$values$finite_mean, mean(H))
  testthat::expect_identical(stored_component$hessian$diagonal, diag(H))
})

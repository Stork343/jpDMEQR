root <- v2_require_root()
source(file.path(root, "scripts", "00_source_v2.R"), local = FALSE)
source_v2_module(root, "profile_v2", envir = .GlobalEnv)
source_v2_module(root, "simulation_v2", envir = .GlobalEnv)
source_v2_module(root, "metrics_v2", envir = .GlobalEnv)
sys.source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"),
           envir = environment())

config_path <- file.path(root, "config", "simulation_pilot.csv")
fixture_commit <- current_commit_v2(root)
fixture_sha <- file_checksum_v2(config_path)

# 构造一个分片运行目录，仅用于核对恢复时的身份检查。
make_shard_fixture <- function(run_id,
                               final = TRUE,
                               allow_missing_benchmarks = FALSE,
                               base_seed = 20260817,
                               max_reps = 1,
                               experiment_ids = "P01",
                               sha = fixture_sha,
                               recorded_commit = fixture_commit,
                               lock_pid = NULL,
                               complete = FALSE,
                               request = TRUE) {
  fixture_root <- tempfile("resume-root-")
  dir.create(fixture_root, recursive = TRUE, showWarnings = FALSE)
  writeLines(fixture_commit, file.path(fixture_root, "implementation_commit.txt"))
  run_dir <- file.path(fixture_root, "results", "raw", "simulation", run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  if (request) {
    write_json_v2(list(
      run_id = run_id,
      final = final,
      allow_missing_benchmarks = allow_missing_benchmarks,
      base_seed = base_seed,
      max_reps = if (is.finite(max_reps)) max_reps else "Inf",
      experiment_ids = experiment_ids
    ), file.path(run_dir, "run_request.json"))
  }
  writeLines(sha, file.path(run_dir, "config_sha256.txt"))
  writeLines(recorded_commit, file.path(run_dir, "implementation_commit.txt"))
  if (!is.null(lock_pid)) {
    writeLines(c(paste0("pid=", lock_pid), "started_utc=fixture"), file.path(run_dir, "RUNNING.lock"))
  }
  if (complete) {
    write_json_v2(list(pass = TRUE), file.path(run_dir, "completed_manifest.json"))
  }
  list(root = fixture_root, run_dir = run_dir)
}

resume_call <- function(fixture, run_id) {
  suppressWarnings(run_registry_v2(
    root = fixture$root,
    config_path = config_path,
    run_id = run_id,
    max_reps = 1,
    experiment_ids = "P01",
    final = TRUE,
    allow_missing_benchmarks = FALSE,
    resume = TRUE
  ))
}

testthat::test_that("process liveness probe distinguishes live and absent pids", {
  testthat::expect_true(pid_running_v2(Sys.getpid()))
  testthat::expect_false(pid_running_v2(999999L))
  testthat::expect_false(pid_running_v2(NA_integer_))
})

testthat::test_that("an existing run directory is refused without the resume flag", {
  fixture <- make_shard_fixture("fixture_no_resume")
  testthat::expect_error(
    run_registry_v2(
      root = fixture$root,
      config_path = config_path,
      run_id = "fixture_no_resume",
      experiment_ids = "P01",
      final = TRUE,
      allow_missing_benchmarks = FALSE
    ),
    "already exists and is immutable"
  )
})

testthat::test_that("resume refuses a directory without a recorded request", {
  fixture <- make_shard_fixture("fixture_no_request", request = FALSE)
  err <- tryCatch(resume_call(fixture, "fixture_no_request"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "Cannot resume without run_request.json")
})

testthat::test_that("resume refuses a completed shard", {
  fixture <- make_shard_fixture("fixture_done", complete = TRUE)
  err <- tryCatch(resume_call(fixture, "fixture_done"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "already complete")
})

testthat::test_that("resume refuses a shard held by a running process", {
  fixture <- make_shard_fixture("fixture_locked", lock_pid = Sys.getpid())
  err <- tryCatch(resume_call(fixture, "fixture_locked"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "held by a running process")
})

testthat::test_that("resume refuses a changed registry checksum", {
  fixture <- make_shard_fixture("fixture_sha", sha = "0000")
  err <- tryCatch(resume_call(fixture, "fixture_sha"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "registry checksum differs")
})

testthat::test_that("resume refuses a changed implementation commit", {
  fixture <- make_shard_fixture("fixture_commit", recorded_commit = "other-commit")
  err <- tryCatch(resume_call(fixture, "fixture_commit"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "implementation commit differs")
})

testthat::test_that("resume refuses a changed replication cap", {
  fixture <- make_shard_fixture("fixture_cap", max_reps = 5)
  err <- tryCatch(resume_call(fixture, "fixture_cap"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "replication cap differs")
})

testthat::test_that("resume refuses a changed experiment set", {
  fixture <- make_shard_fixture("fixture_experiments", experiment_ids = "P02")
  err <- tryCatch(resume_call(fixture, "fixture_experiments"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "experiment set differs")
})

testthat::test_that("resume refuses a changed seed, final flag or benchmark policy", {
  seed_fixture <- make_shard_fixture("fixture_seed", base_seed = 1)
  seed_err <- tryCatch(resume_call(seed_fixture, "fixture_seed"),
                       error = function(e) conditionMessage(e))
  testthat::expect_match(seed_err, "base seed differs")

  final_fixture <- make_shard_fixture("fixture_final", final = FALSE)
  final_err <- tryCatch(resume_call(final_fixture, "fixture_final"),
                        error = function(e) conditionMessage(e))
  testthat::expect_match(final_err, "final flag differs")

  policy_fixture <- make_shard_fixture("fixture_policy", allow_missing_benchmarks = TRUE)
  policy_err <- tryCatch(resume_call(policy_fixture, "fixture_policy"),
                         error = function(e) conditionMessage(e))
  testthat::expect_match(policy_err, "benchmark policy differs")
})

testthat::test_that("the contract schema list is the column requirement used on resume", {
  columns <- contract_replication_columns_v2(root)
  testthat::expect_true(length(columns) >= 55)
  testthat::expect_true(all(c("run_id", "n_clusters_original", "p_original",
                              "selected_names", "m_rule", "screen_dimension",
                              "post_refit_hessian_condition_number") %in% columns))
  testthat::expect_false(any(c("mu", "screen_size", "cluster_size_rule",
                               "hessian_condition") %in% columns))
})

testthat::test_that("every recorded table has a contract section", {
  sections <- lapply(2:5, function(section) contract_schema_columns_v2(root, section))
  testthat::expect_true(all(vapply(sections, function(x) "run_id" %in% x, logical(1))))
  testthat::expect_true(all(c("experiment_id", "replicate", "method_id") %in% sections[[1]]))
  testthat::expect_true(all(c("coordinate", "estimated_se", "dantzig_mu") %in% sections[[2]]))
  testthat::expect_true(all(c("kkt_residual", "precision_row_error") %in% sections[[3]]))
  testthat::expect_true(all(c("screen_engine", "overlap_count") %in% sections[[4]]))
  testthat::expect_identical(sections[[1]], contract_replication_columns_v2(root))
})

testthat::test_that("the screening record carries the run identifier", {
  dat <- generate_profile_qr_data_v2(
    n_clusters = 8, p = 6, s = 2, tau = 0.5, q = 1, m_values = 3:4,
    m_rule = "uniform", sigma_b0 = 1, error_dist = "gaussian", seed = 11
  )
  registry <- read_simulation_registry_v2(config_path)
  config <- config_row_to_list_v2(registry[1, , drop = FALSE])
  prepared <- prepare_analysis_data_v2(dat, config, seed = 11)
  record <- screening_record_v2(list(screening = prepared$screening, config = config,
                                     seed = 11, run_id = "fixture_run"))
  testthat::expect_equal(record$run_id, "fixture_run")
  testthat::expect_equal(record$experiment_id, config$experiment_id)
  testthat::expect_true(all(c("run_id", "experiment_id", "screen_method",
                              "overlap_count") %in% names(record)))
  contract <- contract_schema_columns_v2(root, 5L)
  testthat::expect_identical(names(record), contract)
})

testthat::test_that("the screening record separates disjoint and same-sample designs", {
  # 维度缩小到测试规模，筛选设计沿用注册表中相应单元的设置。
  dat <- generate_profile_qr_data_v2(
    n_clusters = 12, p = 40, s = 3, tau = 0.5, q = 1, m_values = 3:4,
    m_rule = "uniform", sigma_b0 = 1, error_dist = "gaussian", seed = 21
  )
  registry <- read_simulation_registry_v2(
    file.path(root, "config", "simulation_main_b200.csv")
  )
  record_for <- function(experiment_id) {
    row <- registry[registry$experiment_id == experiment_id, , drop = FALSE]
    config <- config_row_to_list_v2(row)
    config$screen_dim <- 20L
    prepared <- prepare_analysis_data_v2(dat, config, seed = 21)
    screening_record_v2(list(screening = prepared$screening, config = config,
                             seed = 21, run_id = "fixture_run"))
  }
  n_clusters <- length(unique(dat$cluster_id))

  split_record <- record_for("F01")
  testthat::expect_equal(split_record$screen_method, "split_quantile_score")
  testthat::expect_equal(split_record$overlap_count, 0L)
  testthat::expect_gt(split_record$screen_cluster_count, 0L)
  testthat::expect_gt(split_record$fit_cluster_count, 0L)
  testthat::expect_equal(split_record$screen_cluster_count +
                           split_record$fit_cluster_count, n_clusters)

  for (experiment_id in c("F05", "F04")) {
    same_sample <- record_for(experiment_id)
    testthat::expect_equal(same_sample$screen_cluster_count, n_clusters)
    testthat::expect_equal(same_sample$fit_cluster_count, n_clusters)
    testthat::expect_equal(same_sample$overlap_count, n_clusters)
  }

  none_record <- record_for("A01")
  testthat::expect_equal(none_record$screen_method, "none")
  testthat::expect_equal(none_record$screen_cluster_count, 0L)
  testthat::expect_equal(none_record$fit_cluster_count, n_clusters)
  testthat::expect_equal(none_record$overlap_count, 0L)
})

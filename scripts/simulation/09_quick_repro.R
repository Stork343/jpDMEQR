#!/usr/bin/env Rscript
# 审稿人尺度的本地快速复现入口（演示层级）：运行 config/simulation_quick_repro.csv
# 里的缩减单元。输出用于验证管线在本地可运行并评估计算时长，永远不进入稿件；
# 稿件数字只来自冻结注册表（config/simulation_main_b200.csv）在匹配冻结清单下的运行。
# 前置条件：同一求解配置下的严格几何门已通过
# （scripts/simulation/00_validate_profile_geometry.R --strict=true）。

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/09_quick_repro.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source_v2_module(root, "benchmark_adapters_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
config <- cli$config %||% file.path(root, "config", "simulation_quick_repro.csv")
jobs <- as.integer(cli$jobs %||% 1L)
max_reps <- as.numeric(cli$`max-reps` %||% Inf)
run_id <- cli$`run-id` %||% paste0("quick_repro_", format(Sys.time(), "%Y%m%dT%H%M%S"))
experiment_ids <- cli$`experiment-id` %||% ""
experiment_ids <- if (nzchar(experiment_ids)) {
  trimws(strsplit(experiment_ids, ",", fixed = TRUE)[[1]])
} else character()
# 与 02_run_main.R 相同的任务内并行设置：mu-CV 折并行与簇并行都保持确定性。
options(jpDMEQR.cv_cores = as.integer(Sys.getenv("JPDMEQR_CV_CORES", unset = "1")))
options(jpDMEQR.cluster_cores = as.integer(Sys.getenv("JPDMEQR_CLUSTER_CORES", unset = as.character(jobs))))

if (!file.exists(config)) stop("Configuration file is missing: ", config)

governance <- validate_governance_files_v2(root)
if (!governance$pass) {
  stop("Governance files missing: ", paste(governance$missing, collapse = ", "))
}
registry <- validate_registry_contract_v2(root, config, "preflight")
if (!registry$pass) {
  stop("Quick reproduction registry validation failed: ",
       paste(registry$problems, collapse = "; "))
}
if (length(registry$warnings)) {
  warning(paste(registry$warnings, collapse = "; "))
}
geometry <- validate_geometry_gate_v2(root, final = FALSE)
if (!geometry$pass) {
  stop("Geometry gate must pass before the quick reproduction run (", geometry$reason,
       "). Run first, under the same solver-configuration environment:",
       " Rscript scripts/simulation/00_validate_profile_geometry.R --strict=true")
}

solver_config <- solver_configuration_v2()
message("QUICK REPRODUCTION TIER: demonstration outputs, never manuscript numbers.")
message("Solver configuration: ",
        paste(names(solver_config), as.character(unlist(solver_config)),
              sep = "=", collapse = " "))

ans <- run_registry_v2(
  root = root,
  config_path = config,
  run_id = run_id,
  jobs = jobs,
  max_reps = max_reps,
  allow_missing_benchmarks = TRUE,
  experiment_ids = experiment_ids,
  final = FALSE
)

cat("\nStatus counts by method\n")
print(table(ans$metrics$method_id, ans$metrics$status))
cat("\nMedian runtime per replication, seconds\n")
print(stats::aggregate(runtime_sec ~ experiment_id + method_id,
                       data = ans$metrics, FUN = stats::median, na.rm = TRUE))
message("Quick reproduction run written to: ", ans$run_dir)
message("Summarise with: Rscript scripts/simulation/03_summarise.R --run-id=", run_id)

bad <- ans$metrics[ans$metrics$status %in% c("failed", "not_implemented"), , drop = FALSE]
if (nrow(bad)) {
  print(bad[, c("experiment_id", "replicate", "method_id", "status",
                "failure_stage", "failure_message"), drop = FALSE])
}
if (nrow(bad) || length(ans$errors)) {
  stop("Quick reproduction run has failures; inspect ", ans$run_dir)
}
message("Quick reproduction run completed without method failures.")

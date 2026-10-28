args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("需要提供新的运行名称和计算进程数量。")
root <- normalizePath(".", mustWork = TRUE)
run_id <- args[[1L]]
workers <- as.integer(args[[2L]])
if (!nzchar(run_id) || grepl("[/\\\\]", run_id) || is.na(workers) || workers < 1L) {
  stop("运行名称或计算进程数量无效。")
}
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, environment())
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
configuration <- solver_configuration_v2()
if (!identical(configuration$cluster_solver, "compiled") ||
    !identical(configuration$dantzig_solver, "dual_warm") ||
    !identical(Sys.getenv("JPDMEQR_RESULT_STORAGE"), "streaming")) {
  stop("需要已通过检查的计算配置和 streaming 结果保存方式。")
}
output <- file.path(root, "results", "raw", "study_pipeline", run_id)
if (dir.exists(output)) stop("运行目录已经存在。")
dir.create(output, recursive = TRUE)
commit <- current_commit_v2(root)
configurations <- setNames(lapply(c("preflight", "pilot", "main"), function(phase) {
  path <- file.path("config", paste0("simulation_pair_", phase, "_20261022.yml"))
  registry <- read_simulation_registry_v2(path)
  check <- validate_registry_contract_v2(root, path, phase)
  if (!check$pass) stop(paste(check$problems, collapse = "; "))
  list(path = path, scenarios = nrow(registry), replications = sum(registry$replications))
}), c("preflight", "pilot", "main"))
expected <- c(preflight = 20L, pilot = 1200L, main = 45850L)
if (!identical(as.numeric(vapply(configurations, function(x) x$replications, numeric(1))),
               as.numeric(expected))) stop("登记的重复数量发生变化。")
request <- list(run_id = run_id, commit_sha = commit, configurations = configurations,
                workers = workers, solver_configuration = configuration,
                storage_mode = "streaming", started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
write_json_v2(request, file.path(output, "request.json"))
records <- list()
run_stage <- function(stage, script, arguments) {
  if (!identical(current_commit_v2(root), commit)) stop("运行期间的 source commit 发生变化。")
  started <- Sys.time()
  command <- processx::run(
    file.path(R.home("bin"), "Rscript"), c(script, arguments), wd = root,
    env = Sys.getenv(), stdout = file.path(output, paste0(stage, ".stdout.log")),
    stderr = file.path(output, paste0(stage, ".stderr.log")),
    error_on_status = FALSE, timeout = Inf
  )
  records[[stage]] <<- list(stage = stage, script = script, arguments = arguments,
                           exit_status = command$status,
                           elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
                           completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  write_json_v2(records, file.path(output, "stage_records.json"))
  if (command$status != 0L) stop("计算阶段未通过：", stage, "，记录位于 ", output)
  if (!identical(current_commit_v2(root), commit)) stop("运行期间的 source commit 发生变化。")
}
preflight_id <- paste0(run_id, "_preflight")
pilot_id <- paste0(run_id, "_pilot")
main_id <- paste0(run_id, "_main")
run_stage("preflight", "scripts/simulation/02_run_main.R",
          c(paste0("--config=", configurations$preflight$path), paste0("--run-id=", preflight_id),
            paste0("--jobs=", workers)))
micro <- read_json_v2(micro_preflight_manifest_path_v2(root, "profile_pair_20261022"))
if (!isTRUE(micro$pass) || !identical(micro$commit_sha, commit)) stop("计算前检查未通过。")
run_stage("pilot", "scripts/simulation/01_run_pilot.R",
          c(paste0("--config=", configurations$pilot$path), paste0("--run-id=", pilot_id),
            paste0("--jobs=", workers)))
run_stage("pilot_summary", "scripts/simulation/03_summarise.R",
          c(paste0("--run-id=", pilot_id), "--pilot-gate=true"))
pilot <- read_json_v2(pilot_gate_manifest_path_v2(root, "profile_pair_20261022"))
if (!isTRUE(pilot$pass) || !identical(pilot$commit_sha, commit) ||
    !isTRUE(pilot$complete_registered_replications)) stop("pilot 的正式检查未通过。")
run_stage("freeze", "scripts/simulation/05_preflight_freeze.R",
          paste0("--config=", configurations$main$path))
verify_freeze_manifest_v2(root, configurations$main$path)
run_stage("main", "scripts/simulation/02_run_main.R",
          c(paste0("--config=", configurations$main$path), paste0("--run-id=", main_id),
            paste0("--jobs=", workers)))
run_stage("main_summary", "scripts/simulation/03_summarise.R", paste0("--run-id=", main_id))
write_json_v2(c(request, list(completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                             stages = records, preflight_run_id = preflight_id,
                             pilot_run_id = pilot_id, main_run_id = main_id)),
              file.path(output, "completed_manifest.json"))

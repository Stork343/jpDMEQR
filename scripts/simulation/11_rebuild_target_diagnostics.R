# 目标资产诊断表的重建脚本：目标资产里保存着每次拟合的收敛、KKT 残差、身份误差、
# 带宽规则与耗时字段，本脚本按这些记录重新写出
# results/intermediate/targets/profile_target_diagnostics.csv。
# 分批构建目标资产时，诊断表由单次构建写出；本脚本使该表始终是全部资产的完整记录。
#
# 用法：Rscript --vanilla scripts/simulation/11_rebuild_target_diagnostics.R

args <- commandArgs(trailingOnly = TRUE)
root <- getOption("jpDMEQR.root", Sys.getenv("JPDMEQR_ROOT", unset = "."))
if (length(args)) {
  for (a in args) {
    if (grepl("^--root=", a)) root <- sub("^--root=", "", a)
  }
}
root <- normalizePath(root, mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

target_dir <- file.path(root, "results", "intermediate", "targets")
assets <- sort(list.files(target_dir, pattern = "\\.rds$", full.names = TRUE))
if (!length(assets)) stop("No target assets found under ", target_dir)

rows <- list()
index <- 0L
for (path in assets) {
  obj <- readRDS(path)
  if (!length(obj$fits)) next
  experiment_id <- obj$experiment_id %||% sub("\\.rds$", "", basename(path))
  for (rr in seq_along(obj$fits)) {
    fit <- obj$fits[[rr]]
    index <- index + 1L
    rows[[index]] <- data.frame(
      experiment_id = experiment_id,
      asset_type = "profile_target",
      repeat_id = rr,
      seed = fit$seed %||% NA_integer_,
      n_population = fit$n_population %||% obj$n_population %||% NA_integer_,
      bandwidth_rule = fit$bandwidth_rule %||% NA_character_,
      h_target = fit$h_target %||% NA_real_,
      h_analysis_rule = fit$h_analysis_rule %||% NA_character_,
      smoothing_target_distance = obj$smoothing_target_distance_by_repeat[rr] %||% NA_real_,
      converged = fit$converged %||% NA,
      kkt_residual = fit$kkt_residual %||% NA_real_,
      profile_identity_error = fit$profile_identity_error %||% NA_real_,
      max_nuisance_gradient = fit$max_nuisance_gradient %||% NA_real_,
      nuisance_converged = fit$nuisance_converged %||% NA,
      runtime_sec = fit$runtime_sec %||% NA_real_,
      cpu_sec = fit$cpu_sec %||% NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
if (!length(rows)) stop("No fit records found in the target assets.")
merged <- do.call(rbind, rows)
merged <- merged[order(merged$asset_type, merged$experiment_id, merged$repeat_id), ,
                 drop = FALSE]
write_atomic_csv_v2(merged, file.path(target_dir, "profile_target_diagnostics.csv"))
message(sprintf("Wrote %d fit records for %d target assets to %s",
                nrow(merged), length(unique(merged$experiment_id)),
                file.path(target_dir, "profile_target_diagnostics.csv")))
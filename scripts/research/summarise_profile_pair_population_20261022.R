root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2")) source_v2_module(root, module)
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
run_id <- cli$`run-id` %||% stop("Supply a population-summary run-id.")
out <- file.path(root, "results", "scratch", "profile_pair_population_final_20261022", run_id)
if (dir.exists(out)) stop("The population summary directory already exists.")
dir.create(out, recursive = TRUE)
paths <- c("config/simulation_pair_pilot_20261022.yml", "config/simulation_pair_main_20261022.yml",
           "config/simulation_pair_preflight_20261022.yml")
configs <- unlist(lapply(paths, function(path) {
  registry <- read_simulation_registry_v2(file.path(root, path))
  lapply(seq_len(nrow(registry)), function(i) config_row_to_list_v2(registry[i, , drop = FALSE]))
}), recursive = FALSE)
commit <- current_commit_v2(root)
target_identity <- list()
rows <- list()
for (config in configs) {
  if (!identical(config$target_mode, "profile_mc")) next
  id <- config$experiment_id
  path <- profile_pair_target_asset_path_v2(root, id)
  asset <- readRDS(path)
  validate_profile_pair_target_asset_v2(asset, config, final = TRUE, expected_commit = commit)
  target_identity[[id]] <- list(beta = asset$beta_star_mc, seeds = asset$seeds,
                                 commit = asset$implementation_commit)
  rows[[length(rows) + 1L]] <- data.frame(
    experiment_id = id, asset_kind = "target", p = length(asset$target_names),
    n_population = asset$n_population, repeats = asset$repeats,
    target_source = asset$target_source, h_target = asset$h_target, h_analysis = asset$h_analysis,
    population_engine = asset$acceptance$population_engine,
    implementation_commit = asset$implementation_commit,
    kkt_role = "population_root", max_kkt_residual = asset$max_kkt_residual,
    max_nuisance_gradient = asset$max_nuisance_gradient,
    profile_identity_error = asset$profile_identity_error,
    covariance_identity_error = max(vapply(asset$repeat_assets, function(record) {
      max(vapply(record[c("target", "refinement", "analysis")], function(fit) {
        fit$moments$covariance_identity_error
      }, numeric(1)))
    }, numeric(1))),
    maximum_coordinate_mcse = asset$target_mc_se,
    maximum_inverse_defect = NA_real_, path = path)
  rm(asset)
  invisible(gc(FALSE))
}
for (config in configs) {
  if (!"PROFILE-DQR-POP-H" %in% config$methods) next
  id <- config$experiment_id
  path <- profile_pair_direction_asset_path_v2(root, id)
  asset <- readRDS(path)
  validate_profile_pair_direction_asset_v2(asset, config, final = TRUE, expected_commit = commit)
  if (identical(config$target_mode, "profile_mc")) {
    target <- target_identity[[id]]
    if (!identical(asset$beta_target, target$beta) ||
        !identical(asset$target_asset_identity$beta_star_mc, target$beta) ||
        !identical(asset$target_asset_identity$seeds, target$seeds) ||
        !identical(asset$target_asset_identity$implementation_commit, target$commit)) {
      stop("The population direction and independent target identities differ: ", id)
    }
  }
  rows[[length(rows) + 1L]] <- data.frame(
    experiment_id = id, asset_kind = "direction", p = length(asset$target_names),
    n_population = asset$n_population, repeats = asset$repeats,
    target_source = asset$target_source, h_target = asset$h_target, h_analysis = asset$h_analysis,
    population_engine = asset$acceptance$population_engine,
    implementation_commit = asset$implementation_commit,
    kkt_role = "not_applicable_direction_evaluation", max_kkt_residual = NA_real_,
    max_nuisance_gradient = asset$max_nuisance_gradient,
    profile_identity_error = asset$profile_identity_error,
    covariance_identity_error = asset$covariance_identity_error,
    maximum_coordinate_mcse = max(asset$sigma_mc_se),
    maximum_inverse_defect = max(vapply(asset$directions, `[[`, numeric(1), "residual")),
    path = path)
  rm(asset)
  invisible(gc(FALSE))
}
assets <- do.call(rbind, rows)
stopifnot(sum(assets$asset_kind == "target") == 8L,
          sum(assets$asset_kind == "direction") == 27L,
          all(assets$n_population >= 100000L), all(assets$repeats >= 4L))
utils::write.csv(assets, file.path(out, "population_assets.csv"), row.names = FALSE)
audit_rows <- list()
for (id in c("E03", "E04", "E05", "E06", "LAMBDA1_BASELINE")) {
  path <- file.path(root, "results", "intermediate", "profile_pair_20261022", "construction_v3",
                    "structural_audits", paste0(id, ".rds"))
  asset <- readRDS(path)
  if (!identical(asset$asset_kind, "paired_structural_target_audit") ||
      !identical(asset$score_definition, "profile_bandwidth_pair_2_half_minus_full")) {
    stop("An incompatible structural audit was encountered: ", id)
  }
  validate_profile_pair_target_asset_v2(asset$paired_target_asset, asset$config,
                                        final = TRUE, expected_commit = commit)
  check <- validate_structural_target_audit_v2(asset, final = TRUE)
  if (!check$valid) stop(paste(check$problems, collapse = "; "))
  audit_rows[[length(audit_rows) + 1L]] <- data.frame(
    experiment_id = id, lambda_gamma = asset$config$lambda_gamma,
    registered_p = asset$registered_p, audit_p = asset$audit_p,
    n_population = asset$n_population, repeats = asset$repeats,
    maximum_displacement = max(abs(asset$target_displacement)),
    maximum_coordinate_mcse = asset$target_mc_se,
    max_kkt_residual = asset$max_kkt_residual,
    max_nuisance_gradient = asset$max_nuisance_gradient,
    profile_identity_error = asset$profile_identity_error,
    implementation_commit = asset$implementation_commit, path = path)
  rm(asset)
  invisible(gc(FALSE))
}
audits <- do.call(rbind, audit_rows)
utils::write.csv(audits, file.path(out, "structural_audits.csv"), row.names = FALSE)
manifest <- list(pass = TRUE, target_count = 8L, direction_count = 27L, structural_audit_count = 5L,
                 construction_commits = unique(c(assets$implementation_commit, audits$implementation_commit)),
                 validation_commit = commit,
                 minimum_population = min(c(assets$n_population, audits$n_population)),
                 minimum_repeats = min(c(assets$repeats, audits$repeats)),
                 max_kkt = max(c(assets$max_kkt_residual, audits$max_kkt_residual), na.rm = TRUE),
                 max_nuisance_gradient = max(c(assets$max_nuisance_gradient, audits$max_nuisance_gradient)),
                 max_profile_identity_error = max(c(assets$profile_identity_error, audits$profile_identity_error)),
                 max_covariance_identity_error = max(assets$covariance_identity_error),
                 max_inverse_defect = max(assets$maximum_inverse_defect, na.rm = TRUE))
jsonlite::write_json(manifest, file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines(c(
  "# Paired population validation",
  "",
  sprintf("%d 个总体目标和 %d 个总体方向通过了逐个回读校验。每个对象使用 %d 个独立 clusters，并保留 %d 次独立重复。",
          manifest$target_count, manifest$direction_count, manifest$minimum_population, manifest$minimum_repeats),
  "",
  sprintf("全部 target roots 的最大 KKT 残差为 %.8g。nuisance gradient 的最大值为 %.8g。",
          manifest$max_kkt, manifest$max_nuisance_gradient),
  sprintf("profile identity 的最大误差为 %.8g，cross-bandwidth covariance identity 的最大误差为 %.8g。总体 inverse defect 的最大值为 %.8g。",
          manifest$max_profile_identity_error, manifest$max_covariance_identity_error, manifest$max_inverse_defect),
  "",
  sprintf("五项 nuisance-penalty 审计保留各自的 %d 维研究设计，原注册维度为 %d。其逐坐标位移与 MCSE 保存在对应坐标文件，structural_audits.csv 保存审计验收数据。",
          unique(audits$audit_p), unique(audits$registered_p)),
  "",
  paste0("总体对象的真实 construction commit 为 ", paste(manifest$construction_commits, collapse = ", "),
         "。本次校验 commit 为 ", manifest$validation_commit, "。"),
  "population_assets.csv 保存全部对象路径、计算位置来源及数值诊断。manifest.json 记录验收结果。"
), file.path(out, "validation_report.md"))
writeLines("COMPLETE", file.path(out, "COMPLETE.txt"))
print(manifest)

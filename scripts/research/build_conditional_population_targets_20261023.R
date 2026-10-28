args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1L] else ".", mustWork = TRUE)
id <- if (length(args) >= 2L) args[2L] else stop("Supply one registered conditional-scale experiment ID.")
mode <- if (length(args) >= 3L) args[3L] else "both"
if (!mode %in% c("both", "profile", "sqr")) stop("Choose profile, sqr or both target types.")
.libPaths(c(file.path(root, "results", "scratch", "score_validity_20261026", "library"), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2")) source_v2_module(root, module, environment())
options(jpDMEQR.source_root = root)
specification <- yaml::read_yaml(file.path(root, "config", "joint_profile_submission_20261027.yml"))
registry_path <- file.path(root, specification$generating_registry)
registry <- utils::read.csv(registry_path, stringsAsFactors = FALSE)
row <- registry[registry$experiment_id == id, , drop = FALSE]
if (nrow(row) != 1L) stop("The declared experiment is absent or repeated.")
config <- config_row_to_list_v2(row)
override <- specification$target_mode_overrides[[id]]
if (!is.null(override)) config$target_mode <- override
conditional_scale_target_conditions_v2(config, marginal = mode != "profile")
commit <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
paths <- c(profile = profile_pair_conditional_scale_target_asset_path_v2(root, id),
  sqr = sqr_marginal_scale_target_asset_path_v2(root, id))
selected <- if (mode == "both") names(paths) else mode
if (any(file.exists(paths[selected]))) stop("A selected immutable target object already exists.")
checkpoint_directory <- file.path(root, "results", "intermediate", "profile_pair_20261022", "construction_v4",
  "checkpoints", paste0(id, "_", mode, "_", format(Sys.time(), "%Y%m%dT%H%M%OS6")))
if (dir.exists(checkpoint_directory)) stop("The conditional target checkpoint directory already exists.")
dir.create(checkpoint_directory, recursive = TRUE)
checkpoint <- function(label, value) {
  path <- file.path(checkpoint_directory, paste0(label, ".rds"))
  if (file.exists(path)) stop("A completed conditional target checkpoint already exists.")
  saveRDS(value, path)
}
for (type in selected) {
  path <- paths[[type]]
  if (type == "profile") {
    asset <- build_profile_pair_conditional_scale_target_v2(config,
      n_population = 100000L, repeats = 4L, implementation_commit = commit, checkpoint = checkpoint)
    validation <- validate_profile_pair_conditional_scale_target_v2(asset, config, final = TRUE, expected_commit = commit)
  } else {
    asset <- build_sqr_marginal_scale_target_v2(config,
      n_population = 100000L, repeats = 4L, implementation_commit = commit, checkpoint = checkpoint)
    validation <- validate_sqr_marginal_scale_target_v2(asset, config, final = TRUE, expected_commit = commit)
  }
  if (!isTRUE(validation$valid)) stop("The actual conditional target validator did not pass.")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(asset, path)
  table <- data.frame(experiment_id = id, target_type = asset$asset_kind, coordinate = asset$target_names,
    target = as.numeric(asset$beta_star_mc), target_mcse = as.numeric(asset$target_mc_se_by_coordinate),
    n_population = asset$n_population, independent_repeats = asset$repeats, stringsAsFactors = FALSE)
  utils::write.csv(table, sub("\\.rds$", ".csv", path), row.names = FALSE)
  jsonlite::write_json(list(asset_kind = asset$asset_kind, population_engine = asset$population_engine,
    construction_version = asset$construction_version, method_version = asset$method_version,
    implementation_commit = commit, source_paths = asset$source_paths, registry_path = registry_path,
    full_registered_p = config$p, complete_signal_support = seq_len(config$s),
    n_population = asset$n_population, independent_repeats = asset$repeats, seeds = asset$seeds,
    passed = validation$valid, actual_object_path = path, checkpoint_directory = checkpoint_directory,
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    sub("\\.rds$", ".json", path), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  cat(type, path, "\n")
  flush.console()
}

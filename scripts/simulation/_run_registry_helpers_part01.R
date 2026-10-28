# Shared CLI, provenance, screening, target-loading, and mapping helpers.

parse_cli_args_v2 <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- pair[1]
    value <- if (length(pair) > 1L) paste(pair[-1], collapse = "=") else "true"
    out[[key]] <- value
  }
  out
}

as_bool_cli_v2 <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}

current_commit_v2 <- function(root = ".") {
  ans <- tryCatch(
    system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(ans)) return(trimws(ans[1]))
  # Git-free fallback (offline compute carrier): the kit stamps the exact
  # frozen commit in implementation_commit.txt at build time.
  stamp <- file.path(root, "implementation_commit.txt")
  if (file.exists(stamp)) {
    s <- trimws(readLines(stamp, warn = FALSE)[1])
    if (nzchar(s)) return(s)
  }
  "unknown"
}

pid_running_v2 <- function(pid) {
  if (!is.finite(pid)) return(FALSE)
  if (identical(as.integer(pid), as.integer(Sys.getpid()))) return(TRUE)
  if (!requireNamespace("ps", quietly = TRUE)) {
    stop("The ps package is required to check whether a shard lock is active.")
  }
  as.integer(pid) %in% ps::ps_pids()
}

current_branch_v2 <- function(root = ".") {
  ans <- tryCatch(
    system2("git", c("-C", shQuote(root), "rev-parse", "--abbrev-ref", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(ans)) return(trimws(ans[1]))
  stamp <- file.path(root, "implementation_branch.txt")
  if (file.exists(stamp)) {
    s <- trimws(readLines(stamp, warn = FALSE)[1])
    if (nzchar(s)) return(s)
  }
  "unknown"
}

file_checksum_v2 <- function(path, algo = "sha256") {
  if (!file.exists(path)) stop("Cannot checksum missing file: ", path)
  if (!requireNamespace("digest", quietly = TRUE)) stop("Package 'digest' is required.")
  if (tolower(tools::file_ext(path)) %in% c("yml", "yaml")) {
    specification <- yaml::read_yaml(path)
    if (identical(specification$registry_kind, "profile_pair_registry")) {
      source_path <- file.path(dirname(path), specification$registry_source)
      if (!file.exists(source_path)) stop("The referenced frozen registry is missing.")
      return(digest::digest(list(
        specification = digest::digest(file = path, algo = algo),
        registry = digest::digest(file = source_path, algo = algo)
      ), algo = algo))
    }
  }
  digest::digest(file = path, algo = algo)
}

object_checksum_v2 <- function(object, algo = "sha256") {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Package 'digest' is required.")
  digest::digest(object, algo = algo)
}

write_json_v2 <- function(object, path, pretty = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required.")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(object, path, pretty = pretty, auto_unbox = TRUE,
                       null = "null", digits = NA)
  invisible(path)
}

read_json_v2 <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' is required.")
  if (!file.exists(path)) stop("Missing JSON file: ", path)
  jsonlite::read_json(path, simplifyVector = TRUE)
}

# 输出模式以 docs/RESULTS_CONTRACT.md 的对应节为准：恢复分片时用它核对
# 既有对象与当前记录层是否使用同一套列。
contract_schema_columns_v2 <- function(root, section) {
  path <- file.path(root, "docs", "RESULTS_CONTRACT.md")
  if (!file.exists(path)) stop("Missing ", path)
  lines <- readLines(path, warn = FALSE)
  start <- grep(paste0("^## ", as.integer(section), "\\."), lines)
  if (length(start) != 1L) {
    stop("Cannot locate section ", section, " in docs/RESULTS_CONTRACT.md.")
  }
  fences <- grep("^```", lines)
  fences <- fences[fences > start]
  if (length(fences) < 2L) {
    stop("The schema block of section ", section, " is not fenced.")
  }
  block <- trimws(lines[(fences[1] + 1L):(fences[2] - 1L)])
  block <- sub("\\s+#.*$", "", block)
  block <- sub("\\s{2,}.*$", "", block)
  block[nzchar(block)]
}

contract_replication_columns_v2 <- function(root) contract_schema_columns_v2(root, 2L)

hardware_manifest_v2 <- function() {
  list(
    sysname = Sys.info()[["sysname"]] %||% NA_character_,
    release = Sys.info()[["release"]] %||% NA_character_,
    machine = Sys.info()[["machine"]] %||% NA_character_,
    nodename = Sys.info()[["nodename"]] %||% NA_character_,
    logical_cores = parallel::detectCores(logical = TRUE),
    physical_cores = parallel::detectCores(logical = FALSE),
    r_version = R.version.string,
    platform = R.version$platform
  )
}

prepare_analysis_data_v2 <- function(dat, config, seed) {
  dat <- fit_working_nuisance_design_v2(dat, config$fit_random_effects)
  # For profile_mc scenarios the inferential target is defined on the frozen
  # target submodel (active coordinates plus designated null controls), not
  # on the full generated design. Restrict the analysis design to those
  # columns so that the profile-MC asset covers every fitted column, as
  # required by docs/SIMULATION_FREEZE_DECISIONS.md (R6). The DGP still
  # generates the full design; only the fitted/inferred columns are the
  # submodel.
  if (identical(config$target_mode, "profile_mc") &&
      !identical(config$method_version, "profile_pair_20261022")) {
    target_columns <- population_target_columns_v2(config)
    if (max(target_columns) > ncol(dat$X)) {
      stop("Target submodel columns exceed the generated design.")
    }
    dat <- subset_features_v2(dat, target_columns)
  }
  full_names <- colnames(dat$X)
  full_active <- dat$active
  target_names_global <- resolve_target_coordinates_v2(config, dat)
  target_idx_global <- match(target_names_global, full_names)
  target_idx_global <- target_idx_global[!is.na(target_idx_global)]

  selected_global <- seq_len(ncol(dat$X))
  screen_clusters <- character()
  fit_clusters <- unique(dat$cluster_id)
  scores <- NULL
  forced_target_names <- character()
  screen_engine <- "none"
  screen_method <- config$screen_method %||% "none"

  if (screen_method == "split_quantile_score") {
    split <- split_clusters_v2(
      dat$cluster_id,
      fraction = config$screen_fraction %||% 0.5,
      seed = seed + 17L
    )
    if (length(intersect(split$a, split$b))) {
      stop("The screen set and the fit set share clusters.")
    }
    screen_dat <- subset_clusters_v2(dat, split$a)
    d <- min(as.integer(config$screen_dim), ncol(screen_dat$X))
    scr <- quantile_score_screen_v2(
      screen_dat$y, screen_dat$X, screen_dat$cluster_id,
      config$tau, d
    )
    selected_global <- scr$selected
    scores <- scr$scores
    screen_engine <- scr$method
    screen_clusters <- split$a
    fit_clusters <- split$b
    dat <- subset_clusters_v2(dat, fit_clusters)
  } else if (screen_method == "variance_filter") {
    d <- min(as.integer(config$screen_dim), ncol(dat$X))
    scr <- cluster_weighted_variance_screen_v2(dat$X, dat$cluster_id, d)
    selected_global <- scr$selected
    scores <- scr$scores
    screen_engine <- scr$method
  } else if (screen_method == "ca_iqr_sis") {
    d <- min(as.integer(config$screen_dim), ncol(dat$X))
    scr <- ca_iqr_sis_screen_v2(
      dat$y, dat$X, dat$cluster_id, config$tau, d,
      quantiles = sort(unique(c(0.25, 0.5, 0.75, config$tau))),
      rounds = 2L
    )
    selected_global <- scr$selected
    scores <- scr$score_history
    screen_engine <- paste(scr$fit_engine, collapse = "|")
  } else if (!screen_method %in% c("none", "")) {
    stop("Unknown screen_method: ", screen_method)
  }

  # 同一样本内的筛选（结果盲的方差筛选与 CA-IQR-SIS）在全部聚类上进行；
  # 记录中把这些聚类写成筛选集合，使重合数量如实反映同一样本设计。
  if (screen_method %in% c("variance_filter", "ca_iqr_sis")) {
    screen_clusters <- fit_clusters
  }

  if (isTRUE(config$force_target_inclusion) && screen_method != "none") {
    missing_targets <- setdiff(target_idx_global, selected_global)
    if (length(missing_targets)) {
      forced_target_names <- full_names[missing_targets]
      selected_global <- sort(unique(c(selected_global, missing_targets)))
    }
  }
  selected_global <- sort(unique(selected_global))
  selected_names <- full_names[selected_global]
  sure_inclusion <- all(full_active %in% selected_global)
  dat <- subset_features_v2(dat, selected_global)

  screening_info <- list(
    method = screen_method,
    engine = screen_engine,
    selected_global = selected_global,
    selected_names = selected_names,
    scores = scores,
    screen_clusters = screen_clusters,
    fit_clusters = fit_clusters,
    overlap_count = length(intersect(screen_clusters, fit_clusters)),
    forced_target_names = forced_target_names,
    target_names_global = target_names_global,
    screen_dimension_requested = config$screen_dim,
    screen_dimension_realised = length(selected_global),
    sure_inclusion = sure_inclusion,
    active_global = full_active,
    full_feature_names = full_names
  )
  list(data = dat, screening = screening_info)
}

load_profile_target_v2 <- function(root, config, dat, final = FALSE) {
  if (identical(config$method_version, "joint_profile_design_20261026")) {
    original <- dat
    original$beta0 <- make_beta_v2(config$p, config$s, config$signal)
    target <- submission_target(config, original)
    return(list(beta = target$beta, target_type = target$type,
                target_mc_se = max(target$mcse),
                target_mc_se_by_coordinate = target$mcse, asset_path = target$path))
  }
  if (identical(config$target_mode, "structural")) {
    return(list(
      beta = dat$beta0,
      target_type = "structural",
      target_mc_se = 0,
      target_mc_se_by_coordinate = setNames(rep(0, length(dat$beta0)), names(dat$beta0)),
      smoothing_target_distance = NA_real_,
      smoothing_target_distance_l1 = NA_real_,
      asset_path = NA_character_
    ))
  }
  if (identical(config$target_mode, "profile_mc")) {
    paired <- identical(config$method_version, "profile_pair_20261022")
    path <- if (paired) profile_pair_target_asset_path_v2(root, config$experiment_id) else
      target_asset_path_v2(root, config$experiment_id)
    if (!file.exists(path)) stop("Profile-MC target is required but missing: ", path)
    obj <- readRDS(path)
    valid <- if (paired) validate_profile_pair_target_asset_v2(
      obj, config, final = final,
      expected_commit = if (isTRUE(final)) current_commit_v2(root) else NULL
    ) else validate_profile_target_asset_v2(
      obj,
      min_population = 100000L,
      min_repeats = 4L,
      expected_commit = if (isTRUE(final)) current_commit_v2(root) else NULL,
      final = final
    )
    if (!valid$valid) stop("Invalid profile target asset: ", paste(valid$problems, collapse = "; "))
    beta <- obj$beta_star_mc
    if (!all(names(dat$beta0) %in% names(beta))) {
      stop("Profile-MC target does not cover the analysis design.")
    }
    se_by <- obj$target_mc_se_by_coordinate %||%
      setNames(rep(obj$target_mc_se %||% NA_real_, length(beta)), names(beta))
    return(list(
      beta = beta[names(dat$beta0)],
      target_type = "profile_mc",
      target_mc_se = max(se_by[names(dat$beta0)], na.rm = TRUE),
      target_mc_se_by_coordinate = se_by[names(dat$beta0)],
      smoothing_target_distance = obj$smoothing_target_distance %||% NA_real_,
      smoothing_target_distance_l1 = obj$smoothing_target_distance_l1 %||% NA_real_,
      asset_path = path,
      asset = obj
    ))
  }
  stop("Unknown target_mode: ", config$target_mode)
}

load_population_direction_asset_v2 <- function(root, config, required = FALSE,
                                               final = FALSE) {
  if (identical(config$method_version, "joint_profile_design_20261026")) {
    if ("JOINT-PROFILE-POP-H" %in% config$methods) return(submission_direction(config))
    return(NULL)
  }
  paired <- identical(config$method_version, "profile_pair_20261022")
  path <- if (paired) profile_pair_direction_asset_path_v2(root, config$experiment_id) else
    population_direction_asset_path_v2(root, config$experiment_id)
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Population-direction asset is missing: ", path)
    return(NULL)
  }
  obj <- readRDS(path)
  if (paired) {
    validate_profile_pair_direction_asset_v2(
      obj, config, final = final,
      expected_commit = if (isTRUE(final)) current_commit_v2(root) else NULL
    )
    return(obj)
  }
  if (!identical(obj$experiment_id, config$experiment_id)) {
    stop("Population-direction asset experiment_id mismatch.")
  }
  if (isTRUE(final) && !identical(obj$implementation_commit, current_commit_v2(root))) {
    stop("Population-direction asset implementation commit is stale.")
  }
  obj
}

common_lqmm_subset_v2 <- function(dat, target_names, tau, max_features = 20L,
                                  active = NULL, seed = 20261022L,
                                  screen_fraction = 0.5) {
  roles <- NULL
  fit_data <- dat
  if (!is.null(active)) {
    keep <- sort(unique(c(as.integer(active), match(target_names, colnames(dat$X)))))
    keep <- keep[!is.na(keep)]
  } else {
    if (!is.finite(screen_fraction) || screen_fraction <= 0 || screen_fraction >= 1) {
      stop("Low-dimensional screening requires a fraction strictly within (0,1).")
    }
    subjects <- unique(dat$cluster_id)
    set.seed(seed)
    screen_count <- max(1L, min(length(subjects) - 1L,
                                floor(length(subjects) * screen_fraction)))
    screen_subjects <- sample(subjects, screen_count)
    fit_subjects <- setdiff(subjects, screen_subjects)
    screen_data <- subset_clusters_v2(dat, screen_subjects)
    fit_data <- subset_clusters_v2(dat, fit_subjects)
    scr <- quantile_score_screen_v2(
      screen_data$y, screen_data$X, screen_data$cluster_id, tau,
      d = min(max_features, ncol(dat$X))
    )
    roles <- list(screen_clusters = as.character(screen_subjects),
                   fit_clusters = as.character(fit_subjects), overlap_count = 0L,
                   screen_fraction = screen_fraction, seed = seed)
    target_idx <- match(target_names, colnames(dat$X))
    target_idx <- target_idx[!is.na(target_idx)]
    keep <- unique(c(target_idx, scr$selected))
  }
  if (length(keep) > max_features) {
    target_idx <- match(target_names, colnames(dat$X))
    target_idx <- target_idx[!is.na(target_idx)]
    non_target <- setdiff(keep, target_idx)
    keep <- unique(c(target_idx, head(non_target, max_features - length(target_idx))))
  }
  list(data = subset_features_v2(fit_data, keep), keep = keep, roles = roles)
}

map_beta_to_full_v2 <- function(beta, full_names) {
  out <- setNames(rep(0, length(full_names)), full_names)
  if (!is.null(beta)) {
    if (is.null(names(beta))) stop("Benchmark beta_hat must be named.")
    common <- intersect(names(beta), full_names)
    out[common] <- beta[common]
  }
  out
}

apply_feature_map_to_test_v2 <- function(test_dat, screening, fit_random_effects) {
  test_dat <- fit_working_nuisance_design_v2(test_dat, fit_random_effects)
  selected <- screening$selected_global
  subset_features_v2(test_dat, selected)
}

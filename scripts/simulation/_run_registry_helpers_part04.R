# Registry validation, evidence gates, and freeze-manifest helpers.

required_governance_files_v2 <- function() {
  c(
    "AGENTS.md",
    "docs/METHOD_SPECIFICATION.md",
    "docs/SIMULATION_FREEZE_DECISIONS.md",
    "docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md",
    "docs/EMPIRICAL_FREEZE_DECISIONS.md",
    "config/benchmark_requirements.csv",
    "docs/RESULTS_CONTRACT.md",
    "docs/SIMULATION_PROTOCOL.md"
  )
}

validate_governance_files_v2 <- function(root) {
  paths <- file.path(root, required_governance_files_v2())
  missing <- paths[!file.exists(paths)]
  list(
    pass = !length(missing),
    missing = sub(paste0("^", normalizePath(root, mustWork = TRUE), "/?"), "", missing)
  )
}

read_benchmark_requirements_v2 <- function(root) {
  path <- file.path(root, "config", "benchmark_requirements.csv")
  req <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("method_id", "registry_policy", "final_gate", "allowed_metrics",
              "implementation_state", "acceptance_test", "reference_identifier")
  missing <- setdiff(needed, names(req))
  if (length(missing)) stop("benchmark_requirements.csv missing: ", paste(missing, collapse = ", "))
  if (anyDuplicated(req$method_id)) stop("Duplicate method_id in benchmark_requirements.csv.")
  req$final_gate <- tolower(as.character(req$final_gate)) %in% c("true", "1", "yes")
  req
}

registry_methods_v2 <- function(cfg_df) {
  sort(unique(unlist(lapply(cfg_df$methods, parse_pipe_character_v2), use.names = FALSE)))
}

registry_type_for_path_v2 <- function(config_path) {
  if (tolower(tools::file_ext(config_path)) %in% c("yml", "yaml")) {
    registry <- read_simulation_registry_v2(config_path)
    phase <- attr(registry, "registry_phase")
    if (is.null(phase)) stop("The registry specification is missing its phase.")
    return(phase)
  }
  # 注册表类型决定校验强度。两份主注册表（原始规模与 B=200 修订）都按主注册表
  # 校验，文件名前缀决定类型，无法判定的文件名直接中止而不静默降级。
  name <- basename(config_path)
  if (grepl("^simulation_main", name)) return("main")
  if (grepl("^simulation_pilot", name)) return("pilot")
  if (grepl("^simulation_preflight", name)) return("preflight")
  stop("Unknown registry file name: ", name)
}

validate_registry_contract_v2 <- function(root, config_path,
                                          registry_type = c("main", "pilot", "preflight")) {
  registry_type <- match.arg(registry_type)
  cfg <- read_simulation_registry_v2(config_path)
  # Normalise column types from the CSV before any numeric/logical comparison.
  numeric_cols <- intersect(c(
    "n_clusters", "p", "s", "tau", "q", "rho_x", "signal", "copula_rho",
    "sigma_b0", "sigma_b1", "rho_b", "x_b_corr", "nonlinear_re_strength",
    "lambda_gamma", "h_multiplier", "screen_dim", "screen_fraction",
    "target_coordinate_count", "workers", "replications"
  ), names(cfg))
  for (nm in numeric_cols) {
    cfg[[nm]] <- suppressWarnings(as.numeric(as.character(cfg[[nm]])))
  }
  bool_cols <- intersect(c("heteroskedastic", "informative_size",
                           "force_target_inclusion"), names(cfg))
  for (nm in bool_cols) {
    cfg[[nm]] <- tolower(as.character(cfg[[nm]])) %in% c("true", "1", "yes", "y")
  }
  problems <- character()
  warnings <- character()
  add_problem <- function(x) problems <<- c(problems, x)
  add_warning <- function(x) warnings <<- c(warnings, x)

  if (registry_type == "main" && nrow(cfg) != 78L) {
    add_problem(paste0("Main registry must contain 78 experiment rows; found ", nrow(cfg), "."))
  }
  if (registry_type == "pilot") {
    ids <- cfg$experiment_id
    if (!all(sprintf("P%02d", 1:6) %in% ids) ||
        any(!ids %in% sprintf("P%02d", 1:7))) {
      add_problem("Pilot registry must contain P01--P06 (P07 allowed as scaling diagnostic).")
    }
  }
  if (any(!cfg$target_mode %in% c("structural", "profile_mc"))) {
    add_problem("Unknown target_mode in registry.")
  }
  if (any(!cfg$screen_method %in% c("none", "variance_filter", "ca_iqr_sis",
                                      "split_quantile_score"))) {
    add_problem("Unknown screen_method in registry.")
  }
  numeric_positive <- c("n_clusters", "p", "s", "q", "h_multiplier", "replications")
  for (nm in numeric_positive) {
    if (any(!is.finite(cfg[[nm]]) | cfg[[nm]] <= 0)) add_problem(paste0(nm, " must be positive."))
  }
  if (any(cfg$s > cfg$p)) add_problem("A registry row has s > p.")
  if (any(cfg$tau <= 0 | cfg$tau >= 1)) add_problem("All tau values must lie in (0,1).")

  req <- read_benchmark_requirements_v2(root)
  methods <- registry_methods_v2(cfg)
  unknown <- setdiff(methods, req$method_id)
  if (length(unknown)) add_problem(paste("Registry contains unknown methods:", paste(unknown, collapse = ", ")))
  excluded <- intersect(methods, req$method_id[req$registry_policy == "excluded"])
  if (length(excluded)) add_problem(paste("Registry contains excluded methods:", paste(excluded, collapse = ", ")))

  if (registry_type == "main") {
    A <- cfg[cfg$module == "A_theory_scaling", , drop = FALSE]
    for (p in c(500, 1000, 2000)) {
      got <- sort(unique(A$n_clusters[A$p == p & A$s == 5]))
      if (!identical(got, c(100, 200, 400))) {
        add_problem(paste0("Module A lacks complete n sequence for p=", p, ", s=5."))
      }
    }
    got_s10 <- sort(unique(A$n_clusters[A$p == 1000 & A$s == 10]))
    if (!identical(got_s10, c(100, 200, 400))) {
      add_problem("Module A lacks complete n sequence for p=1000, s=10.")
    }

    B <- cfg[cfg$module == "B_error_robustness", , drop = FALSE]
    family_key <- ifelse(B$heteroskedastic,
                         paste0("heteroskedastic_", B$error_dist), B$error_dist)
    required_families <- c(
      "gaussian", "t3", "laplace", "skew_chisq3", "contaminated_normal",
      "asymmetric_laplace", "heteroskedastic_gaussian", "heteroskedastic_t3"
    )
    for (fam in required_families) {
      if (sum(family_key == fam) < 2L) add_problem(paste0("Error family appears fewer than twice: ", fam))
    }
    for (fam in c("gaussian", "t3", "contaminated_normal")) {
      if (!all(c(0.25, 0.5, 0.75) %in% B$tau[family_key == fam])) {
        add_problem(paste0("Core error family lacks all three quantiles: ", fam))
      }
    }

    if (!all(c(0.4, 1.1) %in% cfg$signal[cfg$module == "H_signal_strength"])) {
      add_problem("Weak/strong signal rows are incomplete.")
    }
    if (!any(cfg$error_dependence == "gaussian_copula_ar1")) {
      add_problem("Gaussian-copula dependence sensitivity is missing.")
    }
    if (!any(cfg$m_rule == "geometric_imbalance")) {
      add_problem("Truncated-geometric imbalance sensitivity is missing.")
    }
    if (!any(cfg$x_b_corr != 0 & cfg$target_mode == "profile_mc")) {
      add_problem("X--random-effect correlation profile-target row is missing.")
    }
    if (!any(cfg$informative_size & cfg$target_mode == "profile_mc")) {
      add_problem("Informative cluster-size profile-target row is missing.")
    }
    if (!any(cfg$nonlinear_re_strength != 0 & cfg$target_mode == "profile_mc")) {
      add_problem("Nonlinear random-effect profile-target row is missing.")
    }
    screen_set <- unique(cfg$screen_method[cfg$module == "F_screening"])
    if (!all(c("none", "variance_filter", "ca_iqr_sis", "split_quantile_score") %in% screen_set)) {
      add_problem("Module F screening routes are incomplete.")
    }
    if (!all(c(0.25, 0.5, 2, 4) %in% cfg$lambda_gamma[cfg$module == "E_tuning"])) {
      add_problem("Module E nuisance-penalty sensitivity grid is incomplete.")
    }
    G <- cfg[cfg$module == "G_computation", , drop = FALSE]
    if (max(G$target_coordinate_count, na.rm = TRUE) <= min(G$target_coordinate_count, na.rm = TRUE)) {
      add_problem("Module G does not vary target-coordinate count.")
    }
    if (length(unique(G$workers)) < 2L) add_problem("Module G does not vary worker count.")

    E04 <- cfg[cfg$experiment_id == "E06", , drop = FALSE]
    if (nrow(E04) != 1L || E04$lambda_gamma != 4 || E04$target_mode != "structural") {
      add_problem("Frozen strong-ridge structural-target row is missing or mislabeled.")
    }
  }

  if (any(cfg$target_mode == "profile_mc" & cfg$p > 500)) {
    add_warning("One or more profile_mc rows have p>500; population target cost may be prohibitive.")
  }
  list(
    pass = length(problems) == 0L,
    problems = unique(problems),
    warnings = unique(warnings),
    row_count = nrow(cfg),
    experiment_ids = cfg$experiment_id,
    methods = methods,
    config_sha256 = file_checksum_v2(config_path),
    registry_type = registry_type
  )
}

benchmark_acceptance_manifest_path_v2 <- function(root, method_id) {
  file.path(root, "results", "preflight", "benchmarks", method_id, "acceptance.json")
}

validate_benchmark_gate_v2 <- function(root, cfg_df, final = TRUE) {
  req <- read_benchmark_requirements_v2(root)
  methods <- registry_methods_v2(cfg_df)
  rows <- vector("list", length(methods))
  for (ii in seq_along(methods)) {
    method <- methods[ii]
    rr <- req[req$method_id == method, , drop = FALSE]
    if (!nrow(rr)) {
      rows[[ii]] <- data.frame(method_id = method, pass = FALSE,
                               reason = "method missing from benchmark requirements")
      next
    }
    if (rr$registry_policy == "excluded") {
      rows[[ii]] <- data.frame(method_id = method, pass = FALSE,
                               reason = "method is excluded from final registry")
      next
    }
    manifest_path <- benchmark_acceptance_manifest_path_v2(root, method)
    if (!file.exists(manifest_path)) {
      required <- isTRUE(final) && rr$final_gate
      rows[[ii]] <- data.frame(
        method_id = method,
        pass = !required,
        reason = if (required) "acceptance manifest missing" else "optional manifest missing"
      )
      next
    }
    manifest <- tryCatch(read_json_v2(manifest_path), error = function(e) NULL)
    if (is.null(manifest)) {
      rows[[ii]] <- data.frame(method_id = method, pass = FALSE,
                               reason = "acceptance manifest is unreadable")
      next
    }
    checks <- c(
      isTRUE(manifest$unit_test_pass),
      isTRUE(manifest$limiting_case_pass),
      isTRUE(manifest$fidelity_check_pass),
      isTRUE(manifest$schema_check_pass),
      isTRUE(manifest$deterministic_seed_pass),
      identical(as.character(manifest$commit_sha), current_commit_v2(root))
    )
    rows[[ii]] <- data.frame(
      method_id = method,
      pass = all(checks),
      reason = if (all(checks)) "accepted" else "one or more acceptance checks failed or commit is stale"
    )
  }
  tab <- do.call(rbind, rows)
  list(pass = all(tab$pass), table = tab)
}

profile_pair_gate_directory_v2 <- function(root) {
  file.path(root, "results", "preflight", "profile_pair_20261022")
}

registry_method_version_v2 <- function(registry) {
  if (!"method_version" %in% names(registry)) return(NULL)
  versions <- unique(as.character(registry$method_version))
  if (length(versions) != 1L || is.na(versions)) stop("The registry requires one method version.")
  versions
}

geometry_manifest_path_v2 <- function(root, method_version = NULL) {
  if (identical(method_version, "profile_pair_20261022")) {
    return(file.path(profile_pair_gate_directory_v2(root), "geometry_manifest.json"))
  }
  file.path(root, "results", "preflight", "geometry", "manifest.json")
}

# 几何审计身份文件里记录的求解路径必须与当前配置相同：一次在默认配置下取得的
# 通过结论不能为另一次采用加速开关的运行担保。
validate_geometry_record_v2 <- function(x, expected_commit,
                                         current_configuration = solver_configuration_v2(),
                                         final = TRUE, method_version = NULL) {
  same_configuration <- configuration_matches_v2(x$solver_configuration, current_configuration)
  seed_count <- x$seed_count %||% 0
  checks <- c(
    isTRUE(x$pass),
    isTRUE(x$strict),
    isTRUE(length(seed_count) == 1L && is.numeric(seed_count) && is.finite(seed_count) &&
             seed_count >= if (final) 20L else 1L),
    identical(as.character(x$commit_sha), as.character(expected_commit)),
    isTRUE(x$dantzig_checked),
    isTRUE(x$dantzig_pass),
    isTRUE(x$assembly_modes_checked),
    isTRUE(x$assembly_pass),
    is.null(method_version) || identical(x$method_version, method_version),
    same_configuration
  )
  reason <- if (all(checks)) {
    "accepted"
  } else if (!same_configuration) {
    "geometry evidence was measured under a different solver configuration"
  } else {
    "geometry evidence failed or is stale"
  }
  list(pass = all(checks), reason = reason, manifest = x)
}

validate_geometry_gate_v2 <- function(root, final = TRUE, method_version = NULL) {
  path <- geometry_manifest_path_v2(root, method_version)
  if (!file.exists(path)) return(list(pass = FALSE, reason = "strict geometry manifest missing"))
  validate_geometry_record_v2(read_json_v2(path), current_commit_v2(root),
                              solver_configuration_v2(), final, method_version)
}

structural_target_audit_asset_path_v2 <- function(root, experiment_id) {
  file.path(root, "results", "intermediate", "target_audits",
            paste0(experiment_id, ".rds"))
}

validate_structural_target_audit_v2 <- function(obj, expected_commit = NULL,
                                                final = TRUE,
                                                expected_dependency_hash = NULL) {
  problems <- character()
  if (is.null(obj$target_displacement) || any(!is.finite(obj$target_displacement))) {
    problems <- c(problems, "target displacement is missing or non-finite")
  }
  if (isTRUE(final) && (is.null(obj$n_population) || min(obj$n_population) < 100000L)) {
    problems <- c(problems, "population sample is below 100000")
  }
  if (isTRUE(final) && (is.null(obj$repeats) || obj$repeats < 4L)) {
    problems <- c(problems, "repeat count is below four")
  }
  # Convergence acceptance follows the frozen thresholds in
  # docs/METHOD_SPECIFICATION.md (Section 11): KKT residual and nuisance
  # gradient within the acceptance tolerances. The solver-internal
  # convergence flags are stricter (1e-8) and can sit on a numerical plateau
  # that still satisfies the frozen acceptance threshold; they are kept as
  # diagnostics only (same rule as the strict geometry gate).
  if (!isTRUE(obj$all_converged) &&
      !(is.numeric(obj$max_kkt_residual) &&
        max(obj$max_kkt_residual, na.rm = TRUE) <= 1e-5 &&
        is.numeric(obj$max_nuisance_gradient) &&
        max(obj$max_nuisance_gradient, na.rm = TRUE) <= 1e-7)) {
    problems <- c(problems, "one or more audit fits did not converge")
  }
  if (!is.null(expected_dependency_hash)) {
    if (!identical(obj$dependency_hash %||% NULL, expected_dependency_hash)) {
      problems <- c(problems, "structural audit dependency hash is stale")
    }
  } else if (!is.null(expected_commit) &&
             !identical(obj$implementation_commit, expected_commit)) {
    problems <- c(problems, "structural audit implementation commit is stale")
  }
  if (!is.null(obj$target_displacement) && !is.null(obj$target_mc_se_by_coordinate)) {
    displacement <- abs(as.numeric(obj$target_displacement))
    mcse <- as.numeric(obj$target_mc_se_by_coordinate[names(obj$target_displacement)])
    tolerance <- pmax(0.02, 5 * mcse)
    if (any(!is.finite(tolerance)) || any(displacement > tolerance)) {
      problems <- c(problems,
                    "population displacement is not compatible with the frozen structural target tolerance")
    }
  }
  list(valid = length(problems) == 0L, problems = unique(problems))
}

validate_target_gate_v2 <- function(root, cfg_df, final = TRUE, expected_config_sha = NULL) {
  if ("method_version" %in% names(cfg_df) && any(cfg_df$method_version %in% "profile_pair_20261022")) {
    return(validate_profile_pair_target_gate_v2(root, cfg_df, final))
  }
  profile_ids <- cfg_df$experiment_id[cfg_df$target_mode == "profile_mc"]
  pop_h_ids <- cfg_df$experiment_id[vapply(cfg_df$methods, function(x) {
    "PROFILE-DQR-POP-H" %in% parse_pipe_character_v2(x)
  }, logical(1))]
  structural_audit_ids <- cfg_df$experiment_id[
    cfg_df$module == "E_tuning" & cfg_df$target_mode == "structural" &
      cfg_df$lambda_gamma %in% c(0.25, 0.5, 2, 4)
  ]
  rows <- list(); kk <- 0L
  cfg_lookup <- split(cfg_df, cfg_df$experiment_id)
  for (id in profile_ids) {
    kk <- kk + 1L
    path <- target_asset_path_v2(root, id)
    if (!file.exists(path)) {
      rows[[kk]] <- data.frame(experiment_id = id, asset = "profile_target",
                               pass = FALSE, reason = "missing")
    } else {
      obj <- readRDS(path)
      dep_hash <- target_dependency_hash_v2(
        config_row_to_list_v2(cfg_lookup[[id]]),
        n_population = 100000L, repeats = 4L
      )
      val <- validate_profile_target_asset_v2(
        obj, min_population = 100000L, min_repeats = 4L,
        expected_commit = if (final) current_commit_v2(root) else NULL,
        expected_dependency_hash = if (final) dep_hash else NULL,
        final = final
      )
      # Reused assets (authorised by the dependency hash) may carry the
      # pre-change full-registry checksum; the dependency hash is then the
      # authoritative identity and the config checksum comparison is skipped.
      checksum_ok <- is.null(expected_config_sha) ||
        (!final && identical(as.character(obj$config_sha256), as.character(expected_config_sha))) ||
        (isTRUE(final) && !is.null(dep_hash) &&
           identical(obj$dependency_hash %||% NULL, dep_hash))
      pass <- val$valid && checksum_ok
      reason <- c(val$problems, if (!checksum_ok) "config checksum mismatch")
      rows[[kk]] <- data.frame(experiment_id = id, asset = "profile_target",
                               pass = pass,
                               reason = if (pass) "accepted" else paste(reason, collapse = "; "))
    }
  }
  for (id in unique(structural_audit_ids)) {
    kk <- kk + 1L
    path <- structural_target_audit_asset_path_v2(root, id)
    if (!file.exists(path)) {
      rows[[kk]] <- data.frame(experiment_id = id, asset = "structural_target_audit",
                               pass = FALSE, reason = "missing")
    } else {
      obj <- readRDS(path)
      dep_hash <- target_dependency_hash_v2(
        config_row_to_list_v2(cfg_lookup[[id]]),
        n_population = 100000L, repeats = 4L
      )
      val <- validate_structural_target_audit_v2(
        obj, expected_commit = if (final) current_commit_v2(root) else NULL,
        expected_dependency_hash = if (final) dep_hash else NULL,
        final = final
      )
      checksum_ok <- is.null(expected_config_sha) ||
        (!final && identical(as.character(obj$config_sha256), as.character(expected_config_sha))) ||
        (isTRUE(final) && !is.null(dep_hash) &&
           identical(obj$dependency_hash %||% NULL, dep_hash))
      pass <- val$valid && checksum_ok
      reason <- c(val$problems, if (!checksum_ok) "config checksum mismatch")
      rows[[kk]] <- data.frame(
        experiment_id = id, asset = "structural_target_audit",
        pass = pass,
        reason = if (pass) "accepted" else paste(reason, collapse = "; ")
      )
    }
  }
  for (id in unique(pop_h_ids)) {
    kk <- kk + 1L
    path <- population_direction_asset_path_v2(root, id)
    if (!file.exists(path)) {
      rows[[kk]] <- data.frame(experiment_id = id, asset = "population_direction",
                               pass = FALSE, reason = "missing")
    } else {
      obj <- readRDS(path)
      checksum_ok <- is.null(expected_config_sha) ||
        identical(as.character(obj$config_sha256), as.character(expected_config_sha))
      # Nuisance acceptance follows the frozen threshold in
      # docs/METHOD_SPECIFICATION.md (Section 11): maximum nuisance gradient
      # below 1e-7. The solver-internal convergence flag (1e-8) is stricter
      # and is retained as a diagnostic only; small-cluster fits can sit on a
      # numerical plateau between 1e-8 and 1e-7 that still satisfies the
      # frozen acceptance threshold (same rule as the strict geometry gate).
      nuisance_ok <- !is.null(obj$max_nuisance_gradient) &&
        max(obj$max_nuisance_gradient, na.rm = TRUE) <= 1e-7
      cfg_row <- config_row_to_list_v2(cfg_lookup[[id]])
      n_analysis <- as.integer(cfg_row$n_clusters)
      h_analysis <- cfg_row$h_multiplier * n_analysis^(-3 / 10)
      dep_ok <- identical(obj$dependency_hash %||% NULL,
                          pop_h_dependency_hash_v2(cfg_row, 100000L, 4L,
                                                   n_analysis, h_analysis))
      # Round-5 amendment section 12: the dependency hash (h_inf/target/DGP
      # based) is the authoritative reuse identity; a file-level config
      # checksum change (e.g. the methods column) does not invalidate the
      # POP-H asset, mirroring the profile-target branch.
      checksum_ok <- checksum_ok || dep_ok
      pass <- identical(obj$experiment_id, id) && checksum_ok &&
        (!final || (obj$n_population >= 100000L && obj$repeats >= 4L &&
                      identical(obj$implementation_commit, current_commit_v2(root)))) &&
        dep_ok &&
        isTRUE(all.equal(as.numeric(obj$h_analysis), h_analysis, tolerance = 1e-12)) &&
        nuisance_ok &&
        max(vapply(obj$directions, `[[`, numeric(1), "residual")) <= 1e-5
      rows[[kk]] <- data.frame(experiment_id = id, asset = "population_direction",
                               pass = pass,
                               reason = if (pass) "accepted" else "invalid, inaccurate, or stale")
    }
  }
  tab <- if (length(rows)) do.call(rbind, rows) else
    data.frame(experiment_id = character(), asset = character(), pass = logical(), reason = character())
  list(pass = !nrow(tab) || all(tab$pass), table = tab)
}

micro_preflight_manifest_path_v2 <- function(root, method_version = NULL) {
  if (identical(method_version, "profile_pair_20261022")) {
    return(file.path(profile_pair_gate_directory_v2(root), "micro_preflight_manifest.json"))
  }
  file.path(root, "results", "preflight", "micro_preflight_manifest.json")
}

# 手稿产物的预期集合：按登记表包含的模块推断应当写出哪些片段与图形。
# 生成条件写在构建脚本内部，这里给出与登记表内容对应的清单，构建结束时逐项核对文件是否写出。
expected_manuscript_assets_v2 <- function(registry) {
  fragments <- c("tab_simulation_main.tex", "tab_simulation_methods.tex",
                 "tab_simulation_inference.tex", "tab_simulation_theory.tex",
                 "tab_simulation_failures.tex")
  figures <- if (isTRUE(any(registry$target_coordinate_count > 0L))) {
    c("fig_type1_power.pdf", "fig_null_pvalue_uniformity.pdf",
      "fig_studentised_quantiles.pdf")
  } else character(0)
  if (any(registry$module == "A_theory_scaling")) {
    fragments <- c(fragments, "tab_simulation_rates.tex")
    figures <- c(figures, "fig_rate_slopes.pdf", "fig_coverage_clusters.pdf",
                 "fig_se_sd_ratio.pdf")
  }
  if (any(registry$module == "G_computation")) {
    fragments <- c(fragments, "tab_simulation_computation.tex")
    figures <- c(figures, "fig_computation_time.pdf")
  }
  if (any(registry$module == "F_screening")) {
    fragments <- c(fragments, "tab_simulation_selection.tex")
  }
  list(fragments = fragments, figures = figures)
}

pilot_gate_manifest_path_v2 <- function(root, method_version = NULL) {
  if (identical(method_version, "profile_pair_20261022")) {
    return(file.path(profile_pair_gate_directory_v2(root), "pilot_gate.json"))
  }
  file.path(root, "results", "preflight", "pilot_gate.json")
}

validate_named_gate_manifest_v2 <- function(root, path, label,
                                            require_final = FALSE) {
  if (!file.exists(path)) return(list(pass = FALSE, reason = paste(label, "manifest missing")))
  x <- read_json_v2(path)
  pass <- isTRUE(x$pass) && identical(as.character(x$commit_sha), current_commit_v2(root))
  if (pass && isTRUE(require_final) && !isTRUE(x$final)) {
    return(list(
      pass = FALSE,
      reason = paste(label, "manifest was produced without the final flag"),
      manifest = x
    ))
  }
  list(pass = pass, reason = if (pass) "accepted" else paste(label, "gate failed or is stale"), manifest = x)
}

freeze_manifest_path_v2 <- function(root, method_version = NULL) {
  if (identical(method_version, "profile_pair_20261022")) {
    return(file.path(profile_pair_gate_directory_v2(root), "freeze_manifest.json"))
  }
  file.path(root, "results", "preflight", "freeze_manifest.json")
}

freeze_evidence_paths_v2 <- function(root, cfg_df, config_path) {
  version <- registry_method_version_v2(cfg_df)
  if (identical(version, "profile_pair_20261022")) {
    methods <- registry_methods_v2(cfg_df)
    return(rbind(
      data.frame(
        label = c("config", "registry_source", "geometry_manifest",
                  paste0("benchmark:", methods), "micro_preflight_manifest", "pilot_gate_manifest"),
        path = c(normalizePath(config_path), attr(cfg_df, "registry_source_file"),
                 geometry_manifest_path_v2(root, version),
                 vapply(methods, function(m) benchmark_acceptance_manifest_path_v2(root, m), character(1)),
                 micro_preflight_manifest_path_v2(root, version),
                 pilot_gate_manifest_path_v2(root, version)),
        stringsAsFactors = FALSE
      ), profile_pair_freeze_asset_paths_v2(root, cfg_df)
    ))
  }
  # 冻结清单逐项记录门禁依据文件的校验和。运行时重新计算这些校验和，任何
  # 依据文件被重新生成、替换或缺失都会使既有清单失效，因此一次较早的通过
  # 结论不能为之后改动过的资产或门禁担保。
  profile_ids <- cfg_df$experiment_id[cfg_df$target_mode == "profile_mc"]
  pop_h_ids <- cfg_df$experiment_id[vapply(cfg_df$methods, function(x) {
    "PROFILE-DQR-POP-H" %in% parse_pipe_character_v2(x)
  }, logical(1))]
  audit_ids <- cfg_df$experiment_id[
    cfg_df$module == "E_tuning" & cfg_df$target_mode == "structural" &
      cfg_df$lambda_gamma %in% c(0.25, 0.5, 2, 4)
  ]
  methods <- registry_methods_v2(cfg_df)
  audit_unique <- unique(audit_ids)
  # paste0 会把零长度输入折叠成只有前缀的一个标签，sprintf 与输入长度一致，
  # 因此空分组不会在依据清单里留下没有路径的条目。
  labelled <- function(prefix, ids) {
    if (!length(ids)) return(character(0))
    sprintf("%s%s", prefix, ids)
  }
  data.frame(
    label = c(
      "config",
      "geometry_manifest",
      labelled("benchmark:", methods),
      labelled("profile_target:", profile_ids),
      labelled("population_direction:", pop_h_ids),
      labelled("structural_audit:", audit_unique),
      "micro_preflight_manifest",
      "pilot_gate_manifest"
    ),
    path = c(
      normalizePath(config_path, mustWork = TRUE),
      geometry_manifest_path_v2(root),
      vapply(methods, function(m) benchmark_acceptance_manifest_path_v2(root, m),
             character(1)),
      vapply(profile_ids, function(id) target_asset_path_v2(root, id), character(1)),
      vapply(pop_h_ids, function(id) population_direction_asset_path_v2(root, id),
             character(1)),
      vapply(audit_unique, function(id) structural_target_audit_asset_path_v2(root, id),
             character(1)),
      micro_preflight_manifest_path_v2(root),
      pilot_gate_manifest_path_v2(root)
    ),
    stringsAsFactors = FALSE
  )
}

freeze_evidence_v2 <- function(root, cfg_df, config_path) {
  tab <- freeze_evidence_paths_v2(root, cfg_df, config_path)
  tab$sha256 <- vapply(tab$path, function(p) {
    if (file.exists(p)) file_checksum_v2(p) else NA_character_
  }, character(1))
  tab
}

freeze_preflight_v2 <- function(root,
                                config_path = file.path(root, "config", "simulation_main_b200.csv"),
                                final = TRUE,
                                write_manifest = TRUE) {
  governance <- validate_governance_files_v2(root)
  registry <- validate_registry_contract_v2(root, config_path, "main")
  cfg <- read_simulation_registry_v2(config_path)
  method_version <- registry_method_version_v2(cfg)
  benchmark <- validate_benchmark_gate_v2(root, cfg, final = final)
  geometry <- validate_geometry_gate_v2(root, final = final, method_version = method_version)
  targets <- validate_target_gate_v2(
    root, cfg, final = final, expected_config_sha = file_checksum_v2(config_path)
  )
  micro <- validate_named_gate_manifest_v2(
    root, micro_preflight_manifest_path_v2(root, method_version), "micro-preflight",
    require_final = TRUE
  )
  pilot <- validate_named_gate_manifest_v2(
    root, pilot_gate_manifest_path_v2(root, method_version), "pilot", require_final = TRUE
  )
  pass <- governance$pass && registry$pass && benchmark$pass && geometry$pass &&
    targets$pass && micro$pass && pilot$pass
  manifest <- list(
    pass = pass,
    final = final,
    method_version = method_version,
    score_definition = if (identical(method_version, "profile_pair_20261022"))
      "profile_bandwidth_pair_2_half_minus_full" else NULL,
    commit_sha = current_commit_v2(root),
    branch = current_branch_v2(root),
    config_path = normalizePath(config_path, mustWork = TRUE),
    config_sha256 = file_checksum_v2(config_path),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    governance = governance,
    registry = registry,
    benchmark = list(pass = benchmark$pass, table = benchmark$table),
    geometry = geometry,
    targets = list(pass = targets$pass, table = targets$table),
    micro_preflight = micro,
    pilot = pilot,
    evidence = freeze_evidence_v2(root, cfg, config_path),
    hardware = hardware_manifest_v2(),
    solver_configuration = solver_configuration_v2()
  )
  if (isTRUE(write_manifest)) {
    path <- freeze_manifest_path_v2(root, method_version)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write_json_v2(manifest, path)
  }
  manifest
}

requires_freeze_manifest_v2 <- function(config_path) {
  # 微型预检清单是冻结清单的输入之一，微型预检本身不等待冻结清单
  identical(registry_type_for_path_v2(config_path), "main")
}

verify_freeze_manifest_v2 <- function(root, config_path) {
  method_version <- registry_method_version_v2(read_simulation_registry_v2(config_path))
  path <- freeze_manifest_path_v2(root, method_version)
  if (!file.exists(path)) stop("Freeze manifest is missing. Run 05_preflight_freeze.R.")
  x <- read_json_v2(path)
  checks <- c(
    isTRUE(x$pass),
    isTRUE(x$final),
    identical(as.character(x$commit_sha), current_commit_v2(root)),
    identical(as.character(x$config_sha256), file_checksum_v2(config_path))
  )
  if (!all(checks)) {
    stop("Freeze manifest is failed or stale for the current commit/configuration.")
  }
  # 求解路径也必须与清单一致：几何门在冻结求解路径下取得的通过结论不能为
  # 另一条求解路径下的运行担保，因此运行入口在这里再次比对六个字段。
  if (!configuration_matches_v2(x$solver_configuration)) {
    stop("Freeze manifest was created under a different solver configuration.")
  }
  recorded <- x$evidence
  if (is.null(recorded) || !nrow(recorded) || !all(c("label", "path", "sha256") %in% names(recorded))) {
    stop("Freeze manifest records no evidence fingerprints. Regenerate it with 05_preflight_freeze.R.")
  }
  verify_freeze_evidence_v2(root, config_path, recorded)
  invisible(x)
}

validate_profile_pair_final_gate_v2 <- function(root, config_path) {
  path <- freeze_manifest_path_v2(root, "profile_pair_20261022")
  if (!file.exists(path)) stop("PAIRED_FINAL_GATE_INCOMPLETE: final freeze manifest is missing.")
  manifest <- read_json_v2(path)
  if (!identical(manifest$method_version, "profile_pair_20261022") ||
      !identical(manifest$score_definition, "profile_bandwidth_pair_2_half_minus_full")) {
    stop("PAIRED_FINAL_GATE_INCOMPLETE: final checks must certify the paired score definition.")
  }
  verify_freeze_manifest_v2(root, config_path)
}

verify_freeze_evidence_v2 <- function(root, config_path, recorded) {
  current <- freeze_evidence_v2(root, read_simulation_registry_v2(config_path), config_path)
  recorded_sha <- stats::setNames(as.character(recorded$sha256), as.character(recorded$label))
  current_sha <- stats::setNames(as.character(current$sha256), as.character(current$label))
  absent <- setdiff(names(recorded_sha), names(current_sha))
  added <- setdiff(names(current_sha), names(recorded_sha))
  changed <- names(recorded_sha)[!vapply(names(recorded_sha), function(label) {
    label %in% names(current_sha) &&
      identical(unname(recorded_sha[[label]]), unname(current_sha[[label]]))
  }, logical(1))]
  if (length(absent) || length(added) || length(changed)) {
    stop("Freeze evidence changed since the manifest was written: ",
         paste(unique(c(absent, added, changed)), collapse = ", "))
  }
  invisible(TRUE)
}

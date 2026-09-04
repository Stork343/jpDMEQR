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

geometry_manifest_path_v2 <- function(root) {
  file.path(root, "results", "preflight", "geometry", "manifest.json")
}

validate_geometry_gate_v2 <- function(root, final = TRUE) {
  path <- geometry_manifest_path_v2(root)
  if (!file.exists(path)) return(list(pass = FALSE, reason = "strict geometry manifest missing"))
  x <- read_json_v2(path)
  checks <- c(
    isTRUE(x$pass),
    isTRUE(x$strict),
    as.integer(x$seed_count %||% 0) >= if (final) 20L else 1L,
    identical(as.character(x$commit_sha), current_commit_v2(root)),
    isTRUE(x$dantzig_checked),
    isTRUE(x$dantzig_pass)
  )
  list(pass = all(checks), reason = if (all(checks)) "accepted" else "geometry evidence failed or is stale",
       manifest = x)
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

micro_preflight_manifest_path_v2 <- function(root) {
  file.path(root, "results", "preflight", "micro_preflight_manifest.json")
}

pilot_gate_manifest_path_v2 <- function(root) {
  file.path(root, "results", "preflight", "pilot_gate.json")
}

validate_named_gate_manifest_v2 <- function(root, path, label) {
  if (!file.exists(path)) return(list(pass = FALSE, reason = paste(label, "manifest missing")))
  x <- read_json_v2(path)
  pass <- isTRUE(x$pass) && identical(as.character(x$commit_sha), current_commit_v2(root))
  list(pass = pass, reason = if (pass) "accepted" else paste(label, "gate failed or is stale"), manifest = x)
}

freeze_manifest_path_v2 <- function(root) {
  file.path(root, "results", "preflight", "freeze_manifest.json")
}

freeze_preflight_v2 <- function(root,
                                config_path = file.path(root, "config", "simulation_main.csv"),
                                final = TRUE,
                                write_manifest = TRUE) {
  governance <- validate_governance_files_v2(root)
  registry <- validate_registry_contract_v2(root, config_path, "main")
  cfg <- read_simulation_registry_v2(config_path)
  benchmark <- validate_benchmark_gate_v2(root, cfg, final = final)
  geometry <- validate_geometry_gate_v2(root, final = final)
  targets <- validate_target_gate_v2(
    root, cfg, final = final, expected_config_sha = file_checksum_v2(config_path)
  )
  micro <- validate_named_gate_manifest_v2(
    root, micro_preflight_manifest_path_v2(root), "micro-preflight"
  )
  pilot <- validate_named_gate_manifest_v2(
    root, pilot_gate_manifest_path_v2(root), "pilot"
  )
  pass <- governance$pass && registry$pass && benchmark$pass && geometry$pass &&
    targets$pass && micro$pass && pilot$pass
  manifest <- list(
    pass = pass,
    final = final,
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
    hardware = hardware_manifest_v2()
  )
  if (isTRUE(write_manifest)) write_json_v2(manifest, freeze_manifest_path_v2(root))
  manifest
}

verify_freeze_manifest_v2 <- function(root, config_path) {
  path <- freeze_manifest_path_v2(root)
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
  invisible(x)
}

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(".", mustWork = TRUE)
run_id <- if (length(args)) args[[1L]] else stop("An immutable run ID is required.")
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, environment())
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
source(file.path(root, "scripts", "research", "audit_profile_pair_geometry_20261022.R"))
configuration <- profile_pair_audit_configuration_20261022("production")
environment_fields <- c(dantzig_solver = "JPDMEQR_DANTZIG_SOLVER",
                        cluster_solver = "JPDMEQR_CLUSTER_SOLVER", assembly = "JPDMEQR_ASSEMBLY",
                        cv_workers = "JPDMEQR_CV_CORES", hessian_refresh = "JPDMEQR_HESSIAN_REFRESH",
                        shared_first_stage = "JPDMEQR_SHARED_FIRST_STAGE",
                        matrix_backend = "JPDMEQR_MATRIX_BACKEND", identity_detail = "JPDMEQR_IDENTITY_DETAIL",
                        proximal_scaling = "JPDMEQR_PROXIMAL_SCALING",
                        unpenalised_optimizer = "JPDMEQR_UNPENALISED_OPTIMIZER")
do.call(Sys.setenv, setNames(lapply(configuration, as.character), environment_fields))
config_path <- file.path(root, "config", "simulation_application_method_preflight_pair_20261022.csv")
registry <- read_simulation_registry_v2(config_path)
geometry <- validate_geometry_gate_v2(root, final = TRUE, method_version = "profile_pair_20261022")
benchmarks <- validate_benchmark_gate_v2(root, registry, final = TRUE)
targets <- validate_target_gate_v2(root, registry, final = TRUE)
if (!geometry$pass || !benchmarks$pass || !targets$pass) stop("The numerical prerequisites have not passed.")
result <- run_registry_v2(root, config_path, run_id, jobs = 2L,
                          allow_missing_benchmarks = FALSE, final = TRUE)
expected_methods <- sum(vapply(registry$methods, function(x) length(parse_pipe_character_v2(x)), integer(1))) * 2L
expected_coordinates <- expected_methods * 4L
schema_pass <- identical(names(result$metrics), contract_schema_columns_v2(root, 2L)) &&
  identical(names(result$coordinates), contract_schema_columns_v2(root, 3L)) &&
  identical(names(result$theory), contract_schema_columns_v2(root, 4L)) &&
  identical(names(result$screening), contract_schema_columns_v2(root, 5L))
pass <- length(result$errors) == 0L && nrow(result$metrics) == expected_methods &&
  nrow(result$coordinates) == expected_coordinates && schema_pass &&
  all(result$metrics$status %in% c("ok", "warning")) &&
  !anyDuplicated(result$metrics[c("experiment_id", "replicate", "method_id")]) &&
  !anyDuplicated(result$coordinates[c("experiment_id", "replicate", "method_id", "coordinate")])
manifest <- list(
  pass = pass, final = TRUE, scope = "application_paired_method_software_preflight",
  commit_sha = current_commit_v2(root), method_version = "profile_pair_20261022",
  score_definition = "profile_bandwidth_pair_2_half_minus_full",
  config_path = config_path, run_dir = normalizePath(result$run_dir),
  schema_pass = schema_pass, method_rows = nrow(result$metrics),
  coordinate_rows = nrow(result$coordinates), expected_method_rows = expected_methods,
  expected_coordinate_rows = expected_coordinates, runner_failure_count = length(result$errors),
  solver_configuration = solver_configuration_v2(),
  initial_configuration = list(clinical_acceleration = Sys.getenv("JPDMEQR_INITIAL_ACCELERATION", unset = "none")),
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write_json_v2(manifest, file.path(result$run_dir, "method_preflight_manifest.json"))
if (!pass) stop("The paired method software preflight failed: ", result$run_dir)
destination <- file.path(profile_pair_gate_directory_v2(root), "application_method_preflight_manifest.json")
dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
write_json_v2(manifest, destination)
print(manifest)

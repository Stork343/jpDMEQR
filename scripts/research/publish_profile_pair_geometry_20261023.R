args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("需要提供真实 geometry manifest 路径。")
root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
source_path <- normalizePath(args[[1L]], mustWork = TRUE)
record <- read_json_v2(source_path)
check <- validate_geometry_record_v2(
  record, current_commit_v2(root), solver_configuration_v2(), final = TRUE,
  method_version = "profile_pair_20261022"
)
if (!check$pass) stop("geometry 检查未通过：", check$reason)
destination <- geometry_manifest_path_v2(root, "profile_pair_20261022")
if (file.exists(destination)) stop("geometry 发布路径已经存在，需要保留原始记录。")
dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
if (!file.copy(source_path, destination, overwrite = FALSE)) stop("geometry manifest 复制失败。")
if (!identical(read_json_v2(destination), record)) stop("geometry manifest 读取结果不一致。")
message(destination)

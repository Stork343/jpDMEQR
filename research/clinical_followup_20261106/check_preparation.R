library(data.table)
library(yaml)
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
source_config <- read_yaml(args[1L])
reference <- source_config$output_root
out <- tempfile("compare_preparation_", tmpdir = "results/scratch/clinical_followup_20261106")
stopifnot(dir.create(out), length(source_config$quantiles) == 3L)
config <- source_config
config$output_root <- file.path(out, "fits")
config_file <- file.path(out, "config.yml")
write_yaml(config, config_file)
prepare <- "research/clinical_search_20261105/prepare_discovery.R"
stopifnot(system2("Rscript", c(prepare, config_file, "necessary")) == 0L)
stopifnot(system2("Rscript", c(prepare, config_file, "offsets")) == 0L)
checked <- 0L
for (fold in seq_len(config$outer_folds)) {
  label <- paste0("outer_", fold)
  original <- readRDS(file.path(reference, label, "recipes.rds"))
  direct <- readRDS(file.path(config$output_root, label, "recipes.rds"))
  offsets <- readRDS(file.path(config$output_root, label, "offset_recipes.rds"))
  stopifnot(identical(names(direct), names(original)), identical(offsets, original))
  for (name in names(direct)) {
    stopifnot(identical(direct[[name]]$direction_recipe$base_recipe, original[[name]]$direction_recipe$base_recipe))
    checked <- checked + 1L
  }
}
adapter <- "research/clinical_search_20261105/run_discovery.R"
implementation <- file.path(config$output_root, "implementation")
methods <- c("CLINICAL_QR","ENET","PRESS_RIDGE")
for (method in methods) {
  stopifnot(system2("Rscript", c(adapter,"outer_1",method,"all","1",implementation)) == 0L)
  fitted <- fread(file.path(config$output_root,"outer_1","models",paste0(method,"_all"),"prediction.csv"))
  if (method %in% c("CLINICAL_QR","ENET")) {
    original <- fread(file.path(reference,"outer_1","models",paste0(method,"_all"),"prediction.csv"))
    stopifnot(identical(fitted$sample_id, original$sample_id), max(abs(fitted$prediction-original$prediction)) < 1e-10)
  } else {
    stopifnot(nrow(fitted) > 0L, all(is.finite(fitted$prediction)))
  }
}
dest <- "results/preflight/clinical_followup_20261106"
write_json(list(study_id = source_config$study_id, verified_training_designs = checked,
  complete_offset_objects_identical = TRUE, evaluated_adapters = methods, reserved_outcomes_used = FALSE),
  file.path(dest,"preparation_equivalence.json"),auto_unbox = TRUE,pretty = TRUE)

library(yaml)
library(jsonlite)
source("research/clinical_quantile_20261105/engine.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
study <- read_yaml(args[1L])
config <- read_yaml("config/clinical_quantile_20261105.yml")
root <- file.path(config$output_root, study$study_id)
if (dir.exists(root)) stop("计算目录已经存在")
dir.create(root, recursive = TRUE)
stopifnot(file.copy("config/clinical_quantile_20261105.yml", file.path(root, "method_config.yml")),
  file.copy(args[1L], file.path(root, "source_config.yml")),
  file.copy("research/clinical_quantile_20261105/engine.R", file.path(root, "engine.R")),
  file.copy("research/clinical_quantile_20261105/qp_path.R", file.path(root, "qp_path.R")),
  file.copy("research/clinical_quantile_20261105/elastic_path.R", file.path(root, "elastic_path.R")))
write_json(list(status = "running", method = config$method_id,
  source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  formal_validation_used = FALSE, manuscript_adoption = FALSE),
  file.path(root, "status.json"), auto_unbox = TRUE, pretty = TRUE)
panel <- cp_panel(study$input_directory)
records <- diagnostics <- list()
for (fold in seq_len(study$outer_folds)) {
  outer <- file.path(study$output_root, paste0("outer_", fold))
  rows <- fread(file.path(outer, "training.csv"))
  evaluation <- fread(file.path(outer, "evaluation.csv"))
  recipes <- readRDS(file.path(outer, "recipes.rds"))
  stopifnot(all(rows$role == "training"), all(evaluation$role == "discovery_evaluation"))
  for (tau in unlist(config$quantiles)) {
    directory <- file.path(root, paste0("outer_", fold), sprintf("tau%03d", as.integer(100 * tau)))
    write_json(list(status = "running", fold = fold, tau = tau, completed = length(diagnostics)),
      file.path(root, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
    started <- proc.time()[["elapsed"]]
    fitted <- cq_fit(rows, evaluation, panel, recipes, study, config, tau, directory)
    elapsed <- proc.time()[["elapsed"]] - started
    record <- data.table(study_id = study$study_id, fold, method = "CQ_ENET", tau,
      sample_id = evaluation$sample_id, strain_id = evaluation$strain_id, y = evaluation$y,
      prediction = fitted$prediction)
    fwrite(record, file.path(directory, "prediction.csv"))
    records[[length(records) + 1L]] <- record
    diagnostics[[length(diagnostics) + 1L]] <- data.table(study_id = study$study_id, fold, tau,
      alpha = fitted$alpha, lambda = fitted$lambda, selected_index = fitted$grid[fitted$choice]$index,
      nonzero_proteins = sum(abs(fitted$beta) > 1e-7 & startsWith(names(fitted$beta), "P_")),
      clinical_columns = length(fitted$clinical_columns),
      path_objective_difference = fitted$path_objective_difference, gap = fitted$diagnostics$gap,
      stationarity = fitted$diagnostics$stationarity, runtime_sec = elapsed)
    fwrite(rbindlist(diagnostics), file.path(root, "diagnostics.csv"))
    writeLines(trimws(capture.output(warnings()), which = "right"), file.path(directory, "warnings.txt"))
    write_json(list(status = "completed", runtime_sec = elapsed), file.path(directory, "status.json"),
      auto_unbox = TRUE, pretty = TRUE)
    cat(study$study_id, "completed", fold, tau, elapsed, "seconds\n")
    flush.console()
    gc()
  }
}
fwrite(rbindlist(records), file.path(root, "predictions.csv"))
write_json(list(status = "completed", fits = length(diagnostics), formal_validation_used = FALSE,
  manuscript_adoption = FALSE), file.path(root, "status.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(root, "sessionInfo.txt"))

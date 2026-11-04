library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
library(glmnet)
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
source("research/drugmatrix_application_20261101/prediction.R")
source("research/empirical_selection_20261101/prediction_core.R")

config <- read_yaml("config/discovery_drugmatrix_20261101.yml")
panel <- readRDS(config$input_panel)
annotation <- fread(config$input_annotation)
panel$allowed_probes <- intersect(annotation[TYPE == config$probe_type, ID], colnames(panel$X))
rows <- readRDS(config$input_growth)
directory <- "results/scratch/empirical_selection_20261101/core_equivalence"
if (dir.exists(directory)) stop("检查目录已经存在")
dir.create(directory, recursive = TRUE)
reference <- readRDS("results/raw/application/drugmatrix_20261101/discovery/ALT/JOINT_qp_tau075/fit.rds")
fit <- ep_fit_prediction(rows[exploration_role == "training"], rows[exploration_role == "evaluation"],
  panel, "JOINT", 0.75, config, function(stage, details) cat(stage, "\n"), directory,
  function(rows, panel, config) dm_recipe(rows, panel), dm_features, "compound_id")
error <- max(abs(fit$prediction - reference$prediction), abs(fit$beta - reference$beta),
               abs(fit$loss - reference$loss), abs(fit$oof - reference$oof))
stopifnot(error < 1e-8, fit$choice == reference$choice)
output <- "results/preflight/empirical_selection_20261101/prediction_core"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write_json(list(status = "passed", actual_dataset = "DrugMatrix", maximum_error = error,
                selected_index = fit$choice, formal_validation_used = FALSE),
           file.path(output, "equivalence.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("预测流程一致性检查完成\n")

library(data.table)
library(yaml)
library(jsonlite)
source("research/drug_response_application_20261101/data.R")
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")

config <- read_yaml("config/discovery_gdsc_20261101.yml")
data <- readRDS("data/derived/drug_response_20261101/discovery_cases.rds")
rows <- data[DRUG_NAME == "Dabrafenib" & exploration_role == "training" &
               analysis_curve & inner_fold != 1L]
expression <- readRDS("data/derived/drug_response_20261101/expression.rds")
recipe <- drug_recipe(rows, expression, config)
recipe$genes <- recipe$genes[1:3]
recipe$centers <- recipe$centers[1:3]
recipe$scales <- recipe$scales[1:3]
recipe$medians <- recipe$medians[1:3]
features <- drug_features(rows, expression, recipe)
selected <- c(which(features$penalty_factor == 0),
               which(startsWith(colnames(features$X), recipe$genes[1L])))
X <- features$X[, selected, drop = FALSE]
checks <- rbindlist(lapply(c(0.25, 0.5, 0.75), function(tau) {
  prepared <- pc_prepare_initial(features$y, X, features$Z, features$cluster_id,
    tau, matrix(1, 1L, 1L), rep(0, ncol(X)), fit_always_included = FALSE)
  started <- proc.time()[["elapsed"]]
  conic <- pc_orthogonal_path(prepared, tau, 0L, gap_tolerance = 1e-8, solver = jp_joint_qp)
  conic_time <- proc.time()[["elapsed"]] - started
  started <- proc.time()[["elapsed"]]
  reference <- pc_orthogonal_path(prepared, tau, 0L, gap_tolerance = 1e-8,
                                  max_iterations = 300000L, solver = joint_profile_pdal_cpp)
  reference_time <- proc.time()[["elapsed"]] - started
  prediction_difference <- drop(X %*% (conic$coefficients[, 1L] - reference$coefficients[, 1L]))
  distance <- sum(prediction_difference^2 * prepared$inverse_size) / prepared$n
  bound <- 4 * prepared$n * (max(conic$fitted$duality_gap, 0) + max(reference$fitted$duality_gap, 0) + 1e-12)
  stopifnot(distance <= bound,
            abs(conic$fitted$primal - reference$fitted$primal) < 2e-8)
  data.table(tau = tau, clusters = prepared$n, variables = ncol(X),
             conic_gap = conic$fitted$duality_gap, reference_gap = reference$fitted$duality_gap,
             prediction_distance = distance, certificate_bound = bound,
             conic_time_sec = conic_time, reference_time_sec = reference_time)
}))
directory <- "results/preflight/empirical_selection_20261101/selected_qp"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
fwrite(checks, file.path(directory, "checks.csv"))
write_json(list(status = "passed", actual_observations_used = TRUE,
                objective = "complete_joint_profile_with_prediction_ridge",
                solver = "clarabel_selected_qp"),
           file.path(directory, "checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(checks)

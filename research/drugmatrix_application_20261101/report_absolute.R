library(data.table)
library(jsonlite)
library(yaml)
source("research/drugmatrix_application_20261101/prediction.R")
source("research/empirical_selection_20261101/prediction_core.R")

config <- read_yaml("config/discovery_drugmatrix_absolute_20261104.yml")
rows <- readRDS(config$input_growth)[exploration_role == "evaluation"]
stopifnot(all(rows$role == "discovery"))
rows[, y := alt_u_l]
panel <- readRDS(config$input_panel)
paths <- c(JOINT_raw = file.path(config$output_root, "JOINT_tau075"),
           POOLED_raw = file.path(config$output_root, "POOLED_tau075"),
           JOINT_log = "results/raw/application/drugmatrix_20261101/discovery/ALT/JOINT_prediction_tau075",
           POOLED_log = "results/raw/application/drugmatrix_20261101/discovery/ALT/POOLED_tau075",
           ENET_log = "results/raw/application/drugmatrix_20261101/discovery/ALT/ENET_weighted_cv")
records <- checks <- list()
for (name in names(paths)) {
  directory <- paths[[name]]
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed")
  saved <- fread(file.path(directory, "predictions.csv"))[tau == 0.75]
  stopifnot(identical(saved$sample_id, rows$sample_id), identical(saved$compound_id, rows$compound_id))
  fit <- readRDS(file.path(directory, "fit.rds"))
  stopifnot(!length(intersect(fit$recipe$training_compounds, rows$compound_id)))
  if (endsWith(name, "_raw")) {
    reconstructed <- drop(dm_predict(dm_features(rows, panel, fit$recipe), fit$beta)) + fit$calibration[[1L]]["all"]
    error <- max(abs(reconstructed - saved$prediction), abs(saved$y - rows$y))
    prediction <- saved$prediction
  } else {
    error <- max(abs(exp(saved$y) - rows$y))
    prediction <- exp(saved$prediction)
  }
  stopifnot(is.finite(error), error < 1e-7, all(is.finite(prediction)))
  checks[[name]] <- data.table(method = name, maximum_identity_error = error)
  records[[name]] <- data.table(compound_id = rows$compound_id, loss = ep_loss(rows$y, prediction, 0.75))[
    , .(loss = mean(loss)), by = compound_id][, method := name]
}
wide <- dcast(rbindlist(records), compound_id ~ method, value.var = "loss")
stopifnot(!anyNA(wide))
set.seed(20261104L)
n <- nrow(wide)
weights <- rmultinom(10000L, n, rep(1 / n, n)) / n
comparisons <- rbindlist(lapply(setdiff(names(wide), c("compound_id", "JOINT_raw")), function(name) {
  reference <- wide[[name]]
  difference <- reference - wide$JOINT_raw
  bootstrap <- drop(crossprod(weights, difference)) / drop(crossprod(weights, reference))
  data.table(comparator = name, compounds = n, joint_loss = mean(wide$JOINT_raw),
    comparator_loss = mean(reference), relative_gain = mean(difference) / mean(reference),
    lower = quantile(bootstrap, 0.025), upper = quantile(bootstrap, 0.975),
    units = "ALT_U_per_L", phase = "discovery", manuscript_adoption = FALSE)
}))
output <- "results/preflight/empirical_selection_20261101/drugmatrix_absolute"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(comparisons, file.path(output, "comparisons.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
write_json(list(stage = "discovery", formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(comparisons)

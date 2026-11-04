library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")

config <- read_yaml("config/development_ridge_direction_20261105.yml")
base <- read_yaml(config$base_config)
base$conditions <- base$assays$NaCl
panel <- proteome_panel(config$input_directory)
recipes <- readRDS(file.path(config$output_root, "recipes.rds"))
rows <- fread(file.path(config$input_directory, "NaCl", "discovery.csv"))
training <- rows[exploration_role == "training"]
evaluation <- rows[exploration_role == "evaluation"]
directory <- file.path(config$output_root, "NaCl", "ADJUSTMENT_MEAN_cv")
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
started <- proc.time()[["elapsed"]]
oof <- rep(NA_real_, nrow(training))
normal_equation_error <- 0
for (fold in sort(unique(training$inner_fold))) {
  recipe <- recipes[[paste0("inner_", fold)]]$base_recipe
  fitting <- proteome_features(training[inner_fold != fold], panel, recipe)
  validation <- proteome_features(training[inner_fold == fold], panel, recipe)
  low <- which(fitting$penalty_factor == 0)
  model <- lm.wfit(fitting$X[, low, drop = FALSE], fitting$y, fitting$weight)
  stopifnot(model$rank == length(low), all(is.finite(model$coefficients)))
  normal_equation_error <- max(normal_equation_error,
    abs(crossprod(fitting$X[, low, drop = FALSE], model$residuals * fitting$weight)))
  prediction <- drop(validation$X[, low, drop = FALSE] %*% model$coefficients) * recipe$y_scale + recipe$y_center
  oof[which(training$inner_fold == fold)] <- prediction
}
stopifnot(all(is.finite(oof)), normal_equation_error < 1e-8)
recipe <- recipes$training$base_recipe
fitting <- proteome_features(training, panel, recipe)
test <- proteome_features(evaluation, panel, recipe)
low <- which(fitting$penalty_factor == 0)
model <- lm.wfit(fitting$X[, low, drop = FALSE], fitting$y, fitting$weight)
stopifnot(model$rank == length(low), all(is.finite(model$coefficients)))
beta <- setNames(numeric(ncol(fitting$X)), colnames(fitting$X))
beta[low] <- model$coefficients
quantiles <- unlist(config$quantiles)
calibration <- lapply(quantiles, function(value) ep_calibrate(training, oof, value, "strain_id", "condition"))
predicted <- drop(ep_predict(test, beta))
prediction <- vapply(calibration, function(offset) predicted + offset[evaluation$condition], numeric(nrow(evaluation)))
fit <- list(beta = beta, recipe = recipe, prediction = prediction, quantiles = quantiles,
             calibration = calibration, oof = oof, normal_equation_error = normal_equation_error)
saveRDS(fit, file.path(directory, "fit.rds"))
predictions <- rbindlist(lapply(seq_along(quantiles), function(j)
  data.table(sample_id = evaluation$sample_id, strain_id = evaluation$strain_id, condition = evaluation$condition,
    y = evaluation$y, tau = quantiles[j], prediction = prediction[, j],
    loss = ep_loss(evaluation$y, prediction[, j], quantiles[j]))))
fwrite(predictions, file.path(directory, "predictions.csv"))
losses <- predictions[, .(loss = mean(loss)), by = .(strain_id, tau)]
fwrite(losses, file.path(directory, "strain_loss.csv"))
summary <- losses[, .(strains = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = tau]
elapsed <- proc.time()[["elapsed"]] - started
summary[, `:=`(method = "ADJUSTMENT_MEAN", runtime_sec = elapsed, phase = "method_development")]
fwrite(summary, file.path(directory, "summary.csv"))
write_json(list(method = "ADJUSTMENT_MEAN", config = config, code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE)
write_json(list(status = "completed", runtime_sec = elapsed, formal_validation_used = FALSE,
                manuscript_adoption = FALSE), file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(summary)

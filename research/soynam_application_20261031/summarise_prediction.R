library(data.table)
library(yaml)
config <- read_yaml("config/application_soynam_20261031.yml")
arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments)) arguments[1L] else "all"
if (!mode %in% c("all", "joint")) stop("Unknown prediction summary mode.")
input <- readRDS(file.path(config$results_directory, "input.rds"))
root <- file.path(config$results_directory, "prediction", config$execution_version)
destination <- if (mode == "all") "results/tables/soynam_20261031" else "results/preflight/soynam_20261031/joint_prediction"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
methods <- if (mode == "all") unlist(config$prediction$methods) else c("joint_canopy", "joint_genomic")
registry <- CJ(method = methods, tau = unlist(config$quantiles), sorted = FALSE)
registry[, case_id := paste0(method, "_q", sprintf("%02d", round(100 * tau)))]
complete <- file.exists(file.path(root, registry$case_id, "completion.json"))
if (!all(complete)) stop("Incomplete registered predictions: ", paste(registry$case_id[!complete], collapse = ", "))
heldout <- input$line_folds[prediction_fold == config$prediction$heldout_fold, cluster_id]
predictions <- list()
tuning <- list()
for (i in seq_len(nrow(registry))) {
  directory <- file.path(root, registry$case_id[i])
  request <- readRDS(file.path(directory, "request.rds"))
  model <- readRDS(file.path(directory, "model.rds"))
  current <- fread(file.path(directory, "predictions.csv"))
  if (!setequal(heldout, request$test_lines) || !setequal(heldout, current$cluster_id) ||
      length(intersect(request$training_lines, request$test_lines))) stop("Prediction roles differ from their fixed allocation.")
  counts <- current[, .N, by = cluster_id]
  if (any(counts$N != 2L) || anyDuplicated(current, by = c("cluster_id", "Year"))) stop("The heldout pairs are incomplete.")
  model$tuning[, `:=`(method = registry$method[i], tau = registry$tau[i])]
  tuning[[i]] <- model$tuning
  predictions[[i]] <- current
}
predictions <- rbindlist(predictions)
line_loss <- predictions[, .(loss = mean(loss), raw_loss = mean(raw_loss),
  calibration = mean(covered), raw_calibration = mean(raw_covered)), by = .(method, tau, cluster_id)]
metrics <- line_loss[, .(lines = .N, loss = mean(loss), loss_se = sd(loss) / sqrt(.N),
  raw_loss = mean(raw_loss), calibration = mean(calibration), raw_calibration = mean(raw_calibration)), by = .(method, tau)]
year_metrics <- predictions[, .(lines = uniqueN(cluster_id), loss = mean(loss), raw_loss = mean(raw_loss),
  calibration = mean(covered), raw_calibration = mean(raw_covered)), by = .(method, tau, Year)]
critical <- qnorm(1 - (1 - config$confidence) / (2 * config$prediction$comparison_family_size))
comparisons <- list()
comparator_ids <- if (mode == "all") c("joint_genomic", "pooled_quantile_lasso") else "joint_genomic"
for (tau_value in unlist(config$quantiles)) for (comparator in comparator_ids) {
  selected <- line_loss[tau == tau_value & method %in% c("joint_canopy", comparator)]
  wide <- dcast(selected, cluster_id ~ method, value.var = c("loss", "raw_loss"))
  if (nrow(wide) != length(heldout) || anyNA(wide)) stop("The paired prediction losses are incomplete.")
  for (calibrated in c(TRUE, FALSE)) {
    prefix <- if (calibrated) "loss_" else "raw_loss_"
    a <- wide[[paste0(prefix, "joint_canopy")]]
    b <- wide[[paste0(prefix, comparator)]]
    difference <- b - a
    estimate <- mean(difference)
    se <- sd(difference) / sqrt(length(difference))
    relative <- 100 * (1 - mean(a) / mean(b))
    relative_influence <- 100 * (-(a - mean(a)) / mean(b) + mean(a) * (b - mean(b)) / mean(b)^2)
    relative_se <- sd(relative_influence) / sqrt(length(relative_influence))
    comparisons[[length(comparisons) + 1L]] <- data.table(tau = tau_value, comparator, calibrated,
      lines = length(difference), canopy_loss = mean(a), comparator_loss = mean(b),
      improvement = estimate, se, lower = estimate - critical * se, upper = estimate + critical * se,
      relative_percent = relative, relative_se, relative_lower = relative - critical * relative_se,
      relative_upper = relative + critical * relative_se,
      p_bonferroni = pmin(1, config$prediction$comparison_family_size * 2 * pnorm(-abs(estimate / se))))
  }
}
fwrite(metrics, file.path(destination, "prediction_metrics.csv"))
fwrite(year_metrics, file.path(destination, "prediction_year_metrics.csv"))
fwrite(rbindlist(comparisons), file.path(destination, "prediction_comparisons.csv"))
fwrite(rbindlist(tuning, fill = TRUE), file.path(destination, "prediction_tuning.csv"))
print(metrics)
print(rbindlist(comparisons)[calibrated == TRUE])

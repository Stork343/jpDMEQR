library(data.table)
library(jsonlite)
library(yaml)
source("research/empirical_selection_20261101/prediction_core.R")
config <- read_yaml("config/application_charite_hb_20261105.yml")
root <- config$output_root
release <- read_json(file.path(root, "validation_release.json"), simplifyVector = TRUE)
stopifnot(release$status == "frozen", release$primary_quantile == config$decision_quantile)
rows <- fread(file.path(config$input_directory, "validation.csv"))
prediction <- rbindlist(list(fread(file.path(root, "validation", "predictions_R.csv")),
  fread(file.path(root, "validation", "predictions_boost.csv"))))
methods <- c(config$target_method, unlist(config$comparators))
stopifnot(setequal(unique(prediction$method), methods), !anyDuplicated(prediction[, .(sample_id, method, tau)]),
  nrow(prediction) == nrow(rows) * length(methods) * length(config$quantiles))
for (method_name in methods) for (value in unlist(config$quantiles))
  stopifnot(setequal(prediction[method == method_name & tau == value]$sample_id, rows$sample_id))
prediction <- merge(prediction, rows[, .(sample_id, actual_strain_id = strain_id, y)], by = "sample_id", sort = FALSE)
stopifnot(identical(prediction$strain_id, prediction$actual_strain_id), all(is.finite(prediction$prediction)))
prediction[, actual_strain_id := NULL]
raw <- copy(prediction)
raw[, `:=`(scoring = "raw", loss = ep_loss(y, prediction, tau))]
primary <- copy(prediction)
setorder(primary, method, sample_id, tau)
primary[, prediction := sort(pmax(prediction, 0)), by = .(method, sample_id)]
primary[, `:=`(scoring = "primary", loss = ep_loss(y, prediction, tau))]
all_predictions <- rbindlist(list(raw, primary))
losses <- all_predictions[, .(loss = mean(loss)), by = .(scoring, strain_id, method, tau)]
stopifnot(all(losses[, .N, by = .(scoring, method, tau)]$N == uniqueN(rows$strain_id)))
summary <- losses[, .(patients = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = .(scoring, method, tau)]
set.seed(20261105L)
n <- uniqueN(rows$strain_id)
weights <- rmultinom(config$validation$bootstrap_repetitions, n, rep(1 / n, n)) / n
comparisons <- list()
for (score_name in c("primary", "raw")) for (value in unlist(config$quantiles)) {
  wide <- dcast(losses[scoring == score_name & tau == value], strain_id ~ method, value.var = "loss")
  stopifnot(!anyNA(wide), identical(wide$strain_id, sort(unique(rows$strain_id))))
  for (method_name in unlist(config$comparators)) {
    reference <- wide[[method_name]]
    difference <- reference - wide[[config$target_method]]
    bootstrap_difference <- drop(crossprod(weights, difference))
    bootstrap_gain <- bootstrap_difference / drop(crossprod(weights, reference))
    stopifnot(all(is.finite(bootstrap_gain)))
    lower <- unname(quantile(bootstrap_gain, config$validation$prediction_alpha))
    equivalent <- max(abs(difference)) <= 1e-8
    comparisons[[length(comparisons) + 1L]] <- data.table(scoring = score_name, tau = value,
      target_method = config$target_method, comparator = method_name, patients = n,
      proposed_loss = mean(wide[[config$target_method]]), comparator_loss = mean(reference),
      absolute_gain = mean(difference), relative_gain = mean(difference) / mean(reference),
      lower95 = quantile(bootstrap_gain, 0.025), upper95 = quantile(bootstrap_gain, 0.975),
      one_sided_lower = lower, absolute_one_sided_lower = quantile(bootstrap_difference, config$validation$prediction_alpha),
      numerically_equivalent = equivalent, superiority_pass = lower > 0 && !equivalent,
      formal_validation_used = TRUE)
  }
}
comparison <- rbindlist(comparisons)
calibration <- all_predictions[, .(below = mean(y < prediction), at_or_below = mean(y <= prediction)),
  by = .(scoring, method, tau, strain_id)][,
  .(patients = .N, fraction_below = mean(below), fraction_at_or_below = mean(at_or_below)), by = .(scoring, method, tau)]
crossing <- raw[, {
  wide <- dcast(.SD, sample_id ~ tau, value.var = "prediction")
  gap <- cbind(wide[["0.25"]] - wide[["0.5"]], wide[["0.5"]] - wide[["0.75"]])
  list(records = nrow(wide), crossing_records = sum(apply(gap, 1L, max) > 1e-10),
    maximum_crossing = max(0, gap), negative_predictions = sum(prediction < 0))
}, by = method]
output <- config$report_directory
fwrite(summary, file.path(output, "prediction_summary.csv"))
fwrite(comparison, file.path(output, "prediction_comparisons.csv"))
fwrite(calibration, file.path(output, "quantile_calibration.csv"))
fwrite(crossing, file.path(output, "quantile_crossings.csv"))
fwrite(all_predictions, file.path(root, "validation", "scored_predictions.csv"))
decision <- comparison[scoring == "primary" & tau == config$decision_quantile]
stopifnot(nrow(decision) == length(config$comparators))
write_json(list(status = "completed", stage = "formal_validation", patients = n, records = nrow(rows),
  primary_quantile = config$decision_quantile, primary_scoring = config$primary_prediction_rule,
  prediction_alpha = config$validation$prediction_alpha, superiority_pass = all(decision$superiority_pass),
  formal_validation_used = TRUE, manuscript_adoption = FALSE), file.path(output, "prediction_validation_status.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(decision)

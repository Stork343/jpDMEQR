library(yaml)
source("research/clinical_quantile_20261105/engine.R")
source("research/clinical_quantile_20261105/orthogonal.R")
config <- read_yaml("config/clinical_quantile_smooth_20261106.yml")
records <- list()
for (path in unlist(config$source_studies)) {
  study <- read_yaml(path)
  root <- file.path(config$output_root, paste0("orthogonal_", study$study_id))
  if (!dir.exists(root)) next
  panel <- cp_panel(study$input_directory)
  for (fold in seq_len(study$outer_folds)) {
    rows <- fread(file.path(study$output_root, paste0("outer_", fold), "training.csv"))
    central <- readRDS(file.path(config$output_root, study$study_id, paste0("outer_", fold), "tau075", "fit.rds"))
    features <- cq_features(rows, panel, central$recipe, central$columns)
    nuisance <- readRDS(file.path(root, paste0("outer_", fold), "nuisance.rds"))
    density <- cq_density(features, nuisance$density_fits[[1L]]$beta,
      nuisance$density_fits[[2L]]$beta, nuisance$delta, config)
    residual <- rows$y - central$oof
    weight <- features$weight / features$n
    average <- sum(weight * residual)
    deviation <- sqrt(sum(weight * (residual - average)^2))
    quartiles <- vapply(c(0.25, 0.75), function(tau)
      unname(ep_calibrate(rows, central$oof, tau, "strain_id", "condition")[1L]), numeric(1))
    spread <- diff(quartiles) / (qnorm(0.75) - qnorm(0.25))
    stopifnot(deviation > 0, spread > 0)
    bandwidth <- 0.9 * min(deviation, spread) * features$n^(-0.2)
    marginal <- sum(weight * dnorm(residual / bandwidth)) / bandwidth
    conditional_average <- sum(weight * density$value) / central$recipe$y_scale
    records[[length(records) + 1L]] <- data.table(study_id = study$study_id, fold,
      patients = features$n, oof_residual_density = marginal, spacing_density_average = conditional_average,
      ratio = conditional_average / marginal, residual_bandwidth = bandwidth,
      constrained_fraction = (density$lower_bound_count + density$upper_bound_count) / nrow(rows))
  }
}
result <- rbindlist(records)
fwrite(result, file.path(config$report_directory, "density_diagnostics.csv"))
print(result)

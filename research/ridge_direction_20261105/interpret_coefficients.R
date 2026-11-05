library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")

config <- read_yaml("config/development_ridge_offset_series_20261105.yml")
panel <- proteome_panel(config$input_directory)
selected <- fread("results/preflight/empirical_selection_20261101/ridge_offset_series/selected_coefficients.csv")
annotation <- fread("data/raw/yeast_1011_20261104/Teyssonniere_2024/data_upload/CompleteGeneAnnot_strand.csv")
records <- checks <- list()
for (assay_name in c("Benomyl", "Caffeine")) {
  variables <- selected[assay == assay_name & role == "selected" & startsWith(variable, "P_")]$variable
  genes <- unique(sub("_condition1$", "", sub("^P_", "", variables)))
  evaluation <- fread(file.path(config$input_directory, assay_name, "discovery.csv"))[exploration_role == "evaluation"]
  for (value in c(0.25, 0.5, 0.75)) {
    run_id <- sprintf("tau%03d", as.integer(value * 100))
    fit <- readRDS(file.path(config$output_root, assay_name, paste0("JOINT_OFFSET_", run_id), "fit.rds"))
    recipe <- fit$recipe
    direction <- recipe$direction_recipe
    base <- direction$base_recipe
    d <- length(base$low_columns)
    p <- length(base$genes)
    base_dimension <- length(direction$background$beta)
    stopifnot(nrow(base$contrasts) == 2L, ncol(base$contrasts) == 1L,
              identical(base$low_columns, seq_len(d)), d %% 2L == 0L,
              base_dimension == d + 2L * p, length(fit$beta) == base_dimension + 2L)
    low_base <- (d - 2L) / 2L
    swapped <- c(2L, 1L, seq.int(3L + low_base, d), seq.int(3L, 2L + low_base),
                 d + p + seq_len(p), d + seq_len(p))
    stopifnot(identical(sort(as.integer(swapped)), seq_len(base_dimension)))
    scale <- rep(1, length(fit$beta))
    scale[recipe$penalized] <- recipe$scale
    correction <- fit$beta * recipe$y_scale / scale
    mean_beta <- direction$background$beta
    embedded_projection <- matrix(0, base_dimension, 2L)
    embedded_projection[direction$low_columns, ] <- direction$projection
    direction_map <- cbind(mean_beta, mean_beta[swapped]) - embedded_projection
    expanded <- base$y_scale * mean_beta + correction[seq_len(base_dimension)] +
      drop(direction_map %*% (correction[base_dimension + seq_len(2L)] / direction$direction_scale))
    expanded[1L] <- expanded[1L] + base$y_center - sum(correction[recipe$penalized] * recipe$center)
    H <- cbind(mean = 1, base$contrasts)
    expanded[1:2] <- expanded[1:2] + solve(H, fit$calibration[[1L]][rownames(H)])
    features <- proteome_features(evaluation, panel, base)
    reconstruction <- drop(features$X %*% expanded)
    error <- max(abs(reconstruction - fit$prediction[, 1L]))
    stopifnot(error < 1e-9)
    checks[[length(checks) + 1L]] <- data.table(assay = assay_name, tau = value,
      reconstructed_predictions = length(reconstruction), maximum_error = error)
    for (gene in genes) {
      j <- match(gene, base$genes)
      stopifnot(!is.na(j))
      indices <- c(d + j, d + p + j)
      total <- drop(H %*% expanded[indices]) / base$scale[j]
      location <- drop(H %*% mean_beta[indices]) * base$y_scale / base$scale[j]
      info <- annotation[systematic_name == gene]
      stopifnot(nrow(info) == 1L)
      records[[length(records) + 1L]] <- data.table(assay = assay_name, protein_id = gene,
        gene_symbol = info$Gene, condition = rownames(H), tau = value,
        total_predictive_association = total, mean_background_association = location,
        quantile_correction = total - location,
        covariate_unit = "one_log2_unit_relative_to_sample_median",
        response_unit = "growth_ratio_to_YPD", phase = "method_development")
    }
  }
}
output <- "results/preflight/empirical_selection_20261101/ridge_offset_series"
fwrite(rbindlist(checks), file.path(output, "coefficient_expansion_checks.csv"))
fwrite(rbindlist(records), file.path(output, "predictive_associations.csv"))
write_json(list(status = "passed", method = "algebraic_expansion_of_ridge_offset_and_correction",
                formal_validation_used = FALSE, confidence_intervals_computed = FALSE),
           file.path(output, "coefficient_interpretation_status.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(records)[tau == 0.25, .(assay, gene_symbol, condition, total_predictive_association,
                                      mean_background_association, quantile_correction)])

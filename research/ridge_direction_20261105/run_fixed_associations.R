library(data.table)
library(jsonlite)
library(yaml)
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/fixed_inference.R")
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")

config_path <- "config/development_proteome_associations_20261105.yml"
config <- read_yaml(config_path)
root <- config$output_root
if (dir.exists(root)) stop("运行目录已经存在")
dir.create(root, recursive = TRUE)
stopifnot(file.copy(config_path, file.path(root, "config.yml")))
source_paths <- c("research/ridge_direction_20261105/run_fixed_associations.R",
  "research/joint_profile_manuscript_20261026/fixed_inference.R",
  "research/score_validity_20261026/small_profile.R")
dir.create(file.path(root, "implementation"))
stopifnot(all(file.copy(source_paths, file.path(root, "implementation"))))
write_json(list(phase = config$phase, code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  formal_validation_used = FALSE, manuscript_adoption = FALSE, config = config),
  file.path(root, "request.json"), auto_unbox = TRUE, pretty = TRUE)
panel <- proteome_panel(config$input_directory)
selection <- fread(config$selection_table)
annotation <- fread("data/raw/yeast_1011_20261104/Teyssonniere_2024/data_upload/CompleteGeneAnnot_strand.csv")
records <- diagnostics <- checks <- designs <- list()

for (assay_name in names(config$assays)) {
  genes <- unlist(config$assays[[assay_name]])
  selected <- selection[assay == assay_name & role == "selected" & startsWith(variable, "P_")]$variable
  selected <- unique(sub("_condition1$", "", sub("^P_", "", selected)))
  stopifnot(setequal(genes, selected))
  rows <- fread(file.path(config$input_directory, assay_name, "discovery.csv"))
  training <- rows[exploration_role == "training"]
  estimation <- rows[exploration_role == config$estimation_role]
  stopifnot(all(rows$role == "discovery"), !length(intersect(training$strain_id, estimation$strain_id)))
  recipes <- readRDS(file.path(config$recipe_directory, assay_name, "recipes.rds"))
  recipe <- recipes$training$direction_recipe$base_recipe
  stopifnot(identical(recipe$training_strains, sort(unique(training$strain_id))))
  features <- proteome_features(estimation, panel, recipe)
  protein_columns <- unlist(lapply(genes, function(gene) c(paste0("P_", gene), paste0("P_", gene, "_condition1"))))
  columns <- c(colnames(features$X)[features$penalty_factor == 0], protein_columns)
  X <- features$X[, columns, drop = FALSE]
  data <- sb_data(features$y, X, features$Z, features$cluster_id)
  stopifnot(nrow(recipe$contrasts) == 2L, ncol(recipe$contrasts) == 1L, ncol(X) < data$n)
  gram <- crossprod(X, X * data$weight) / data$n
  mean_beta <- drop(qr.solve(X * sqrt(data$weight), data$y * sqrt(data$weight)))
  mean_residual <- data$y - drop(X %*% mean_beta)
  mean_scores <- rowsum(X * (data$weight * mean_residual), data$group, reorder = FALSE)
  mean_influence <- mean_scores %*% solve(gram)
  mean_kkt <- max(abs(crossprod(X, data$weight * mean_residual) / data$n))
  stopifnot(mean_kkt < 1e-9)
  designs[[length(designs) + 1L]] <- data.table(assay = assay_name, clusters = data$n,
    observations = nrow(X), dimension = ncol(X), selected_proteins = length(genes),
    zero_responses = sum(estimation$y == 0), response_scale = recipe$y_scale)

  for (multiplier in unlist(config$bandwidth_multipliers)) {
    directory <- file.path(root, assay_name, sprintf("bandwidth_%03d", as.integer(100 * multiplier)))
    dir.create(directory, recursive = TRUE)
    fits <- list()
    for (value in unlist(config$quantiles)) {
      cat(format(Sys.time(), tz = "UTC", usetz = TRUE), assay_name, multiplier, value, "\n")
      flush.console()
      fit <- jp_fit_fixed(data$y, X, data$Z, data$cluster_id, value, Lambda = matrix(config$working_lambda),
        h = data$n^(-0.3) * multiplier)
      stopifnot(max(fit$diagnostics$full_kkt, fit$diagnostics$half_kkt) < 1e-8,
        max(fit$diagnostics$full_nuisance_gradient, fit$diagnostics$half_nuisance_gradient) < 1e-8)
      coordinate <- match(protein_columns[1L], colnames(X))
      order <- c(coordinate, setdiff(seq_len(ncol(X)), coordinate))
      reference_data <- sb_data(data$y, X[, order, drop = FALSE], data$Z, data$cluster_id)
      reference <- ds_pair_small(reference_data, fit$full$beta[order], value,
        data$n^(-0.3) * multiplier, matrix(config$working_lambda))
      reference_se <- sqrt(sum((reference$influence - mean(reference$influence))^2)) / data$n
      estimate_error <- abs(reference$estimate - fit$beta[coordinate])
      se_error <- abs(reference_se - fit$table$se[coordinate])
      covariance_error <- max(abs(crossprod(scale(fit$influence, center = TRUE, scale = FALSE)) / data$n^2 - fit$covariance))
      stopifnot(max(estimate_error, se_error) < 2e-6, covariance_error < 1e-12)
      checks[[length(checks) + 1L]] <- data.table(assay = assay_name, multiplier, tau = value,
        estimate_error, se_error, covariance_error)
      diagnostics[[length(diagnostics) + 1L]] <- cbind(data.table(assay = assay_name,
        multiplier, tau = value), as.data.table(fit$diagnostics))
      fits[[sprintf("tau%03d", as.integer(value * 100))]] <- fit
    }
    model_names <- c(names(fits), "mean")
    coefficients <- unlist(c(lapply(fits, function(fit) unname(fit$beta)), list(mean_beta)), use.names = FALSE)
    influence <- do.call(cbind, c(lapply(fits, `[[`, "influence"), list(mean_influence)))
    centered <- scale(influence, center = TRUE, scale = FALSE)
    covariance <- crossprod(centered) / data$n^2
    contrasts <- labels <- list()
    add_contrast <- function(gene, condition, label, model_weights, condition_weights) {
      vector <- numeric(length(coefficients))
      locations <- match(c(paste0("P_", gene), paste0("P_", gene, "_condition1")), columns)
      for (model_name in names(model_weights)) {
        indices <- (match(model_name, model_names) - 1L) * ncol(X) + locations
        vector[indices] <- model_weights[[model_name]] * condition_weights * recipe$y_scale / recipe$scale[gene]
      }
      contrasts[[length(contrasts) + 1L]] <<- vector
      labels[[length(labels) + 1L]] <<- data.table(assay = assay_name, protein_id = gene,
        gene_symbol = annotation[systematic_name == gene]$Gene, condition, contrast = label, multiplier)
    }
    H <- cbind(1, recipe$contrasts)
    for (gene in genes) {
      for (condition in rownames(H)) {
        for (model_name in model_names) add_contrast(gene, condition, model_name,
          setNames(1, model_name), H[condition, ])
        add_contrast(gene, condition, "tau025_minus_tau050", c(tau025 = 1, tau050 = -1), H[condition, ])
        add_contrast(gene, condition, "tau025_minus_tau075", c(tau025 = 1, tau075 = -1), H[condition, ])
        add_contrast(gene, condition, "tau025_minus_mean", c(tau025 = 1, mean = -1), H[condition, ])
      }
      condition_difference <- H[2L, ] - H[1L, ]
      for (model_name in model_names) add_contrast(gene, "second_minus_first", model_name,
        setNames(1, model_name), condition_difference)
      add_contrast(gene, "second_minus_first", "tau025_minus_tau050", c(tau025 = 1, tau050 = -1), condition_difference)
      add_contrast(gene, "second_minus_first", "tau025_minus_tau075", c(tau025 = 1, tau075 = -1), condition_difference)
    }
    C <- do.call(rbind, contrasts)
    result <- rbindlist(labels)
    result[, estimate := drop(C %*% coefficients)]
    result[, se := sqrt(diag(C %*% covariance %*% t(C))) ]
    direct_se <- sqrt(colSums((centered %*% t(C))^2)) / data$n
    stopifnot(nrow(result) == 20L * length(genes), max(abs(direct_se - result$se)) < 1e-10,
      all(is.finite(result$se)), all(result$se > 0))
    critical <- qnorm((1 + config$confidence) / 2)
    family_critical <- qnorm(1 - (1 - config$confidence) / (2 * config$contrast_family_size_per_bandwidth))
    result[, `:=`(lower = estimate - critical * se, upper = estimate + critical * se,
      family_lower = estimate - family_critical * se, family_upper = estimate + family_critical * se,
      p_value = 2 * pnorm(-abs(estimate / se)), clusters = data$n,
      phase = config$phase, formal_validation_used = FALSE, manuscript_adoption = FALSE)]
    saveRDS(list(fits = fits, mean_beta = mean_beta, mean_influence = mean_influence,
      covariance = covariance, influence = influence, contrasts = C, labels = result,
      cluster_ids = names(data$rows), recipe = recipe, columns = columns, config = config),
      file.path(directory, "fit.rds"))
    fwrite(result, file.path(directory, "associations.csv"))
    write_json(list(status = "completed", formal_validation_used = FALSE, manuscript_adoption = FALSE),
      file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
    records[[length(records) + 1L]] <- result
  }
}

result <- rbindlist(records)
stopifnot(all(result[, .N, by = multiplier]$N == config$contrast_family_size_per_bandwidth))
dir.create(config$report_directory, recursive = TRUE, showWarnings = FALSE)
fwrite(result, file.path(config$report_directory, "associations.csv"))
fwrite(rbindlist(diagnostics), file.path(config$report_directory, "diagnostics.csv"))
fwrite(rbindlist(checks), file.path(config$report_directory, "reference_checks.csv"))
fwrite(rbindlist(designs), file.path(config$report_directory, "designs.csv"))
capture.output(sessionInfo(), file = file.path(root, "session.txt"))
write_json(list(status = "completed", phase = config$phase, formal_validation_used = FALSE,
  manuscript_adoption = FALSE, family_size_per_bandwidth = config$contrast_family_size_per_bandwidth),
  file.path(config$report_directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(result[multiplier == config$reference_bandwidth_multiplier & (family_lower > 0 | family_upper < 0)])

library(data.table)
library(readxl)
data <- as.data.table(read_excel(
  "data/raw/candidate_datasets_20261030/Database_MIR_A2_milk.xlsx",
  .name_repair = "minimal"))
eligible <- fread("data/interim/candidate_datasets_20261030/milk_eligible_identifiers.csv")
data <- data[match(eligible$ID_Sample, ID_Sample)]
if (anyNA(data$ID_Sample) || !identical(data$ID_Cow, as.numeric(eligible$ID_Cow)))
  stop("The eligible identifiers do not match the source records.")
if (any(data$Milk_typology != ifelse(data$Betacasein == "A2A2", "A2", "A1")))
  stop("The recorded milk type differs from the documented genotype rule.")
data[, `:=`(A2_milk = as.integer(Milk_typology == "A2"),
  farm_control = interaction(Farm, Milk_Control, drop = TRUE))]
design <- model.matrix(~ 0 + farm_control + A2_milk, data = data)
design_rank <- qr(design)$rank
if (design_rank != ncol(design)) stop("The proposed fixed covariates are not linearly independent.")
cows <- unique(data[, .(ID_Cow, Farm, Milk_typology)])
support <- cows[, .(cows = .N), by = .(Farm, Milk_typology)][order(Farm, Milk_typology)]
if (any(support[, .N, by = Farm]$N != 2L))
  stop("A farm does not contain both recorded milk types.")
bands <- as.character(240:1299)
spectra <- as.matrix(data[, ..bands])
centered <- sweep(spectra, 2L, colMeans(spectra), "-")
norms <- sqrt(colSums(centered^2))
adjacent_correlation <- colSums(centered[, -ncol(centered), drop = FALSE] *
  centered[, -1L, drop = FALSE]) / (norms[-length(norms)] * norms[-1L])
timing <- data[, .(visits = .N, interval_days = as.numeric(difftime(max(date), min(date), units = "days"))),
  by = ID_Cow]
score <- log2(data[["SCC (cell/microL)"]] / 100) + 3
frequency <- table(score)
report <- list(target_variable = "A2_milk", target_definition = "A2A2 versus A1A1 or A1A2",
  always_included_columns = colnames(design), always_included_rank = design_rank,
  spectral_variables = ncol(spectra), total_design_columns = ncol(design) + ncol(spectra),
  spectral_adjacent_correlation_quantiles = as.list(quantile(adjacent_correlation, c(0, .25, .5, .75, 1))),
  cow_counts_by_milk_type = cows[, .(cows = .N), by = Milk_typology],
  eligible_SCC_units = "cells per microliter", score_formula = "log2(SCC/100)+3",
  eligible_score_quantiles = as.list(quantile(score, c(0, .25, .5, .75, 1))),
  score_distinct_values = length(frequency), score_modal_fraction = max(frequency) / length(score),
  repeated_interval_days_quantiles = as.list(quantile(timing[visits == 2L, interval_days], c(0, .25, .5, .75, 1))),
  sampling_dates = as.character(range(as.Date(data$date))),
  model_fitting_executed = FALSE,
  audited = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
output <- "results/preflight/dataset_search_20261030"
fwrite(support, file.path(output, "milk_target_support.csv"))
jsonlite::write_json(report, file.path(output, "milk_design_assessment.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(report)

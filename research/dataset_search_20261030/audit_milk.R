library(data.table)
library(readxl)
source_file <- "data/raw/candidate_datasets_20261030/Database_MIR_A2_milk.xlsx"
data <- as.data.table(read_excel(source_file, .name_repair = "minimal"))
required <- c("ID_Sample", "ID_Cow", "Farm", "Betacasein", "Kappacasein", "Betalactoglobulin", "SCC (cell/microL)", "date", "Milk_Control")
if (!all(required %in% names(data))) stop("The observed data do not contain the documented identifiers and measurements.")
bands <- names(data)[grepl("^[0-9]+$", names(data))]
if (length(bands) != 1060L || !setequal(as.integer(bands), 240:1299))
  stop("The observed spectral grid differs from the deposited description.")
X <- as.matrix(data[, ..bands])
SCC <- data[["SCC (cell/microL)"]]
if (!is.numeric(X) || !is.numeric(SCC)) stop("The recorded spectra or SCC are not numerical measurements.")
metadata_complete <- complete.cases(data[, ..required])
data[, analysis_complete := metadata_complete & is.finite(SCC) & SCC > 0 &
  rowSums(!is.finite(X)) == 0L]
sizes <- data[, .(observations = .N, farms = uniqueN(Farm), controls = uniqueN(Milk_Control),
  dates = uniqueN(date), complete_observations = sum(analysis_complete),
  beta_values = uniqueN(Betacasein), kappa_values = uniqueN(Kappacasein),
  lacto_values = uniqueN(Betalactoglobulin)), by = ID_Cow]
sizes[, identifier_consistent := observations <= 2L & farms == 1L &
  controls == observations & dates == observations & beta_values == 1L &
  kappa_values == 1L & lacto_values == 1L]
genotypes <- rbindlist(lapply(c("Betacasein", "Kappacasein", "Betalactoglobulin"), function(field) {
  records <- unique(data[, c("ID_Cow", field), with = FALSE])
  counts <- records[, .N, by = ID_Cow]
  data.table(variable = field, cows_with_multiple_values = sum(counts$N > 1L),
    recorded_missing = sum(is.na(data[[field]])))
}))
eligible_ids <- sizes[identifier_consistent & complete_observations == observations, ID_Cow]
eligible <- data[ID_Cow %in% eligible_ids]
paired_ids <- sizes[ID_Cow %in% eligible_ids & observations == 2L, ID_Cow]
paired <- data[analysis_complete & ID_Cow %in% paired_ids]
paired_SCC <- paired[["SCC (cell/microL)"]]
phenotype <- log2(paired_SCC / 100) + 3
if (any(!is.finite(phenotype))) stop("The recorded SCC does not yield finite SCS values.")
paired_sizes <- paired[, .(observations = .N), by = ID_Cow]
output <- "results/preflight/dataset_search_20261030"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
distribution <- sizes[, .(cows = .N), by = observations][order(observations)]
cohorts <- list(eligible = eligible, paired = paired)
farms <- rbindlist(lapply(names(cohorts), function(cohort) {
  rows <- cohorts[[cohort]][, .(cows = uniqueN(ID_Cow), observations = .N), by = Farm]
  rows[, cohort := cohort]
}))
frequencies <- rbindlist(lapply(names(cohorts), function(cohort) {
  cows <- unique(cohorts[[cohort]][, .(ID_Cow, Farm, Betacasein, Kappacasein, Betalactoglobulin)])
  rbindlist(lapply(c("Betacasein", "Kappacasein", "Betalactoglobulin"), function(field) {
    counts <- cows[, .(cows = .N), by = c("Farm", field)]
    setnames(counts, field, "genotype")
    counts[, `:=`(cohort = cohort, variable = field)]
    counts
  }))
}))
flow <- data.table(stage = c("deposited", "identifier_consistent", "complete_records", "two_controls"),
  cows = c(nrow(sizes), sum(sizes$identifier_consistent), length(eligible_ids), length(paired_ids)),
  observations = c(nrow(data), sizes[identifier_consistent == TRUE, sum(observations)], nrow(eligible), nrow(paired)))
fwrite(distribution, file.path(output, "milk_cluster_sizes.csv"))
fwrite(farms, file.path(output, "milk_farm_counts.csv"))
fwrite(genotypes, file.path(output, "milk_genotype_consistency.csv"))
fwrite(frequencies, file.path(output, "milk_genotype_counts.csv"))
fwrite(flow, file.path(output, "milk_cohort_flow.csv"))
audit_directory <- "data/interim/candidate_datasets_20261030"
dir.create(audit_directory, recursive = TRUE, showWarnings = FALSE)
fwrite(sizes[identifier_consistent == FALSE], file.path(audit_directory, "milk_identifier_review.csv"))
fwrite(eligible[, .(ID_Sample, ID_Cow, Farm, Milk_Control)],
  file.path(audit_directory, "milk_eligible_identifiers.csv"))
fwrite(paired[, .(ID_Sample, ID_Cow, Farm, Milk_Control)],
  file.path(audit_directory, "milk_paired_identifiers.csv"))
report <- list(dataset_doi = "10.34810/data2338", source_file = source_file,
  source_bytes = file.info(source_file)$size, observations = nrow(data), cows = nrow(sizes),
  farms = uniqueN(data$Farm), spectral_variables = length(bands),
  duplicate_sample_identifiers = sum(duplicated(data$ID_Sample)),
  duplicate_cow_control = sum(duplicated(data[, .(ID_Cow, Milk_Control)])),
  missing_required_metadata = lapply(data[, ..required], function(x) sum(is.na(x))),
  cows_recorded_in_multiple_farms = sum(sizes$farms > 1L),
  repeated_cows_without_distinct_dates = sum(sizes$observations > 1L & sizes$dates < 2L),
  nonfinite_spectral_cells = sum(!is.finite(X)), nonpositive_SCC = sum(SCC <= 0, na.rm = TRUE),
  missing_SCC = sum(!is.finite(SCC)), constant_spectral_variables = sum(apply(X, 2L, function(x) length(unique(x)) == 1L)),
  identifier_excluded_cows = sum(!sizes$identifier_consistent),
  identifier_excluded_observations = sizes[identifier_consistent == FALSE, sum(observations)],
  eligible_cows = length(eligible_ids), eligible_observations = nrow(eligible),
  eligible_cluster_size_range = range(sizes[ID_Cow %in% eligible_ids, observations]),
  complete_repeated_cows = nrow(paired_sizes), complete_repeated_observations = nrow(paired),
  repeated_cluster_size_range = range(paired_sizes$observations),
  repeated_SCC_quantiles = as.list(quantile(paired_SCC, c(0, .25, .5, .75, 1))),
  repeated_SCS_quantiles = as.list(quantile(phenotype, c(0, .25, .5, .75, 1))),
  genotype_consistency = genotypes, license = "CC BY-NC-ND 4.0",
  audited = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
jsonlite::write_json(report, file.path(output, "milk_audit.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(output, "milk_session.txt"))
print(report)
print(frequencies[, .(cows = sum(cows)), by = .(cohort, variable, genotype)])

library(data.table)
library(readxl)
library(yaml)
config <- read_yaml("config/application_milk_20261031.yml")
source("research/score_validity_20261026/load_joint.R")
output <- config$results_directory
if (dir.exists(output)) stop("The application directory already exists.")
dir.create(output, recursive = TRUE)
data <- as.data.table(read_excel(config$source_file, .name_repair = "minimal"))
eligible <- fread(config$eligible_identifiers)
paired <- fread(config$paired_identifiers)
data <- data[match(eligible$ID_Sample, ID_Sample)]
if (anyNA(data$ID_Sample) || anyDuplicated(data$ID_Sample) ||
    !identical(data$ID_Cow, as.numeric(eligible$ID_Cow))) stop("Invalid source identifiers.")
if (any(data$Milk_typology != ifelse(data$Betacasein == "A2A2", "A2", "A1")))
  stop("The recorded milk type conflicts with the target definition.")
set.seed(config$split_seed)
ordered_cows <- sample(sort(unique(data$ID_Cow)))
data[, cluster_id := sprintf("C%05d", match(ID_Cow, ordered_cows))]
data[, farm_control := interaction(Farm, Milk_Control, drop = TRUE)]
bands <- as.character(seq.int(config$spectral_columns[1], config$spectral_columns[2]))
fixed <- model.matrix(~ 0 + farm_control, data = data)
X <- cbind(A2_milk = as.integer(data$Milk_typology == "A2"), fixed, as.matrix(data[, ..bands]))
y <- log2(data[[config$outcome_field]] / 100) + 3
Z <- matrix(1, nrow(data), 1L, dimnames = list(NULL, "cow_intercept"))
penalty_factor <- c(rep(0, 1L + ncol(fixed)), rep(1, length(bands)))
if (any(!is.finite(X)) || any(!is.finite(y))) stop("Nonfinite application measurements.")
fold_tables <- list()
cohorts <- list()
checks <- list()
for (cohort in c("eligible", "paired")) {
  index <- if (cohort == "eligible") seq_len(nrow(data)) else which(data$ID_Sample %in% paired$ID_Sample)
  current <- sb_data(y[index], X[index, , drop = FALSE], Z[index, , drop = FALSE], data$cluster_id[index])
  assignment <- pc_folds(current$rows, 2L)
  subjects <- unique(data[index, .(ID_Cow, cluster_id, Farm, Milk_typology)])
  subjects[, outer_fold := assignment[match(cluster_id, names(current$rows))]]
  fold_tables[[cohort]] <- subjects
  for (fold in 1:2) {
    training <- sb_subset(current, unlist(current$rows[assignment != fold], use.names = FALSE))
    validation <- sb_subset(current, unlist(current$rows[assignment == fold], use.names = FALSE))
    if (length(intersect(names(training$rows), names(validation$rows)))) stop("Overlapping cow roles.")
    inner <- pc_folds(training$rows, 3L)
    for (inner_fold in 1:3) {
      rows <- unlist(training$rows[inner != inner_fold], use.names = FALSE)
      if (qr(training$X[rows, penalty_factor == 0, drop = FALSE])$rank != sum(penalty_factor == 0))
        stop("The fixed covariates lack training-fold rank.")
    }
    checks[[length(checks) + 1L]] <- data.table(cohort = cohort, fold = fold,
      training_cows = training$n, estimation_cows = validation$n,
      p = ncol(X), always_included = sum(penalty_factor == 0), overlap = 0L)
  }
  cohorts[[cohort]] <- list(y = current$y, X = current$X, Z = current$Z,
    cluster_id = current$cluster_id, penalty_factor = penalty_factor, n = current$n,
    sample_id = data$ID_Sample[index], cow_id = data$ID_Cow[index])
}
registry <- rbindlist(lapply(config$settings, function(setting) {
  data.table(case_id = paste0(setting$setting, "_q", sprintf("%02d", 100 * unlist(config$quantiles))),
    setting = setting$setting, cohort = setting$cohort, tau = unlist(config$quantiles),
    lambda = setting$lambda, h_multiplier = setting$h_multiplier)
}))
registry[, `:=`(n = vapply(cohort, function(name) cohorts[[name]]$n, integer(1)),
  coordinates = config$coordinate, p = ncol(X))]
registry[, h := h_multiplier * n^(-0.3)]
saveRDS(cohorts, file.path(output, "cohorts.rds"), compress = FALSE)
saveRDS(fold_tables, file.path(output, "cow_folds.rds"))
fwrite(registry, file.path(output, "registry.csv"))
file.copy("config/application_milk_20261031.yml", file.path(output, "configuration.yml"), overwrite = FALSE)
preflight <- "results/preflight/milk_scs_20261031"
dir.create(preflight, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(checks), file.path(preflight, "input_checks.csv"))
fwrite(registry, file.path(preflight, "registry.csv"))
manifest <- list(application_id = config$application_id, method_version = config$method_version,
  n = lapply(cohorts, function(cohort) cohort$n), p = ncol(X), cases = nrow(registry),
  split_seed = config$split_seed, source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  created = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
jsonlite::write_json(manifest, file.path(preflight, "input_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
dir.create(config$backup_directory, recursive = TRUE, showWarnings = FALSE)
files <- list.files(output, full.names = TRUE)
if (!all(file.copy(files, config$backup_directory, overwrite = FALSE))) stop("Input backup failed.")
print(registry)

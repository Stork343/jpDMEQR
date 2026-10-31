library(data.table)
library(yaml)
library(mgcv)
library(splines)
config <- read_yaml("config/application_soynam_20261031.yml")
root <- config$results_directory
if (dir.exists(root)) stop("The registered application directory already exists.")
source("research/score_validity_20261026/load_joint.R")
audited <- readRDS(config$audited_panel)
observed <- as.data.table(audited$observations)
setorder(observed, Family, ID, Year)
observed[, Family := factor(Family, levels = as.character(sort(as.integer(unique(Family)))))]
observed[, Block := factor(Block)]
line_table <- unique(observed[, .(ID, Family)])
set.seed(config$split_seed)
ordered_lines <- unlist(lapply(levels(observed$Family), function(family) {
  sample(line_table[Family == family, ID])
}), use.names = FALSE)
observed[, cluster_id := sprintf("L%05d", match(ID, ordered_lines))]
line_table[, `:=`(cluster_id = sprintf("L%05d", match(ID, ordered_lines)),
  inference_fold = (match(ID, ordered_lines) - 1L) %% 2L + 1L,
  prediction_fold = (match(ID, ordered_lines) - 1L) %% config$prediction$folds + 1L)]
reference <- new.env()
load(config$reference_file, envir = reference)
checks <- as.data.table(reference$data.check.in)[, .(Year = as.integer(as.character(year)),
  Block = as.integer(BLOCK), Row = as.numeric(ROW), Col = as.numeric(COL),
  Strain = as.character(strain), check_y = yield / config$response_divisor)]
checks <- checks[Year %in% unlist(config$years) & complete.cases(checks)]
if (anyDuplicated(checks[, .(Year, Block, Row, Col)])) stop("Check plot coordinates are not unique.")
plot_keys <- observed[, .(Year, Block = as.integer(as.character(Block)), Row, Col)]
if (nrow(merge(plot_keys, checks, by = c("Year", "Block", "Row", "Col"))))
  stop("Calibration check plots overlap the analysed line plots.")
checks[, `:=`(Block = factor(Block, levels = levels(observed$Block)),
  StrainYear = factor(paste(Strain, Year, sep = "_")))]
calibration <- mgcv::gam(check_y ~ StrainYear + Block +
  te(Row, Col, by = Block, k = c(5, 5), bs = c("cr", "cr")),
  data = checks, method = "REML", select = TRUE)
prediction_data <- data.frame(Row = observed$Row, Col = observed$Col, Block = observed$Block,
  StrainYear = factor(paste("IA3023", observed$Year, sep = "_"), levels = levels(checks$StrainYear)))
if (anyNA(prediction_data$StrainYear)) stop("The reference check is absent from an observed year.")
terms <- predict(calibration, newdata = prediction_data, type = "terms")
smooth_columns <- grepl("^te\\(", colnames(terms))
if (sum(smooth_columns) != length(levels(observed$Block))) stop("The check calibration lacks a field surface.")
surface <- rowSums(terms[, smooth_columns, drop = FALSE])
if (any(!is.finite(surface)) || sd(surface) <= 1e-8) stop("The check-plot surface is not numerically identified.")
canopy <- sapply(unlist(config$years), function(year) observed$ACC / config$canopy_divisor * (observed$Year == year))
colnames(canopy) <- unlist(config$coordinates)
block <- model.matrix(~0 + Block, observed)
family <- model.matrix(~Family, observed)[, -1L, drop = FALSE]
maturity <- splines::ns(observed$DTM, knots = unlist(config$maturity_knots),
  Boundary.knots = unlist(config$maturity_boundary))
maturity_year <- do.call(cbind, lapply(unlist(config$years), function(year) {
  current <- maturity * as.numeric(observed$Year == year)
  colnames(current) <- paste0("DTM_", year, "_", seq_len(ncol(current)))
  current
}))
fixed <- cbind(canopy, block, family, maturity_year, check_surface = surface)
if (qr(fixed, tol = 1e-10)$rank != ncol(fixed)) stop("The fixed design is not linearly independent.")
spatial <- do.call(cbind, lapply(levels(observed$Block), function(field) {
  x <- (observed$Row - 10.5) / 10
  y <- (observed$Col - 40) / 40
  current <- cbind(x, y, x^2, x * y, y^2) * as.numeric(observed$Block == field)
  colnames(current) <- paste0("field_", field, "_", c("row", "column", "row2", "row_column", "column2"))
  current
}))
if (qr(cbind(fixed, spatial), tol = 1e-10)$rank != ncol(fixed) + ncol(spatial))
  stop("The spatial sensitivity design is not linearly independent.")
input <- list(observations = observed, genotypes = audited$genotypes, fixed_design = fixed,
  spatial_polynomial = spatial, years = unlist(config$years), response_divisor = config$response_divisor,
  line_folds = line_table, calibration = calibration)
dir.create(root, recursive = TRUE)
saveRDS(input, file.path(root, "input.rds"), compress = FALSE)
fwrite(line_table, file.path(root, "line_folds.csv"))
registry <- rbindlist(lapply(config$settings, function(setting) {
  data.table(case_id = paste0(setting$setting, "_q", sprintf("%02d", 100 * unlist(config$quantiles))),
    setting = setting$setting, cohort = setting$cohort, lambda = setting$lambda,
    h_multiplier = setting$h_multiplier, tau = unlist(config$quantiles))
}))
fwrite(registry, file.path(root, "registry.csv"))
if (!file.copy("config/application_soynam_20261031.yml", file.path(root, "configuration.yml")))
  stop("The configuration snapshot was not saved.")
output <- "results/preflight/soynam_20261031"
calibration_summary <- summary(calibration)
manifest <- list(lines = nrow(line_table), records = nrow(observed), genomic_markers = ncol(audited$genotypes),
  marker_by_year_columns = 2L * ncol(audited$genotypes), fixed_columns = ncol(fixed),
  total_columns = ncol(fixed) + 2L * ncol(audited$genotypes),
  check_records = nrow(checks), check_strain_years = nlevels(checks$StrainYear),
  check_model_rank = calibration$rank, check_model_edf = sum(calibration$edf),
  check_model_deviance_explained = calibration_summary$dev.expl,
  check_surface_range = range(surface), registered_inference_cases = nrow(registry),
  inference_family_size = config$family_size, response_unit = config$response_unit,
  response_divisor = config$response_divisor,
  prepared = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
jsonlite::write_json(manifest, file.path(output, "input_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
fwrite(registry, file.path(output, "registry.csv"))
fwrite(line_table[, .(lines = .N), by = .(Family, inference_fold, prediction_fold)],
  file.path(output, "fold_counts.csv"))
dir.create(config$backup_directory, recursive = TRUE, showWarnings = FALSE)
if (!all(file.copy(list.files(root, full.names = TRUE), config$backup_directory, overwrite = FALSE)))
  stop("The prepared input backup failed.")
print(manifest)

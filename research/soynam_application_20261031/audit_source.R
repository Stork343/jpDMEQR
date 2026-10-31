library(data.table)
met <- new.env()
load("data/raw/soynam_20261031/met.RData", envir = met)
reference <- new.env()
load("data/raw/soynam_20261031/soyin_reference.RData", envir = reference)
observed <- as.data.table(met$Obs)
observed[, `:=`(ID = as.character(ID), Year = as.integer(as.character(Year)))]
source <- as.data.table(reference$data.line.in)[, .(ID = as.character(strain), Year = as.integer(as.character(year)),
  Block = BLOCK, Row = ROW, Col = COL, Family = as.character(family), source_yield = yield,
  source_maturity = R8, source_canopy = AvgCC)]
keys <- c("ID", "Year", "Block", "Row", "Col")
if (anyDuplicated(observed[, ..keys]) || anyDuplicated(source[, ..keys])) stop("Field plot identifiers are not unique.")
source_panel <- source[Year %in% c(2013L, 2014L)]
observed_panel <- observed[Year %in% c(2013L, 2014L)]
plot_keys <- c("Year", "Block", "Row", "Col")
if (anyDuplicated(source_panel[, ..plot_keys]) || anyDuplicated(observed_panel[, ..plot_keys]))
  stop("The two-year physical plot identifiers are not unique.")
units <- merge(observed_panel, source_panel[, .(source_ID = ID, Year, Block, Row, Col, Family,
  source_yield, source_maturity, source_canopy)], by = plot_keys, all.x = TRUE, sort = FALSE)
if (anyNA(units$source_yield) || anyNA(units$Family) || nrow(units) != nrow(observed_panel))
  stop("Yield provenance matching is incomplete.")
family_map <- unique(units[, .(ID, Family)])
if (anyDuplicated(family_map$ID)) stop("A line has inconsistent family assignments across years.")
joined <- units
genotypes <- as.matrix(met$Gen)
if (!is.numeric(genotypes) || anyNA(genotypes) || any(!genotypes %in% c(0, 1, 2)))
  stop("Invalid deposited genotype dosage.")
if (anyNA(match(observed$ID, rownames(genotypes)))) stop("A phenotype lacks its genotype row.")
record <- list(observations = nrow(observed), lines = uniqueN(observed$ID), markers = ncol(genotypes),
  source_matched_records = nrow(units),
  yield_maximum_difference = max(abs(units$YLD - units$source_yield), na.rm = TRUE),
  yield_within_source_rounding = mean(abs(units$YLD - units$source_yield) <= 0.051, na.rm = TRUE),
  maturity_difference_records = sum(units$DTM != units$source_maturity, na.rm = TRUE),
  canopy_maximum_difference = max(abs(units$ACC - units$source_canopy), na.rm = TRUE),
  source_identifier_alias_records = sum(units$ID != units$source_ID),
  missing_genotype_calls = sum(is.na(genotypes)),
  extreme_maturity_records = sum(joined$DTM > 180),
  source_version = "NAM 1.8.0", source_commit = "e2980d8fc7be457e0a5e65abead7b96a1ce01111")
year_counts <- observed[, .(observations = .N, lines = uniqueN(ID)), by = Year]
base <- joined[Year %in% c(2013L, 2014L)]
counts <- base[, .(observations = .N, years = uniqueN(Year)), by = ID]
paired_ids <- counts[observations == 2L & years == 2L, ID]
paired <- base[ID %in% paired_ids]
if (anyDuplicated(paired[, .(ID, Year)])) stop("The two-year panel has duplicate genotype-year records.")
genotypes <- genotypes[paired_ids, , drop = FALSE]
constant <- apply(genotypes, 2L, function(column) length(unique(column)) == 1L)
duplicate <- duplicated(genotypes, MARGIN = 2L)
canonical <- genotypes - 1
orientation <- vapply(seq_len(ncol(canonical)), function(j) {
  positions <- which(canonical[, j] != 0)
  if (!length(positions)) return(1)
  sign(canonical[positions[1L], j])
}, numeric(1))
canonical <- sweep(canonical, 2L, orientation, "*")
affine_duplicate <- duplicated(canonical, MARGIN = 2L)
record$two_year_lines <- length(paired_ids)
record$two_year_observations <- nrow(paired)
record$families <- uniqueN(paired$Family)
record$constant_markers <- sum(constant)
record$duplicate_markers <- sum(duplicate)
record$complementary_marker_aliases <- sum(affine_duplicate) - sum(duplicate)
record$valid_markers <- sum(!constant & !affine_duplicate)
record$year_range <- range(paired$Year)
output <- "results/preflight/soynam_20261031"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create("data/interim/soynam_20261031", recursive = TRUE, showWarnings = FALSE)
fwrite(year_counts, file.path(output, "source_year_counts.csv"))
fwrite(data.table(marker = colnames(genotypes), constant = constant, exact_duplicate = duplicate,
  affine_duplicate = affine_duplicate),
  file.path(output, "marker_audit.csv"))
fwrite(units[DTM > 180, .(Year, DTM, source_maturity)], file.path(output, "maturity_audit.csv"))
fwrite(unique(units[ID != source_ID, .(ID, source_ID, Family)]),
  "data/interim/soynam_20261031/source_identifier_aliases.csv")
saveRDS(list(observations = paired, genotypes = genotypes[, !constant & !affine_duplicate, drop = FALSE]),
  "data/interim/soynam_20261031/audited_two_year_panel.rds", compress = FALSE)
jsonlite::write_json(record, file.path(output, "data_audit.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(record)
print(year_counts)

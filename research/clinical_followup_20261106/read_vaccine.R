library(Biobase)
library(data.table)
library(matrixStats)
library(yaml)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, args[1L] %in% c("inventory","discovery"))
folder <- "data/raw/vaccine_20261106"
array <- readRDS(file.path(folder, "all_noNorm_withResponse_eset.rds"))
serum <- as.data.table(readRDS(file.path(folder, "immdata_all.rds"))$elisa)
studies <- c("SDY67","SDY269")
metadata <- as.data.table(pData(array))[study_accession %in% studies]
metadata[, array_id := as.character(uid)]
metadata[, day := as.numeric(as.character(study_time_collected))]
metadata <- metadata[study_time_collected_unit == "Days" & day %in% c(0,3,7,28,75)]
stopifnot(!anyDuplicated(metadata$array_id), all(metadata$array_id %in% colnames(exprs(array))))
assays <- serum[study_accession %in% studies & study_time_collected_unit == "Days" & analyte == "CXL10_HUMAN"]
assays[, day := as.numeric(as.character(study_time_collected))]
baseline <- unique(assays[day==0, .(participant_id,study_accession,vaccine,gender,age_reported)])
baseline <- merge(baseline, unique(metadata[day==0, .(participant_id,study_accession)]),
                  by=c("participant_id","study_accession"))
stopifnot(!anyDuplicated(baseline$participant_id))
if (args[1L] == "inventory") {
  fwrite(baseline[order(study_accession,participant_id)], file.path(folder,"roster.csv"))
  shared <- merge(unique(metadata[, .(participant_id,study_accession,day)]),
                  unique(assays[, .(participant_id,study_accession,day)]), by=c("participant_id","study_accession","day"))
  print(shared[, .(participants=uniqueN(participant_id)), by=.(study_accession,day)][order(study_accession,day)])
  quit(save="no")
}
study <- read_yaml("results/raw/clinical_followup_settings_20261106/vaccine_cxcl10.yml")
discovery_ids <- fread(file.path(study$input_directory,"roles.csv"))[role=="discovery",patient]
metadata <- metadata[participant_id %in% discovery_ids]
assays <- serum[participant_id %in% discovery_ids & study_accession %in% studies &
                  study_time_collected_unit=="Days" & unit_reported=="pg/ml" &
                  analyte %in% c("CXL10_HUMAN","IFNG_HUMAN","IFNA2_HUMAN","CCL2_HUMAN","CCL5_HUMAN")]
assays[, day := as.numeric(as.character(study_time_collected))]
assays <- assays[day %in% c(0,3,7,28,75) & (day==0 | analyte=="CXL10_HUMAN")]
fwrite(assays[,.(patient=participant_id,study_accession,day,analyte,value_reported,value_preferred,unit_reported)],
       file.path(study$input_directory,"source_assays.csv"))
stopifnot(identical(featureNames(array),rownames(exprs(array))), !anyDuplicated(featureNames(array)))
data <- exprs(array)[,metadata$array_id,drop=FALSE]
key <- paste(metadata$participant_id,metadata$day,sep="|")
groups <- split(seq_len(ncol(data)),key)
ranks <- vapply(groups, function(ids) {
  v <- if (length(ids)==1L) data[,ids] else rowMedians(data[,ids,drop=FALSE],na.rm=TRUE)
  rank(v, ties.method="average", na.last="keep")/sum(is.finite(v))
}, numeric(nrow(data)))
rownames(ranks) <- rownames(data)
stopifnot(all(c("IFIT1","IFIT3","ISG15","MX1","OAS1","CXCL10") %in% rownames(ranks)))
columns <- data.table(sample_key=colnames(ranks), t(ranks))
fwrite(columns,file.path(study$input_directory,"ranked_expression.csv.gz"))

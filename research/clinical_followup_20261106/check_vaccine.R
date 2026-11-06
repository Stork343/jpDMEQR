library(Biobase)
library(data.table)
library(matrixStats)
library(jsonlite)
library(yaml)
study <- read_yaml("results/raw/clinical_followup_settings_20261106/vaccine_cxcl10.yml")
inputs <- study$input_directory
src <- study$source_directory
eligible <- fread(file.path(inputs,"roles.csv"))[role=="discovery",patient]
array <- readRDS(file.path(src,"all_noNorm_withResponse_eset.rds"))
p <- as.data.table(pData(array))
p[, id := paste(participant_id,study_time_collected,study_time_collected_unit,biosample_accession,sep="_")]
p <- p[participant_id %in% eligible & study_time_collected_unit=="Days",
       .(patient=participant_id,day=as.numeric(as.character(study_time_collected)),id)]
p <- p[day %in% c(0,3,7,28,75)]
stopifnot(all(p$id %in% colnames(exprs(array))),!anyDuplicated(p$id))
program <- c("MX1","ISG15","IFIT1","OAS1","IFIT3")
scores <- p[, {
  x <- exprs(array)[,id,drop=FALSE]
  v <- if (ncol(x)==1L) x[,1L] else rowMedians(x,na.rm=TRUE)
  values <- v[is.finite(v)]
  q <- v[program]
  percentile <- vapply(q, function(threshold) {
    if (!is.finite(threshold)) return(NA_real_)
    (sum(values < threshold)+(1+sum(values==threshold))/2)/length(values)
  }, numeric(1))
  list(score=10*mean(percentile))
}, by=.(patient,day)]
cytokines <- as.data.table(readRDS(file.path(src,"immdata_all.rds"))$elisa)
target <- cytokines[participant_id %in% eligible & analyte=="CXL10_HUMAN" & unit_reported=="pg/ml" &
                      study_time_collected_unit=="Days",.(patient=participant_id,day=as.integer(as.character(study_time_collected)),
                        concentration=suppressWarnings(as.numeric(as.character(value_reported))))]
target <- target[day %in% c(0,3,7,28,75) & is.finite(concentration) & concentration>=0,
                 .(concentration=median(concentration)),by=.(patient,day)]
series <- merge(target,scores,by=c("patient","day"))
baseline <- series[day==0,.(patient,baseline_score=score,baseline_c=concentration)]
ref <- merge(series[day>0],baseline,by="patient")
ref[, sample_id := paste("FLU",patient,day,sep="_")]
saved <- fread(file.path(inputs,"discovery.csv"))
features <- fread(file.path(inputs,"features.csv.gz"),select=c("sample_id","C_isg_change","C_baseline_cxcl10"))
ref <- merge(ref,merge(saved,features,by="sample_id"),by="sample_id")
stopifnot(nrow(ref)==nrow(saved), nrow(ref)==sum(series$day>0 & series$patient %in% baseline$patient),
          !anyDuplicated(ref$sample_id), all(ref$role=="discovery"))
endpoint_error <- max(abs(ref$y-log1p(ref$concentration)))
baseline_error <- max(abs(ref$C_baseline_cxcl10-log1p(ref$baseline_c)))
expected_change <- ref$score-ref$baseline_score
stopifnot(identical(is.na(expected_change),is.na(ref$C_isg_change)))
target_error <- max(abs(expected_change-ref$C_isg_change),na.rm=TRUE)
stopifnot(max(endpoint_error,baseline_error,target_error)<1e-10)
write_json(list(source="separate assay read and direct counting of all genes below each interferon gene in the published ExpressionSet",
  patients=uniqueN(ref$patient),records=nrow(ref),maximum_endpoint_difference=endpoint_error,
  maximum_contrast_difference=target_error, maximum_baseline_difference=baseline_error, reserved_outcomes_used=FALSE),
  file.path(study$report_directory,"endpoint_verification.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)

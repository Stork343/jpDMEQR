library(data.table)
library(GEOquery)
folder <- "data/raw/trialnet_20261107"
arrays <- lapply(c("metadata1.gz","metadata2.gz"),function(file)
  getGEO(filename=file.path(folder,file),getGPL=FALSE))
original <- rbindlist(lapply(arrays,function(object) {
  data <- as.data.table(pData(object))
  data.table(libid=sub(" .*","",data$title),patient=data[["subject:ch1"]],
             day=as.numeric(data[["visit day:ch1"]]),batch=data[["library batch:ch1"]])
}))
original[,study:=sub("_.*","",patient)]
quality <- fread(file.path(folder,"quality.csv"),colClasses="character")
quality[,libid:=make.unique(sub("^(lib[0-9]+).*","\\1",libId),sep="_")]
numeric_column <- function(x) {
  as.numeric(sub("%","",x,fixed=TRUE))/ifelse(grepl("%",x,fixed=TRUE),100,1)
}
quality[,reads:=numeric_column(PF_HQ_ALIGNED_READS)]
quality[,aligned:=numeric_column(mapped_reads_w_dups)]
quality[,coverage:=numeric_column(MEDIAN_CV_COVERAGE)]
thresholds <- data.table(study=c("ABATE","START","T1DAL","TN02","TN05","TN09"),
  reads_cut=c(1e6,1e6,1e6,1.5e6,5e6,1.5e6),aligned_cut=c(.9,.8,.8,.8,.8,.8),
  coverage_cut=c(1,1,1,.7,1,.8))
duplicate_libraries <- c("lib1473","lib1476","lib1477","lib1478","lib4899","lib4976")
controls <- merge(original[!libid %in% duplicate_libraries],
                  quality[,.(libid,reads,aligned,coverage)],by="libid")
controls <- merge(controls,thresholds,by="study")
controls <- controls[reads>=reads_cut & aligned>=aligned_cut & coverage<=coverage_cut]
stopifnot(nrow(controls)==471L,!anyDuplicated(controls$libid))
controls <- controls[day>=-30 & day<=30,.(patient,study,libid,day,batch)]
controls[,treatment:="placebo"]
controls[,source:="controls"]
clinical <- fread(file.path(folder,"patients.csv"),select=c("study","participant_id","age_years","sex"))
setnames(clinical,"participant_id","patient")
controls <- merge(controls,clinical,by=c("patient","study"),all.x=TRUE)

rituximab <- fread(file.path(folder,"tn05_rna.csv"))
rituximab <- rituximab[treatment=="rituximab" & rnaseq_baseline_visit=="Y"]
clinical_rtx <- fread(file.path(folder,"tn05_patients.csv"),select=c("participant_id","age_years","sex"))
rituximab <- merge(rituximab,clinical_rtx,by="participant_id")
rituximab <- rituximab[,.(patient=participant_id,study="TN05",libid,day=7*rnaseq_visit_weeks,
                         batch="TN05",treatment,source="rituximab",age_years,sex)]
abatacept <- as.data.table(read.delim(file.path(folder,"tn09_rna.txt"),row.names=1,
                                     check.names=FALSE,colClasses=c(id="character")))
abatacept <- abatacept[treatment=="active" & Day>=-30 & Day<=30]
clinical_aba <- fread(file.path(folder,"tn09_clinical.txt"),select=c("Random_ID","Gender"),
                      colClasses=c(Random_ID="character"))
stopifnot(!anyDuplicated(clinical_aba$Random_ID))
abatacept <- merge(abatacept,clinical_aba,by.x="id",by.y="Random_ID",all.x=TRUE)
abatacept <- abatacept[,.(patient=paste0("TN09_",id),study="TN09",libid=libID,day=Day,
                         batch="TN09",treatment="abatacept",source="abatacept",age_years=age,sex=Gender)]
all <- rbindlist(list(controls,rituximab,abatacept),use.names=TRUE)
stopifnot(all(is.finite(all$day)),!anyDuplicated(all$libid),
          all[,all(uniqueN(treatment)==1L),by=patient]$V1)
fwrite(all,file.path(folder,"baseline_rna.csv"))
roster <- unique(all[,.(patient,study,treatment,source,age_years,sex)])
stopifnot(!anyDuplicated(roster$patient))
fwrite(roster[order(patient)],file.path(folder,"roster.csv"))
print(roster[,.(patients=.N),by=.(study,treatment)])

library(data.table)
library(yaml)
library(matrixStats)
source("research/empirical_selection_20261101/prediction_core.R")
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>0L)
records <- shifts <- list()
for (name in args) {
  folder <- "results/raw/clinical_followup_settings_20261106"
  study <- read_yaml(file.path(folder,paste0(name,".yml")))
  method <- read_yaml(file.path(folder,paste0("cq_",name,".yml")))
  for (fold in seq_len(study$outer_folds)) {
    outer <- file.path(study$output_root,paste0("outer_",fold))
    training <- fread(file.path(outer,"training.csv"))
    held <- fread(file.path(outer,"evaluation.csv"))
    stopifnot(all(training$role=="training"),all(held$role=="discovery_evaluation"),
              !length(intersect(training$strain_id,held$strain_id)))
    quantiles <- as.numeric(unlist(method$quantiles))
    original <- adjusted <- matrix(NA_real_,nrow(held),length(quantiles))
    for (j in seq_along(quantiles)) {
      tau <- quantiles[j]
      path <- file.path(method$output_root,study$study_id,paste0("outer_",fold),sprintf("tau%03d",round(tau*100)))
      fit <- readRDS(file.path(path,"fit.rds"))
      predicted <- fread(file.path(path,"prediction.csv"))
      stopifnot(identical(predicted$sample_id,held$sample_id),length(fit$oof)==nrow(training),
                identical(fit$recipe$training_strains,sort(unique(training$strain_id))))
      correction <- ep_calibrate(training,fit$oof,tau,"strain_id","condition")
      original[,j] <- predicted$prediction
      adjusted[,j] <- original[,j]+correction[held$condition]
      shifts[[length(shifts)+1L]] <- data.table(study_id=study$study_id,fold,tau,correction=unname(correction))
    }
    for (version in c("registered","training_oof_intercept")) {
      matrix <- if (version=="registered") original else adjusted
      matrix <- t(apply(matrix,1L,sort))
      matrix[matrix<0] <- 0
      for (j in seq_along(quantiles)) {
        tau <- quantiles[j]
        records[[length(records)+1L]] <- data.table(study_id=study$study_id,version,fold,tau,
          patient=held$strain_id,loss=ep_loss(held$y,matrix[,j],tau),
          below=as.numeric(held$y<matrix[,j]),at_or_below=as.numeric(held$y<=matrix[,j]))
      }
    }
  }
}
patient <- rbindlist(records)[,lapply(.SD,mean),by=.(study_id,version,tau,patient),.SDcols=c("loss","below","at_or_below")]
summary <- patient[,.(patients=.N,loss=mean(loss),below=mean(below),at_or_below=mean(at_or_below)),by=.(study_id,version,tau)]
summary[,c("stage","formal_validation_used","manuscript_adoption"):=list("discovery_method_diagnostic",FALSE,FALSE)]
output <- "results/scratch/clinical_followup_20261106"
fwrite(summary,file.path(output,"calibration_diagnostic.csv"))
fwrite(rbindlist(shifts),file.path(output,"calibration_shifts.csv"))
print(summary)

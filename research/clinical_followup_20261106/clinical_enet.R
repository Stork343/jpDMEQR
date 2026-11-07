library(data.table)
library(yaml)
library(jsonlite)
library(glmnet)
library(matrixStats)
source("research/clinical_search_20261105/panel.R")
source("research/empirical_selection_20261101/prediction_core.R")
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
study <- read_yaml(args[1L])
stopifnot(isTRUE(read_json(file.path(study$report_directory,"decision.json"))$necessary_pass))
panel <- cp_panel(study$input_directory)
panel$X <- panel$X[,startsWith(colnames(panel$X),"C_"),drop=FALSE]
tau <- study$decision_quantile
predictions <- list()
for (fold in seq_len(study$outer_folds)) {
  outer <- file.path(study$output_root,paste0("outer_",fold))
  rows <- fread(file.path(outer,"training.csv"))
  held <- fread(file.path(outer,"evaluation.csv"))
  stopifnot(all(rows$role=="training"),all(held$role=="discovery_evaluation"),
            !length(intersect(rows$strain_id,held$strain_id)))
  directory <- file.path(outer,"models","CLINICAL_ENET_all")
  dir.create(directory,showWarnings=FALSE)
  fit_file <- file.path(directory,"fit.rds")
  if (!file.exists(fit_file)) {
    feature_names <- colnames(panel$X)
    stopifnot(all(startsWith(feature_names,"C_")))
    progress <- function(stage, details=list()) cat(fold,stage,"\n")
    fit <- ep_fit_prediction(rows,held,panel,"ENET",NA_real_,study,progress,directory,
              cp_recipe,cp_features,"strain_id","condition",ep_predict,"oof")
    saveRDS(fit,paste0(fit_file,".partial"))
    stopifnot(file.rename(paste0(fit_file,".partial"),fit_file))
  }
  fit <- readRDS(fit_file)
  stopifnot(identical(sort(unique(rows$strain_id)),fit$recipe$training_strains))
  x <- cp_features(held,panel,fit$recipe)
  base <- drop(ep_predict(x,fit$beta))
  rebuilt <- vapply(fit$calibration,function(offset) base+offset[held$condition],numeric(nrow(held)))
  stopifnot(max(abs(rebuilt-fit$prediction))<1e-9)
  stored <- rbindlist(lapply(seq_along(fit$quantiles),function(j)
    data.table(sample_id=held$sample_id,strain_id=held$strain_id,method="CLINICAL_ENET",
               tau=fit$quantiles[j],prediction=rebuilt[,j])))
  fwrite(stored,file.path(directory,"prediction.csv"))
  write_json(list(status="completed",stage="verified_discovery_prediction",formal_validation_used=FALSE),
             file.path(directory,"status.json"),auto_unbox=TRUE,pretty=TRUE)
  selected <- pmax(0,apply(rebuilt,1,sort)[match(tau,fit$quantiles),])
  predictions[[fold]] <- data.table(patient=held$strain_id,loss=ep_loss(held$y,selected,tau))
}
patient <- rbindlist(predictions)[,.(loss=mean(loss)),by=patient]
method <- read_yaml(file.path(dirname(args[1L]),paste0("cq_",basename(args[1L]))))
saved <- fread(file.path(method$report_directory,"discovery_comparisons.csv"))[decision_quantile==TRUE]
loss <- unique(saved$current_loss)
stopifnot(length(loss)==1L,nrow(patient)==read_json(file.path(study$report_directory,"decision.json"))$patients)
gain <- 1-loss/mean(patient$loss)
write_json(list(method="CLINICAL_ENET",patients=nrow(patient),tau=tau,loss=mean(patient$loss),cq_loss=loss,
                gain=gain,passed=gain>=study$screening_minimum_relative_gain,formal_validation_used=FALSE),
           file.path(study$report_directory,"clinical_enet.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)

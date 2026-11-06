library(data.table)
library(Matrix)
library(yaml)
source("research/clinical_search_20261105/panel.R")
source("research/empirical_selection_20261101/prediction_core.R")
settings <- "results/raw/clinical_followup_settings_20261106"
studies <- c("ispy_diameter","pw_ispy2_volume_20261106","pw_ra_westchina_20261105","physionet_ph")
records <- list()
for (name in studies) {
  cfg <- read_yaml(file.path(settings,paste0(name,".yml")))
  outer <- file.path(cfg$output_root,"outer_1")
  original <- readRDS(file.path(outer,"models","ENET_all","inner_1.rds"))
  source_rows <- fread(file.path(outer,"training.csv"))
  stopifnot(all(source_rows$role=="training"))
  panel <- cp_panel(cfg$input_directory)
  design <- cp_features(source_rows[inner_fold!=1],panel,original$recipe)
  test <- cp_features(source_rows[inner_fold==1],panel,original$recipe)
  started <- proc.time()[["elapsed"]]
  fit <- ep_mean_path(design,0,cfg)
  elapsed <- proc.time()[["elapsed"]]-started
  indexes <- which(original$path$grid$alpha==0)
  stopifnot(length(indexes)==cfg$elastic_net_path_length)
  beta_error <- max(abs(as.matrix(fit$beta)-as.matrix(original$path$beta[,indexes])))
  predictions <- as.matrix(ep_predict(test,fit$beta))
  prediction_error <- max(abs(predictions-original$prediction[,indexes]))
  grid_error <- max(abs(fit$lambda-original$path$diagnostics[[1L]]$lambda))
  stopifnot(is.finite(grid_error),grid_error<1e-12,beta_error<1e-7,prediction_error<1e-7)
  records[[name]] <- data.table(study_id=cfg$study_id, rows=nrow(design$X),
    penalized_variables=sum(design$penalty_factor!=0),lambda_values=length(indexes),
    system_dimension=fit$diagnostics$system_dimension,solver=fit$diagnostics$solver,
    coefficient_difference=beta_error,prediction_difference=prediction_error,lambda_difference=grid_error,
    maximum_system_residual=max(fit$diagnostics$linear_system_error),seconds=elapsed,reserved_outcomes_used=FALSE)
}
out <- rbindlist(records)
stopifnot(uniqueN(out$solver)==2L)
fwrite(out,"results/preflight/clinical_followup_20261106/ridge_equivalence.csv")
print(out)

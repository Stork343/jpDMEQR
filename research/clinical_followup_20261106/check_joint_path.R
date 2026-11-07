library(data.table)
library(yaml)
library(jsonlite)
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
source("research/clinical_search_20261105/panel.R")
Rcpp::sourceCpp("research/clinical_followup_20261106/intercept_qp.cpp",
               cacheDir="results/scratch/clinical_followup_20261106/intercept_cache",rebuild=FALSE)
cfg <- read_yaml("results/raw/clinical_followup_settings_20261106/physionet_ph.yml")
outer <- file.path(cfg$output_root,"outer_1")
rows <- fread(file.path(outer,"training.csv"))
recipe <- readRDS(file.path(outer,"recipes.rds"))$inner_1$direction_recipe$base_recipe
data <- cp_panel(cfg$input_directory)
fitting <- cp_features(rows[inner_fold!=1],data,recipe)
evaluation <- cp_features(rows[inner_fold==1],data,recipe)
rank <- qr(fitting$X*sqrt(fitting$weight/fitting$n),tol=1e-5)$rank
maximum <- min(recipe$search_dimension,floor(fitting$n/log(fitting$n)),rank-sum(fitting$penalty_factor==0))
tau <- cfg$decision_quantile
Lambda <- matrix(min(unlist(cfg$working_penalties)))
prepared <- pc_prepare_initial(fitting$y,fitting$X,fitting$Z,fitting$cluster_id,tau,Lambda,
                               fitting$penalty_factor,fit_always_included=FALSE)
solver <- function(X,Z,outcome,cluster,inverse_size,Lambda,coordinate_penalty,quantile,
                    theta,gamma,dual,max_iterations=30000L,gap_tolerance=1e-8) {
  clinical_intercept_qp(X,Z,outcome,cluster,inverse_size,Lambda,coordinate_penalty,quantile,
                       theta,gamma,dual,max_iterations,min(gap_tolerance,1e-10))
}
started <- proc.time()[["elapsed"]]
fast <- pc_orthogonal_path(prepared,tau,0:maximum,solver=solver,gap_tolerance=1e-8)
fast_seconds <- proc.time()[["elapsed"]]-started
saveRDS(fast,"results/scratch/clinical_followup_20261106/joint_fast_path.rds")
cat("Certified path",fast_seconds,"seconds",maximum,"extensions\n")
flush.console()
started <- proc.time()[["elapsed"]]
original <- pc_orthogonal_path(prepared,tau,0:maximum,solver=jp_joint_qp,gap_tolerance=1e-8)
original_seconds <- proc.time()[["elapsed"]]-started
saveRDS(original,"results/scratch/clinical_followup_20261106/joint_reference_path.rds")
error <- max(abs(evaluation$X%*%(fast$coefficients-original$coefficients)))*recipe$y_scale
same <- identical(fast$selected_path,original$selected_path)
checks <- list(study_id=cfg$study_id,training_patients=fitting$n,evaluation_patients=evaluation$n,
  tau=tau,lambda=Lambda[1],columns=ncol(fitting$X),design_rank=rank,extensions=maximum,
  identical_selected_paths=same,maximum_evaluation_prediction_difference=error,
  fast_seconds=fast_seconds,reference_seconds=original_seconds,reserved_outcomes_used=FALSE)
write_json(checks,"results/preflight/clinical_followup_20261106/joint_path_equivalence.json",
           auto_unbox=TRUE,pretty=TRUE,digits=NA)
stopifnot(same,error<1e-4)

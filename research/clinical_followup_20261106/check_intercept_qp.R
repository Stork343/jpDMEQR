library(data.table)
library(yaml)
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
source("research/clinical_search_20261105/panel.R")
cache <- "results/scratch/clinical_followup_20261106/intercept_cache"
dir.create(cache,showWarnings=FALSE)
Rcpp::sourceCpp("research/clinical_followup_20261106/intercept_qp.cpp",cacheDir=cache,rebuild=FALSE)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
cfg <- read_yaml(args[1L])
outer <- file.path(cfg$output_root,"outer_1")
rows <- fread(file.path(outer,"training.csv"))[inner_fold!=1]
recipe <- readRDS(file.path(outer,"recipes.rds"))$inner_1$direction_recipe$base_recipe
features <- cp_features(rows,cp_panel(cfg$input_directory),recipe)
stopifnot(all(rows$role=="training"),identical(sort(unique(rows$strain_id)),recipe$training_strains))
group <- match(rows$strain_id,sort(unique(rows$strain_id)))
patients <- split(seq_len(nrow(rows)),group)
square_weight <- sqrt(features$weight/features$n)
weighted <- features$X*square_weight
full <- qr(weighted,tol=1e-8)
dimensions <- unique(pmin(c(12L,80L,160L),full$rank))
cases <- CJ(dimension=dimensions,tau=unlist(cfg$quantiles),lambda=unlist(cfg$working_penalties))
records <- list()
out <- "results/preflight/clinical_followup_20261106/intercept_qp_checks.csv"
for (j in seq_len(nrow(cases))) {
  case <- cases[j]
  columns <- full$pivot[seq_len(case$dimension)]
  basis <- qr.Q(qr(weighted[,columns,drop=FALSE]))
  design <- basis/square_weight
  initial <- drop(crossprod(basis,features$y*square_weight))
  nuisance <- pc_exact_profile(features$y-drop(design%*%initial),features$Z,patients,matrix(case$lambda),case$tau)
  arguments <- list(design,features$Z,features$y,group-1L,features$weight,matrix(case$lambda),
                    numeric(ncol(design)),case$tau,initial,nuisance$gamma,nuisance$dual,
                    max_iterations=100000L,gap_tolerance=1e-10)
  clock <- proc.time()[["elapsed"]]
  conic <- do.call(jp_joint_qp,arguments)
  conic_seconds <- proc.time()[["elapsed"]]-clock
  clock <- proc.time()[["elapsed"]]
  fitted <- do.call(clinical_intercept_qp,arguments)
  seconds <- proc.time()[["elapsed"]]-clock
  objective_error <- abs(fitted$primal-conic$primal)
  prediction_error <- max(abs(design%*%(fitted$theta-conic$theta)))
  stopifnot(objective_error<1e-8,prediction_error<1e-4,
            fitted$duality_gap<1e-10*max(1,abs(fitted$primal)))
  records[[j]] <- data.table(study_id=cfg$study_id,patients=features$n,records=nrow(rows),
    minimum_visits=min(lengths(patients)),maximum_visits=max(lengths(patients)),dimension=case$dimension,
    tau=case$tau,lambda=case$lambda,objective_error,prediction_error,conic_seconds,seconds,
    gap=fitted$duality_gap,iterations=fitted$iterations,polished=fitted$polished,reserved_outcomes_used=FALSE)
  fwrite(rbindlist(records),out)
  print(records[[j]])
  flush.console()
}

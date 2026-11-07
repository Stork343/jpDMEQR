library(data.table)
library(yaml)
source("research/clinical_search_20261105/panel.R")
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
study <- read_yaml(args[1L])
panel <- cp_panel(study$input_directory)
records <- list()
for (fold in seq_len(study$outer_folds)) {
  outer <- file.path(study$output_root,paste0("outer_",fold))
  rows <- fread(file.path(outer,"training.csv"))
  recipes <- readRDS(file.path(outer,"recipes.rds"))
  for (key in names(recipes)) {
    current <- if (key=="training") rows else rows[inner_fold!=as.integer(sub("inner_","",key))]
    recipe <- recipes[[key]]$direction_recipe$base_recipe
    features <- cp_features(current,panel,recipe)
    weighted <- features$X*sqrt(features$weight/features$n)
    rank <- qr(weighted,tol=1e-5)$rank
    singular <- svd(weighted,nu=0,nv=0)$d
    stopifnot(rank==sum(singular>1e-5),rank>=sum(features$penalty_factor==0))
    records[[length(records)+1L]] <- data.table(fold,panel=key,patients=features$n,
        columns=ncol(features$X),independent_directions=rank,
        estimable_extensions=rank-sum(features$penalty_factor==0),
        declared_dimension=recipe$search_dimension,registered_size=floor(features$n/log(features$n)))
  }
}
result <- rbindlist(records)
result[,search_limit:=min(estimable_extensions,declared_dimension,registered_size),by=fold]
fwrite(result,file.path(study$report_directory,"search_rank.csv"))
print(result[,.(search_limit=unique(search_limit)),by=fold])

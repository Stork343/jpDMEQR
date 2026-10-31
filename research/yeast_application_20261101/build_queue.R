library(data.table)
library(yaml)

config <- read_yaml("config/discovery_yeast_20261101.yml")
data <- readRDS("data/derived/yeast_20261101/growth_discovery.rds")
traits <- setdiff(sort(unique(data$trait)), config$excluded_control_traits)
traits <- c(config$initial_trait, setdiff(traits, config$initial_trait))
jobs <- rbindlist(lapply(traits, function(trait_name) {
  current <- data[trait == trait_name]
  stopifnot(uniqueN(current[exploration_role == "training", strain_id]) >= 100L,
            uniqueN(current[exploration_role == "evaluation", strain_id]) >= 50L)
  data.table(trait = trait_name, method = c(rep("JOINT", 3L), "ENET"),
             tau = c(unlist(config$quantiles), NA_real_),
             run_id = c("qp_tau025", "qp_tau050", "qp_tau075", "cv"))
}))
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path("results/raw/application/yeast_20261101/discovery", trait, paste0(method, "_", run_id))]
jobs[, status := vapply(directory, function(path) if (file.exists(file.path(path, "status.json"))) "verify_existing" else "pending", character(1))]
directory <- "results/raw/application/yeast_20261101/screening_20261101"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(directory, "jobs.csv"))) stop("计算清单已经存在")
fwrite(jobs, file.path(directory, "jobs.csv"))
cat(length(traits), "conditions", nrow(jobs), "tasks\n")

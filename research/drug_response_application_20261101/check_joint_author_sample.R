library(jsonlite)
source("research/score_validity_20261026/load_joint.R")

output_dir <- "results/preflight/empirical_selection_20261101/joint_author_example"
raw_dir <- "results/raw/application/drug_response_20261101/joint_author_example"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(raw_dir, "fit.rds"))) stop("结果文件已经存在")
objects <- new.env()
load("results/scratch/empirical_selection_20261101/fbioselect_source/test_gdsc_data.RData", objects)
fd <- objects$fd
dose <- splines::bs(fd$DosageSQ, df = 4L, degree = 3L, intercept = FALSE)
colnames(dose) <- paste0("dose", seq_len(ncol(dose)))
low <- cbind(model.matrix(~ DRUG_NAME, fd), dose)
genes <- as.matrix(fd[, -(1:5)])
X <- cbind(low, genes)
request <- list(
  role = "numerical_feasibility_only",
  method = "joint_profile_design_20261026",
  tau = 0.5, Lambda = 1, seed = 20261102L,
  experimental_ids = length(unique(fd$ID)), rows = nrow(fd), columns = ncol(X),
  coordinates = colnames(genes)[1L], manuscript_adoption = FALSE
)
write_json(request, file.path(output_dir, "request.json"), auto_unbox = TRUE, pretty = TRUE)
set.seed(request$seed)
fit <- ds_fit_profile(fd$Response, X, matrix(1, nrow(fd), 1L), fd$ID,
                      tau = request$tau, coordinates = request$coordinates,
                      Lambda = matrix(request$Lambda, 1L, 1L),
                      penalty_factor = c(rep(0, ncol(low)), rep(1, ncol(genes))))
stopifnot(all(is.finite(fit$table$beta_tilde)), all(is.finite(fit$table$se)),
          all(fit$table$max_kkt < 1e-8))
saveRDS(fit, file.path(raw_dir, "fit.rds"))
report <- c(request, list(status = "completed", timing = as.list(fit$timing),
                         max_kkt = max(fit$table$max_kkt),
                         validation_comparison_completed = FALSE))
write_json(report, file.path(output_dir, "execution.json"), auto_unbox = TRUE, pretty = TRUE)
print(report[c("status", "experimental_ids", "rows", "columns", "timing", "max_kkt")])

source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/fixed_inference.R")
args <- commandArgs(trailingOnly = TRUE)
job_id <- as.integer(args[1L])
input <- "results/scratch/joint_profile_manuscript_20261026/module_inputs"
jobs <- read.csv(file.path(input, "jobs.csv"), stringsAsFactors = FALSE)
if (!job_id %in% jobs$job_id) stop("Invalid module job.")
job <- jobs[jobs$job_id == job_id, ]
output <- file.path("results/raw/application/joint_profile_20261026", sprintf("job_%03d", job_id))
if (dir.exists(output)) stop("This module job already exists.")
dir.create(output, recursive = TRUE)
saveRDS(list(job = job, input = input, method_version = "joint_profile_fixed_design_20261026",
             started = Sys.time(), session = sessionInfo()), file.path(output, "request.rds"))
globalCallingHandlers(error = function(condition) saveRDS(list(job = job, condition = condition), file.path(output, "failure.rds")))
stages <- readRDS(file.path(input, paste0("stages_", job$model, ".rds")))
n <- length(stages) - 1L
h <- job$h_multiplier * n^(-0.3)
estimates <- matrix(NA_real_, n, ncol(stages[[1L]]$X), dimnames = list(NULL, colnames(stages[[1L]]$X)))
diagnostics <- vector("list", length(stages))
full <- NULL
for (i in seq_along(stages)) {
  stage <- stages[[i]]
  Z <- if (job$random_slope) cbind(intercept = 1, centered_time = stage$time) else stage$Z
  result <- jp_fit_fixed(stage$y, stage$X, Z, stage$cluster_id, job$tau,
                          Lambda = diag(job$lambda_gamma, ncol(Z)), h = h,
                          initial = if (is.null(full)) NULL else full$beta)
  if (i == 1L) full <- result else estimates[i - 1L, ] <- result$beta
  diagnostics[[i]] <- cbind(data.frame(omitted_subject = stage$omitted_subject), result$diagnostics)
}
if (any(!is.finite(estimates))) stop("Incomplete subject deletion estimates.")
centered <- sweep(estimates, 2L, colMeans(estimates), "-")
workflow_se <- sqrt((n - 1) / n * colSums(centered^2))
table <- full$table
table$model <- job$model
table$tau <- job$tau
table$specification <- job$specification
table$lambda_gamma <- job$lambda_gamma
table$h_multiplier <- job$h_multiplier
table$random_slope <- job$random_slope
table$workflow_se <- workflow_se
table$workflow_lower <- table$beta_hat - qnorm(.975) * table$workflow_se
table$workflow_upper <- table$beta_hat + qnorm(.975) * table$workflow_se
family_size <- if (job$tau == .1) 10L else 30L
critical <- qnorm(1 - .05 / (2 * family_size))
table$family_size <- family_size
table$simultaneous_lower <- table$beta_hat - critical * table$se
table$simultaneous_upper <- table$beta_hat + critical * table$se
table$workflow_simultaneous_lower <- table$beta_hat - critical * table$workflow_se
table$workflow_simultaneous_upper <- table$beta_hat + critical * table$workflow_se
saveRDS(list(job = job, table = table, full = full, leave_one_estimates = estimates,
             diagnostics = do.call(rbind, diagnostics), completed = Sys.time()), file.path(output, "result.rds"))
write.csv(table, file.path(output, "estimates.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(output, "diagnostics.csv"), row.names = FALSE)
cat("Completed module job", job_id, "\n")

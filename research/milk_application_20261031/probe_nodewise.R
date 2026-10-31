source("research/score_validity_20261026/load_joint.R")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("Provide the Gaussian solver type and iteration limit.")
cohort <- readRDS("results/raw/application/milk_scs_20261031/cohorts.rds")$eligible
data <- sb_data(cohort$y, cohort$X, cohort$Z, cohort$cluster_id)
outer <- pc_folds(data$rows, 2L)
training <- sb_subset(data, unlist(data$rows[outer != 1L], use.names = FALSE))
factor <- training$X * sqrt(training$weight)
scale <- sqrt(colSums(factor^2) / training$n)
standardised <- sweep(factor, 2L, scale, "/")
folds <- pc_folds(training$rows, 5L)
indices <- unlist(training$rows[folds != 1L], use.names = FALSE)
clusters <- sum(folds != 1L)
design <- standardised[indices, -1L, drop = FALSE]
response <- standardised[indices, 1L]
ratios <- pc_ratio_path(training$n, 24L)
penalties <- max(abs(drop(crossprod(design, response)))) / clusters * ratios
started <- proc.time()[[3L]]
fit <- glmnet::glmnet(design * sqrt(length(indices) / clusters),
  response * sqrt(length(indices) / clusters), family = "gaussian", alpha = 1,
  lambda = penalties, intercept = FALSE, standardize = FALSE,
  thresh = 1e-13, maxit = as.integer(arguments[2L]), type.gaussian = arguments[1L])
elapsed <- proc.time()[[3L]] - started
coefficients <- as.matrix(fit$beta)
gradient <- -crossprod(design, response - design %*% coefficients) / clusters
defects <- vapply(seq_len(ncol(coefficients)), function(j) {
  beta <- coefficients[, j]
  max(ifelse(abs(beta) > 1e-10, abs(gradient[, j] + penalties[j] * sign(beta)),
    pmax(abs(gradient[, j]) - penalties[j], 0)))
}, numeric(1))
print(list(type = arguments[1L], seconds = elapsed, status = fit$jerr,
  fitted_lambdas = ncol(coefficients), max_kkt = max(defects),
  df = fit$df, defects = defects))
saveRDS(list(fit = fit, elapsed = elapsed, defects = defects),
  paste0("results/scratch/dataset_search_20261030/nodewise_", arguments[1L], ".rds"))
if (fit$jerr != 0L || ncol(coefficients) != length(penalties) || max(defects) > 1e-6)
  stop("The actual design path did not pass its convergence check.")

library(jsonlite)
library(fbioSelectr)
source("research/drug_response_application_20261101/fbio_adapter.R")

objects <- new.env()
load("results/scratch/empirical_selection_20261101/fbioselect_source/test_gdsc_data.RData", objects)
fd <- objects$fd
genes <- as.matrix(fd[, -(1:5)])
Z <- model.matrix(~ DRUG_NAME, fd)
threshold <- floor(length(unique(fd$ID)) / log(length(unique(fd$ID))))
model <- drug_fbio_fit(fd$Response, genes, fd$DosageSQ, Z, as.character(fd$ID),
                       weights = fd$weights, threshold = threshold)
prediction <- drug_fbio_predict(model, genes, fd$DosageSQ, Z)
chosen <- as.character(head(model$fit$gene.ranking.marginal.screening$hdcov, threshold))
author_basis <- fbioSelectr::BasisFuncApproximation(cbind(Z, genes[, chosen, drop = FALSE]),
  fd$DosageSQ, degrFreedom = 5L, degree = 3L, intercept = TRUE, knots = NULL)
author_prediction <- drop(predict(model$fit$pen.regression.fit,
  X = author_basis[, -1L, drop = FALSE], type = "response"))
prediction_error <- max(abs(prediction - author_prediction))
group <- c(rep(attr(Z, "assign") + 1L, each = 5L),
           rep(seq.int(ncol(Z) + 1L, ncol(Z) + threshold), each = 5L))
reference_cv <- grpreg::cv.grpreg(author_basis[, -1L, drop = FALSE], fd$Response,
  group[-1L], nfolds = 5L, fold = model$fold, penalty = "grSCAD", weights = fd$weights)
lambda_error <- abs(model$lambda - reference_cv$lambda.min)
fold_integrity <- all(vapply(split(model$fold, as.character(fd$ID)), function(value) length(unique(value)) == 1L, logical(1)))
repeat_model <- drug_fbio_fit(fd$Response, genes, fd$DosageSQ, Z, as.character(fd$ID),
                              weights = fd$weights, threshold = threshold)
determinism_error <- max(abs(prediction - drug_fbio_predict(repeat_model, genes, fd$DosageSQ, Z)))

# 使用实际观测检查零惩罚下的最小二乘结果。
small <- cbind(model$basis[, -1L, drop = FALSE], genes[, 1:2, drop = FALSE])
zero <- grpreg::grpreg(small, fd$Response, group = seq_len(ncol(small)),
                        penalty = "grSCAD", lambda = c(1e-4, 0), eps = 1e-12)
zero_prediction <- drop(predict(zero, X = small, type = "response", lambda = 0))
least_squares <- lm.fit(cbind(1, small), fd$Response)$fitted.values
limit_error <- max(abs(zero_prediction - least_squares))
stopifnot(prediction_error < 1e-10, lambda_error < 1e-12,
          determinism_error < 1e-10, fold_integrity, limit_error < 1e-6)
report <- list(
  status = "passed", role = "adapter_fidelity",
  prediction_identity_error = prediction_error,
  penalty_identity_error = lambda_error,
  deterministic_error = determinism_error,
  zero_penalty_prediction_error = limit_error,
  experimental_id_fold_integrity = fold_integrity,
  native_sample_rows = nrow(fd), native_sample_genes = ncol(genes),
  formal_comparison_completed = FALSE,
  package_version = as.character(packageVersion("fbioSelectr")),
  grpreg_version = as.character(packageVersion("grpreg"))
)
output_dir <- "results/preflight/empirical_selection_20261101/fbio_adapter"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_json(report, file.path(output_dir, "checks.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(report)

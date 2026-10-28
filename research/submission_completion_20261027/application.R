submission_application_base <- function(processing, repetition, fold, inner_fold, metadata) {
  suffix <- if (inner_fold == 0L) "outer" else sprintf("inner%02d", inner_fold)
  key <- sprintf("prediction_r%02d_f%02d_%s", repetition, fold, suffix)
  path <- if (processing == "reference") {
    file.path("results/scratch/submission_application_20261027/preprocessing", key, "prepared_base.rds")
  } else file.path("data/derived/GSE65391/processing_20261023_v1", processing, "preprocessing_cache", key, "prepared_base.rds")
  base <- readRDS(path)
  if (!identical(base$sample_ids, rownames(metadata))) stop("The preprocessing sample identifiers differ.")
  base
}

submission_application_design <- function(base, metadata, model, tau, screen_dimension) {
  clinical <- base$clinical[[model]]$X
  training <- base$training_rows
  columns <- which(vapply(seq_len(ncol(clinical)), function(j) length(unique(clinical[training, j])) > 3L, logical(1)))
  centers <- colMeans(clinical[training, columns, drop = FALSE])
  scales <- apply(clinical[training, columns, drop = FALSE], 2L, sd)
  if (any(!is.finite(scales)) || any(scales <= 0)) stop("Invalid training clinical scales.")
  clinical[, columns] <- sweep(sweep(clinical[, columns, drop = FALSE], 2L, centers, "-"), 2L, scales, "/")
  modules <- base$module_scores
  colnames(modules) <- paste0("module__", sub("^module__", "", colnames(modules)))
  X <- cbind(clinical, modules)
  penalty <- numeric(ncol(X))
  ranking <- NULL
  if (screen_dimension > 0L) {
    genes <- base$genes[, seq_len(5000L), drop = FALSE]
    colnames(genes) <- paste0("gene__", sub("^gene__", "", colnames(genes)))
    cluster <- as.character(metadata$subject_id[training])
    score <- tau - (metadata$log_c3[training] <= quantile(metadata$log_c3[training], tau, names = FALSE, type = 7))
    contribution <- rowsum(genes[training, , drop = FALSE] * score, cluster, reorder = FALSE)
    contribution <- contribution / as.numeric(table(factor(cluster, levels = rownames(contribution))))
    ranking <- abs(colMeans(contribution))
    selected <- order(ranking, decreasing = TRUE)[seq_len(screen_dimension)]
    X <- cbind(X, genes[, selected, drop = FALSE])
    penalty <- c(penalty, rep(1, length(selected)))
  }
  list(X = X, penalty_factor = penalty, clinical_transform = list(columns = columns, center = centers, scale = scales), ranking = ranking)
}

submission_qgee <- function(metadata, design, training_rows, tau, seed,
                            kernel = submission_qgee_kernel(), allow_failed_penalties = TRUE) {
  set.seed(seed)
  rows <- which(training_rows)
  rows <- rows[order(as.character(metadata$subject_id[rows]), metadata$cumulative_time[rows])]
  X <- design$X[rows, , drop = FALSE]
  y <- metadata$log_c3[rows]
  ids <- as.character(metadata$subject_id[rows])
  if (colnames(X)[1L] != "(Intercept)" || any(X[, 1L] != 1)) stop("Invalid QGEE intercept.")
  attr(X, "assign") <- c(0L, seq_len(ncol(X) - 1L))
  nobs <- as.integer(table(factor(ids, levels = unique(ids))))
  initial <- as.numeric(quantreg::rq.fit.br(X, y, tau = tau)$coefficients)
  ratio <- if (length(nobs) > ncol(X)) .001 else .01
  penalties <- exp(seq(log(10), log(10 * ratio), length.out = 30L))
  one_penalty <- function(lambda) kernel(X, y, nobs, tau, "exchangeable", lambda,
    design$penalty_factor[-1L], TRUE, initial,
    geeVerse::qpgeeControl(epsilon = 1e-4, maxit = 1000L, standardize = FALSE))
  fits <- lapply(penalties, function(lambda) {
    if (!allow_failed_penalties) return(one_penalty(lambda))
    tryCatch(one_penalty(lambda), error = function(condition) list(hbic = Inf, converged = FALSE,
      iterations = NA_integer_, failure_message = conditionMessage(condition), condition = condition))
  })
  hbic <- vapply(fits, `[[`, numeric(1), "hbic")
  if ((!allow_failed_penalties && any(!is.finite(hbic))) || !any(is.finite(hbic))) stop("No eligible QGEE HBIC path.")
  eligible <- which(is.finite(hbic))
  selected <- eligible[which.min(hbic[eligible])]
  fit <- fits[[selected]]
  beta <- setNames(as.numeric(fit$coefficients), colnames(X))
  if (any(!is.finite(beta))) stop("Nonfinite QGEE coefficients.")
  if (max(abs(drop(X %*% beta) - drop(fit$fitted.values))) > 1e-10) stop("QGEE prediction coordinates differ.")
  list(beta_prediction = beta, converged = isTRUE(fit$converged), status = if (isTRUE(fit$converged)) "ok" else "nonconverged",
       selected_penalty = penalties[selected], lambda_path = data.frame(lambda = penalties, hbic = hbic,
         converged = vapply(fits, `[[`, logical(1), "converged"), iterations = vapply(fits, function(x) as.integer(x$iterations), integer(1)),
         failure_message = vapply(fits, function(x) if (is.null(x$failure_message)) "" else x$failure_message, character(1))),
       penalty_failure_policy = if (allow_failed_penalties) "source_HBIC_individual_penalties" else "complete_finite_path",
       penalty_failures = lapply(which(vapply(fits, function(x) !is.null(x$condition), logical(1))), function(j)
         list(lambda = penalties[j], condition = fits[[j]]$condition)),
       method = "QGEE-SCAD", method_version = "geeVerse_0.3.1_HBIC_training_clinical_units",
       numerical_backend = "cluster_batched_products_LAPACK_SVD", package_version = as.character(packageVersion("geeVerse")))
}

submission_lqmm <- function(metadata, design, training_rows, tau, seed) {
  set.seed(seed)
  rows <- which(training_rows)
  rows <- rows[order(as.character(metadata$subject_id[rows]), metadata$cumulative_time[rows])]
  X <- design$X[rows, , drop = FALSE]
  fixed_names <- sprintf("covariate_%03d", seq_len(ncol(X)))
  data <- data.frame(response = metadata$log_c3[rows], subject = factor(metadata$subject_id[rows]), X, check.names = FALSE)
  names(data)[seq_len(ncol(X)) + 2L] <- fixed_names
  formula <- reformulate(fixed_names, response = "response", intercept = FALSE)
  fit <- lqmm::lqmm(fixed = formula, random = ~1, group = "subject", tau = tau, data = data,
    covariance = "pdSymm", nK = 7L,
    control = lqmm::lqmmControl(LP_tol_ll = 1e-6, LP_max_iter = 2000L, UP_max_iter = 100L))
  beta <- setNames(as.numeric(coef(fit)), colnames(X))
  if (any(!is.finite(beta))) stop("Nonfinite LQMM coefficients.")
  converged <- all(unlist(fit$opt) >= 0)
  list(beta_prediction = beta, converged = converged, status = if (converged) "ok" else "nonconverged",
       method = "LQMM", method_version = "lqmm_official_training_clinical_units", fit_object = fit,
       convergence_information = fit$opt, package_version = as.character(packageVersion("lqmm")))
}

source("research/submission_completion_20261027/qgee_lapack.R")
source("research/submission_completion_20261027/qgee_algebra.R")
submission_compile_qgee <- function() {
  lock <- filelock::lock("results/scratch/submission_completion_20261027/qgee_compile.lock", timeout = 60000)
  if (is.null(lock)) stop("Unable to acquire the QGEE compilation lock.")
  on.exit(filelock::unlock(lock), add = TRUE)
  Rcpp::sourceCpp("research/submission_completion_20261027/qgee_algebra.cpp",
    cacheDir = "results/scratch/submission_completion_20261027/qgee_cache", rebuild = FALSE,
    showOutput = FALSE, env = .GlobalEnv)
}
submission_compile_qgee()
rm(submission_compile_qgee)

drug_fbio_fit <- function(y, genes, dose, Z, cluster_id, weights = NULL,
                          threshold = floor(2 * length(unique(cluster_id)) /
                                              log(length(unique(cluster_id)))),
                          nfolds = 5L, seed = 20261103L, fold = NULL) {
  stopifnot(nrow(genes) == length(y), nrow(Z) == length(y),
            length(dose) == length(y), length(cluster_id) == length(y),
            !anyNA(cluster_id), all(is.finite(y)), all(is.finite(genes)),
            all(is.finite(Z)), all(is.finite(dose)), all(Z[, 1L] == 1),
            !is.null(attr(Z, "assign")), ncol(Z) >= 2L,
            threshold >= 1L, threshold <= ncol(genes))
  if (is.null(weights)) weights <- 1 / as.numeric(table(cluster_id)[cluster_id])
  stopifnot(length(weights) == length(y), all(is.finite(weights)), all(weights > 0))
  if (is.null(fold)) {
    clusters <- sort(unique(as.character(cluster_id)))
    stopifnot(length(clusters) >= nfolds)
    set.seed(seed)
    allocation <- sample(rep(seq_len(nfolds), length.out = length(clusters)))
    fold <- allocation[match(as.character(cluster_id), clusters)]
  }
  stopifnot(length(fold) == length(y), !anyNA(fold), length(unique(fold)) == nfolds)
  stopifnot(all(vapply(split(fold, as.character(cluster_id)),
                      function(value) length(unique(value)) == 1L, logical(1))))
  environment <- new.env(parent = asNamespace("fbioSelectr"))
  environment$fbioSelect <- fbioSelectr::fbioSelect
  environment(environment$fbioSelect) <- environment
  environment$cv.grpreg <- function(X, y, group, ..., nfolds) {
    grpreg::cv.grpreg(X, y, group, ..., nfolds = nfolds, fold = fold)
  }
  started <- proc.time()[["elapsed"]]
  fit <- environment$fbioSelect(
    y = y, id = cluster_id, t.var = dose, Z = Z, X = genes,
    low.dim.cov.number = ncol(Z) - 1L, threshold.size = threshold,
    type.of.penalty = "grSCAD", degree = 3L, intercept = TRUE,
    degrFreedom = 5L, knots = NULL, weights = weights,
    nfolds = nfolds, seed = seed, ncores = 1L
  )
  selected <- as.character(head(fit$gene.ranking.marginal.screening$hdcov, threshold))
  basis <- splines::bs(dose, df = 5L, degree = 3L, intercept = TRUE)
  model <- list(
    fit = fit, selected = selected, basis = basis, Z_columns = colnames(Z),
    fold = fold, cluster_id = as.character(cluster_id),
    lambda = fit$pen.regression.fit$lambda,
    runtime_sec = proc.time()[["elapsed"]] - started,
    source = "fbioSelectr_0.1.0_with_explicit_cv_folds"
  )
  stopifnot(all(is.finite(fit$pen.regression.fit$beta)),
            all(fit$pen.regression.fit$iter < 10000L))
  model
}

drug_fbio_predict <- function(model, genes, dose, Z) {
  stopifnot(identical(colnames(Z), model$Z_columns),
            all(model$selected %in% colnames(genes)), nrow(genes) == nrow(Z),
            length(dose) == nrow(Z), all(Z[, 1L] == 1))
  basis <- predict(model$basis, newx = dose)
  design <- cbind(Z, genes[, model$selected, drop = FALSE])
  expanded <- do.call(cbind, lapply(seq_len(ncol(design)), function(j) basis * design[, j]))
  prediction <- drop(predict(model$fit$pen.regression.fit, X = expanded[, -1L, drop = FALSE],
                             type = "response"))
  stopifnot(length(prediction) == nrow(Z), all(is.finite(prediction)))
  prediction
}

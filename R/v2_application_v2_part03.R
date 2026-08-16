filter_top_training_mad_v2 <- function(expression, training_samples, d = 5000L) {
  expression <- as.matrix(expression)
  train_idx <- match(training_samples, colnames(expression))
  if (anyNA(train_idx)) stop("Unknown training samples.")
  mad_values <- apply(expression[, train_idx, drop = FALSE], 1, stats::mad,
                      constant = 1, na.rm = TRUE)
  d <- min(as.integer(d), nrow(expression))
  keep <- order(mad_values, decreasing = TRUE)[seq_len(d)]
  list(expression = expression[keep, , drop = FALSE],
       feature_ids = rownames(expression)[keep], mad = mad_values[keep])
}

scale_from_training_v2 <- function(X, training_rows, min_sd = 1e-8) {
  X <- as.matrix(X)
  mu <- colMeans(X[training_rows, , drop = FALSE])
  s <- apply(X[training_rows, , drop = FALSE], 2, stats::sd)
  keep <- is.finite(s) & s > min_sd
  Xs <- sweep(sweep(X[, keep, drop = FALSE], 2, mu[keep], "-"), 2, s[keep], "/")
  list(X = Xs, center = mu[keep], scale = s[keep], keep = keep)
}

build_clinical_design_v2 <- function(metadata, include_sledai = FALSE) {
  md <- metadata
  md$centered_time <- md$cumulative_time - mean(md$cumulative_time, na.rm = TRUE)
  vars <- c("centered_time", "age", "gender", "race", "treatment", "batch")
  if (include_sledai) vars <- c(vars, "sledai")
  missing <- setdiff(vars, names(md))
  if (length(missing)) stop("Clinical fields missing: ", paste(missing, collapse = ", "))
  formula <- stats::as.formula(paste("~", paste(vars, collapse = "+")))
  X <- stats::model.matrix(formula, data = md, na.action = stats::na.pass)
  list(X = X, variables = vars, formula = formula)
}

prepare_clinical_design_split_v2 <- function(metadata,
                                             training_rows,
                                             include_sledai = FALSE,
                                             rare_level_min = 5L) {
  md <- as.data.frame(metadata, stringsAsFactors = FALSE)
  n <- nrow(md)
  if (is.logical(training_rows)) {
    if (length(training_rows) != n) stop("training_rows logical vector has wrong length.")
    train <- which(training_rows)
  } else {
    train <- as.integer(training_rows)
  }
  if (!length(train) || any(train < 1L | train > n)) stop("Invalid training_rows.")

  required <- c("cumulative_time", "age", "gender", "race", "treatment", "batch")
  if (include_sledai) required <- c(required, "sledai")
  missing <- setdiff(required, names(md))
  if (length(missing)) stop("Clinical fields missing: ", paste(missing, collapse = ", "))

  time_center <- mean(as.numeric(md$cumulative_time[train]), na.rm = TRUE)
  if (!is.finite(time_center)) stop("Training cumulative_time has no finite values.")
  md$centered_time <- as.numeric(md$cumulative_time) - time_center

  numeric_vars <- c("centered_time", "age")
  if (include_sledai) numeric_vars <- c(numeric_vars, "sledai")
  numeric_imputation <- setNames(numeric(length(numeric_vars)), numeric_vars)
  for (v in numeric_vars) {
    x <- suppressWarnings(as.numeric(md[[v]]))
    med <- stats::median(x[train], na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    x[!is.finite(x)] <- med
    md[[v]] <- x
    numeric_imputation[v] <- med
  }

  factor_vars <- c("gender", "race", "treatment", "batch")
  factor_levels <- vector("list", length(factor_vars))
  names(factor_levels) <- factor_vars
  for (v in factor_vars) {
    x <- standardize_geo_missing_v2(md[[v]])
    x[is.na(x)] <- "Missing"
    train_tab <- table(x[train])
    retained <- names(train_tab)[train_tab >= rare_level_min]
    retained <- union(retained, "Missing")
    # Levels absent or rare in training are pooled before model.matrix is formed.
    x[!x %in% retained] <- "Other"
    levels_v <- sort(unique(c(retained, "Other")))
    md[[v]] <- factor(x, levels = levels_v)
    factor_levels[[v]] <- levels_v
  }

  rhs <- c("centered_time", "age", factor_vars)
  if (include_sledai) rhs <- c(rhs, "sledai")
  formula <- stats::as.formula(paste("~", paste(rhs, collapse = " + ")))
  X <- stats::model.matrix(formula, data = md, na.action = stats::na.pass)
  if (any(!is.finite(X))) stop("Clinical model matrix contains non-finite entries.")

  list(
    X = X,
    formula = formula,
    time_center = time_center,
    numeric_imputation = numeric_imputation,
    factor_levels = factor_levels,
    variables = rhs
  )
}

score_frozen_modules_v2 <- function(expression_gene,
                                    module_signatures,
                                    minimum_fraction = 0.60) {
  expression_gene <- assert_numeric_matrix_v2(expression_gene, "expression_gene")
  module_signatures <- as.data.frame(module_signatures, stringsAsFactors = FALSE)
  required <- c("module_id", "gene_symbol")
  missing <- setdiff(required, names(module_signatures))
  if (length(missing)) stop("Module signature fields missing: ", paste(missing, collapse = ", "))
  module_signatures <- module_signatures[
    nzchar(module_signatures$module_id) & nzchar(module_signatures$gene_symbol), , drop = FALSE
  ]
  if (!nrow(module_signatures)) stop("Module signature file is empty; confirmatory analysis is blocked.")

  modules <- split(module_signatures$gene_symbol, module_signatures$module_id)
  rows <- vector("list", length(modules))
  score <- matrix(NA_real_, nrow = ncol(expression_gene), ncol = length(modules),
                  dimnames = list(colnames(expression_gene), names(modules)))
  for (ii in seq_along(modules)) {
    genes <- unique(as.character(modules[[ii]]))
    observed <- intersect(genes, rownames(expression_gene))
    fraction <- length(observed) / max(length(genes), 1L)
    eligible <- length(observed) > 0L && fraction >= minimum_fraction
    if (eligible) {
      score[, ii] <- colMeans(expression_gene[observed, , drop = FALSE])
    }
    rows[[ii]] <- data.frame(
      module_id = names(modules)[ii],
      listed_genes = length(genes),
      observed_genes = length(observed),
      observed_fraction = fraction,
      eligible = eligible,
      stringsAsFactors = FALSE
    )
  }
  list(scores = score, audit = do.call(rbind, rows))
}

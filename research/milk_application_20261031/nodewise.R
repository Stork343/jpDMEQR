milk_nodewise_compute <- function(factor, rows, coordinate) {
  n <- length(rows)
  p <- ncol(factor)
  N <- nrow(factor)
  diagonal <- sqrt(colSums(factor^2) / n)
  if (any(!is.finite(diagonal)) || any(diagonal <= 1e-10)) stop("A design coordinate is unidentified.")
  standardised <- sweep(factor, 2L, diagonal, "/")
  other <- setdiff(seq_len(p), coordinate)
  response <- standardised[, coordinate]
  design <- standardised[, other, drop = FALSE]
  solve_path <- function(indices, clusters, ratios) {
    penalty_max <- max(abs(drop(crossprod(design[indices, , drop = FALSE], response[indices]))) / clusters)
    penalties <- penalty_max * ratios
    if (length(other) == 1L) {
      coefficient <- vapply(penalties, function(penalty) {
        soft_threshold_weighted_v2(sum(design[indices, 1L] * response[indices]) / clusters, penalty) /
          (sum(design[indices, 1L]^2) / clusters)
      }, numeric(1))
      return(list(coefficient = matrix(coefficient, 1L), penalties = penalties, kkt = 0))
    }
    fit <- glmnet::glmnet(design[indices, , drop = FALSE] * sqrt(length(indices) / clusters),
      response[indices] * sqrt(length(indices) / clusters), family = "gaussian", alpha = 1,
      lambda = penalties, intercept = FALSE, standardize = FALSE, thresh = 1e-13,
      maxit = 10000000L, type.gaussian = "covariance")
    if (fit$jerr != 0L || ncol(fit$beta) != length(penalties)) stop("The design path did not converge.")
    coefficient <- as.matrix(fit$beta)
    residual <- response[indices] - design[indices, , drop = FALSE] %*% coefficient
    gradient <- -crossprod(design[indices, , drop = FALSE], residual) / clusters
    defect <- vapply(seq_along(penalties), function(j) {
      max(ifelse(abs(coefficient[, j]) > 1e-10,
        abs(gradient[, j] + penalties[j] * sign(coefficient[, j])),
        pmax(abs(gradient[, j]) - penalties[j], 0)))
    }, numeric(1))
    if (max(defect) > 1e-6) stop("The design path failed its KKT check.")
    list(coefficient = coefficient, penalties = penalties, kkt = defect)
  }
  ratios <- pc_ratio_path(n, 24L)
  folds <- pc_folds(rows, 5L)
  loss <- matrix(0, 5L, length(ratios))
  path_kkt <- matrix(0, 5L, length(ratios))
  for (f in 1:5) {
    training <- unlist(rows[folds != f], use.names = FALSE)
    validation <- unlist(rows[folds == f], use.names = FALSE)
    path <- solve_path(training, sum(folds != f), ratios)
    residual <- response[validation] - design[validation, , drop = FALSE] %*% path$coefficient
    loss[f, ] <- colSums(residual^2)
    path_kkt[f, ] <- path$kkt
  }
  selected <- which.min(colSums(loss))
  path <- solve_path(seq_len(N), n, ratios[selected])
  coefficient <- as.numeric(path$coefficient[, 1L])
  penalty <- rep(path$penalties[1L], length(other))
  residual <- response - drop(design %*% coefficient)
  gradient <- -drop(crossprod(design, residual)) / n
  defect <- ifelse(abs(coefficient) > 1e-10,
    abs(gradient + penalty * sign(coefficient)), pmax(abs(gradient) - penalty, 0))
  if (max(defect) > 1e-6) stop("The selected design path failed its KKT check.")
  residual_tau2 <- sum(residual^2) / n + sum(penalty * abs(coefficient))
  tau2 <- sum(response * residual) / n
  identity_error <- abs(residual_tau2 - tau2 - sum(coefficient * (gradient + penalty * sign(coefficient))))
  if (identity_error > 1e-8 || !is.finite(tau2) || tau2 <= 1e-10) stop("Invalid design curvature identity.")
  contrast <- numeric(p)
  contrast[coordinate] <- 1
  contrast[other] <- -coefficient
  omega <- contrast / diagonal / (diagonal[coordinate] * tau2)
  product <- drop(crossprod(factor, factor %*% omega)) / n
  normalisation_error <- abs(product[coordinate] - 1)
  if (normalisation_error > 1e-6) stop("The design direction normalisation failed.")
  list(omega = omega, coefficient = coefficient, penalty = penalty, kkt = max(defect),
    normalisation_error = normalisation_error, inverse_defect = max(abs(product - replace(numeric(p), coordinate, 1))),
    tau2 = tau2, residual_tau2 = residual_tau2, identity_error = identity_error, diagonal = diagonal,
    tuning = list(rule = "cluster_curvature_cv", ratios = ratios, fold_loss_sum = loss,
      selected = selected, selected_ratio = ratios[selected], path_kkt = path_kkt))
}

milk_nodewise_cached <- function(factor, rows, coordinate) {
  index <- which(vapply(milk_design_cache, function(entry) {
    identical(entry$subjects, names(rows)) && identical(entry$coordinate, coordinate)
  }, logical(1)))
  if (length(index) != 1L || !identical(milk_design_cache[[index]]$factor, factor))
    stop("The requested design is absent from the verified nodewise cache.")
  milk_design_cache[[index]]$fit
}

milk_curvature_supports <- function(factor, n, coordinate, ratios) {
  p <- ncol(factor)
  other <- setdiff(seq_len(p), coordinate)
  diagonal <- sqrt(colSums(factor^2) / n)
  if (any(!is.finite(diagonal)) || any(diagonal <= 1e-10)) stop("Unidentified profile coordinate.")
  if (!length(other)) return(rep(list(integer()), length(ratios)))
  standardised <- sweep(factor, 2L, diagonal, "/")
  design <- standardised[, other, drop = FALSE] * sqrt(nrow(factor) / n)
  response <- standardised[, coordinate] * sqrt(nrow(factor) / n)
  maximum <- max(abs(drop(crossprod(design, response)))) / nrow(factor)
  if (maximum == 0) return(rep(list(integer()), length(ratios)))
  if (length(other) == 1L) {
    coefficient <- matrix(vapply(maximum * ratios, function(lambda) {
      soft_threshold_weighted_v2(sum(design[, 1L] * response) / nrow(factor), lambda) /
        mean(design[, 1L]^2)
    }, numeric(1)), nrow = 1L)
  } else {
    path <- glmnet::glmnet(design, response, lambda = maximum * ratios,
      family = "gaussian", intercept = FALSE, standardize = FALSE, alpha = 1,
      thresh = 1e-13, maxit = 10000000L, type.gaussian = "covariance")
    if (path$jerr != 0L || ncol(path$beta) != length(ratios)) stop("Incomplete curvature path.")
    coefficient <- as.matrix(path$beta)
  }
  gradient <- -crossprod(design, response - design %*% coefficient) / nrow(factor)
  defects <- vapply(seq_along(ratios), function(j) {
    max(ifelse(abs(coefficient[, j]) > 1e-10,
      abs(gradient[, j] + maximum * ratios[j] * sign(coefficient[, j])),
      pmax(abs(gradient[, j]) - maximum * ratios[j], 0)))
  }, numeric(1))
  if (max(defects) > 1e-6) stop("The curvature path failed its KKT check.")
  lapply(seq_along(ratios), function(j) other[abs(coefficient[, j]) > 1e-10])
}

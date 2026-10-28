sb_data <- function(y, X, Z, cluster_id) {
  rows <- cluster_index_v2(cluster_id)$rows
  group <- integer(length(y))
  weight <- numeric(length(y))
  for (i in seq_along(rows)) {
    group[rows[[i]]] <- i
    weight[rows[[i]]] <- 1 / length(rows[[i]])
  }
  list(y = y, X = X, Z = Z, cluster_id = cluster_id, rows = rows, n = length(rows),
       group = group, weight = weight, order = as.integer(unlist(rows) - 1L),
       pointers = as.integer(c(0L, cumsum(lengths(rows)))))
}

sb_subset <- function(data, index) {
  sb_data(data$y[index], data$X[index, , drop = FALSE],
          data$Z[index, , drop = FALSE], data$cluster_id[index])
}

sb_budget <- function(n, p, included = 0L, multiplier = 1) {
  if (length(multiplier) != 1L || !is.finite(multiplier) || multiplier <= 0) stop("Invalid search budget multiplier.")
  capacity <- floor(multiplier * min(floor(n / log(n)), ceiling(log(n)^3)))
  as.integer(min(p - included, max(0L, capacity - included)))
}

sb_profile <- function(data, residual, tau, h, Lambda, gamma = NULL) {
  full <- compiled_cluster_profiles_v2(residual, data$Z, data$rows, tau, h, Lambda,
                                       gamma_init = gamma,
                                       control = list(grad_tol = 1e-10, newton_maxit = 150L))
  half <- compiled_cluster_profiles_v2(residual, data$Z, data$rows, tau, h / 2, Lambda,
                                       gamma_init = full$gamma,
                                       control = list(grad_tol = 1e-10, newton_maxit = 150L))
  list(full = full, half = half)
}

sb_geometry <- function(data, phi, Lambda, coordinates = seq_len(ncol(data$X)), factor = FALSE) {
  sb_geometry_cpp(data$X[, coordinates, drop = FALSE], data$Z, phi,
                   data$order, data$pointers, Lambda, factor)
}

sb_initial <- function(data, tau, Lambda, penalty_factor, budget_multiplier = 1) {
  folds <- pc_folds(data$rows, 3L)
  min_training <- min(vapply(1:3, function(f) sum(folds != f), integer(1)))
  sizes <- 0:sb_budget(min_training, ncol(data$X), sum(penalty_factor == 0), budget_multiplier)
  paths <- training <- validation <- vector("list", 3L)
  loss <- matrix(0, 3L, length(sizes))
  for (f in 1:3) {
    training[[f]] <- sb_subset(data, unlist(data$rows[folds != f], use.names = FALSE))
    validation[[f]] <- sb_subset(data, unlist(data$rows[folds == f], use.names = FALSE))
    train <- training[[f]]
    valid <- validation[[f]]
    prepared <- pc_prepare_initial(train$y, train$X, train$Z, train$cluster_id,
                                    tau, Lambda, penalty_factor)
    path <- pc_orthogonal_path(prepared, tau, sizes, gap_tolerance = min(1e-8, train$n^(-3)),
                               max_iterations = 100000L)
    for (j in seq_along(sizes)) {
      profile <- pc_exact_profile(valid$y - drop(valid$X %*% path$coefficients[, j]),
                                   valid$Z, valid$rows, Lambda, tau)
      loss[f, j] <- sum(profile$loss)
    }
    paths[[f]] <- path
  }
  selected <- which.min(colSums(loss))
  list(training = training, validation = validation, paths = paths,
       selected = selected, size = sizes[selected], loss = loss,
       supports = lapply(paths, function(path) path$selected_path[[selected]]),
       coefficients = lapply(paths, function(path) path$coefficients[, selected]),
       fold_subjects = lapply(validation, function(data) names(data$rows)))
}

sb_lasso_supports <- function(factor, n, coordinate, ratios) {
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
    path <- glmnet::glmnet(design, response, lambda = maximum * ratios, family = "gaussian",
                           intercept = FALSE, standardize = FALSE, alpha = 1,
                           thresh = 1e-13, maxit = 1000000L)
    if (path$jerr != 0L || ncol(path$beta) != length(ratios)) stop("Incomplete nodewise path.")
    coefficient <- as.matrix(path$beta)
  }
  lapply(seq_along(ratios), function(j) other[abs(coefficient[, j]) > 1e-10])
}

sb_contrast <- function(H, reference, eta) {
  d <- ncol(H)
  if (d == 1L) return(1)
  nuisance <- 2:d
  system <- H[nuisance, nuisance, drop = FALSE] + diag(eta * reference[nuisance], d - 1L)
  factor <- chol(system)
  coefficient <- backsolve(factor, forwardsolve(t(factor), H[nuisance, 1L]))
  c(1, -drop(coefficient))
}

sb_graph <- function(initial, coordinates, tau, h, Lambda, budget_multiplier = 1) {
  min_training <- min(vapply(initial$training, `[[`, integer(1), "n"))
  ratios <- 2^(-(0:ceiling(2 * log2(min_training))) / 4)
  p <- ncol(initial$training[[1L]]$X)
  factors <- vector("list", 3L)
  for (f in 1:3) {
    train <- initial$training[[f]]
    valid <- initial$validation[[f]]
    beta <- initial$coefficients[[f]]
    training_profile <- sb_profile(train, train$y - drop(train$X %*% beta), tau, h, Lambda)
    validation_profile <- sb_profile(valid, valid$y - drop(valid$X %*% beta), tau, h, Lambda)
    factors[[f]] <- list(
      train = lapply(training_profile, function(profile) sb_geometry(train, profile$phi, Lambda, factor = TRUE)),
      valid = lapply(validation_profile, function(profile) sb_geometry(valid, profile$phi, Lambda, factor = TRUE)))
  }
  results <- vector("list", length(coordinates))
  for (j in seq_along(coordinates)) {
    k <- coordinates[j]
    losses <- matrix(0, 3L, length(ratios))
    supports <- vector("list", 3L)
    for (f in 1:3) {
      ntrain <- initial$training[[f]]$n
      F <- factors[[f]]$train
      V <- factors[[f]]$valid
      full <- sb_lasso_supports(F$full, ntrain, k, ratios)
      half <- sb_lasso_supports(F$half, ntrain, k, ratios)
      supports[[f]] <- vector("list", length(ratios))
      for (r in seq_along(ratios)) {
        graph <- union(full[[r]], half[[r]])
        chosen <- c(k, setdiff(union(initial$supports[[f]], graph), k))
        supports[[f]][[r]] <- chosen
        if (length(graph) > sb_budget(ntrain, p, multiplier = budget_multiplier)) {
          losses[f, r] <- Inf
          next
        }
        Hfull <- crossprod(F$full[, chosen, drop = FALSE]) / ntrain
        Hhalf <- crossprod(F$half[, chosen, drop = FALSE]) / ntrain
        reference <- (diag(Hfull) + diag(Hhalf)) / 2
        cfull <- sb_contrast(Hfull, reference, 1 / ntrain)
        chalf <- sb_contrast(Hhalf, reference, 1 / ntrain)
        losses[f, r] <- sum((V$full[, chosen, drop = FALSE] %*% cfull)^2) +
          sum((V$half[, chosen, drop = FALSE] %*% chalf)^2)
      }
    }
    selected <- which.min(colSums(losses))
    if (!is.finite(sum(losses[, selected]))) stop("No admissible graph candidate.")
    support <- sort(unique(unlist(lapply(supports, `[[`, selected))))
    results[[j]] <- list(coordinate = k, support = c(k, setdiff(support, k)),
                         ratios = ratios, loss = losses, selected = selected,
                         inner_supports = lapply(supports, `[[`, selected))
  }
  results
}

sb_prepare_refit <- function(data, support, beta, tau, h, Lambda, gap_tolerance = 1e-12) {
  small <- data
  small$X <- data$X[, support, drop = FALSE]
  square_weight <- sqrt(data$weight / data$n)
  decomposition <- qr(small$X * square_weight, tol = 1e-10)
  if (decomposition$rank != length(support)) stop("Selected design is rank deficient.")
  pivot <- decomposition$pivot
  basis <- qr.Q(decomposition) / square_weight
  triangular <- qr.R(decomposition)
  theta <- drop(triangular %*% beta[support[pivot]])
  exact <- pc_exact_profile(data$y - drop(small$X %*% beta[support]), data$Z, data$rows, Lambda, tau)
  fitted <- joint_profile_pdal_cpp(
    basis, data$Z, data$y, data$group - 1L, data$weight, Lambda,
    numeric(length(support)), tau, theta, exact$gamma, exact$dual,
    100000L, min(gap_tolerance, data$n^(-4)))
  coefficient <- numeric(length(support))
  coefficient[pivot] <- drop(backsolve(triangular, fitted$theta))
  anchor <- sb_profile(small, data$y - drop(small$X %*% coefficient), tau, h, Lambda, fitted$gamma)
  Hfull <- sb_geometry(small, anchor$full$phi, Lambda)
  Hhalf <- sb_geometry(small, anchor$half$phi, Lambda)
  reference <- (diag(Hfull) + diag(Hhalf)) / 2
  if (any(reference <= 0) || any(!is.finite(reference))) stop("Nonpositive reference curvature.")
  cfull <- sb_contrast(Hfull, reference, 1 / data$n)
  chalf <- sb_contrast(Hhalf, reference, 1 / data$n)
  list(support = support, coefficient = coefficient, reference = reference,
       full_contrast = cfull, half_contrast = chalf, contrast = 2 * chalf - cfull,
       gap = fitted$duality_gap, iterations = fitted$iterations,
       full_curvature = Hfull, half_curvature = Hhalf)
}

sb_evaluator <- function(folds, coordinate_index, tau, bandwidth, Lambda) {
  cache <- new.env(parent = emptyenv())
  evaluations <- 0L
  validation <- lapply(folds, function(fold) {
    fitted <- fold$refits[[coordinate_index]]
    data <- fold$validation
    X <- data$X[, fitted$support, drop = FALSE]
    projection <- drop(X %*% fitted$contrast)
    data$base <- data$y - drop(X %*% fitted$coefficient) + projection * fitted$coefficient[1L]
    data$X <- matrix(projection, ncol = 1L)
    data
  })
  gamma <- vector("list", length(folds))
  n <- sum(vapply(validation, `[[`, integer(1), "n"))
  function(value) {
    key <- sprintf("%.17g", value)
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache, inherits = FALSE))
    values <- vector("list", length(validation))
    derivative <- 0
    for (f in seq_along(validation)) {
      data <- validation[[f]]
      profile <- compiled_cluster_profiles_v2(data$base - data$X[, 1L] * value,
                                               data$Z, data$rows, tau, bandwidth, Lambda,
                                               gamma_init = gamma[[f]],
                                               control = list(grad_tol = 1e-10, newton_maxit = 150L))
      gamma[[f]] <<- profile$gamma
      values[[f]] <- -as.numeric(rowsum(data$X[, 1L] * profile$psi * data$weight,
                                        data$group, reorder = FALSE))
      derivative <- derivative + data$n * sb_geometry(data, profile$phi, Lambda)[1L, 1L] / n
    }
    score <- unlist(values, use.names = FALSE)
    average <- mean(score)
    se <- sqrt(mean((score - average)^2) / length(score))
    if (!is.finite(se) || se <= 0) stop("Nonpositive validation score variance.")
    evaluations <<- evaluations + 1L
    answer <- list(z = average / se, mean = average, se = se, score = score,
                    derivative = derivative, evaluations = evaluations)
    assign(key, answer, envir = cache)
    answer
  }
}

sb_root <- function(evaluate, center, initial_scale) {
  width <- initial_scale
  bracketed <- FALSE
  for (expansion in 1:24) {
    lower <- center - width
    upper <- center + width
    if (evaluate(lower)$mean < 0 && evaluate(upper)$mean > 0) {
      bracketed <- TRUE
      break
    }
    width <- 2 * width
  }
  if (!bracketed) stop("Unable to bracket the target score.")
  root <- uniroot(function(t) evaluate(t)$mean, c(lower, upper), tol = 1e-8)$root
  at_root <- evaluate(root)
  if (!is.finite(at_root$derivative) || at_root$derivative <= 1e-10) stop("Nonpositive target path curvature.")
  list(estimate = root, derivative = at_root$derivative, root_z = at_root$z,
       score = at_root$score, evaluations = at_root$evaluations,
       bracket = c(lower, upper))
}

sb_prepare <- function(y, X, Z, cluster_id, tau, coordinates, Lambda, penalty_factor, h, gap_tolerance = 1e-12) {
  data <- sb_data(y, X, Z, cluster_id)
  assignment <- pc_folds(data$rows, 5L)
  folds <- vector("list", 5L)
  for (f in 1:5) {
    training <- sb_subset(data, unlist(data$rows[assignment != f], use.names = FALSE))
    validation <- sb_subset(data, unlist(data$rows[assignment == f], use.names = FALSE))
    initial <- sb_initial(training, tau, Lambda, penalty_factor)
    graphs <- sb_graph(initial, coordinates, tau, h, Lambda)
    beta <- Reduce(`+`, initial$coefficients) / 3
    refits <- lapply(graphs, function(graph) sb_prepare_refit(training, graph$support, beta, tau, h, Lambda, gap_tolerance))
    initial$training <- initial$validation <- initial$paths <- NULL
    folds[[f]] <- list(training_subjects = names(training$rows), validation = validation,
                       initial = initial, graphs = graphs, refits = refits, beta = beta)
  }
  list(folds = folds, centers = Reduce(`+`, lapply(folds, `[[`, "beta"))[coordinates] / 5,
       coordinates = coordinates, n = data$n)
}

sb_pair_roots <- function(prepared, coordinate_index, tau, h, Lambda, initial_scale, confidence) {
  center <- prepared$centers[coordinate_index]
  full <- sb_root(sb_evaluator(prepared$folds, coordinate_index, tau, h, Lambda), center, initial_scale)
  half <- sb_root(sb_evaluator(prepared$folds, coordinate_index, tau, h / 2, Lambda), center, initial_scale)
  estimate <- 2 * half$estimate - full$estimate
  influence <- -2 * half$score / half$derivative + full$score / full$derivative
  se <- sqrt(mean((influence - mean(influence))^2) / prepared$n)
  if (!is.finite(se) || se <= 0) stop("Nonpositive paired root variance.")
  critical <- qnorm((1 + confidence) / 2)
  list(estimate = estimate, se = se, lower = estimate - critical * se, upper = estimate + critical * se,
       full = full, half = half, influence = influence,
       max_root_z = max(abs(c(full$root_z, half$root_z))),
       min_derivative = min(full$derivative, half$derivative))
}

sb_fit <- function(y, X, Z, cluster_id, tau, coordinates,
                    Lambda = diag(ncol(Z)), penalty_factor = rep(1, ncol(X)),
                    h = length(unique(cluster_id))^(-0.3), confidence = 0.95,
                    gap_tolerance = 1e-12, keep_prepared = FALSE) {
  started <- proc.time()[[3L]]
  validate_tau_h_v2(tau, h)
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  if (length(y) != nrow(X) || length(cluster_id) != nrow(X) || nrow(Z) != nrow(X) ||
      any(!is.finite(y)) || anyNA(cluster_id)) stop("Invalid observations.")
  if (is.null(colnames(X)) || anyNA(colnames(X)) || any(!nzchar(colnames(X))) || anyDuplicated(colnames(X))) {
    stop("X needs unique nonempty column names.")
  }
  if (length(penalty_factor) != ncol(X) || any(!penalty_factor %in% c(0, 1))) stop("Invalid inclusion indicators.")
  if (is.character(coordinates)) coordinates <- match(coordinates, colnames(X))
  if (!length(coordinates) || anyNA(coordinates) || anyDuplicated(coordinates) ||
      any(coordinates < 1 | coordinates > ncol(X))) stop("Invalid target coordinates.")
  if (length(confidence) != 1L || !is.finite(confidence) || confidence <= 0 || confidence >= 1) stop("Invalid confidence level.")
  Lambda <- make_lambda_matrix_v2(Lambda, ncol(Z))
  prepared <- sb_prepare(y, X, Z, cluster_id, tau, coordinates, Lambda, penalty_factor, h, gap_tolerance)
  training_seconds <- proc.time()[[3L]] - started
  inference_started <- proc.time()[[3L]]
  details <- lapply(seq_along(coordinates), function(j) {
    scale <- stats::sd(y) / sqrt(prepared$n) / sqrt(mean(X[, coordinates[j]]^2))
    sb_pair_roots(prepared, j, tau, h, Lambda, scale, confidence)
  })
  table <- data.frame(coordinate = colnames(X)[coordinates],
                       beta_tilde = vapply(details, `[[`, numeric(1), "estimate"),
                       se = vapply(details, `[[`, numeric(1), "se"),
                       lower = vapply(details, `[[`, numeric(1), "lower"),
                       upper = vapply(details, `[[`, numeric(1), "upper"),
                       min_derivative = vapply(details, `[[`, numeric(1), "min_derivative"),
                       max_root_z = vapply(details, `[[`, numeric(1), "max_root_z"))
  record <- lapply(prepared$folds, function(fold) list(
    training_subjects = fold$training_subjects, validation_subjects = names(fold$validation$rows),
    inner_subjects = fold$initial$fold_subjects, initial_size = fold$initial$size,
    initial_loss = fold$initial$loss, initial_supports = fold$initial$supports,
    initial_coefficients = fold$initial$coefficients, graphs = fold$graphs,
    refits = fold$refits))
  list(table = table, details = details, training = record, h = h, confidence = confidence,
       timing = c(training = training_seconds, inference = proc.time()[[3L]] - inference_started,
                  total = proc.time()[[3L]] - started),
       method_version = "honest_paired_roots_20261026",
       prepared = if (keep_prepared) prepared else NULL)
}

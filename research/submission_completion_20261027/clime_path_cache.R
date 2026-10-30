submission_clime_path_guard <- function(H, gamma, cache) {
  if (!is.environment(cache)) stop("The CLIME path cache must belong to one fixed Hessian.")
  if (is.null(cache$hessian)) cache$hessian <- H else
    if (!identical(cache$hessian, H)) stop("The CLIME path cache belongs to a different Hessian.")
  p <- nrow(H)
  maximum <- 4L * p
  value <- cache$first_direction
  if (is.null(value) || (!value$exhausted && !value$infeasible && gamma < value$minimum_knot)) {
    invisible(loadNamespace("fastclime"))
    if (packageVersion("fastclime") != "1.4.1.1") stop("The CLIME path requires fastclime 1.4.1.1.")
    canonical <- numeric(p)
    canonical[1L] <- 1
    started <- proc.time()[[3L]]
    path <- .C("dantzig", as.double(H), as.double(canonical), beta = numeric(p * maximum),
      as.integer(p), as.double(gamma), as.integer(maximum), gamma = rep(NA_real_, maximum),
      NAOK = TRUE, PACKAGE = "fastclime")
    indices <- which(!is.na(path$gamma))
    knots <- path$gamma[indices]
    if (!length(knots) || !any(is.finite(knots))) stop("The native CLIME path returned no finite knot.")
    value <- list(coordinate = 1L, requested_gamma = gamma,
      minimum_knot = min(knots[is.finite(knots)]), final_knot = tail(knots, 1L),
      knots = length(indices), maximum_pivots = maximum,
      elapsed = proc.time()[[3L]] - started,
      exhausted = length(indices) == maximum && !any(is.finite(knots) & knots <= gamma),
      infeasible = length(indices) < maximum && all(is.finite(knots)) && tail(knots, 1L) > gamma)
    cache$first_direction <- value
  }
  if (value$exhausted && gamma < value$minimum_knot) {
    condition <- simpleError("The shared CLIME path exhausted its declared pivot budget.")
    condition$failure_stage <- "precision_solver"
    condition$diagnostics <- c(value, list(gamma = gamma, status = "parametric_iteration_limit",
      shared_path = TRUE))
    stop(condition)
  }
  value
}

soy_compile_pooled <- function() {
  prefix <- normalizePath("results/scratch/dependencies/highs-native", mustWork = TRUE)
  cache <- "results/scratch/soynam_20261031/pooled_lp_cache"
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(file.path(cache, "compile.lock"), timeout = 60000)
  if (is.null(lock)) stop("The pooled LP compiler is already active.")
  on.exit(filelock::unlock(lock), add = TRUE)
  previous <- Sys.getenv(c("PKG_CPPFLAGS", "PKG_LIBS"), unset = NA_character_)
  on.exit({
    for (name in names(previous)) {
      if (is.na(previous[[name]])) Sys.unsetenv(name) else
        do.call(Sys.setenv, setNames(list(previous[[name]]), name))
    }
  }, add = TRUE)
  Sys.setenv(PKG_CPPFLAGS = paste0("-I", shQuote(file.path(prefix, "include/highs"))),
    PKG_LIBS = paste0("-L", shQuote(file.path(prefix, "lib")), " -lhighs -Wl,-rpath,", shQuote(file.path(prefix, "lib"))))
  Rcpp::sourceCpp("research/soynam_application_20261031/pooled_lp.cpp", cacheDir = cache,
    rebuild = FALSE, showOutput = FALSE, env = .GlobalEnv)
}
soy_compile_pooled()
rm(soy_compile_pooled)

soy_pooled_compiled_path <- function(data, tau, nlambda = 40L, lambda_min = .01, lambda = NULL, progress = FALSE) {
  p <- ncol(data$X)
  penalised <- which(data$penalty_factor > 0)
  means <- numeric(p)
  means[penalised] <- colMeans(data$X[, penalised, drop = FALSE])
  X <- sweep(data$X, 2L, means, "-")
  scale <- sqrt(colMeans(X^2))
  if (any(!is.finite(scale)) || any(scale <= 1e-10)) stop("A pooled training column has zero scale.")
  X <- sweep(X, 2L, scale, "/")
  fitted <- soy_pooled_lp_cpp(X, data$y, data$penalty_factor, tau, nlambda, lambda_min,
    if (is.null(lambda)) numeric() else lambda, progress)
  beta <- fitted$beta / scale
  intercept <- -drop(means %*% beta)
  blocks <- grep("^Block[0-9]+$", colnames(data$X))
  beta[blocks, ] <- sweep(beta[blocks, , drop = FALSE], 2L, intercept, "+")
  identity <- max(abs(data$X %*% beta - X %*% fitted$beta))
  if (identity > 1e-9) stop("The compiled pooled coordinate mapping changed its predictions.")
  rownames(beta) <- colnames(data$X)
  list(beta = beta, lambda = fitted$lambda,
    diagnostics = data.table(point = seq_along(fitted$lambda), lambda = fitted$lambda,
      primal_dual_gap = fitted$primal_dual_gap, dual_violation = fitted$dual_violation,
      active_constraints = fitted$active_constraints, selected_variables = fitted$selected_variables,
      generations = fitted$generations, iterations = fitted$iterations, mapping_error = identity),
    solver = paste0("HiGHS_", fitted$highs_version, "_certified_pooled_path"))
}

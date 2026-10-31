library(data.table)
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
source("research/soynam_application_20261031/data.R")
Rcpp::sourceCpp("research/soynam_application_20261031/pair_qp.cpp",
  cacheDir = "results/scratch/soynam_20261031/pair_qp_cache", rebuild = FALSE)
input <- readRDS("results/raw/application/soynam_yield_20261031/input.rds")
data <- soy_design(input)
rows <- cluster_index_v2(data$cluster_id)$rows
outer <- pc_folds(rows, 2L)
training <- soy_subset(data, unlist(rows[outer != 1L], use.names = FALSE))
inner <- pc_folds(training$rows, 3L)
training <- sb_subset(training, unlist(training$rows[inner != 1L], use.names = FALSE))
square_weight <- sqrt(training$weight / training$n)
grid <- rbind(expand.grid(dimension = 56L, tau = c(.25, .5, .75), lambda = c(.5, 1, 2)),
  data.frame(dimension = c(120L, 245L, 419L), tau = .5, lambda = 1))
records <- list()
for (i in seq_len(nrow(grid))) {
  case <- grid[i, ]
  selected <- c(seq_len(data$fixed_dimension),
    if (case$dimension > data$fixed_dimension) as.integer(seq(data$fixed_dimension + 1L, data$p,
      length.out = case$dimension - data$fixed_dimension)) else integer())
  decomposition <- qr(training$X[, selected, drop = FALSE] * square_weight, tol = 1e-10)
  if (decomposition$rank != length(selected)) stop("The real-data verification design is rank deficient.")
  Q <- qr.Q(decomposition)
  X <- Q / square_weight
  initial <- drop(crossprod(Q, training$y * square_weight))
  nuisance <- pc_exact_profile(training$y - drop(X %*% initial), training$Z,
    training$rows, matrix(case$lambda), case$tau)
  arguments <- list(X, training$Z, training$y, training$group - 1L, training$weight,
    matrix(case$lambda), numeric(ncol(X)), case$tau, initial, nuisance$gamma, nuisance$dual,
    max_iterations = 100000L, gap_tolerance = 1e-10)
  started <- proc.time()[[3L]]
  conic <- do.call(jp_joint_qp, arguments)
  conic_seconds <- proc.time()[[3L]] - started
  started <- proc.time()[[3L]]
  fitted <- do.call(soy_pair_qp_cpp, arguments)
  elapsed <- proc.time()[[3L]] - started
  row <- data.table(clusters = training$n, columns = ncol(X), tau = case$tau, lambda = case$lambda,
    conic_seconds = conic_seconds, structured_seconds = elapsed, iterations = fitted$iterations,
    face_checks = fitted$face_checks, polished = fitted$polished,
    conic_gap = conic$duality_gap, structured_gap = fitted$duality_gap,
    objective_difference = abs(conic$primal - fitted$primal),
    prediction_difference = max(abs(X %*% (conic$theta - fitted$theta))))
  if (row$objective_difference > 1e-8 || row$prediction_difference > 1e-4)
    stop("The structured quadratic solver disagrees with the reference.")
  records[[length(records) + 1L]] <- row
  fwrite(rbindlist(records), "results/preflight/soynam_20261031/pair_qp_checks.csv")
  print(row)
  flush.console()
}

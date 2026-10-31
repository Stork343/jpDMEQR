library(data.table)
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
source("research/soynam_application_20261031/data.R")
input <- readRDS("results/raw/application/soynam_yield_20261031/input.rds")
data <- soy_design(input)
rows <- cluster_index_v2(data$cluster_id)$rows
outer <- pc_folds(rows, 2L)
training <- soy_subset(data, unlist(rows[outer != 1L], use.names = FALSE))
inner <- pc_folds(training$rows, 3L)
training <- sb_subset(training, unlist(training$rows[inner != 1L], use.names = FALSE))
square_weight <- sqrt(training$weight / training$n)
records <- list()
for (dimension in c(56L, 120L, 245L)) {
  selected <- c(seq_len(data$fixed_dimension),
    if (dimension > data$fixed_dimension) as.integer(seq(data$fixed_dimension + 1L, data$p,
      length.out = dimension - data$fixed_dimension)) else integer())
  decomposition <- qr(training$X[, selected, drop = FALSE] * square_weight, tol = 1e-10)
  if (decomposition$rank != length(selected)) stop("The real-data benchmark design is rank deficient.")
  Q <- qr.Q(decomposition)
  X <- Q / square_weight
  initial <- drop(crossprod(Q, training$y * square_weight))
  nuisance <- pc_exact_profile(training$y - drop(X %*% initial), training$Z,
    training$rows, matrix(1, 1L, 1L), 0.5)
  arguments <- list(X, training$Z, training$y, training$group - 1L, training$weight,
    matrix(1, 1L, 1L), numeric(ncol(X)), 0.5, initial, nuisance$gamma, nuisance$dual,
    max_iterations = 100000L, gap_tolerance = 1e-10)
  started <- proc.time()[[3L]]
  conic <- do.call(jp_joint_qp, arguments)
  conic_seconds <- proc.time()[[3L]] - started
  record <- data.table(clusters = training$n, columns = ncol(X), conic_seconds = conic_seconds,
    conic_gap = conic$duality_gap, conic_iterations = conic$iterations)
  print(record)
  flush.console()
  if (dimension == 56L) {
    started <- proc.time()[[3L]]
    primal_dual <- do.call(joint_profile_pdal_cpp, arguments)
    record[, `:=`(reference_seconds = proc.time()[[3L]] - started,
      reference_gap = primal_dual$duality_gap,
      objective_difference = abs(primal_dual$primal - conic$primal),
      prediction_difference = max(abs(X %*% (primal_dual$theta - conic$theta))))]
    if (record$objective_difference > 1e-8 || record$prediction_difference > 1e-4)
      stop("The two actual-data solvers disagree.")
  }
  records[[length(records) + 1L]] <- record
  fwrite(rbindlist(records, fill = TRUE), "results/preflight/soynam_20261031/selected_solver_benchmark.csv")
}

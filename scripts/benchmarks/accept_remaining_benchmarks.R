root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "benchmarks", "accept_all_benchmarks_20261022.R"))
run_benchmark_acceptance_20261022(root, methods = c(
  "BIAS-ADJ-LQMM", "DOUBLE-PEN-QLMM", "QGEE-SCAD", "QIF-SEE"
))

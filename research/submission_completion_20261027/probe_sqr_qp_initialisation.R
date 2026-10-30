source("research/submission_completion_20261027/load_comparators.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Provide an input directory, candidate number and coordinate.")
root <- args[1L]
candidate <- as.integer(args[2L])
coordinate <- as.integer(args[3L])
input <- readRDS(file.path(root, "inputs.rds"))
gamma <- input$cells$gamma[candidate]
started <- proc.time()[[3L]]
value <- submission_clime_qp_candidate(input$H, coordinate, gamma)
elapsed <- proc.time()[[3L]] - started
certificate <- if (value$status == "candidate")
  submission_clime_certificate(input$H, coordinate, gamma, value$omega) else NULL
saveRDS(list(candidate = candidate, coordinate = coordinate, gamma = gamma,
  elapsed = elapsed, value = value, certificate = certificate),
  file.path(root, sprintf("qp_%02d_%04d.rds", candidate, coordinate)))
print(list(candidate = candidate, coordinate = coordinate, gamma = gamma, elapsed = elapsed,
  status = value$status, certified = certificate$verified,
  qp_dimensions = vapply(value$qp_history, function(x) length(x$selected), integer(1))))

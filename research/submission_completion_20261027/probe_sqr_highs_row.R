source("research/submission_completion_20261027/load_comparators.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Provide an input directory, candidate and coordinate.")
root <- args[1L]
candidate <- as.integer(args[2L])
coordinate <- as.integer(args[3L])
input <- readRDS(file.path(root, "inputs.rds"))
H <- input$H
p <- nrow(H)
gamma <- input$cells$gamma[candidate]
canonical <- numeric(p)
canonical[coordinate] <- 1
started <- proc.time()[[3L]]
model <- highs::highs_model(L = rep(1, 2L * p), lower = rep(0, 2L * p),
  upper = rep(Inf, 2L * p), A = cbind(H, -H), lhs = canonical - gamma, rhs = canonical + gamma)
control <- highs::highs_control(threads = 1L, solver = "simplex", presolve = "on",
  simplex_strategy = 1L, simplex_scale_strategy = 2L, log_to_console = FALSE,
  primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8)
solver <- highs::highs_solver(model, control)
solver$solve(log_to_console = FALSE)
status <- highs_model_status_label_v2(solver$status())
elapsed <- proc.time()[[3L]] - started
if (!status %in% c("optimal", "infeasible")) stop("The complete row programme did not reach an accepted status.")
saveRDS(list(candidate = candidate, coordinate = coordinate, gamma = gamma, status = status,
  elapsed = elapsed, info = solver$info(), solution = solver$solution()),
  file.path(root, sprintf("highs_%02d_%04d.rds", candidate, coordinate)))
print(list(candidate = candidate, coordinate = coordinate, gamma = gamma,
  status = status, elapsed = elapsed, info = solver$info()))

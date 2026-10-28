args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the actual full-matrix precision.rds record.")
path <- normalizePath(args[1], mustWork = TRUE)
precision <- readRDS(path)
grid_index <- match(precision$gamma, precision$raw_grid_records[[1L]]$mu)
records <- lapply(seq_along(precision$raw_grid_records), function(k) {
  grid <- precision$raw_grid_records[[k]]
  record <- grid$solver_records[[grid_index]]
  primal <- sum(abs(grid$omega[[grid_index]]))
  dual_minimum <- record$info$objective_function_value
  data.frame(coordinate = k, gamma = precision$gamma,
               original_primal_residual = grid$residual[grid_index],
               primal_objective = primal, dual_minimum_objective = dual_minimum,
               relative_duality_gap = abs(primal + dual_minimum) / max(1, primal),
               max_solver_primal_infeasibility = record$info$max_primal_infeasibility,
               max_solver_dual_infeasibility = record$info$max_dual_infeasibility,
               solver_primal_feasible = identical(record$info$primal_solution_status, "Feasible"),
               solver_dual_feasible = identical(record$info$dual_solution_status, "Feasible"),
               basis_valid = isTRUE(record$basis$valid), stringsAsFactors = FALSE)
})
table <- do.call(rbind, records)
table$passed <- with(table,
                     is.finite(relative_duality_gap) & relative_duality_gap <= 1e-6 &
                       original_primal_residual <= gamma + 1e-6 &
                       max_solver_primal_infeasibility <= 1e-8 &
                       max_solver_dual_infeasibility <= 1e-8 &
                       solver_primal_feasible & solver_dual_feasible & basis_valid)
utils::write.csv(table, file.path(dirname(path), "LP_KKT_records.csv"), row.names = FALSE)
jsonlite::write_json(list(source_file = path, actual_LP_rows = nrow(table),
                          passed = all(table$passed), max_relative_duality_gap = max(table$relative_duality_gap),
                          max_original_residual_excess = max(table$original_primal_residual - table$gamma)),
                     file.path(dirname(path), "LP_KKT_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(table[which.max(table$relative_duality_gap), ])
if (!all(table$passed)) stop("An actual LP record failed its original feasibility or KKT check.")

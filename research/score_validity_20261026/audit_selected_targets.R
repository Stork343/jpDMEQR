args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args)) args[1L] else "joint_independent"
root <- file.path("results/scratch/score_validity_20261026", run_id)
files <- sort(list.files(root, pattern = "^result[.]rds$", recursive = TRUE, full.names = TRUE))
records <- list()
for (path in files) {
  fit <- readRDS(path)
  cfg <- fit$request$config
  beta <- cfg$signal * (-1)^(seq_len(cfg$s) + 1L)
  for (f in seq_along(fit$training)) {
    for (j in seq_len(nrow(fit$table))) {
      support <- fit$training[[f]]$graphs[[j]]$support
      target <- support[1L]
      omitted <- setdiff(seq_len(cfg$s), support)
      displacement <- 0
      if (length(omitted)) {
        covariance <- cfg$rho_x^abs(outer(support, support, "-"))
        cross <- cfg$rho_x^abs(outer(support, omitted, "-"))
        displacement <- drop(solve(covariance, cross %*% beta[omitted]))[1L]
      }
      records[[length(records) + 1L]] <- data.frame(case_id = cfg$case_id, replicate = fit$request$replicate,
                                                   coordinate = fit$table$coordinate[j], fold = f,
                                                   nonzero = target <= cfg$s, target_shift = displacement,
                                                   support_size = length(support),
                                                   covered = fit$table$covered[j])
    }
  }
}
table <- do.call(rbind, records)
table$material_shift <- abs(table$target_shift) > 1e-8
table$absolute_shift <- abs(table$target_shift)
directory <- file.path("results/tables/score_validity_20261026", run_id)
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
write.csv(table, file.path(directory, "selected_target_audit.csv"), row.names = FALSE)
print(aggregate(cbind(material_shift, absolute_shift) ~ case_id + nonzero, table, mean), row.names = FALSE)

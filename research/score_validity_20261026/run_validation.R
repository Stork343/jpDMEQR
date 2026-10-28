library(processx)
args <- commandArgs(trailingOnly = TRUE)
validation_run <- if (length(args)) args[1L] else paste0("joint_", format(Sys.time(), "%Y%m%dT%H%M%S"))
registry_path <- if (length(args) >= 2L) args[2L] else "research/score_validity_20261026/evaluation.csv"
if (!grepl("^[A-Za-z0-9_-]+$", validation_run)) stop("Use a simple run identifier.")
root <- file.path("results/scratch/score_validity_20261026", validation_run)
if (dir.exists(root)) stop("This validation run already exists.")
dir.create(root, recursive = TRUE)
launch <- function() {
  processes <- list()
  on.exit({for (process in processes) if (process$is_alive()) process$kill_tree()}, add = TRUE)
  for (worker in 1:4) {
    processes[[worker]] <- processx::process$new(
      file.path(R.home("bin"), "Rscript"),
      c("research/score_validity_20261026/run_joint_evaluation.R", as.character(worker), "4",
        validation_run, "1800000000", "new", registry_path),
      stdout = file.path(root, paste0("worker_", worker, ".log")), stderr = "2>&1")
  }
  repeat {
    alive <- vapply(processes, function(process) process$is_alive(), logical(1))
    for (worker in which(!alive)) {
      status <- processes[[worker]]$get_exit_status()
      if (!is.null(status) && status != 0L) stop(paste("Validation worker", worker, "failed; read its log in", root))
    }
    if (!any(alive)) break
    Sys.sleep(0.5)
  }
}
launch()
result <- processx::run(file.path(R.home("bin"), "Rscript"),
                        c("research/score_validity_20261026/summarise_evaluation.R", validation_run))
cat(result$stdout)

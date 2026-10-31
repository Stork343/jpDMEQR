library(jsonlite)

root <- "results/raw/application/drug_response_20261101/discovery/Dabrafenib"
attempts <- data.frame(run_id = c("JOINT_tau050", "JOINT_tau075"),
                       gap = c(1.55727e-8, 1.82648e-8))
for (j in seq_len(nrow(attempts))) {
  directory <- file.path(root, attempts$run_id[j])
  stopifnot(file.exists(file.path(directory, "request.json")),
            !file.exists(file.path(directory, "fit.rds")),
            !file.exists(file.path(directory, "status.json")))
  write_json(list(status = "numerical_failure", failure_stage = "selection_fit",
                   gap = attempts$gap[j], requested_gap = 1e-8,
                   solver = "joint_profile_pdal_cpp", exit_code = 1L,
                   manuscript_adoption = FALSE),
             file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
}

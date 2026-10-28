root <- normalizePath(".", mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[1]) else 2L
run_id <- if (length(args) > 1L) args[2] else paste0("matched_targets_", format(Sys.time(), "%Y%m%dT%H%M%S"))
if (!is.finite(workers) || workers < 1L) stop("Supply a positive worker count.")
config <- yaml::read_yaml(file.path(root, "config", "application_matched_prospective_20261022.yml"))
jobs <- expand.grid(model = config$models, q = config$random_effect_dimensions,
                      error = config$error_families, stringsAsFactors = FALSE)
jobs$case_id <- sprintf("AM_%s_q%d_%s", jobs$model, jobs$q, jobs$error)
jobs$run_id <- paste0(run_id, "_", jobs$case_id)
output <- file.path(root, "results", "raw", "application_matched", run_id)
if (dir.exists(output)) stop("The matched target launcher directory already exists.")
dir.create(output, recursive = TRUE)
write.csv(jobs, file.path(output, "scheduled_target_jobs.csv"), row.names = FALSE)
cluster <- parallel::makePSOCKcluster(workers)
parallel::clusterExport(cluster, c("root", "output"), envir = environment())
result <- parallel::parLapplyLB(cluster, split(jobs, seq_len(nrow(jobs))), function(job) {
  setwd(root)
  log <- file.path(output, paste0(job$case_id, "_console.log"))
  arguments <- c("scripts/application/11_build_matched_targets_GSE65391.R",
                 paste0("--model=", job$model), paste0("--q=", job$q),
                 paste0("--error=", job$error), paste0("--run-id=", job$run_id), "--engine=streamed")
  status <- system2(file.path(R.home("bin"), "Rscript"), arguments, stdout = log, stderr = log)
  data.frame(case_id = job$case_id, run_id = job$run_id, exit_status = status,
               target_asset = file.path(root, "results", "raw", "application_matched", job$run_id,
                                         "matched_target_asset.rds"), stringsAsFactors = FALSE)
})
parallel::stopCluster(cluster)
result <- do.call(rbind, result)
write.csv(result, file.path(output, "target_job_results.csv"), row.names = FALSE)
if (any(result$exit_status != 0L)) stop("A matched target job terminated; read its actual console and checkpoints.")
if (any(!file.exists(result$target_asset))) stop("A completed target asset is missing.")
jsonlite::write_json(list(run_id = run_id, completed_target_cases = nrow(result),
                          completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                     file.path(output, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
print(result)

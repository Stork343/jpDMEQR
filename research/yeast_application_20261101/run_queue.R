library(data.table)
library(jsonlite)
library(processx)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
implementation_directory <- normalizePath(args[1L], mustWork = TRUE)
workers <- as.integer(args[2L])
stopifnot(is.finite(workers), workers >= 1L, workers <= 4L)
directory <- "results/raw/application/yeast_20261101/screening_20261101"
jobs <- fread(file.path(directory, "jobs.csv"))
log_directory <- "results/logs/yeast_20261101/screening"
dir.create(log_directory, recursive = TRUE, showWarnings = FALSE)
jobs[, `:=`(pid = NA_integer_, exit_code = NA_integer_)]
for (j in seq_len(nrow(jobs))) {
  status_file <- file.path(jobs$directory[j], "status.json")
  if (file.exists(status_file)) {
    status <- read_json(status_file, simplifyVector = TRUE)
    stopifnot(status$status == "completed", file.exists(file.path(jobs$directory[j], "fit.rds")))
    jobs[j, status := "completed"]
  } else if (dir.exists(jobs$directory[j])) stop("未完成的运行需要核查", jobs$directory[j])
}
active <- list()
write_state <- function() {
  fwrite(jobs, file.path(directory, "progress.csv"))
  summary <- jobs[, .N, by = status]
  write_json(list(updated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                   workers = workers, counts = summary, manuscript_adoption = FALSE),
             file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
}
write_state()
repeat {
  changed <- FALSE
  completed <- character()
  for (key in names(active)) {
    child <- active[[key]]
    if (!child$is_alive()) {
      j <- as.integer(key)
      code <- child$get_exit_status()
      result_exists <- file.exists(file.path(jobs$directory[j], "status.json"))
      if (code == 0L && result_exists) {
        saved <- read_json(file.path(jobs$directory[j], "status.json"), simplifyVector = TRUE)
        jobs[j, status := if (saved$status == "completed") "completed" else "failed"]
      } else jobs[j, status := "failed"]
      jobs[j, exit_code := code]
      cat(format(Sys.time(), tz = "UTC", usetz = TRUE), jobs$trait[j], jobs$method[j], jobs$status[j], "\n")
      flush.console()
      completed <- c(completed, key)
      changed <- TRUE
    }
  }
  active[completed] <- NULL
  pending <- which(jobs$status == "pending")
  stop_new <- any(jobs$status == "failed")
  count <- if (stop_new) 0L else min(workers - length(active), length(pending))
  if (count > 0L) {
    for (j in head(pending, count)) {
      log_path <- file.path(log_directory, sprintf("job_%03d.log", jobs$job_id[j]))
      child <- processx::process$new(file.path(R.home("bin"), "Rscript"),
        args = c(file.path(implementation_directory, "run_discovery.R"),
                   jobs$trait[j], jobs$method[j], jobs$run_id[j],
                   as.character(jobs$tau[j]), implementation_directory),
        stdout = log_path, stderr = "2>&1", cleanup = FALSE)
      active[[as.character(j)]] <- child
      jobs[j, `:=`(status = "running", pid = child$get_pid())]
      changed <- TRUE
      cat(format(Sys.time(), tz = "UTC", usetz = TRUE), "started", jobs$job_id[j], jobs$trait[j], jobs$method[j], "\n")
      flush.console()
    }
  }
  if (changed) write_state()
  if (!length(active) && (stop_new || !any(jobs$status == "pending"))) break
  Sys.sleep(2)
}
if (any(jobs$status == "failed")) stop("计算清单中存在失败记录")
cat("计算清单执行完成\n")

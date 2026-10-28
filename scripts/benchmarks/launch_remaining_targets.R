#!/usr/bin/env Rscript
# Batch launcher for the remaining final-scale assets.
# Waits for free CPU cores and launches the next scenario in the queue.
#
# Queue (each needs 10^5 clusters x 4 repeats):
#   profile_mc targets + directions:  (already launched: C06 C07 C14 C15 C16)
#   A-module POP-H directions:         A02 A03 A04 A05 A06 A07 A08 A09 A10 A11 A12
#   D-module POP-H directions:         D01 D02 D03 D04 D05 D06
#   E-module structural audits:        E04 E05 E06
#
# Usage: Rscript scripts/benchmarks/launch_remaining_targets.R
#   --max-cores=<n>    maximum simultaneous Rscript target jobs (default 7,
#                      leaving one core for the pilot/other work)
#   --ids=<csv>        override the queue (e.g. "A02,A03")

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/benchmarks/launch_remaining_targets.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
max_cores <- as.integer(get_arg("max-cores", "7"))
queue <- if (!is.null(get_arg("ids"))) {
  trimws(strsplit(get_arg("ids"), ",", fixed = TRUE)[[1]])
} else {
  c("A02","A03","A04","A05","A06","A07","A08","A09","A10","A11","A12",
    "D01","D02","D03","D04","D05","D06",
    "E04","E05","E06")
}

config_path <- file.path(root, "config", "simulation_main_b200.csv")
n_pop <- 100000L
repeats <- 4L

launch_one <- function(id) {
  cmd <- sprintf(
    "Rscript scripts/simulation/04_build_profile_targets.R --config=%s --n-population=%d --repeats=%d --direction-repeats=%d --experiment-id=%s",
    shQuote(config_path), n_pop, repeats, repeats, id
  )
  cat("LAUNCH", id, "\n")
  system2("Rscript", c(
    "scripts/simulation/04_build_profile_targets.R",
    paste0("--config=", config_path),
    paste0("--n-population=", n_pop),
    paste0("--repeats=", repeats),
    paste0("--direction-repeats=", repeats),
    paste0("--experiment-id=", id)
  ), wait = FALSE)
  Sys.sleep(3)  # let the process start before counting cores
}

count_rscript <- function() {
  out <- tryCatch(system2("powershell", c("-Command",
    "(Get-Process Rscript -ErrorAction SilentlyContinue).Count"),
    stdout = TRUE), error = function(e) "0")
  as.integer(trimws(out[length(out)]))
}

# Launch loop: keep filling up to max_cores until the queue is empty.
remaining <- queue
while (length(remaining)) {
  n <- count_rscript()
  cat(sprintf("running Rscript processes: %d, queue left: %d\n", n, length(remaining)))
  if (n < max_cores) {
    launch_one(remaining[1])
    remaining <- remaining[-1]
  } else {
    Sys.sleep(60)  # wait for a slot
  }
}
cat("Queue complete.\n")

library(data.table)
library(yaml)
library(jsonlite)
library(processx)
library(filelock)
config <- read_yaml("config/application_soynam_20261031.yml")
if (!isTRUE(config$manuscript_adoption)) stop("This application has not been adopted for the manuscript.")
root <- normalizePath(".")
progress_directory <- file.path(config$results_directory, "finalisation")
dir.create(progress_directory, recursive = TRUE, showWarnings = FALSE)
lock <- filelock::lock(file.path(progress_directory, "finalisation.lock"), timeout = 0)
if (is.null(lock)) stop("The application finalisation process is already active.")
progress <- function(stage, detail = list()) {
  write_json(c(list(stage = stage, time = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)), detail),
    file.path(progress_directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), stage, "\n")
  flush.console()
}
options(error = function() {
  write_json(list(status = "failed", error = geterrmessage(), time = as.character(Sys.time())),
    file.path(progress_directory, "failure.json"), auto_unbox = TRUE, pretty = TRUE)
  quit(save = "no", status = 1L, runLast = FALSE)
})
prediction_root <- file.path(config$results_directory, "prediction", config$execution_version)
jobs <- CJ(method = unlist(config$prediction$methods), tau = unlist(config$quantiles))
jobs[, case_id := paste0(method, "_q", sprintf("%02d", round(100 * tau)))]
previous_complete <- -1L
repeat {
  complete <- file.exists(file.path(prediction_root, jobs$case_id, "completion.json"))
  if (sum(complete) != previous_complete) {
    progress("waiting_for_registered_predictions", list(completed = sum(complete), required = nrow(jobs)))
    previous_complete <- sum(complete)
  }
  if (all(complete)) break
  if (any(file.exists(file.path(prediction_root, jobs$case_id[!complete], "failure.rds"))))
    stop("A registered prediction requires numerical inspection.")
  status_file <- file.path(prediction_root, "status.csv")
  if (!file.exists(status_file) || as.numeric(difftime(Sys.time(), file.info(status_file)$mtime, units = "secs")) > 90)
    stop("The prediction queue is not updating its execution record.")
  Sys.sleep(15)
}
run <- function(command, args, directory = root) {
  result <- processx::run(command, args, wd = directory, error_on_status = FALSE,
    echo = FALSE, echo_cmd = FALSE, timeout = Inf)
  if (nzchar(result$stdout)) cat(result$stdout)
  if (nzchar(result$stderr)) cat(result$stderr)
  if (result$status != 0L) stop("A finalisation command failed: ", command, " ", paste(args, collapse = " "))
  invisible(result)
}
for (script in c("summarise_inference.R", "summarise_prediction.R", "descriptive.R",
                 "export_manuscript.R", "plot_results.R", "report_results.R", "check_delivery.R")) {
  progress(script)
  run(file.path(R.home("bin"), "Rscript"), file.path("research/soynam_application_20261031", script))
}
for (document in c("sjs_profile_quantile_main.tex", "sjs_profile_quantile_supplement.tex")) {
  progress("compile_manuscript", list(document = document))
  run("latexmk", c("-xelatex", "-silent", "-interaction=nonstopmode", "-halt-on-error", document), file.path(root, "manuscript"))
}
progress("verify_pdf")
Sys.setenv(PYTHONPATH = file.path(root, "results/scratch/python_dependencies"),
  PYTHONPYCACHEPREFIX = file.path(root, "results/scratch/runtime/pycache"))
run("python3", "research/soynam_application_20261031/check_pdf.py")
progress("backup_delivery")
run(file.path(R.home("bin"), "Rscript"), "research/soynam_application_20261031/backup_delivery.R")
generated <- list.files("manuscript/generated/soynam_profile", full.names = TRUE)
tables <- list.files("results/tables/soynam_20261031", pattern = "[.]csv$", full.names = TRUE)
preflight <- file.path("results/preflight/soynam_20261031",
  c("delivery_checks.csv", "delivery_manifest.json", "pdf_checks.json", "session.txt"))
owned <- c(generated, tables, preflight, "research/soynam_application_20261031/RESULTS_zh.md")
if (!all(file.exists(owned))) stop("A final publication artifact is missing.")
progress("save_git_results")
branch <- trimws(run("git", c("branch", "--show-current"))$stdout)
if (branch != "recovery/research-20261028") stop("The current branch differs from the registered research branch.")
run("git", c("-c", "core.whitespace=blank-at-eol,blank-at-eof,space-before-tab,cr-at-eol", "diff", "--check", "--", owned))
run("git", c("add", "-f", "--", owned))
run("git", c("commit", "--only", "-m", "完成全部大豆预测比较与论文交付核验", "--", owned))
progress("push_github")
run("git", c("-c", "core.sshCommand=ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=15",
  "-c", "url.ssh://git@ssh.github.com:443/.insteadOf=https://github.com/", "push", "origin", "recovery/research-20261028"))
commit <- run("git", c("rev-parse", "HEAD"))$stdout
write_json(list(status = "complete", inference_cases = 21L, prediction_cases = 9L,
  commit = trimws(commit), completed = as.character(Sys.time())),
  file.path(progress_directory, "completion.json"), auto_unbox = TRUE, pretty = TRUE)
progress("complete", list(commit = trimws(commit)))
filelock::unlock(lock)

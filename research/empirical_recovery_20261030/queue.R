empirical_queue_root <- "results/raw/application/empirical_recovery_20261030"
empirical_comparator_root <- "results/raw/application_comparison/empirical_recovery_20261030"

empirical_registry <- function() {
  routes <- c("reference", "processed_median_across_mapped_probes",
    "official_R1_R2_normexp_training_quantile")
  jobs <- list()
  add <- function(phase, processing, script, args, output, results, sequence) {
    jobs[[length(jobs) + 1L]] <<- list(phase = phase, processing = processing,
      script = script, arguments = args, output = output,
      expected_results = as.integer(results), sequence = sequence)
  }
  script_root <- "research/joint_profile_manuscript_20261026"
  roles <- expand.grid(repetition = 1:5, fold = 1:5)
  for (route in routes) {
    for (i in seq_len(nrow(roles))) {
      r <- roles$repetition[i]
      f <- roles$fold[i]
      name <- sprintf("r%02d_f%02d", r, f)
      mode <- if (route == "reference") "intercept_only" else "standard"
      add("prediction", route, file.path(script_root, "run_prediction_fold.R"),
        c(route, r, f, mode, "resume"),
        file.path("results/raw/application_prediction/joint_profile_20261026", route, name),
        if (route == "reference") 18L else 6L, i)
      if (route == "reference") add("slope_prediction", route,
        file.path(script_root, "run_prediction_fold.R"), c(route, r, f, "slope_only", "resume"),
        file.path("results/raw/application_prediction/joint_profile_20261026/reference_conic_slope", name), 6L, i)
    }
    for (i in 1:20) {
      role <- sprintf("split%02d_swap%d", (i + 1L) %/% 2L, (i - 1L) %% 2L)
      add("exploratory", route, file.path(script_root, "run_exploratory_role.R"), c(route, role),
        file.path("results/raw/application_exploratory/joint_profile_20261026", route, role), 6L, i)
    }
    if (route != "reference") {
      cases <- expand.grid(model = c("A", "B"), tau = c(.25, .5, .75))
      for (i in seq_len(nrow(cases))) add("processing_modules", route,
        file.path(script_root, "run_processing_module.R"),
        c(route, cases$model[i], cases$tau[i]),
        file.path("results/raw/application/joint_profile_processing_20261026", route,
          paste0(cases$model[i], "_tau", round(100 * cases$tau[i]))), 1L, i)
    }
  }
  comparator <- expand.grid(method = c("LQMM", "POOLED-QR-LASSO", "QGEE-SCAD"),
    model = c("A", "B"), tau = c(.25, .5, .75), processing = routes,
    repetition = 1:5, fold = 1:5, stringsAsFactors = FALSE)
  comparator$job_id <- seq_len(nrow(comparator))
  for (i in seq_len(nrow(comparator))) add("comparators", comparator$processing[i],
    "research/submission_completion_20261027/run_application_comparator.R",
    c(empirical_comparator_root, i), file.path(empirical_comparator_root, sprintf("job_%04d", i)), 1L, i)
  add("comparator_checks", "reference", "research/submission_completion_20261027/check_application.R",
    character(), "results/scratch/submission_application_20261027/checks", 0L, 0L)
  table <- do.call(rbind, lapply(jobs, function(x) as.data.frame(x[setdiff(names(x), "arguments")], stringsAsFactors = FALSE)))
  table$arguments <- I(lapply(jobs, `[[`, "arguments"))
  phase <- match(table$phase, c("exploratory", "prediction", "comparator_checks", "processing_modules", "comparators", "slope_prediction"))
  table <- table[order(table$sequence, phase, match(table$processing, routes)), ]
  rownames(table) <- NULL
  table$job_id <- seq_len(nrow(table))
  if (sum(table$expected_results[table$phase %in% c("prediction", "slope_prediction")]) != 900L ||
      sum(table$expected_results[table$phase == "exploratory"]) != 360L ||
      sum(table$expected_results[table$phase == "processing_modules"]) != 12L || nrow(comparator) != 1350L)
    stop("The empirical registry does not match the complete analysis design.")
  list(jobs = table, comparator = comparator)
}

empirical_job_ready <- function(job) {
  manifest <- if (job$processing == "reference") "input_manifest.json" else
    paste0(job$processing, "_manifest.json")
  path <- file.path("results/preflight/empirical_recovery_20261030", manifest)
  if (!file.exists(path)) return(FALSE)
  if (!isTRUE(jsonlite::read_json(path)$passed)) stop("The processing checks have failed: ", job$processing)
  if (job$phase == "comparator_checks") return(file.exists(
    "results/raw/application_prediction/joint_profile_20261026/reference/r01_f01/case_01/result.rds"))
  if (job$phase == "comparators") return(file.exists(
    "results/scratch/submission_application_20261027/checks/clinical_units.csv"))
  TRUE
}

empirical_job_complete <- function(job) {
  if (job$phase == "comparator_checks") return(file.exists(file.path(job$output, "clinical_units.csv")))
  if (job$phase %in% c("comparators", "processing_modules")) return(file.exists(file.path(job$output, "result.rds")))
  if (!file.exists(file.path(job$output, "completed.rds"))) return(FALSE)
  completion <- readRDS(file.path(job$output, "completed.rds"))
  length(list.files(job$output, "^result[.]rds$", recursive = TRUE)) == job$expected_results &&
    completion$completed == job$expected_results
}

empirical_backup <- function(jobs) {
  base <- file.path(path.expand("~"), "jpDMEQR-backups", "empirical_recovery_20261030")
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  files <- unique(unlist(lapply(jobs$output, function(path) {
    if (!dir.exists(path)) return(character())
    list.files(path, "^(result|fit_[0-9]+|request|completed)[.]rds$", full.names = TRUE, recursive = TRUE)
  })))
  files <- c(files, file.path(empirical_queue_root, c("request.rds", "status.csv")),
    file.path(empirical_comparator_root, "request.rds"))
  files <- files[file.exists(files)]
  manifest <- file.path(empirical_queue_root, "backup_files.txt")
  writeLines(sort(unique(files)), manifest)
  status <- system2("rsync", c("-a", shQuote(paste0("--files-from=", normalizePath(manifest))),
    shQuote(paste0(normalizePath("."), "/")), shQuote(paste0(normalizePath(base), "/"))))
  if (status != 0L) stop("The empirical output backup failed.")
  submission_atomic_csv(data.frame(files = length(files), saved = as.character(Sys.time())),
    file.path(base, "backup_status.csv"))
}

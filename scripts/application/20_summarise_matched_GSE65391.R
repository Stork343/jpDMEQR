root <- normalizePath(".", mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L) || !args[3] %in% c("progress", "complete")) {
  stop("Supply the matched run ID, a new summary ID, progress or complete, and an optional source-retry run ID.")
}
input <- file.path(root, "results", "raw", "application_matched", args[1])
output <- file.path(root, "results", "tables", "application", args[2])
if (dir.exists(output)) stop("The immutable matched summary directory already exists.")
scheduled <- read.csv(file.path(input, "scheduled_replicates.csv"), stringsAsFactors = FALSE)
if (nrow(scheduled) != 4800L || anyDuplicated(scheduled[c("case_id", "replicate")])) {
  stop("The complete 24 by 200 matched schedule is required.")
}
registry <- unique(scheduled[c("case_id", "model", "q", "error_family", "tau")])
if (nrow(registry) != 24L) stop("The complete matched case registry is required.")
process_path <- file.path(input, "job_results.csv")
processes <- if (file.exists(process_path)) read.csv(process_path, stringsAsFactors = FALSE) else data.frame()
retry_input <- if (length(args) == 4L) file.path(root, "results", "raw", "application_matched", args[4]) else NULL
retry_jobs <- retry_processes <- data.frame()
if (!is.null(retry_input)) {
  retry_manifest <- jsonlite::read_json(file.path(retry_input, "run_manifest.json"), simplifyVector = TRUE)
  if (!identical(retry_manifest$original_run, args[1])) stop("The source-retry run references another original study.")
  retry_jobs <- read.csv(file.path(retry_input, "scheduled_replicates.csv"), stringsAsFactors = FALSE)
  path <- file.path(retry_input, "job_results.csv")
  if (file.exists(path)) retry_processes <- read.csv(path, stringsAsFactors = FALSE)
}
bind_tables <- function(tables) {
  if (!length(tables)) return(data.frame())
  columns <- unique(unlist(lapply(tables, names), use.names = FALSE))
  tables <- lapply(tables, function(table) {
    for (name in setdiff(columns, names(table))) table[[name]] <- NA
    table[columns]
  })
  do.call(rbind, tables)
}
accounting <- scheduled[c("case_id", "model", "q", "error_family", "tau", "replicate")]
accounting$status <- "not_started"
accounting$failure_stage <- NA_character_
accounting$failure_message <- NA_character_
accounting$result_run <- args[1]
accounting$original_source_failure <- FALSE
coordinate <- list()
for (i in seq_len(nrow(scheduled))) {
  job <- scheduled[i, ]
  directory <- file.path(input, job$case_id, sprintf("replicate_%04d", job$replicate))
  result_input <- input
  result_processes <- processes
  original_finished <- if (nrow(processes)) processes[processes$case_id == job$case_id & processes$replicate == job$replicate, ] else data.frame()
  before_workflow <- nrow(original_finished) == 1L && original_finished$exit_status != 0L &&
    !file.exists(file.path(directory, "run_request.rds"))
  accounting$original_source_failure[i] <- before_workflow
  if (before_workflow && nrow(retry_jobs)) {
    retry <- retry_jobs[retry_jobs$case_id == job$case_id & retry_jobs$replicate == job$replicate, ]
    if (nrow(retry) > 1L) stop("A source retry was registered more than once.")
    if (nrow(retry) == 1L) {
      if (!identical(retry$target_asset, job$target_asset)) stop("A source retry changed its population target asset.")
      result_input <- retry_input
      result_processes <- retry_processes
      directory <- file.path(retry_input, job$case_id, sprintf("replicate_%04d", job$replicate))
      accounting$result_run[i] <- args[4]
    }
  }
  complete <- file.exists(file.path(directory, "completed_manifest.json"))
  failure <- file.exists(file.path(directory, "coordinate_metrics_failed.csv")) &&
    file.exists(file.path(directory, "failure_condition.rds"))
  if (complete && failure) stop("A replicate contains both completion and failure records.")
  if (complete) {
    table <- read.csv(file.path(directory, "coordinate_metrics.csv"), stringsAsFactors = FALSE)
    expected <- if (job$model == "A") 16L else 17L
    if (nrow(table) != expected || anyDuplicated(table$coordinate) ||
        any(table$status != "ok") || any(table$case_id != job$case_id) ||
        any(table$replicate != job$replicate) || sum(grepl("^module__", table$coordinate)) != 5L ||
        any(!is.finite(table$beta_hat) | !is.finite(table$sandwich_se) | !is.finite(table$jackknife_se)) ||
        any(table$sandwich_se <= 0 | table$jackknife_se <= 0)) {
      stop("A completed matched coordinate table is invalid.")
    }
    accounting$status[i] <- "ok"
    coordinate[[length(coordinate) + 1L]] <- table
  } else if (failure) {
    table <- read.csv(file.path(directory, "coordinate_metrics_failed.csv"), stringsAsFactors = FALSE)
    accounting$status[i] <- "failed"
    accounting$failure_message[i] <- unique(table$failure_message)
    accounting$failure_stage[i] <- if ("failure_stage" %in% names(table)) unique(table$failure_stage) else "unclassified_saved_condition"
  } else if (file.exists(file.path(directory, "run_request.rds"))) {
    accounting$status[i] <- "started_incomplete"
  }
  if (nrow(result_processes)) {
    finished <- result_processes[result_processes$case_id == job$case_id & result_processes$replicate == job$replicate, ]
    if (nrow(finished) > 1L) stop("A scheduled process was recorded more than once.")
    if (nrow(finished) == 1L && !accounting$status[i] %in% c("ok", "failed")) {
      if (finished$exit_status == 0L) stop("A successful process has no completed workflow artifact.")
      accounting$status[i] <- "source_failure_pending_retry"
      accounting$failure_stage[i] <- "before_complete_workflow_record"
      log <- file.path(result_input, "console", sprintf("%s_rep%04d.log", job$case_id, job$replicate))
      accounting$failure_message[i] <- paste(tail(readLines(log, warn = FALSE), 2L), collapse = " ")
    }
  }
}
if (args[3] == "complete" &&
    (!file.exists(file.path(input, "completed_manifest.json")) ||
     (!is.null(retry_input) && !file.exists(file.path(retry_input, "completed_manifest.json"))) ||
     any(!accounting$status %in% c("ok", "failed")))) {
  stop("A complete summary requires every scheduled replicate to have an actual result.")
}
coordinates <- bind_tables(coordinate)
dir.create(output, recursive = TRUE)
write.csv(accounting, file.path(output, "matched_replication_accounting.csv"), row.names = FALSE)
attempts <- list()
for (directory in c(input, retry_input)) {
  path <- file.path(directory, "job_results.csv")
  if (!file.exists(path)) next
  table <- read.csv(path, stringsAsFactors = FALSE)
  table$attempt_run <- basename(directory)
  table$workflow_request_present <- file.exists(file.path(directory, table$case_id,
                                                          sprintf("replicate_%04d", table$replicate), "run_request.rds"))
  attempts[[length(attempts) + 1L]] <- table
}
write.csv(bind_tables(attempts), file.path(output, "matched_all_process_attempts.csv"), row.names = FALSE)
write.csv(coordinates, file.path(output, "matched_coordinate_results.csv"), row.names = FALSE)
fraction <- function(covered) {
  n <- length(covered)
  if (!n) return(c(estimate = NA_real_, mcse = NA_real_, lower = NA_real_, upper = NA_real_))
  p <- mean(covered)
  interval <- stats::binom.test(sum(covered), n, conf.level = .95)$conf.int
  c(estimate = p, mcse = sqrt(p * (1 - p) / n), lower = interval[1], upper = interval[2])
}
target_tables <- list()
for (asset in unique(scheduled$target_asset)) {
  table <- read.csv(file.path(dirname(asset), "population_targets.csv"), stringsAsFactors = FALSE)
  job <- scheduled[match(asset, scheduled$target_asset), ]
  table$model <- job$model
  table$q <- job$q
  table$error_family <- job$error_family
  target_tables[[length(target_tables) + 1L]] <- table
}
targets <- do.call(rbind, target_tables)
write.csv(targets, file.path(output, "matched_population_targets.csv"), row.names = FALSE)
summary <- list()
scenario <- list()
for (i in seq_len(nrow(registry))) {
  case <- registry[i, ]
  counts <- accounting[accounting$case_id == case$case_id, ]
  target <- targets[targets$model == case$model & targets$q == case$q &
                      targets$error_family == case$error_family & targets$tau == case$tau, ]
  for (j in seq_len(nrow(target))) {
    for (source in c("sandwich", "jackknife")) {
      values <- if (nrow(coordinates)) {
        coordinates[coordinates$case_id == case$case_id & coordinates$coordinate == target$coordinate[j], ]
      } else data.frame()
      n <- nrow(values)
      error <- if (n) values$beta_hat - target$target[j] else numeric()
      se <- if (n) values[[paste0(source, "_se")]] else numeric()
      covered <- if (n) values[[paste0(source, "_covered")]] else logical()
      probability <- fraction(covered)
      available <- if (all(counts$status %in% c("ok", "failed"))) {
        fraction(c(covered, rep(FALSE, sum(counts$status == "failed"))))
      } else fraction(logical())
      empirical_sd <- if (n > 1L) stats::sd(values$beta_hat) else NA_real_
      rmse <- if (n) sqrt(mean(error^2)) else NA_real_
      simultaneous <- if (n && grepl("^module__", target$coordinate[j])) {
        fraction(values[[paste0(source, "_simultaneous_covered")]])
      } else fraction(logical())
      row <- data.frame(
        case, coordinate = target$coordinate[j], interval = source, target = target$target[j],
        target_mcse = target$target_mcse[j], target_source = target$target_source[j],
        scheduled = 200L, completed = sum(counts$status == "ok"), failed = sum(counts$status == "failed"),
        pending = sum(!counts$status %in% c("ok", "failed")), n_valid = n,
        bias = if (n) mean(error) else NA_real_, bias_mcse = if (n > 1L) stats::sd(error) / sqrt(n) else NA_real_,
        rmse = rmse, rmse_mcse = if (n > 1L && rmse > 0) stats::sd(error^2) / (2 * rmse * sqrt(n)) else NA_real_,
        empirical_sd = empirical_sd,
        mean_se = if (n) mean(se) else NA_real_,
        mean_se_mcse = if (n > 1L) stats::sd(se) / sqrt(n) else NA_real_,
        se_sd_ratio = if (is.finite(empirical_sd) && empirical_sd > 0) mean(se) / empirical_sd else NA_real_,
        target_mcse_over_mean_se = if (n) target$target_mcse[j] / mean(se) else NA_real_,
        conditional_coverage = probability["estimate"], conditional_coverage_mcse = probability["mcse"],
        conditional_coverage_lower = probability["lower"], conditional_coverage_upper = probability["upper"],
        availability_coverage = available["estimate"], availability_coverage_mcse = available["mcse"],
        availability_coverage_lower = available["lower"], availability_coverage_upper = available["upper"],
        simultaneous_family_size = if (grepl("^module__", target$coordinate[j])) 30L else NA_integer_,
        simultaneous_conditional_coverage = simultaneous["estimate"],
        simultaneous_conditional_coverage_mcse = simultaneous["mcse"],
        simultaneous_conditional_coverage_lower = simultaneous["lower"],
        simultaneous_conditional_coverage_upper = simultaneous["upper"],
        changed_probe_replications = if (n) sum(values$changed_gene_probes > 0) else 0L,
        mean_changed_gene_probes = if (n) mean(values$changed_gene_probes) else NA_real_,
        mean_changed_gene_probes_mcse = if (n > 1L) stats::sd(values$changed_gene_probes) / sqrt(n) else NA_real_,
        stringsAsFactors = FALSE
      )
      summary[[length(summary) + 1L]] <- row
    }
  }
  if (nrow(coordinates)) {
    module <- coordinates[coordinates$case_id == case$case_id & grepl("^module__", coordinates$coordinate), ]
    if (nrow(module)) for (source in c("sandwich", "jackknife")) {
      replicate_mean <- tapply(module[[paste0(source, "_covered")]], module$replicate, mean)
      scenario[[length(scenario) + 1L]] <- data.frame(
        case, interval = source, n_complete_replications = length(replicate_mean),
        module_coverage = mean(replicate_mean),
        module_coverage_mcse = if (length(replicate_mean) > 1L) stats::sd(replicate_mean) / sqrt(length(replicate_mean)) else NA_real_,
        scheduled = 200L, failed = sum(counts$status == "failed"),
        availability_module_coverage = if (all(counts$status %in% c("ok", "failed"))) sum(replicate_mean) / 200 else NA_real_
      )
    }
  }
}
write.csv(bind_tables(summary), file.path(output, "matched_coordinate_summary.csv"), row.names = FALSE)
write.csv(bind_tables(scenario), file.path(output, "matched_module_scenario_summary.csv"), row.names = FALSE)
families <- list()
for (q in c(1L, 2L)) for (error_family in c("normal", "t3")) {
  jobs <- accounting[accounting$q == q & accounting$error_family == error_family, ]
  values <- if (nrow(coordinates)) {
    coordinates[coordinates$q == q & coordinates$error_family == error_family & grepl("^module__", coordinates$coordinate), ]
  } else data.frame()
  complete_replicates <- if (nrow(values)) {
    counts <- table(values$replicate)
    as.integer(names(counts)[counts == 30L])
  } else integer()
  for (source in c("sandwich", "jackknife")) {
    events <- vapply(complete_replicates, function(replication) {
      all(values[[paste0(source, "_simultaneous_covered")]][values$replicate == replication])
    }, logical(1))
    probability <- fraction(events)
    failed_families <- length(unique(jobs$replicate[jobs$status == "failed"]))
    availability <- if (all(jobs$status %in% c("ok", "failed"))) {
      fraction(c(events, rep(FALSE, 200L - length(events))))
    } else fraction(logical())
    families[[length(families) + 1L]] <- data.frame(
      q = q, error_family = error_family, interval = source, family_size = 30L,
      scheduled_families = 200L, complete_families = length(events), failed_families = failed_families,
      conditional_coverage = probability["estimate"], conditional_coverage_mcse = probability["mcse"],
      conditional_coverage_lower = probability["lower"], conditional_coverage_upper = probability["upper"],
      availability_coverage = availability["estimate"], availability_coverage_mcse = availability["mcse"],
      availability_coverage_lower = availability["lower"], availability_coverage_upper = availability["upper"]
    )
  }
}
write.csv(bind_tables(families), file.path(output, "matched_module_family30_summary.csv"), row.names = FALSE)
failure_groups <- accounting[c("case_id", "status", "failure_stage")]
failure_groups$failure_stage[is.na(failure_groups$failure_stage)] <- "not_applicable"
failure <- aggregate(list(replications = accounting$replicate),
                       by = failure_groups, FUN = length)
write.csv(failure, file.path(output, "matched_failure_summary.csv"), row.names = FALSE)
jsonlite::write_json(list(
  input_run = args[1], summary_id = args[2], stage = args[3], scheduled = 4800L,
  source_retry_run = if (length(args) == 4L) args[4] else NA_character_,
  original_source_failures = sum(accounting$original_source_failure),
  successful = sum(accounting$status == "ok"), failed = sum(accounting$status == "failed"),
  incomplete = sum(!accounting$status %in% c("ok", "failed")),
  source_manifest = normalizePath(file.path(input, "run_manifest.json")),
  summary_source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "summary_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

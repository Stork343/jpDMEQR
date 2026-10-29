library(data.table)

root <- "results/raw/application/joint_profile_20261026"
request <- readRDS(file.path(root, "dispatcher.rds"))
completed <- readRDS(file.path(root, "completed.rds"))
jobs <- request$jobs
if (nrow(jobs) != 40L || completed$successful != 40L || completed$failed != 0L)
  stop("The reference-processing module workflows are incomplete.")
modules <- c("M1.2", "M3.4", "M4.11", "M5.12", "M5.15")
records <- checks <- vector("list", nrow(jobs))
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, , drop = FALSE]
  path <- file.path(root, sprintf("job_%03d/result.rds", job$job_id))
  result <- readRDS(path)
  if (!isTRUE(all.equal(result$job[, names(job), drop = FALSE], job, check.attributes = FALSE)))
    stop("A module result differs from its registered job.")
  diagnostics <- result$diagnostics
  deleted <- result$leave_one_estimates
  n <- nrow(deleted)
  if (n != 129L || nrow(diagnostics) != n + 1L ||
      sum(is.na(diagnostics$omitted_subject)) != 1L ||
      anyDuplicated(diagnostics$omitted_subject[!is.na(diagnostics$omitted_subject)]) ||
      diagnostics$n_clusters[1L] != n || any(diagnostics$n_clusters[-1L] != n - 1L) ||
      any(abs(diagnostics$h - diagnostics$h[1L]) > 1e-12))
    stop("The complete subject-deletion design differs.")
  table <- as.data.table(result$table)
  position <- match(table$coordinate, colnames(deleted))
  if (anyNA(position) || any(!is.finite(deleted))) stop("The deletion coefficients are incomplete.")
  centered <- sweep(deleted, 2L, colMeans(deleted))
  computed_se <- sqrt((n - 1) / n * colSums(centered^2))[position]
  error <- max(abs(computed_se - table$workflow_se))
  if (!is.finite(error) || error > 1e-10) stop("The stored workflow covariance differs from its deletion estimates.")
  table <- table[startsWith(coordinate, "module__")]
  if (nrow(table) != 5L || !setequal(table$coordinate, paste0("module__", modules)))
    stop("The declared immune-module family differs.")
  family_size <- if (job$tau == .1) 10L else 30L
  critical <- qnorm(1 - .05 / (2 * family_size))
  score_lower <- table$beta_hat - critical * table$se
  score_upper <- table$beta_hat + critical * table$se
  workflow_lower <- table$beta_hat - critical * table$workflow_se
  workflow_upper <- table$beta_hat + critical * table$workflow_se
  interval_error <- max(abs(c(score_lower - table$simultaneous_lower,
    score_upper - table$simultaneous_upper,
    workflow_lower - table$workflow_simultaneous_lower,
    workflow_upper - table$workflow_simultaneous_upper)))
  if (any(table$family_size != family_size) || !is.finite(interval_error) || interval_error > 1e-10)
    stop("The simultaneous intervals differ from the declared family allocation.")
  table[, `:=`(module = sub("^module__", "", coordinate), source_file = path,
    score_excludes_zero = score_lower * score_upper > 0,
    workflow_excludes_zero = workflow_lower * workflow_upper > 0)]
  records[[i]] <- table
  checks[[i]] <- data.table(job_id = job$job_id, fitted_subjects = n,
    fitted_workflows = nrow(diagnostics), module_intervals = nrow(table),
    workflow_se_difference = error, simultaneous_interval_difference = interval_error)
}
records <- rbindlist(records)
checks <- rbindlist(checks)
reference <- records[specification == "reference"]
setorder(reference, model, tau, module)
if (nrow(reference) != 30L || anyDuplicated(reference[, .(model, tau, module)]))
  stop("The reference module family is incomplete or duplicated.")
sensitivity <- records[specification != "reference", .(workflows = uniqueN(source_file),
  intervals = .N, score_exclusions = sum(score_excludes_zero),
  workflow_exclusions = sum(workflow_excludes_zero)),
  by = .(specification, lambda_gamma, h_multiplier, random_slope)]
output <- "results/tables/submission_20261029/reference_modules"
generated <- "manuscript/generated/joint_profile"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(records, file.path(output, "module_intervals.csv"))
fwrite(reference, file.path(output, "reference_family.csv"))
fwrite(sensitivity, file.path(output, "module_sensitivity_summary.csv"))
fwrite(checks, file.path(output, "checks.csv"))
number <- function(value) formatC(ifelse(abs(value) < .0005, 0, value), digits = 3L, format = "f")
interval <- function(lower, upper) paste0("[", number(lower), ", ", number(upper), "]")
table <- data.frame(Module = reference$module, Model = reference$model, Quantile = reference$tau,
  Estimate = number(reference$beta_hat),
  Score = interval(reference$simultaneous_lower, reference$simultaneous_upper),
  Workflow = interval(reference$workflow_simultaneous_lower, reference$workflow_simultaneous_upper))
names(table) <- c("Module", "Model", "$\\tau$", "Estimate", "Score interval", "Workflow interval")
content <- knitr::kable(table, format = "latex", booktabs = TRUE,
  longtable = TRUE, escape = FALSE, linesep = "",
  caption = "The complete reference module family", label = "joint-modules")
content <- kableExtra::kable_styling(content, latex_options = "repeat_header")
writeLines(as.character(content),
  file.path(generated, "module_intervals.tex"))
jsonlite::write_json(list(source = root, workflows = nrow(jobs), fits = sum(checks$fitted_workflows),
  module_intervals = nrow(records), reference_family_size = nrow(reference),
  maximum_workflow_se_difference = max(checks$workflow_se_difference),
  maximum_simultaneous_interval_difference = max(checks$simultaneous_interval_difference),
  generated = as.character(Sys.time())), file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Verified", nrow(jobs), "module workflows and", sum(checks$fitted_workflows), "fits\n")

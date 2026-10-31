library(data.table)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
drug <- args[1L]
joint_runs <- strsplit(args[2L], ",", fixed = TRUE)[[1L]]
reference_run <- args[3L]
root <- file.path("results/raw/application/drug_response_20261101/discovery", drug)
read_loss <- function(run) {
  directory <- file.path(root, run)
  verification <- read_json(file.path(directory, "artifact_checks.json"), simplifyVector = TRUE)
  stopifnot(verification$status == "prediction_artifacts_verified", verification$training_overlap == 0L)
  fread(file.path(directory, "patient_loss.csv"))
}
joint <- rbindlist(lapply(joint_runs, read_loss))
reference <- read_loss(reference_run)
stopifnot(!anyDuplicated(joint[, .(patient_id, tau)]))
paired <- merge(joint, reference, by = c("patient_id", "tau"), suffixes = c("_joint", "_reference"))
stopifnot(nrow(paired) == nrow(joint), nrow(paired) == nrow(reference))
set.seed(20261106L)
summarise <- function(values, label) {
  n <- nrow(values)
  weights <- rmultinom(10000L, n, rep(1 / n, n)) / n
  difference <- values$loss_reference - values$loss_joint
  bootstrap <- drop(crossprod(weights, difference)) /
    drop(crossprod(weights, values$loss_reference))
  data.table(metric = label, patients = n, joint_loss = mean(values$loss_joint),
    reference_loss = mean(values$loss_reference),
    relative_gain = mean(difference) / mean(values$loss_reference),
    lower = quantile(bootstrap, 0.025), upper = quantile(bootstrap, 0.975))
}
table <- rbindlist(lapply(sort(unique(paired$tau)), function(tau_value)
  summarise(paired[tau == tau_value], paste0("tau_", tau_value))))
overall <- paired[, .(loss_joint = mean(loss_joint), loss_reference = mean(loss_reference)), by = patient_id]
table <- rbind(table, summarise(overall, "average_pinball"))
directory <- file.path("results/preflight/empirical_selection_20261101/comparisons", drug)
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
fwrite(table, file.path(directory, "paired_prediction.csv"))
write_json(list(
  drug = drug, phase = "discovery", reference = reference_run, joint_runs = joint_runs,
  comparison_unit = "patient_id", bootstrap_repetitions = 10000L,
  prediction_advantage_demonstrated = table[metric == "average_pinball", lower] > 0,
  scientific_increment_validated = FALSE, formal_validation_used = FALSE,
  manuscript_adoption = FALSE
), file.path(directory, "decision.json"), auto_unbox = TRUE, pretty = TRUE)
print(table)

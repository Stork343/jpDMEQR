library(data.table)
source("research/submission_completion_20261027/dispatch_state.R")
output <- "results/preflight/sqr_acceleration_20261030"
binding <- jsonlite::read_json(file.path(output, "release_binding.json"))
code_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
if (!isTRUE(binding$passed) || binding$implementation_commit != code_commit)
  stop("The committed production entry point has not passed its numerical check.")
validated_paths <- c("R", "src", "config", "research/submission_completion_20261027",
  "scripts/simulation", "scripts/00_install_dependencies.R", "DESCRIPTION")
if (system2("git", c("diff", "--quiet", code_commit, "--", validated_paths)) != 0L)
  stop("The numerical implementation has uncommitted changes.")
roots <- file.path("results/scratch/sqr_acceleration_20261030",
  c("full_A01_1_certified_bounds", "full_A04_556_certified_bounds",
    "full_A04_779_certified_bounds", "full_A05_737_certified_bounds"))
paired <- rbindlist(lapply(roots, function(root) {
  x <- fread(file.path(root, "checks.csv"))
  x[, source := root]
  x
}))
if (nrow(paired) != 4L || any(paired$interval_difference > 2e-5) ||
    any(paired$coefficient_difference > 2e-5) || any(paired$tuning_difference > 1e-12))
  stop("The complete SQR procedure validations are incomplete.")
registered <- fread("results/scratch/sqr_acceleration_20261030/registered_A07_1/checks.csv")
if (nrow(registered) != 1L || registered$p != 2000L ||
    registered$interval_difference > 2e-5 || registered$maximum_validation_loss_difference > 2e-6)
  stop("The large-design validation is incomplete.")
clime <- fread(file.path(output, "compiled_clime_checks.csv"))
box <- fread(file.path(output, "box_loss_checks.csv"))
if (nrow(clime) != 6L || nrow(box) != 4L ||
    any(clime$coefficient_error > 2e-5, na.rm = TRUE) ||
    any(box$certified_lower_bound > box$independent_minimum + 1e-8))
  stop("The solver or validation-bound checks are incomplete.")
families <- c("mixed", "pooled", "sqr", "qgee", "bias_adjusted")
states <- expand.grid(rep(list(c(FALSE, TRUE)), length(families)))
allocation_checks <- vector("list", nrow(states))
for (i in seq_len(nrow(states))) {
  completed <- setNames(as.logical(states[i, ]), families)
  allocation <- submission_comparator_allocations(completed)
  if (sum(allocation) > 12L || any(allocation[completed] != 0L)) stop("The process allocation exceeds its budget.")
  for (j in which(!completed)) {
    later <- completed
    later[j] <- TRUE
    changed <- submission_comparator_allocations(later)
    if (any(changed[!later] < allocation[!later])) stop("A completion reduces another queue's allocation.")
  }
  allocation_checks[[i]] <- as.data.table(as.list(c(setNames(completed, paste0(families, "_complete")), allocation)))
}
fwrite(paired, file.path(output, "procedure_checks.csv"))
fwrite(registered, file.path(output, "large_design_checks.csv"))
fwrite(rbindlist(allocation_checks), file.path(output, "allocation_checks.csv"))
jsonlite::write_json(list(passed = TRUE, numerical_specification = "submission_comparator_numerics_20261030",
  validated_code_commit = code_commit, validated_paths = validated_paths,
  complete_recorded_datasets = nrow(paired), largest_verified_dimension = registered$p,
  maximum_interval_difference = max(paired$interval_difference, registered$interval_difference),
  maximum_tuning_difference = max(paired$tuning_difference),
  clime_comparisons = nrow(clime), independent_box_comparisons = nrow(box),
  allocation_states = nrow(states), generated = as.character(Sys.time())),
  file.path(output, "acceptance.json"), auto_unbox = TRUE, pretty = TRUE)
print(paired)
print(registered)

library(data.table)
library(jsonlite)
library(yaml)
source("research/clinical_search_20261105/panel.R")
source("research/charite_application_20261105/features.R")
id <- commandArgs(trailingOnly = TRUE)
stopifnot(length(id) == 1L)
root <- file.path("results/raw/application", id)
charite <- startsWith(id, "charite_")
cfg <- read_yaml(file.path(root, if (charite) "linear_implementation" else "implementation", "config.yml"))
panel <- if (charite) ch_panel(cfg$input_directory) else cp_panel(cfg$input_directory)
extract <- if (charite) ch_recipe else cp_recipe
fields <- if (charite) c("genes", "missing_genes", "bmi_median", "impute", "columns", "center", "scale") else
  c("features", "missing_features", "impute", "columns", "center", "scale")
new_feature_key <- if (charite) "genes" else "features"
folders <- if (charite) root else file.path(root, paste0("outer_", seq_len(cfg$outer_folds)))
records <- list()
for (fold in seq_along(folders)) {
  folder <- folders[fold]
  training <- fread(file.path(if (charite) cfg$input_directory else folder, "training.csv"))
  stopifnot(all(training$role == "training"))
  recipes <- readRDS(file.path(folder, "recipes.rds"))
  for (key in names(recipes)) {
    fitting <- if (key == "training") training else training[inner_fold != as.integer(sub("inner_", "", key))]
    archived <- recipes[[key]]$direction_recipe$base_recipe
    fresh <- extract(fitting, panel, cfg)
    stopifnot(setequal(fresh$training_strains, archived$training_strains))
    reference <- if (charite) panel$protein[fitting$sample_id,,drop=FALSE] else panel$X[fitting$sample_id,,drop=FALSE]
    weight <- 1 / as.numeric(table(fitting$strain_id)[fitting$strain_id])
    independent <- drop(crossprod(weight, !is.na(reference))) / uniqueN(fitting$strain_id)
    observed_error <- max(abs(independent[names(fresh$observed)] - fresh$observed))
    stopifnot(is.finite(observed_error), observed_error < 1e-10, max(fresh$observed) <= 1+1e-12)
    old <- archived[[new_feature_key]]
    new <- fresh[[new_feature_key]]
    matched <- all(vapply(fields, function(field) isTRUE(all.equal(fresh[[field]], archived[[field]],
      tolerance=1e-11, check.attributes=TRUE)), logical(1)))
    records[[length(records)+1L]] <- data.table(fold, panel=key, patients=uniqueN(fitting$strain_id),
      original_features=length(old), corrected_features=length(new), removed=length(setdiff(old,new)),
      added=length(setdiff(new,old)), numerical_design_equivalent=matched,
      maximum_observation_fraction_error=observed_error)
  }
}
out <- "results/preflight/clinical_patient_weights_20261106"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
checks <- rbindlist(records)
fwrite(checks, file.path(out, paste0(id, ".csv")))
state <- list(study_id=id, panels=nrow(checks), changed_panels=sum(!checks$numerical_design_equivalent),
  largest_removed=max(checks$removed), largest_added=max(checks$added), reserved_outcomes_used=FALSE)
write_json(state, file.path(out, paste0(id, ".json")), auto_unbox=TRUE, pretty=TRUE)
print(state)

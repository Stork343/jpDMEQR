root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())
output <- file.path(root, "data", "derived", "GSE65391", "module_measurements_20261022_v1")
if (dir.exists(output)) stop("The frozen module measurement object already exists.")
obj <- readRDS(file.path(root, "data", "derived", "GSE65391", "subject_design_20261022_v1", "analysis_inputs.rds"))
rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
expression <- readRDS(obj$expression_file)[, rownames(obj$eligible_metadata), drop = FALSE]
expression <- subset_expression_to_frozen_modules_v2(expression, annotation, members)
dir.create(output, recursive = TRUE)
saveRDS(expression, file.path(output, "observed_module_probe_expression.rds"))
jsonlite::write_json(list(source_expression = obj$expression_file, n_probes = nrow(expression),
                          n_visits = ncol(expression), n_subjects = 129L,
                          scope = "all_unique_gene_probes_and_all_original_probe_members_for_five_modules",
                          source_membership_file = rule$score_membership_file,
                          created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                     file.path(output, "measurement_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

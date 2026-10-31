library(data.table)
library(yaml)
library(jsonlite)
source("research/drug_response_application_20261101/data.R")

config <- read_yaml("config/discovery_gdsc_20261101.yml")
data <- readRDS("data/derived/drug_response_20261101/discovery_cases.rds")
training <- data[DRUG_NAME == "Dabrafenib" & exploration_role == "training" & analysis_curve]
evaluation <- data[DRUG_NAME == "Dabrafenib" & exploration_role == "evaluation" & analysis_curve]
expression <- readRDS("data/derived/drug_response_20261101/expression.rds")
recipe <- drug_recipe(training, expression, config)
rows <- training[1:10]
features <- drug_features(rows, expression, recipe)
components <- drug_features(rows, expression, recipe, expand = FALSE)
native <- fbioSelectr::BasisFuncApproximation(cbind(components$clinical, components$genes),
  components$dose, degrFreedom = 5L, degree = 3L, intercept = TRUE,
  knots = attr(recipe$basis, "knots"))
low_count <- ncol(components$clinical) * 5L
gene_count <- ncol(components$genes)
permutation <- c(seq_len(low_count),
                 as.vector(t(matrix(seq_len(gene_count * 5L) + low_count, nrow = 5L))))
identity_error <- max(abs(features$X - native[, permutation, drop = FALSE]))
zero_overlap <- !length(intersect(recipe$training_patients, evaluation$patient_id))
stopifnot(identity_error < 1e-12, zero_overlap,
          identical(recipe$training_patients, sort(unique(training$patient_id))),
          all(is.finite(features$X)), range(components$dose)[1L] == 0,
          range(components$dose)[2L] == 1)
report <- list(status = "passed", native_basis_identity_error = identity_error,
                training_evaluation_overlap = 0L, genes = length(recipe$genes),
                expanded_variables = ncol(features$X), actual_observations_used = TRUE)
directory <- "results/preflight/empirical_selection_20261101/gdsc_features"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
write_json(report, file.path(directory, "checks.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(report)

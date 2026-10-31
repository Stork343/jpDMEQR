library(data.table)
library(jsonlite)
library(grf)
source("research/yeast_application_20261101/prediction.R")
source("research/yeast_application_20261101/comparators.R")
data <- readRDS("data/derived/yeast_genomic_20261101/growth_discovery.rds")
rows <- data[trait == "Magnesium_Sulfate" & exploration_role == "training"]
rows <- rows[strain_id %in% head(unique(rows$strain_id), 60L)]
panel <- readRDS("data/derived/yeast_genomic_20261101/panel.rds")
features <- yeast_features(rows, panel, yeast_recipe(rows, panel))
settings <- data.frame(regression_splitting = FALSE, mtry_multiplier = 1, min_node_size = 5L, honesty = TRUE)
fit <- yeast_forest_fit(features, settings, 20261201L)
predicted <- yeast_forest_predict(fit, features)
repeated <- yeast_forest_predict(yeast_forest_fit(features, settings, 20261201L), features)
stopifnot(max(abs(predicted - repeated)) < 1e-12,
          all(predicted[, 1L] <= predicted[, 2L]), all(predicted[, 2L] <= predicted[, 3L]))
directory <- "results/preflight/empirical_selection_20261101/yeast_forest"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
write_json(list(status = "passed", actual_observations_used = TRUE, strains = features$n,
                equal_cluster_weights = TRUE, deterministic_error = max(abs(predicted - repeated)),
                quantiles_ordered = TRUE, package_version = as.character(packageVersion("grf"))),
           file.path(directory, "checks.json"), auto_unbox = TRUE, pretty = TRUE)
cat("分位数森林检查通过\n")

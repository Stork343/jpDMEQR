root <- normalizePath(".", mustWork = TRUE)
preprocessing <- readRDS(file.path(root, "results", "raw", "application", "workflow_preflight_20261022_v3",
                                    "workflow_checkpoints", "full_sample_preprocessing.rds"))
metadata <- read.csv(file.path(root, "data", "interim", "GSE65391", "phenotype_parsed.csv"), stringsAsFactors = FALSE)
metadata <- metadata[match(preprocessing$sample_ids, metadata$sample_id), , drop = FALSE]
X <- preprocessing$clinical$X[, c("(Intercept)", "centered_time", "age"), drop = FALSE]
order <- order(as.character(metadata$subject_id), metadata$cumulative_time)
X <- X[order, , drop = FALSE]
attr(X, "assign") <- c(0L, 1L, 2L)
ids <- as.character(metadata$subject_id[order])
output <- file.path(root, "results", "raw", "application_interface_audit", "geeVerse_20261022")
if (dir.exists(output)) stop("The versioned public-interface audit already exists.")
dir.create(output, recursive = TRUE)
withCallingHandlers(
  geeVerse::qpgee(X, log(metadata$c3[order]), nobs = as.integer(table(factor(ids, levels = unique(ids)))),
                   tau = .5, corstr = "exchangeable", method = "HBIC",
                   penalty_factor = c(0, 0), control = geeVerse::qpgeeControl()),
  error = function(condition) {
    saveRDS(condition, file.path(output, "public_api_condition.rds"))
    jsonlite::write_json(list(
      status = "observed_public_api_failure", package = "geeVerse",
      version = as.character(utils::packageVersion("geeVerse")),
      message = conditionMessage(condition), n_subjects = length(unique(ids)),
      n_visits = nrow(X), fixed_dimension = ncol(X),
      source = "official_qpgee_default_with_named_penalty_factor",
      recorded_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ), file.path(output, "public_interface_record.json"), auto_unbox = TRUE, pretty = TRUE)
  }
)

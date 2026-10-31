library(fbioSelectr)
library(jsonlite)

source_dir <- "results/scratch/empirical_selection_20261101/fbioselect_source"
output_dir <- "results/preflight/empirical_selection_20261101/fbioselect_author_example"
raw_dir <- "results/raw/application/drug_response_20261101/author_example"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(raw_dir, "fit.rds")
if (file.exists(output_file)) stop("结果文件已经存在，请使用独立运行目录")
objects <- new.env()
load(file.path(source_dir, "test_gdsc_data.RData"), objects)
fd <- objects$fd
stopifnot(is.data.frame(fd), all(is.finite(as.matrix(fd[, -(1:5)]))),
          all(is.finite(fd$Response)), all(is.finite(fd$weights)))
started <- proc.time()[["elapsed"]]
fit <- fbioSelect(
  y = fd$Response, id = fd$ID, t.var = fd$DosageSQ,
  Z = model.matrix(~ fd$DRUG_NAME), X = fd[, -(1:5)],
  low.dim.cov.number = 1,
  threshold.size = floor(length(unique(fd$ID)) / log(length(unique(fd$ID)))),
  type.of.penalty = "grSCAD", degree = 3, intercept = TRUE,
  degrFreedom = 5, knots = NULL, weights = fd$weights,
  nfolds = 10, seed = 1000, ncores = 1
)
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, output_file)
stopifnot(all(is.finite(fit$pen.regression.fit$beta)))
report <- list(
  status = "author_example_completed",
  role = "implementation_feasibility_only",
  rows = nrow(fd), experimental_ids = length(unique(fd$ID)),
  genes = ncol(fd) - 5L, drugs = sort(unique(fd$DRUG_NAME)),
  observations_per_id = as.list(table(table(fd$ID))),
  runtime_sec = elapsed,
  selected_genes = as.character(fit$active.hdcov.pen.regression),
  native_cv_unit = "row",
  validation_comparison_completed = FALSE,
  manuscript_adoption = FALSE,
  source_repository = "https://github.com/koukoulEv/fbioSelect",
  package_version = as.character(packageVersion("fbioSelectr")),
  source_commit = system2("git", c("-C", source_dir, "rev-parse", "HEAD"), stdout = TRUE)
)
write_json(report, file.path(output_dir, "execution.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(output_dir, "session.txt"))
print(report[c("status", "rows", "experimental_ids", "genes", "runtime_sec")])

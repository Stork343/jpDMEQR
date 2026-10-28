library(data.table)
source <- "results/tables/submission_20261028/gaussian_diagnostics"
registry <- fread("results/raw/simulation/submission_20261027_main_joint/config_snapshot.csv")[module == "A_theory_scaling"]
summary <- merge(fread(file.path(source, "coordinate_summary.csv")),
  registry[, .(experiment_id, n_clusters, p, s)], by = "experiment_id")
figure <- fread(file.path(source, "figure_source.csv"))
accounting <- fread(file.path(source, "accounting.csv"))
population <- fread("results/scratch/submission_completion_20261027/gaussian_curvature_checks/population_check.csv")
if (nrow(accounting) != 10500L || any(!accounting$available) || nrow(summary) != 24L)
  stop("The complete Gaussian scaling calculations are required.")
number <- function(value) formatC(value, format = "f", digits = 3L)
get_value <- function(id) {
  row <- figure[experiment_id == id & metric == "Scaled influence remainder"]
  if (nrow(row) != 1L) stop("A Gaussian result value is unavailable.")
  number(row$value)
}
output <- "manuscript/generated/submission"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
values <- c(GaussianDiagnosticDatasets = as.character(nrow(accounting)),
  GaussianRemainderFiveSmall = get_value("A04"), GaussianRemainderFiveLarge = get_value("A06"),
  GaussianRemainderTenSmall = get_value("A10"), GaussianRemainderTenLarge = get_value("A12"),
  GaussianCurvatureIntegral = formatC(population$quadrature, format = "f", digits = 7L),
  GaussianCurvaturePopulation = formatC(population$population_mean, format = "f", digits = 7L),
  GaussianCurvatureMCSE = formatC(population$population_mcse, format = "f", digits = 7L),
  GaussianCurvatureClusters = as.character(population$clusters))
writeLines(paste0("\\newcommand{\\", names(values), "}{", values, "}"), file.path(output, "gaussian_values.tex"))
setorder(summary, s, p, n_clusters, coordinate_type)
table <- data.frame(Setting = summary$experiment_id, n = summary$n_clusters, p = summary$p, s = summary$s,
  Coordinates = ifelse(summary$coordinate_type == "active", "Active", "Null"),
  Replications = summary$available_replications, RMS = number(summary$remainder_rms),
  MCSE = number(summary$remainder_rms_mcse), Relative = number(summary$relative_remainder_rms))
names(table)[2:4] <- c("$n$", "$p$", "$s$")
content <- knitr::kable(table, format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
  caption = "Gaussian influence remainder by requested coordinate group", label = "gaussian-remainder",
  linesep = "", align = c("l", "r", "r", "r", "l", rep("r", 4L)))
writeLines(as.character(content), file.path(output, "gaussian_remainder.tex"))
cat("Wrote Gaussian numerical values and coordinate table\n")

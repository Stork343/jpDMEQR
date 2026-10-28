args <- base::commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
audit_dir <- file.path(root, "results", "scratch", "profile_score_audit_20261022")
research_dir <- file.path(dirname(root), "research", "profile-score-theory-20261022")
dir.create(research_dir, recursive = TRUE, showWarnings = FALSE)
raw <- utils::read.csv(file.path(audit_dir, "population_targets.csv"))
summary <- raw
summary$numerical_zero_threshold <-
  128 * .Machine$double.eps * pmax(1, abs(summary$beta_star))
summary$numerical_zero_h <- abs(summary$error_h) <= summary$numerical_zero_threshold
summary$numerical_zero_extrapolated <-
  abs(summary$error_extrapolated) <= summary$numerical_zero_threshold
summary$status_h <- ifelse(summary$numerical_zero_h, "numerical_zero", "resolved")
summary$status_extrapolated <-
  ifelse(summary$numerical_zero_extrapolated, "numerical_zero", "resolved")
summary$order_h <- NA_real_
summary$order_extrapolated <- NA_real_
for (ix in split(seq_len(nrow(summary)), interaction(summary$tau, summary$lambda))) {
  ix <- ix[order(summary$h[ix], decreasing = TRUE)]
  for (j in seq_len(length(ix) - 1L)) {
    a <- ix[j]
    b <- ix[j + 1L]
    denominator <- log(summary$h[a] / summary$h[b])
    if (!summary$numerical_zero_h[a] && !summary$numerical_zero_h[b]) {
      summary$order_h[b] <- log(abs(summary$error_h[a] / summary$error_h[b])) / denominator
    }
    if (!summary$numerical_zero_extrapolated[a] &&
        !summary$numerical_zero_extrapolated[b]) {
      summary$order_extrapolated[b] <-
        log(abs(summary$error_extrapolated[a] / summary$error_extrapolated[b])) / denominator
    }
  }
}
stopifnot(!any(is.infinite(summary$order_h)),
          !any(is.infinite(summary$order_extrapolated)),
          all(is.na(summary$order_h[summary$numerical_zero_h])),
          all(is.na(summary$order_extrapolated[summary$numerical_zero_extrapolated])))
utils::write.csv(summary, file.path(audit_dir, "population_target_summary.csv"),
                 row.names = FALSE, na = "NA")

selected <- subset(summary, tau == 0.25 & lambda == 1 &
                     h %in% c(0.05, 0.025, 0.0125, 0.00625))
selected <- selected[order(selected$h, decreasing = TRUE), ]
stopifnot(nrow(selected) == 4L)
selected$order_h[1] <- NA_real_
selected$order_extrapolated[1] <- NA_real_
format_order <- function(x) if (is.na(x)) "NA" else sprintf("%.3f", x)
lines <- c(
  "\\begin{tabular}{@{}rrrrr@{}}",
  "\\toprule",
  "$h$ & $|\\beta_h^\\star-\\beta^\\star|$ & Order & $|\\beta_{\\mathrm{pair},h}^\\star-\\beta^\\star|$ & Order \\\\",
  "\\midrule",
  vapply(seq_len(nrow(selected)), function(j) {
    sprintf("%.5f & %.9f & %s & %.9f & %s \\\\",
            selected$h[j], abs(selected$error_h[j]), format_order(selected$order_h[j]),
            abs(selected$error_extrapolated[j]), format_order(selected$order_extrapolated[j]))
  }, character(1)),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(lines, file.path(research_dir, "population_target_table.tex"))
utils::write.csv(selected, file.path(research_dir, "population_target_table_source.csv"),
                 row.names = FALSE, na = "NA")
print(summary[summary$tau == 0.5,
              c("tau", "lambda", "h", "status_h", "status_extrapolated",
                "order_h", "order_extrapolated")])
cat(paste(lines, collapse = "\n"), "\n")

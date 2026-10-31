library(data.table)
source("research/score_validity_20261026/load_joint.R")
input <- readRDS("results/raw/application/soynam_yield_20261031/input.rds")
ids <- input$line_folds[, head(ID, 10L), by = Family]$V1
rows <- which(input$observations$ID %in% ids)
observed <- input$observations[rows]
fixed_columns <- c("ACC10_2013", "ACC10_2014", paste0("Block", 1:8), "check_surface")
fixed <- input$fixed_design[rows, fixed_columns, drop = FALSE]
G <- input$genotypes[match(observed$ID, rownames(input$genotypes)), 1:8, drop = FALSE]
X <- cbind(fixed, G * as.numeric(observed$Year == 2013), G * as.numeric(observed$Year == 2014))
colnames(X) <- make.unique(colnames(X))
data <- sb_data(observed$YLD / 10, X, matrix(1, nrow(X), 1L), observed$cluster_id)
penalty <- c(rep(0, ncol(fixed)), rep(1, ncol(X) - ncol(fixed)))
checks <- list()
for (tau in c(0.25, 0.5, 0.75)) {
  started <- proc.time()[[3L]]
  prepared <- pc_prepare_initial(data$y, data$X, data$Z, data$cluster_id, tau, matrix(1), penalty)
  reference <- pc_orthogonal_path(prepared, tau, 0:8, gap_tolerance = 1e-9, max_iterations = 100000L)
  reference_seconds <- proc.time()[[3L]] - started
  started <- proc.time()[[3L]]
  prepared <- pc_prepare_initial(data$y, data$X, data$Z, data$cluster_id, tau, matrix(1), penalty,
    fit_always_included = FALSE)
  current <- pc_orthogonal_path(prepared, tau, 0:8, gap_tolerance = 1e-9, max_iterations = 100000L)
  current_seconds <- proc.time()[[3L]] - started
  error <- max(abs(data$X %*% (reference$coefficients - current$coefficients)))
  if (error > 1e-5 || !identical(reference$selected_path, current$selected_path))
    stop("The actual-data orthogonal initialisation check failed.")
  checks[[length(checks) + 1L]] <- data.table(tau = tau, actual_lines = data$n,
    columns = ncol(X), prediction_error = error, same_selected_path = TRUE,
    reference_seconds = reference_seconds, direct_basis_seconds = current_seconds,
    reference_gap = reference$fitted$duality_gap, direct_basis_gap = current$fitted$duality_gap)
}
fwrite(rbindlist(checks), "results/preflight/soynam_20261031/initial_basis_checks.csv")
print(rbindlist(checks))

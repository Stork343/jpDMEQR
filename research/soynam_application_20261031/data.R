soy_design <- function(input, cohort = "full", prediction = FALSE, canopy = TRUE) {
  observed <- input$observations
  fixed <- input$fixed_design
  if (cohort == "maturity_records") {
    omitted <- unique(observed$ID[observed$DTM > 180])
    index <- which(!observed$ID %in% omitted)
    observed <- observed[index]
    fixed <- fixed[index, , drop = FALSE]
  }
  if (cohort == "spatial_polynomial") {
    fixed <- cbind(fixed, input$spatial_polynomial)
  }
  if (prediction) {
    fixed <- fixed[, !grepl("^DTM_", colnames(fixed)), drop = FALSE]
    if (!canopy) fixed <- fixed[, !grepl("^ACC10_", colnames(fixed)), drop = FALSE]
  }
  row_index <- match(observed$ID, rownames(input$genotypes))
  if (anyNA(row_index)) stop("A field observation lacks its genomic row.")
  genetic <- lapply(input$years, function(year) {
    matrix <- input$genotypes[row_index, , drop = FALSE] * as.numeric(observed$Year == year)
    colnames(matrix) <- paste0("SNP_", year, "_", colnames(input$genotypes))
    matrix
  })
  X <- do.call(cbind, c(list(fixed), genetic))
  if (any(!is.finite(X))) stop("The application design contains nonfinite values.")
  list(y = observed$YLD / input$response_divisor, X = X,
    Z = matrix(1, nrow(observed), 1L, dimnames = list(NULL, "line_intercept")),
    cluster_id = observed$cluster_id, penalty_factor = c(rep(0, ncol(fixed)), rep(1, ncol(X) - ncol(fixed))),
    observations = observed, n = uniqueN(observed$ID), p = ncol(X), fixed_dimension = ncol(fixed),
    cohort = cohort, prediction = prediction, canopy = canopy)
}

soy_subset <- function(data, index) {
  selected <- sb_data(data$y[index], data$X[index, , drop = FALSE], data$Z[index, , drop = FALSE],
    data$cluster_id[index])
  selected$penalty_factor <- data$penalty_factor
  selected$original_rows <- index
  selected
}

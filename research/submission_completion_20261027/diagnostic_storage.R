submission_compact_diagnostics <- function(result) {
  compact <- function(value) {
    if (!is.list(value) || is.data.frame(value)) return(value)
    if (is.data.frame(value$history) && nrow(value$history) > 2L &&
        all(c("iteration", "objective", "beta_change", "kkt_residual", "accepted") %in% names(value$history))) {
      history <- value$history
      value$optimizer_history_summary <- list(iteration_records = nrow(history),
        initial = history[1L, , drop = FALSE], final = history[nrow(history), , drop = FALSE],
        accepted_steps = sum(history$accepted), minimum_objective = min(history$objective),
        maximum_kkt_residual = max(history$kkt_residual), maximum_beta_change = max(history$beta_change))
      value$history <- NULL
      if (is.numeric(value$score)) {
        value$optimizer_gradient_summary <- list(coordinates = length(value$score),
          maximum_absolute_score = max(abs(value$score)), l2_norm = sqrt(sum(value$score^2)))
        value$score <- NULL
      }
    }
    for (i in seq_along(value)) value[i] <- list(compact(value[[i]]))
    value
  }
  if (!is.null(result$method_diagnostics[["SQR-DEBIASED-IID"]])) {
    result$method_diagnostics[["SQR-DEBIASED-IID"]] <- compact(result$method_diagnostics[["SQR-DEBIASED-IID"]])
    result$diagnostic_storage <- "SQR_unrounded_end_conditions_and_candidate_coefficients"
  }
  result
}

ep_report_prediction <- function(config, queue, panel, feature_function, output,
                                 target_method = "JOINT", stage = "discovery",
                                 feature_functions = NULL, prediction_functions = NULL) {
  stopifnot(nrow(queue) > 0L, all(queue$status == "completed"), !anyDuplicated(queue$directory))
  records <- checks <- selections <- settings <- list()
  for (j in seq_len(nrow(queue))) {
    directory <- queue$directory[j]
    assay_name <- queue$trait[j]
    evaluation <- data.table::fread(file.path(config$input_directory, assay_name, "discovery.csv"))[exploration_role == "evaluation"]
    request <- jsonlite::read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
    status <- jsonlite::read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
    stopifnot(status$status == "completed", !request$formal_validation_used, all(evaluation$role == "discovery"))
    fit <- readRDS(file.path(directory, "fit.rds"))
    stopifnot(!length(intersect(fit$recipe$training_strains, evaluation$strain_id)))
    current_feature <- if (is.null(feature_functions)) feature_function else feature_functions[[request$method]]
    current_prediction <- if (is.null(prediction_functions)) ep_predict else prediction_functions[[request$method]]
    stopifnot(is.function(current_feature), is.function(current_prediction))
    features <- current_feature(evaluation, panel, fit$recipe)
    base <- drop(current_prediction(features, fit$beta))
    saved <- data.table::fread(file.path(directory, "predictions.csv"))
    error <- 0
    for (k in seq_along(fit$quantiles)) {
      value <- fit$quantiles[k]
      current <- saved[tau == value]
      stopifnot(identical(current$sample_id, evaluation$sample_id), identical(current$strain_id, evaluation$strain_id))
      prediction <- base + fit$calibration[[k]][evaluation$condition]
      error <- max(error, abs(current$y - evaluation$y), abs(current$prediction - prediction),
                   abs(current$loss - ep_loss(evaluation$y, prediction, value)))
    }
    stopifnot(is.finite(error), error < 1e-9)
    checks[[j]] <- data.table::data.table(assay = assay_name, method = request$method, run = basename(directory), maximum_error = error)
    losses <- saved[, .(loss = mean(loss)), by = .(strain_id, tau)]
    losses[, `:=`(assay = assay_name, method = request$method)]
    records[[j]] <- losses
    if (request$method == target_method) {
      active <- which(abs(fit$beta) > 1e-10)
      selections[[length(selections) + 1L]] <- data.table::data.table(assay = assay_name, tau = request$tau,
        variable = names(fit$beta)[active], coefficient_scaled = fit$beta[active],
        coefficient_response_scale = fit$beta[active] * fit$recipe$y_scale,
        role = ifelse(features$penalty_factor[active] == 0, "adjustment", "selected"))
    }
    settings[[j]] <- cbind(data.table::data.table(assay = assay_name, method = request$method,
      tau = request$tau, runtime_sec = status$runtime_sec), data.table::as.data.table(as.data.frame(fit$grid)[fit$choice, , drop = FALSE]))
    rm(features, fit)
    gc()
  }
  losses <- data.table::rbindlist(records)
  stopifnot(!anyDuplicated(losses[, .(assay, strain_id, tau, method)]))
  set.seed(20261105L)
  comparisons <- list()
  for (assay_name in unique(queue$trait)) {
    for (value in sort(unique(losses[assay == assay_name & method == target_method]$tau))) {
      wide <- data.table::dcast(losses[assay == assay_name & tau == value], strain_id ~ method, value.var = "loss")
      stopifnot(!anyNA(wide), nrow(wide) > 1L)
      weights <- rmultinom(10000L, nrow(wide), rep(1 / nrow(wide), nrow(wide))) / nrow(wide)
      for (method_name in setdiff(names(wide), c("strain_id", target_method))) {
        reference <- wide[[method_name]]
        difference <- reference - wide[[target_method]]
        bootstrap <- drop(crossprod(weights, difference)) / drop(crossprod(weights, reference))
        comparisons[[length(comparisons) + 1L]] <- data.table::data.table(assay = assay_name, tau = value,
          target_method = target_method, comparator = method_name, strains = nrow(wide),
          joint_loss = mean(wide[[target_method]]), comparator_loss = mean(reference),
          relative_gain = mean(difference) / mean(reference), lower = quantile(bootstrap, 0.025),
          upper = quantile(bootstrap, 0.975), phase = stage, manuscript_adoption = FALSE)
      }
    }
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  comparison <- data.table::rbindlist(comparisons)
  data.table::fwrite(comparison, file.path(output, "comparisons.csv"))
  data.table::fwrite(data.table::rbindlist(checks), file.path(output, "artifact_checks.csv"))
  data.table::fwrite(data.table::rbindlist(selections), file.path(output, "selected_coefficients.csv"))
  data.table::fwrite(data.table::rbindlist(settings, fill = TRUE), file.path(output, "settings.csv"))
  summary <- losses[, .(strains = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = .(assay, method, tau)]
  data.table::fwrite(summary, file.path(output, "summary.csv"))
  jsonlite::write_json(list(stage = stage, target_method = target_method, completed_runs = nrow(queue), formal_validation_used = FALSE,
                           manuscript_adoption = FALSE), file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
  comparison
}

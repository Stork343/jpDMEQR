drug_calibration <- function(rows, oof, tau) {
  weight <- 1 / as.numeric(table(rows$patient_id)[rows$patient_id])
  residual <- rows$normalized_intensity - oof
  levels <- sort(unique(rows$CONC))
  design <- vapply(levels, function(value) as.numeric(rows$CONC == value), numeric(nrow(rows)))
  fit <- quantreg::rq.fit(design * weight, residual * weight, tau = tau, method = "br")
  stopifnot(all(is.finite(fit$coefficients)))
  list(concentration = levels, offset = as.numeric(fit$coefficients), tau = tau)
}

drug_apply_calibration <- function(calibration, concentration) {
  approx(log(calibration$concentration), calibration$offset,
         xout = log(concentration), rule = 2L, ties = "ordered")$y
}

drug_joint_prediction <- function(rows, evaluation, expression, config, tau, progress, directory) {
  folds <- sort(unique(rows$inner_fold))
  maximum <- min(vapply(folds, function(fold) {
    current <- rows[inner_fold != fold]
    n <- uniqueN(current$patient_id)
    d0 <- qr(drug_clinical(current))$rank * config$dose_basis$degrees_of_freedom
    max(0L, floor(n / log(n)) - d0)
  }, numeric(1)))
  sizes <- 0:as.integer(maximum)
  lambdas <- unlist(config$joint_prediction$Lambda)
  loss <- matrix(0, length(lambdas), length(sizes))
  oof_paths <- lapply(lambdas, function(value) matrix(NA_real_, nrow(rows), length(sizes)))
  certificates <- list()
  for (fold in folds) {
    training_rows <- rows[inner_fold != fold]
    validation_rows <- rows[inner_fold == fold]
    recipe <- drug_recipe(training_rows, expression, config)
    training <- drug_features(training_rows, expression, recipe)
    validation <- drug_features(validation_rows, expression, recipe)
    for (l in seq_along(lambdas)) {
      progress("joint_inner_fit", list(fold = fold, Lambda = lambdas[l],
               training_patients = training$n, variables = ncol(training$X), additions = maximum))
      prepared <- pc_prepare_initial(training$y, training$X, training$Z,
        training$cluster_id, tau, matrix(lambdas[l], 1L, 1L), training$penalty_factor,
        fit_always_included = FALSE)
      path <- pc_orthogonal_path(prepared, tau, sizes,
        gap_tolerance = config$joint_prediction$gap_tolerance,
        max_iterations = as.integer(config$joint_prediction$maximum_iterations), solver = jp_joint_qp)
      predicted <- validation$X %*% path$coefficients
      loss[l, ] <- loss[l, ] + colSums(drug_loss(validation$y, predicted, tau) * validation$weight)
      oof_paths[[l]][which(rows$inner_fold == fold), ] <- predicted
      certificates[[length(certificates) + 1L]] <- list(fold = fold, Lambda = lambdas[l],
        diagnostics = path$diagnostics, final_gap = path$fitted$duality_gap)
      saveRDS(list(recipe = recipe, path = path, prediction = predicted,
                   validation_patients = validation_rows$patient_id, tau = tau, Lambda = lambdas[l]),
              file.path(directory, sprintf("inner_%d_lambda_%g.rds", fold, lambdas[l])))
      rm(prepared, path, predicted)
      gc()
    }
    rm(training, validation)
    gc()
  }
  choice <- arrayInd(which.min(loss), dim(loss))
  selected_lambda <- lambdas[choice[1L]]
  selected_size <- sizes[choice[2L]]
  oof <- oof_paths[[choice[1L]]][, choice[2L]]
  stopifnot(all(is.finite(oof)))
  calibration <- drug_calibration(rows, oof, tau)
  recipe <- drug_recipe(rows, expression, config)
  training <- drug_features(rows, expression, recipe)
  progress("joint_training_fit", list(Lambda = selected_lambda, additions = selected_size,
                                      variables = ncol(training$X), patients = training$n))
  prepared <- pc_prepare_initial(training$y, training$X, training$Z, training$cluster_id,
    tau, matrix(selected_lambda, 1L, 1L), training$penalty_factor, fit_always_included = FALSE)
  fitted <- pc_orthogonal_path(prepared, tau, 0:selected_size,
    gap_tolerance = config$joint_prediction$gap_tolerance,
    max_iterations = as.integer(config$joint_prediction$maximum_iterations), solver = jp_joint_qp)
  beta <- fitted$coefficients[, selected_size + 1L]
  names(beta) <- colnames(training$X)
  final_certificate <- list(gap = fitted$fitted$duality_gap, diagnostics = fitted$diagnostics)
  rm(training, prepared, fitted)
  gc()
  test <- drug_features(evaluation, expression, recipe)
  prediction <- drop(test$X %*% beta)
  list(beta = beta, recipe = recipe, oof = oof, calibration = calibration,
       prediction = prediction + drug_apply_calibration(calibration, evaluation$CONC),
       uncalibrated_prediction = prediction, loss = loss, sizes = sizes, lambdas = lambdas,
       selected_lambda = selected_lambda, selected_size = selected_size,
       certificates = certificates, final_certificate = final_certificate,
       retained_genes = length(recipe$genes), design_variables = length(beta))
}

drug_fbio_prediction <- function(rows, evaluation, expression, config, progress, directory) {
  folds <- sort(unique(rows$inner_fold))
  oof <- rep(NA_real_, nrow(rows))
  diagnostics <- list()
  for (fold in folds) {
    training_rows <- rows[inner_fold != fold]
    validation_rows <- rows[inner_fold == fold]
    recipe <- drug_recipe(training_rows, expression, config)
    training <- drug_features(training_rows, expression, recipe, expand = FALSE)
    validation <- drug_features(validation_rows, expression, recipe, expand = FALSE)
    progress("fbio_inner_fit", list(fold = fold, genes = ncol(training$genes), patients = training$n))
    fit <- drug_fbio_fit(training$y, training$genes, training$dose,
      training$clinical, training$cluster_id, seed = 20261110L + fold)
    oof[which(rows$inner_fold == fold)] <- drug_fbio_predict(fit, validation$genes,
                                                          validation$dose, validation$clinical)
    diagnostics[[length(diagnostics) + 1L]] <- list(fold = fold, lambda = fit$lambda,
                                                   selected = fit$selected, runtime_sec = fit$runtime_sec)
    saveRDS(list(recipe = recipe, fit = fit, validation_patients = validation_rows$patient_id,
                 prediction = oof[which(rows$inner_fold == fold)]),
            file.path(directory, paste0("inner_", fold, ".rds")))
    rm(training, validation, fit)
    gc()
  }
  stopifnot(all(is.finite(oof)))
  recipe <- drug_recipe(rows, expression, config)
  training <- drug_features(rows, expression, recipe, expand = FALSE)
  progress("fbio_training_fit", list(genes = ncol(training$genes), patients = training$n))
  fit <- drug_fbio_fit(training$y, training$genes, training$dose, training$clinical,
                       training$cluster_id, seed = 20261120L)
  rm(training)
  gc()
  test <- drug_features(evaluation, expression, recipe, expand = FALSE)
  mean_prediction <- drug_fbio_predict(fit, test$genes, test$dose, test$clinical)
  calibrations <- lapply(config$quantiles, function(tau) drug_calibration(rows, oof, tau))
  predictions <- vapply(calibrations, function(calibration)
    mean_prediction + drug_apply_calibration(calibration, evaluation$CONC), numeric(nrow(evaluation)))
  list(model = fit, recipe = recipe, mean_oof = oof, mean_prediction = mean_prediction,
       calibration = calibrations, prediction = predictions, diagnostics = diagnostics,
       retained_genes = length(recipe$genes))
}

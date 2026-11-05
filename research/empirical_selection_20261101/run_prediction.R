ep_run_application_prediction <- function(training, evaluation, panel, recipes, config, method, tau,
                                           directory, base_feature_function) {
  stopifnot(method %in% c("JOINT_OFFSET", "POOLED_OFFSET", "JOINT_DIRECT", "ENET", "PRESS_RIDGE"),
    !length(intersect(training$strain_id, evaluation$strain_id)),
    all(training$role == "discovery"), all(evaluation$role == "discovery"))
  if (dir.exists(directory)) stop("运行目录已经存在")
  dir.create(directory, recursive = TRUE)
  progress <- function(stage, details = list()) {
    record <- c(list(stage = stage, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
    jsonlite::write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
    cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
    flush.console()
  }
  jsonlite::write_json(list(method = method, tau = tau, config = config, phase = config$phase,
    code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    formal_validation_used = FALSE, manuscript_adoption = FALSE),
    file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
  recipe_function <- function(current, panel, config) {
    keys <- names(recipes)[vapply(recipes,
      function(recipe) identical(recipe$training_strains, sort(unique(current$strain_id))), logical(1))]
    stopifnot(length(keys) == 1L)
    recipe <- recipes[[keys]]
    if (method %in% c("ENET", "JOINT_DIRECT")) recipe$direction_recipe$base_recipe else recipe
  }
  offset_features <- function(current, panel, recipe) ep_offset_features(current, panel, recipe, base_feature_function)
  started <- proc.time()[["elapsed"]]
  if (method == "PRESS_RIDGE") {
    direction <- recipes$training$direction_recipe
    background <- direction$background
    quantiles <- unlist(config$quantiles)
    oof <- background$leave_out_prediction[training$sample_id] * direction$y_scale + direction$y_center
    calibration <- lapply(quantiles, function(value) ep_calibrate(training, oof, value, "strain_id", "condition"))
    test <- base_feature_function(evaluation, panel, direction$base_recipe)
    base <- drop(ep_predict(test, background$beta))
    prediction <- vapply(calibration, function(value) base + value[evaluation$condition], numeric(nrow(evaluation)))
    fit <- list(beta = background$beta, recipe = direction$base_recipe, prediction = prediction,
      quantiles = quantiles, calibration = calibration, oof = oof, choice = background$choice,
      grid = data.table::data.table(lambda = background$lambda), loss = background$loss)
  } else {
    offset <- method %in% c("JOINT_OFFSET", "POOLED_OFFSET")
    fit <- ep_fit_prediction(training, evaluation, panel,
      if (method == "ENET") "ENET" else if (method == "POOLED_OFFSET") "POOLED" else "JOINT",
      tau, config, progress, directory, recipe_function, if (offset) offset_features else base_feature_function,
      "strain_id", "condition", if (offset) ep_predict_offset else ep_predict,
      if (method == "ENET") "oof" else "training")
  }
  elapsed <- proc.time()[["elapsed"]] - started
  saveRDS(fit, file.path(directory, "fit.rds"))
  predictions <- data.table::rbindlist(lapply(seq_along(fit$quantiles), function(k)
    data.table::data.table(sample_id = evaluation$sample_id, strain_id = evaluation$strain_id,
      condition = evaluation$condition, y = evaluation$y, tau = fit$quantiles[k],
      prediction = fit$prediction[, k], loss = ep_loss(evaluation$y, fit$prediction[, k], fit$quantiles[k]))))
  data.table::fwrite(predictions, file.path(directory, "predictions.csv"))
  summary <- predictions[, .(loss = mean(loss)), by = .(strain_id, tau)][,
    .(clusters = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = tau]
  summary[, `:=`(method = method, runtime_sec = elapsed)]
  data.table::fwrite(summary, file.path(directory, "summary.csv"))
  capture.output(sessionInfo(), file = file.path(directory, "session.txt"))
  capture.output(warnings(), file = file.path(directory, "warnings.txt"))
  jsonlite::write_json(list(status = "completed", runtime_sec = elapsed, formal_validation_used = FALSE,
    manuscript_adoption = FALSE), file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
  progress("completed", list(runtime_sec = elapsed))
  invisible(fit)
}

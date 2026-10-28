source("research/submission_completion_20261027/application.R")
args <- commandArgs(trailingOnly = TRUE)
path <- args[1L]
input <- readRDS(file.path(path, "input.rds"))
globalCallingHandlers(error = function(condition) saveRDS(condition, file.path(path, "failure.rds")))
fit <- submission_qgee(input$metadata, input$design, input$training, input$tau, input$seed)
saveRDS(fit, file.path(path, "result.rds"))

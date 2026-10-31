library(data.table)
library(yaml)
config <- read_yaml("config/application_soynam_20261031.yml")
input <- readRDS(file.path(config$results_directory, "input.rds"))
observed <- copy(input$observations)
destination <- "results/tables/soynam_20261031"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
summary <- rbindlist(lapply(c("YLD", "ACC", "DTM"), function(variable) {
  observed[, .(variable = variable, records = .N, mean = mean(get(variable)),
    sd = sd(get(variable)), minimum = min(get(variable)),
    q25 = quantile(get(variable), .25), median = median(get(variable)),
    q75 = quantile(get(variable), .75), maximum = max(get(variable))), by = Year]
}))
fwrite(summary, file.path(destination, "descriptive.csv"))
folds <- input$line_folds[, .(lines = .N), by = .(Family, inference_fold, prediction_fold)]
fwrite(folds, file.path(destination, "family_fold_counts.csv"))
fwrite(input$line_folds[, .(lines = .N), by = Family], file.path(destination, "family_counts.csv"))
print(summary)

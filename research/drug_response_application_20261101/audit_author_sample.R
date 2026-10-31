library(data.table)
library(jsonlite)
objects <- new.env()
load("results/scratch/empirical_selection_20261101/fbioselect_source/test_gdsc_data.RData", objects)
for (name in ls(objects)) {
  current <- objects[[name]]
  cat("Object:", name, "class:", class(current), "dimensions:", dim(current), "length:", length(current), "\n")
  cat("Names:\n")
  print(head(names(current), 30L))
  if (is.list(current) && !is.data.frame(current)) {
    for (element in head(names(current), 20L)) {
      value <- current[[element]]
      cat(element, "class:", class(value), "dimensions:", dim(value), "length:", length(value), "\n")
      if (is.matrix(value) || is.data.frame(value)) print(head(colnames(value), 12L))
    }
  }
}

library(readxl)
path <- "data/raw/candidate_datasets_20261030/Database_MIR_A2_milk.xlsx"
print(excel_sheets(path))
data <- as.data.frame(read_excel(path, .name_repair = "minimal"))
cat("Dimensions", nrow(data), ncol(data), "\n")
print(names(data)[1:20])
print(tail(names(data), 8L))
print(head(data[, 1:13], 8L))
print(lapply(data[, 1:13], class))

library(data.table)
library(readxl)
data <- as.data.table(read_excel(
  "data/raw/candidate_datasets_20261030/Database_MIR_A2_milk.xlsx",
  .name_repair = "minimal"))
checks <- data[, .(observations = .N, farms = uniqueN(Farm),
  controls = uniqueN(Milk_Control), dates = uniqueN(date),
  beta_values = uniqueN(Betacasein), kappa_values = uniqueN(Kappacasein),
  lacto_values = uniqueN(Betalactoglobulin)), by = ID_Cow]
ambiguous <- checks[observations > 2L | farms > 1L | controls != observations |
  dates != observations | beta_values > 1L | kappa_values > 1L | lacto_values > 1L]
print(ambiguous)
print(data[ID_Cow %in% ambiguous$ID_Cow, 1:13][order(ID_Cow, date)])
print(checks[, .(cows = .N), by = observations][order(observations)])
print(data[, .(observations = .N), by = .(Farm, Milk_Control, date)][order(Farm, Milk_Control, date)])

root <- "."
source("scripts/00_source_v2.R")
source_v2_module(root, "profile_v2")
source_v2_module(root, "simulation_v2")
source("scripts/simulation/_run_registry_helpers.R")
for (phase in c("main", "pilot", "preflight")) {
  path <- file.path("config", paste0("simulation_pair_", phase, "_20261022.yml"))
  current <- read_simulation_registry_v2(path)
  frozen <- read.csv(file.path("config", paste0("simulation_", phase, ".csv")),
                     stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(identical(current[names(frozen)], frozen),
            identical(attr(current, "registry_phase"), phase),
            identical(registry_type_for_path_v2(path), phase),
            all(current$method_version == "profile_pair_20261022"),
            all(current$mean_score_multiplier == 1))
  check <- validate_registry_contract_v2(root, path, phase)
  if (!check$pass) stop(paste(check$problems, collapse = "; "))
  cat(phase, nrow(current), "frozen rows preserved\n")
}

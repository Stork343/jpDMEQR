profile_pair_population_source_paths_v2 <- function(engine = "reference") {
  if (!engine %in% c("reference", "streamed", "conditional_gaussian")) {
    stop("Unknown population source dependency engine.")
  }
  c("R/v2_profile_v2_part01.R", "R/v2_profile_v2_part12.R",
    if (identical(engine, "reference")) c("R/v2_profile_v2_part02.R", "R/v2_profile_v2_part10.R",
                                          "R/v2_profile_v2_part13.R"),
    paste0("R/v2_simulation_v2_part", sprintf("%02d", 1:6), ".R"),
    "src/cluster_profiles.cpp", "src/population_products.cpp",
    "src/Makevars", "src/Makevars.win")
}

profile_pair_population_commit_matches_v2 <- function(asset_commit, expected_commit,
                                                       engine = "reference",
                                                       root = getOption("jpDMEQR.source_root", ".")) {
  if (identical(asset_commit, expected_commit)) return(TRUE)
  if (!all(grepl("^[0-9a-fA-F]{40}$", c(asset_commit, expected_commit)))) return(FALSE)
  ancestor <- system2("git", c("-C", shQuote(root), "merge-base", "--is-ancestor",
                               asset_commit, expected_commit), stdout = FALSE, stderr = FALSE)
  if (!identical(ancestor, 0L)) return(FALSE)
  unchanged <- system2("git", c("-C", shQuote(root), "diff", "--quiet",
                                asset_commit, expected_commit, "--",
                                shQuote(profile_pair_population_source_paths_v2(engine))),
                       stdout = FALSE, stderr = FALSE)
  identical(unchanged, 0L)
}

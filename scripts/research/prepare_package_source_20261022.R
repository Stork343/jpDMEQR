prepare_package_source_20261022 <- function(root, destination) {
  root <- normalizePath(root, mustWork = TRUE)
  if (dir.exists(destination)) stop("The package staging directory already exists.")
  dir.create(destination, recursive = TRUE)
  ignored <- getFromNamespace("inRbuildignore", "tools")
  vcs <- getFromNamespace(".vc_dir_names", "tools")
  hidden <- getFromNamespace(".hidden_file_exclusions", "tools")
  inventory <- character()
  visit <- function(relative = "") {
    names <- list.files(file.path(root, relative), all.files = TRUE, no.. = TRUE)
    paths <- if (nzchar(relative)) file.path(relative, names) else names
    is_directory <- dir.exists(file.path(root, paths))
    excluded <- ignored(paths, root) | names %in% hidden | startsWith(names, "._") |
      names %in% c("Read-and-delete-me", "GNUMakefile") |
      (is_directory & (names %in% c("check", "chm", vcs) |
                         grepl("([Oo]ld|\\.Rcheck)$", names))) |
      (is_directory & grepl("^src.*/[.]deps$", paths)) |
      paths %in% c("src/jpDMEQR_res.rc", "inst/doc/.Rinstignore", "inst/doc/.build.timestamp",
                    "vignettes/.Rinstignore") | grepl("^.Rbuildindex[.]", paths) |
      grepl("^jpDMEQR_[0-9.-]+\\.(tar\\.gz|tar|tar\\.bz2|tar\\.xz|tgz|zip)$", paths)
    for (i in which(!excluded)) {
      if (is_directory[[i]]) visit(paths[[i]]) else inventory <<- c(inventory, paths[[i]])
    }
  }
  visit()
  if (!all(c("DESCRIPTION", "NAMESPACE") %in% inventory)) stop("Package metadata is missing.")
  required_tests <- list.files(file.path(root, "tests", "testthat"), pattern = "^test-.*\\.R$",
                               full.names = FALSE)
  included_tests <- basename(inventory[startsWith(inventory, "tests/testthat/test-")])
  if (!setequal(required_tests, included_tests)) stop("Package staging must retain every test file.")
  destinations <- file.path(destination, inventory)
  invisible(lapply(unique(dirname(destinations)), dir.create, recursive = TRUE, showWarnings = FALSE))
  copied <- file.copy(file.path(root, inventory), destinations, copy.mode = TRUE, copy.date = TRUE)
  if (!all(copied)) stop("Package source copying failed: ", paste(inventory[!copied], collapse = ", "))
  list(path = normalizePath(destination), files = inventory,
       test_file_count = length(required_tests), rules = "R_tools_inRbuildignore_and_build_exclusions")
}

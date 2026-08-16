# Loader for modular simulation-registry helpers.
helper_files <- sort(list.files(file.path(root, "scripts", "simulation"),
                                pattern = "^_run_registry_helpers_part[0-9]+\\.R$",
                                full.names = TRUE))
if (!length(helper_files)) stop("Simulation helper parts are missing.")
invisible(lapply(helper_files, sys.source, envir = environment()))

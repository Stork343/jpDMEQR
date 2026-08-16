# Source modular v2 implementation files when running repository scripts directly.
source_v2_module <- function(root, module, envir = parent.frame()) {
  pattern <- paste0("^v2_", module, "_part[0-9]+\\.R$")
  files <- sort(list.files(file.path(root, "R"), pattern = pattern, full.names = TRUE))
  if (length(files) == 0L) stop("No modular v2 files found for module: ", module)
  invisible(lapply(files, sys.source, envir = envir))
}

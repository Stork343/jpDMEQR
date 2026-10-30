version <- "1.15.1"
directory <- normalizePath("results/scratch/dependencies", mustWork = TRUE)
prefix <- file.path(directory, "highs-native")
header <- file.path(prefix, "include/highs/Highs.h")
libraries <- Sys.glob(file.path(prefix, "lib", "libhighs.*"))
if (!file.exists(header) || !length(libraries)) {
  cmake <- Sys.which("cmake")
  if (!nzchar(cmake)) stop("CMake is required to build the compiled CLIME dependency.")
  archive <- file.path(directory, paste0("HiGHS-", version, ".tar.gz"))
  url <- paste0("https://codeload.github.com/ERGO-Code/HiGHS/tar.gz/refs/tags/v", version)
  if (!file.exists(archive)) {
    options(timeout = max(600, getOption("timeout")))
    download.file(url, archive, mode = "wb")
  }
  source_parent <- file.path(directory, "highs-source")
  source <- file.path(source_parent, paste0("HiGHS-", version))
  if (!dir.exists(source)) {
    dir.create(source_parent, recursive = TRUE, showWarnings = FALSE)
    untar(archive, exdir = source_parent)
  }
  if (!file.exists(file.path(source, "CMakeLists.txt"))) stop("The HiGHS source archive is incomplete.")
  build <- file.path(directory, "highs-native-build-official")
  flags <- c("-S", source, "-B", build, "-DCMAKE_BUILD_TYPE=Release",
    paste0("-DCMAKE_INSTALL_PREFIX=", prefix), "-DBUILD_SHARED_LIBS=ON", "-DBUILD_TESTING=OFF",
    "-DBUILD_EXAMPLES=OFF", "-DBUILD_CXX_EXE=OFF", "-DBUILD_SHARED_EXTRAS_LIB=OFF", "-DZLIB=OFF", "-DHIPO=ON")
  if (system2(cmake, shQuote(flags)) != 0L) stop("HiGHS configuration failed.")
  if (system2(cmake, shQuote(c("--build", build, "--parallel", "2"))) != 0L) stop("HiGHS compilation failed.")
  if (system2(cmake, shQuote(c("--install", build))) != 0L) stop("HiGHS installation failed.")
}
if (!file.exists(header) || !length(Sys.glob(file.path(prefix, "lib", "libhighs.*"))))
  stop("The compiled HiGHS headers and shared library are unavailable.")
pkg_config <- Sys.which("pkg-config")
if (!nzchar(pkg_config)) stop("pkg-config is required to verify the compiled HiGHS version.")
installed_version <- system2(pkg_config, c("--modversion", "highs"), stdout = TRUE, stderr = TRUE,
  env = paste0("PKG_CONFIG_PATH=", shQuote(file.path(prefix, "lib/pkgconfig"))))
if (!is.null(attr(installed_version, "status")) || !identical(unname(installed_version), version))
  stop("The compiled HiGHS version differs from its declared source release.")
jsonlite::write_json(list(version = version, prefix = prefix, source = "ERGO-Code/HiGHS",
  source_tag = paste0("v", version), interface = "C++", build = "Release", worker_threads = 1L),
  file.path(prefix, "jpdmeqr_dependency.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Compiled CLIME dependency available at", prefix, "\n")

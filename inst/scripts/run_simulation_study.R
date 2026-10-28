#!/usr/bin/env Rscript
# v2 管线的模拟复现入口（SJS 重建版）。冻结管线的注册表、运行脚本与治理文件
# 存放在仓库检出里，不随二进制包安装，所以本脚本先定位仓库根目录，再把命令行
# 参数原样转发给快速复现入口 scripts/simulation/09_quick_repro.R。
# 完整协议见 inst/REPRODUCIBILITY.md。
#
# 用法示例（在仓库检出内任意目录执行）：
#   Rscript inst/scripts/run_simulation_study.R
#   Rscript inst/scripts/run_simulation_study.R --run-id=qr_demo --jobs=2
#   Rscript inst/scripts/run_simulation_study.R --experiment-id=QR01 --max-reps=1

find_repository_root <- function() {
  starts <- getwd()
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]),
                                        mustWork = FALSE))
    starts <- c(script_dir, starts)
  }
  for (start in starts) {
    dir <- start
    for (depth in 1:8) {
      if (file.exists(file.path(dir, "scripts", "simulation", "09_quick_repro.R"))) {
        return(normalizePath(dir, mustWork = TRUE))
      }
      parent <- dirname(dir)
      if (identical(parent, dir)) break
      dir <- parent
    }
  }
  NULL
}

root <- find_repository_root()
if (is.null(root)) {
  stop(
    "Repository checkout not found. The v2 simulation pipeline runs from a ",
    "checkout of this repository; it is not part of the installed binary ",
    "package. Clone the repository, install the dependencies with ",
    "scripts/00_install_dependencies.R, and run this script from inside the ",
    "checkout. See inst/REPRODUCIBILITY.md for the full protocol.",
    call. = FALSE
  )
}

entry <- file.path(root, "scripts", "simulation", "09_quick_repro.R")
forwarded <- commandArgs(trailingOnly = TRUE)
message("Repository root: ", root)
message("Forwarding to the v2 quick reproduction entry: ", entry)
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(shQuote(entry), forwarded),
  stdout = "", stderr = ""
)
quit(status = status)

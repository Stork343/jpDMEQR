args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
setwd(root)
source(file.path("scripts", "00_source_v2.R"))
source_v2_module(".", "profile_v2", environment())
out_dir <- file.path("results", "scratch", "profile_chart_audit_20261022")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

integrate_interval <- function(fun) {
  stats::integrate(fun, -1, 1, rel.tol = 1e-9, abs.tol = 1e-11,
                   subdivisions = 400L)$value
}

# 使用两次独立 Normal 观察值的联合密度，逐个区域计算积分。
normal_chart_mean <- function(beta, tau, h, lambda = 1,
                              quantity = "score") {
  stopifnot(quantity %in% c("mass", "score", "hessian"), h >= 0)
  a <- tau - 1
  b <- tau
  C <- 1 / (2 * lambda)
  weight <- function(s1, s2, curvature) {
    if (quantity == "mass") return(rep(1, length(s1 + s2)))
    if (quantity == "score") return((s1 + s2) / 2)
    curvature
  }
  tail_probability <- function(threshold, sign) {
    stats::pnorm(threshold, lower.tail = sign < 0)
  }
  total <- 0

  # 没有 interior dual 坐标时，两个端点与两个残差方向确定区域。
  for (sign1 in c(-1, 1)) for (sign2 in c(-1, 1)) {
    s1 <- if (sign1 < 0) a else b
    s2 <- if (sign2 < 0) a else b
    gamma <- (s1 + s2) / (2 * lambda)
    total <- total + weight(s1, s2, 0) *
      tail_probability(beta + gamma + sign1 * h, sign1) *
      tail_probability(beta + gamma + sign2 * h, sign2)
  }

  # 一个 interior 坐标有两个位置；以 v 代换消除 inverse CDF 导数。
  for (sign in c(-1, 1)) {
    endpoint <- if (sign < 0) a else b
    total <- total + 2 * integrate_interval(function(v) {
      K <- kernel_epan_pdf_v2(v)
      s <- a + kernel_epan_cdf_v2(v)
      gamma <- (s + endpoint) / (2 * lambda)
      curvature <- if (h == 0) rep(lambda, length(v)) else
        lambda * K / (K + 2 * lambda * h)
      weight(s, endpoint, curvature) * stats::dnorm(beta + gamma + h * v) *
        tail_probability(beta + gamma + sign * h, sign) * (C * K + h)
    })
  }

  # 两个 interior 坐标的 Gram 矩阵秩为一，行列式在 h=0 时为零。
  if (h > 0) {
    total <- total + integrate_interval(function(v1) {
      vapply(v1, function(first) integrate_interval(function(v2) {
        K1 <- kernel_epan_pdf_v2(first)
        K2 <- kernel_epan_pdf_v2(v2)
        s1 <- a + kernel_epan_cdf_v2(first)
        s2 <- a + kernel_epan_cdf_v2(v2)
        gamma <- (s1 + s2) / (2 * lambda)
        curvature <- lambda * (K1 + K2) / (K1 + K2 + 2 * lambda * h)
        weight(s1, s2, curvature) *
          stats::dnorm(beta + gamma + h * first) *
          stats::dnorm(beta + gamma + h * v2) *
          (h * C * (K1 + K2) + h^2)
      }), numeric(1))
    })
  }
  total
}

rows <- list()
for (tau in c(0.25, 0.5, 0.75)) {
  beta_star <- stats::uniroot(
    function(beta) normal_chart_mean(beta, tau, 0),
    c(-2, 2), tol = 1e-10
  )$root
  for (h in c(0.2, 0.1, 0.05, 0.025)) {
    ordinary <- stats::uniroot(
      function(beta) normal_chart_mean(beta, tau, h),
      c(-2, 2), tol = 1e-10
    )$root
    paired <- stats::uniroot(function(beta) {
      2 * normal_chart_mean(beta, tau, h / 2) -
        normal_chart_mean(beta, tau, h)
    }, c(-2, 2), tol = 1e-10)$root
    mass <- normal_chart_mean(beta_star, tau, h, quantity = "mass")
    step <- 1e-4
    derivative <- -(normal_chart_mean(beta_star + step, tau, h) -
                      normal_chart_mean(beta_star - step, tau, h)) / (2 * step)
    curvature <- normal_chart_mean(beta_star, tau, h, quantity = "hessian")
    stopifnot(abs(mass - 1) < 1e-7, abs(derivative - curvature) < 1e-6)
    rows[[length(rows) + 1L]] <- data.frame(
      tau = tau, h = h, beta_star = beta_star, beta_h = ordinary,
      beta_paired = paired, ordinary_error = ordinary - beta_star,
      paired_error = paired - beta_star, probability_mass = mass,
      integrated_curvature = curvature, score_derivative = derivative,
      mass_error = abs(mass - 1), derivative_error = abs(derivative - curvature)
    )
    utils::write.csv(do.call(rbind, rows), file.path(out_dir, "normal_joint_charts.csv"),
                     row.names = FALSE)
  }
}
result <- do.call(rbind, rows)
utils::write.csv(result, file.path(out_dir, "normal_joint_charts.csv"), row.names = FALSE)
display <- result[result$tau == 0.25, , drop = FALSE]
table_rows <- vapply(seq_len(nrow(display)), function(i) {
  sprintf("%.3f & %.9f & %.9f \\\\", display$h[i],
          abs(display$ordinary_error[i]), abs(display$paired_error[i]))
}, character(1))
writeLines(c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\caption{Deterministic chart quadrature for two independent Normal outcomes with $X=Z=1$, $\\tau=0.25$ and $\\Lambda=1$.}",
  "\\begin{tabular}{rrr}",
  "\\toprule",
  "$h$ & $|\\beta_h^\\star-\\beta^\\star|$ & $|\\beta_h^{R,\\star}-\\beta^\\star|$ \\\\",
  "\\midrule", table_rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"
), file.path("..", "research", "profile-score-theory-20261022", "joint_chart_table.tex"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Normal joint-density chart checks completed: ", nrow(result))

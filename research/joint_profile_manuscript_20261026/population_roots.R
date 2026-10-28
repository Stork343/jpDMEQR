library(ggplot2)
library(statmod)
output <- "results/tables/joint_profile_manuscript_20261026"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
inverse_kernel <- function(u) 2 * sin(asin(2 * u - 1) / 3)
u <- seq(0, 1, length.out = 1001L)
q <- inverse_kernel(u)
if (max(abs(0.5 + .75 * (q - q^3 / 3) - u)) > 1e-12) stop("The inverse kernel identity failed.")
population_score <- function(beta, tau, lambda, bandwidth) {
  tau - integrate(function(u) pnorm(beta + (tau - 1 + u) / lambda + bandwidth * inverse_kernel(u)),
                    0, 1, subdivisions = 1000L, rel.tol = 1e-12, abs.tol = 1e-12)$value
}
population_root <- function(tau, lambda, bandwidth) {
  uniroot(function(beta) population_score(beta, tau, lambda, bandwidth), c(-8, 8), tol = 1e-12)$root
}
quadrature <- statmod::gauss.quad(512L, kind = "legendre")
nodes <- (quadrature$nodes + 1) / 2
weights <- quadrature$weights / 2
records <- list()
for (tau in c(.25, .5, .75)) {
  for (lambda in c(.5, 1, 2)) {
    target <- population_root(tau, lambda, 0)
    for (h in .4 * 2^(-(0:6))) {
      full <- population_root(tau, lambda, h)
      half <- population_root(tau, lambda, h / 2)
      paired <- 2 * half - full
      checked <- tau - sum(weights * pnorm(full + (tau - 1 + nodes) / lambda + h * inverse_kernel(nodes)))
      if (abs(checked) > 1e-8 || abs(population_score(full, tau, lambda, h)) > 1e-10) stop("The independent population integral check failed.")
      records[[length(records) + 1L]] <- data.frame(tau = tau, lambda = lambda, h = h,
                                                   target = target, single = full, paired = paired,
                                                   single_error = abs(full - target), paired_error = abs(paired - target),
                                                   quadrature_residual = abs(checked))
    }
  }
}
data <- do.call(rbind, records)
write.csv(data, file.path(output, "population_roots.csv"), row.names = FALSE)
lower <- data[data$tau == .25, ]
plot_data <- rbind(data.frame(lower[c("lambda", "h")], error = lower$single_error, method = "Single bandwidth"),
                    data.frame(lower[c("lambda", "h")], error = lower$paired_error, method = "Paired coefficients"))
if (any(plot_data$error <= 0)) stop("The logarithmic population plot requires positive errors.")
plot_data$penalty <- factor(paste0("Lambda = ", plot_data$lambda), levels = paste0("Lambda = ", c(.5, 1, 2)))
plot <- ggplot(plot_data, aes(h, error, colour = method, shape = method, linetype = method)) +
  geom_line(linewidth = .45) + geom_point(size = 1.7) +
  scale_x_log10() + scale_y_log10() + facet_wrap(~ penalty, nrow = 1) +
  scale_colour_manual(values = c("Paired coefficients" = "#0072B2", "Single bandwidth" = "#555555")) +
  scale_shape_manual(values = c("Paired coefficients" = 17, "Single bandwidth" = 16)) +
  scale_linetype_manual(values = c("Paired coefficients" = "solid", "Single bandwidth" = "dashed")) +
  labs(x = "Bandwidth", y = "Absolute population target displacement", colour = NULL, shape = NULL, linetype = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                                 strip.background = element_rect(fill = "white"))
dir.create("manuscript/generated/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("manuscript/generated/figures/fig_joint_smoothing.pdf", plot, width = 6.6, height = 2.7,
         units = "in", device = grDevices::cairo_pdf)
slopes <- do.call(rbind, lapply(split(lower, lower$lambda), function(group) {
  group <- group[order(group$h), ]
  data.frame(lambda = group$lambda[1L], single_order = log(group$single_error[2L] / group$single_error[1L]) / log(2),
               paired_order = log(group$paired_error[2L] / group$paired_error[1L]) / log(2))
}))
write.csv(slopes, file.path(output, "population_orders.csv"), row.names = FALSE)
print(slopes)

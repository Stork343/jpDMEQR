library(data.table)
library(ggplot2)
library(grid)
source <- "results/tables/soynam_20261031"
arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments)) arguments[1L] else "all"
if (!mode %in% c("all", "joint")) stop("Unknown plotting mode.")
prediction_source <- if (mode == "all") source else "results/preflight/soynam_20261031/joint_prediction"
destination <- "manuscript/generated/soynam_profile"
effects <- fread(file.path(source, "inference.csv"))[setting == "reference"]
comparison <- fread(file.path(prediction_source, "prediction_comparisons.csv"))[calibrated == TRUE & comparator == "joint_genomic"]
if (nrow(effects) != 6L || nrow(comparison) != 3L ||
    any(!is.finite(effects$simultaneous_lower)) || any(!is.finite(comparison$relative_lower)))
  stop("The plotted result family is incomplete.")
effects[, year := factor(sub("ACC10_", "", coordinate), levels = c("2013", "2014"))]
style <- theme_classic(base_size = 10.5, base_family = "sans") +
  theme(legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size = 9.5), plot.title = element_text(face = "bold", size = 11),
    axis.title = element_text(size = 10.5), plot.margin = margin(8, 8, 6, 8),
    legend.key.width = unit(10, "pt"))
dodge <- position_dodge(width = .07)
left <- ggplot(effects, aes(tau, estimate_bu_ac, colour = year, shape = year)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = .35) +
  geom_errorbar(aes(ymin = simultaneous_lower, ymax = simultaneous_upper),
    position = dodge, width = .025, linewidth = .55) +
  geom_point(position = dodge, size = 2.3) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(breaks = c(.25, .5, .75), labels = c("0.25", "0.50", "0.75"), limits = c(.15, .85)) +
  scale_y_continuous(limits = c(0, max(effects$simultaneous_upper) * 1.08), expand = expansion(mult = c(0, .01))) +
  labs(title = "A  Canopy association", x = "Quantile", y = "Yield association (bu/ac)\nper 10 percentage points") +
  style + guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE))
right <- ggplot(comparison, aes(tau, relative_percent)) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = .35) +
  geom_errorbar(aes(ymin = relative_lower, ymax = relative_upper), width = .025, linewidth = .55, colour = "#00795D") +
  geom_point(size = 2.3, shape = 15, colour = "#00795D") +
  scale_x_continuous(breaks = c(.25, .5, .75), labels = c("0.25", "0.50", "0.75"), limits = c(.15, .85)) +
  expand_limits(y = 0) + labs(title = "B  Held-out prediction", x = "Quantile", y = "Reduction in pinball loss (%)") +
  style
if (nrow(ggplot_build(left)$data[[3L]]) != 6L || nrow(ggplot_build(right)$data[[3L]]) != 3L)
  stop("The plotted points differ from their source rows.")
draw <- function() {
  first <- ggplotGrob(left)
  second <- ggplotGrob(right)
  if (length(first$heights) != length(second$heights)) stop("The panel layouts have different structures.")
  first$heights <- second$heights <- unit.pmax(first$heights, second$heights)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 2)))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(first)
  popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid.draw(second)
  popViewport()
  popViewport()
}
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
cairo_pdf(file.path(destination, "results.pdf"), width = 7.3, height = 3.7, family = "sans")
draw()
dev.off()
dir.create("results/figures/soynam_20261031", recursive = TRUE, showWarnings = FALSE)
png("results/figures/soynam_20261031/results.png", width = 2190, height = 1110, res = 300, type = "cairo")
draw()
dev.off()
saveRDS(list(left = left, right = right), "results/figures/soynam_20261031/plots.rds")
cat("Exported the two-panel empirical result figure.\n")

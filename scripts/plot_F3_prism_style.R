# =============================================================================
# plot_F3_prism_style.R - Fig F3 (BCAN vs NCAN in the same three-protein condensates, paired)
# drawn in the same GraphPad Prism style as the other panels of the figure:
#   Arial 10 pt for every text, bars #FF8080 (BCAN) / #DE8BF9 (NCAN), light-grey points with black edge,
#   mean +/- SEM with caps, x axis drawn at y = 0, y axis -2 to 6, bracket + asterisks.
#   Lines join the two values of the same preparation (paired design).
# Plot-area size is fixed (PANEL_W x PANEL_H mm) so it matches the companion Prism panels; change it here if needed.
#
# Run in: RStudio (open the LLPS_analysis folder as working directory, then Source)
#   or cmd:  Rscript scripts\plot_F3_prism_style.R
# Input : results/position/figures/FigF3_data_per_replicate_paired.csv
# Output: results/position/figures/FigF3_paired_prism_style.pdf / .svg / .png + _legend.txt
#   (the .svg keeps the text as real Arial text - best for Illustrator)
# Author: Bin Luo, 2026. Written with AI assistance (Claude, Anthropic).
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork)})

PANEL_W <- 15.5   # mm, width of the plotting area (matches the companion Prism panels)
PANEL_H <- 34.3   # mm, height of the plotting area (y from -2 to 6)
FONT    <- "Arial"
COLS    <- c(BCAN = "#FF8080", NCAN = "#DE8BF9")
fig_dir <- file.path("results", "position", "figures")

d <- read_csv(file.path(fig_dir, "FigF3_data_per_replicate_paired.csv"), show_col_types = FALSE)
long <- d %>% select(replicate, BCAN = BCAN_shift_pct, NCAN = NCAN_shift_pct) %>%
  pivot_longer(c(BCAN, NCAN), names_to = "client", values_to = "y") %>%
  mutate(x = ifelse(client == "BCAN", 1, 2))
summ <- long %>% group_by(client, x) %>% summarise(m = mean(y), se = sd(y) / sqrt(n()), .groups = "drop")
p <- t.test(d$NCAN_radial_ratio, d$BCAN_radial_ratio, paired = TRUE)$p.value     # same P as on shift_pct
stars <- if (p < 0.0001) "****" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"

lw_axis <- 0.45    # ~1 pt (2 px at the Prism export scale)
lw_bar  <- 0.38    # error bars, bracket
bw      <- 0.77    # bar width (x units)
g <- ggplot() +
  geom_col(data = summ, aes(x, m, fill = client), width = bw) +
  geom_hline(yintercept = 0, linewidth = lw_axis) +
  geom_errorbar(data = summ, aes(x, ymin = m - se, ymax = m + se), width = 0.37, linewidth = lw_bar) +
  geom_line(data = long, aes(x, y, group = replicate), colour = "grey25", linewidth = 0.35) +
  geom_point(data = long, aes(x, y), shape = 21, fill = "#D4D4D4", colour = "black", size = 1.2, stroke = 0.35) +
  annotate("segment", x = c(1, 1, 2), xend = c(2, 1, 2), y = c(6.37, 6.37, 6.37), yend = c(6.37, 4.68, 5.8),
           linewidth = lw_bar) +
  annotate("text", x = 1.5, y = 6.47, label = stars, vjust = 0, family = FONT, size = 10 / .pt) +
  scale_fill_manual(values = COLS, guide = "none") +
  scale_x_continuous(breaks = c(1, 2), labels = c("BCAN", "NCAN"), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-2, 6, 2), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.202, 2.827), ylim = c(-2, 6), expand = FALSE, clip = "off") +   # bracket may sit above 6
  labs(x = "ACAN + BCAN + NCAN", y = "Peripheral shift\nrelative to ACAN (%)") +
  theme_classic(base_size = 10, base_family = FONT) +
  theme(axis.line.y = element_line(linewidth = lw_axis, lineend = "square"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(linewidth = lw_axis), axis.ticks.length.y = unit(1.3, "mm"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 10, colour = "black"),
        axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black"),
        plot.margin = margin(6, 9, 2, 2, "mm"))
pp <- g + plot_layout(widths = unit(PANEL_W, "mm"), heights = unit(PANEL_H, "mm"))

W <- 54; H <- 66   # mm, whole canvas (plot area fixed above; the rest is text/margins)
ggsave(file.path(fig_dir, "FigF3_paired_prism_style.pdf"), pp, width = W, height = H, units = "mm", device = cairo_pdf)
ggsave(file.path(fig_dir, "FigF3_paired_prism_style.png"), pp, width = W, height = H, units = "mm", dpi = 600)
if (requireNamespace("svglite", quietly = TRUE))
  svglite::svglite(file.path(fig_dir, "FigF3_paired_prism_style.svg"), width = W / 25.4, height = H / 25.4) else
  svg(file.path(fig_dir, "FigF3_paired_prism_style.svg"), width = W / 25.4, height = H / 25.4)
print(pp); invisible(dev.off())
# on Linux the metric-identical Liberation Sans stands in for Arial; name the font Arial in the SVG
svgf <- file.path(fig_dir, "FigF3_paired_prism_style.svg")
writeLines(gsub("\"Liberation Sans\"|'Liberation Sans'|Liberation Sans", "Arial", readLines(svgf, warn = FALSE)), svgf)

leg <- paste0("BCAN and NCAN compared within the same three-protein (ACAN + BCAN + NCAN) condensates. Peripheral shift relative to ACAN = ",
  "(radial ratio - 1) x 100, where the radial ratio is the client's intensity-weighted mean radial position divided by that of ACAN in the same condensate ",
  "(> 0, client sits further towards the rim than ACAN). Lines join the BCAN and NCAN values of the same independent preparation; bars, mean ± SEM; ",
  "n = ", nrow(d), " independent preparations. Two-tailed paired t-test on the radial ratio: P = ", signif(p, 2), " (", stars, ", P < 0.001). ",
  "Asterisks: * P < 0.05, ** P < 0.01, *** P < 0.001, **** P < 0.0001.")
writeLines(c("Figure F3, Prism-style version (FigF3_paired_prism_style)", "", leg),
           file.path(fig_dir, "FigF3_paired_prism_style_legend.txt"), useBytes = TRUE)
cat("P =", signif(p, 3), stars, "\nWritten:", file.path(fig_dir, "FigF3_paired_prism_style.pdf/.svg/.png"), "\n")

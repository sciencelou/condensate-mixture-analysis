# =============================================================================
# plot_driver.R
# Called automatically by make_figures.R (or open in RStudio and click "Source").
#
# Does BCAN raise NCAN incorporation only by making ACAN condensates bigger,
# and does NCAN follow ACAN?
#   Fig E6  NCAN coverage of ACAN condensates at matched condensate size, 2 vs 3 proteins
#   Fig E7  within each image: correlation between NCAN and ACAN amount across condensates
# Each dot = one field of view (FOV); bar = mean +/- SEM; P values on the graphs.
# =============================================================================
if (!exists("res_dir")) res_dir <- "results/three_protein"
for (p in c("ggplot2", "dplyr", "readr", "tidyr", "officer")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr); library(tidyr) })
source(file.path(if (exists("scripts")) scripts else "scripts", "style_prism.R"), local = TRUE)  # shared Prism-like style
fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE)
d <- read_csv(file.path(res_dir, "driver_analysis_objects.csv"), show_col_types = FALSE) %>%
  filter(condition %in% c("ACAN+NCAN", "Triple|NCAN")) %>%
  mutate(mix = factor(ifelse(condition == "ACAN+NCAN", "ACAN + NCAN", "ACAN + BCAN + NCAN"),
                      levels = c("ACAN + NCAN", "ACAN + BCAN + NCAN")))
n_rep <- d %>% group_by(mix) %>% summarise(n = n_distinct(replicate), .groups = "drop")
unit <- if (all(n_rep$n >= 3)) "replicate" else "fov"
dot  <- if (unit == "replicate") "Each dot is one independent preparation (mean of its FOVs)" else "Each dot is one FOV"
fmt_p <- function(p) ifelse(p < 1e-4, "P < 0.0001", paste0("P = ", formatC(signif(p, 2), format = "fg")))
mwu <- function(a, b) if (unit == "replicate") student_p(a, b) else suppressWarnings(wilcox.test(a, b, exact = TRUE))$p.value
test_name <- if (unit == "replicate") TEST_T else "two-sided Mann-Whitney U test"
save_both <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".svg")), p, width = w, height = h, units = "mm", device = svg)
}
bar_layers <- prism_layers()

# ---- Fig E6: matched condensate size ------------------------------------------------
bins <- c(50, 150, 400)
e6 <- d %>% filter(area >= bins[1], area < bins[3]) %>%
  mutate(size = factor(ifelse(area < bins[2], area_range(bins[1], bins[2]), area_range(bins[2], bins[3])),
                       levels = c(area_range(bins[1], bins[2]), area_range(bins[2], bins[3])))) %>%
  group_by(mix, replicate, fov, size) %>% summarise(y = 100 * mean(client_cover), n = n(), .groups = "drop")
if (unit == "replicate") e6 <- e6 %>% group_by(mix, replicate, size) %>% summarise(y = mean(y), .groups = "drop")
st6 <- e6 %>% group_by(size) %>% summarise(p = mwu(y[mix == levels(mix)[1]], y[mix == levels(mix)[2]]), .groups = "drop") %>%
  mutate(p_adj = p, lab = fmt_p(p_adj), yy = 108)
write_csv(e6, file.path(fig_dir, paste0("FigE6_data_per_", unit, ".csv")))
pE6 <- ggplot(e6, aes(mix, y, fill = mix)) + bar_layers +
  geom_segment(data = st6, aes(x = 1.05, xend = 1.95, y = yy, yend = yy), inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = st6, aes(x = 1.5, y = yy, label = lab), inherit.aes = FALSE, vjust = -0.45, size = 2.3) +
  facet_wrap(~size, labeller = labeller(size = function(x) paste0("Condensates\n", x))) +
  scale_fill_manual(values = fill_cols) + scale_y_continuous(limits = c(0, 125), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "NCAN-positive area of\nACAN condensates (%)") + theme_pub
save_both(pE6, "FigE6_NCAN_in_size_matched_condensates", 72, 64)

# ---- Fig E7: within-image NCAN vs ACAN correlation ---------------------------------
e7 <- d %>% filter(area >= 50) %>% group_by(mix, replicate, fov) %>% filter(n() >= 8) %>%
  summarise(y = cor(client, acan, method = "spearman"), n = n(), .groups = "drop")
if (unit == "replicate") e7 <- e7 %>% group_by(mix, replicate) %>% summarise(y = mean(y), .groups = "drop")
p7 <- mwu(e7$y[e7$mix == levels(e7$mix)[1]], e7$y[e7$mix == levels(e7$mix)[2]])
write_csv(e7, file.path(fig_dir, paste0("FigE7_data_per_", unit, ".csv")))
pE7 <- ggplot(e7, aes(mix, y, fill = mix)) + bar_layers + zero_axis() +
  annotate("segment", x = 1.05, xend = 1.95, y = 0.8, yend = 0.8, linewidth = 0.3) +
  annotate("text", x = 1.5, y = 0.8, label = fmt_p(p7), vjust = -0.45, size = 2.3) +
  scale_fill_manual(values = fill_cols) + scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5), expand = c(0, 0)) +
  labs(x = NULL, y = "Spearman's \u03c1\n(NCAN vs ACAN)") + theme_pub
save_both(pE7, "FigE7_NCAN_ACAN_correlation", 48, 66)

# ---- legends -------------------------------------------------------------------
nf <- d %>% distinct(mix, fov) %>% count(mix)
nr <- d %>% distinct(mix, replicate) %>% count(mix)
nfov <- if (unit == "replicate") sprintf("n = %d (ACAN + NCAN) and %d (ACAN + BCAN + NCAN) independent preparations", nr$n[1], nr$n[2]) else
  sprintf("n = %d FOVs (ACAN + NCAN) and %d FOVs (ACAN + BCAN + NCAN) [single preparation per condition; FOVs are technical replicates]", nf$n[1], nf$n[2])
common <- "Same protein concentrations and imaging settings in both mixes ([concentrations, 6% PEG, buffer, time after mixing]). Image analysis is described in Methods."
legE6 <- paste0("NCAN incorporation into ACAN condensates of the same size, without and with BCAN. For ACAN condensates of ", area_range(bins[1], bins[2]), " or ", area_range(bins[2], bins[3]), ", ",
  "the percentage of the condensate area covered by NCAN signal was averaged per field of view (FOV) and then per preparation. Comparing condensates of the same size tests whether BCAN ",
  "increases NCAN incorporation only by making ACAN condensates larger. ", dot, "; bars, mean \u00b1 SEM. ", nfov,
  ". Statistics: ", test_name, " for each size class. ", common)
legE7 <- paste0("Relationship between the amounts of NCAN and ACAN in individual condensates. Within each image, the mean NCAN and mean ACAN intensities of all ACAN condensates ",
  "of at least ", area_txt(50), " were rank-correlated (Spearman correlation); positive values mean that condensates with more ACAN also contain more NCAN. ", dot, "; bars, mean \u00b1 SEM. ",
  nfov, ". Statistics: ", test_name, ". ", if (TRI_SAT) "ACAN was partly saturated in the three-protein images, which compresses ACAN differences between condensates and is expected to weaken, not create, a positive correlation. " else "",
  common)
doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "Figure legends - which protein drives NCAN incorporation (auto-generated by plot_driver.R)", style = "heading 1")
for (L in list(list("Figure E6", "FigE6_NCAN_in_size_matched_condensates", legE6, pE6, 60),
               list("Figure E7", "FigE7_NCAN_ACAN_correlation", legE7, pE7, 62))) {
  writeLines(c(paste0(L[[1]], " (", L[[2]], ")"), "", L[[3]]), file.path(fig_dir, paste0(L[[2]], "_legend.txt")), useBytes = TRUE)
  doc <- officer::body_add_par(doc, paste0(L[[1]], " - ", L[[2]]), style = "heading 2")
  doc <- officer::body_add_par(doc, L[[3]], style = "Normal")
  cap <- paste(strwrap(L[[3]], width = 95), collapse = "\n")
  ggsave(file.path(fig_dir, paste0(L[[2]], "_with_legend.pdf")),
         L[[4]] + labs(caption = cap) + theme(plot.caption = element_text(hjust = 0, size = 6.5), plot.caption.position = "plot"),
         width = 180, height = L[[5]] + 12 + 3.2 * length(strsplit(cap, "\n")[[1]]), units = "mm", device = cairo_pdf)
}
print(doc, target = file.path(fig_dir, "Figure_legends_driver.docx"))
sink(file.path(fig_dir, "stats_driver.txt"))
cat("Unit of analysis:", ifelse(unit == "replicate", "independent preparation (mean of its FOVs)", "FOV"), "| test:", test_name, "\n\n")
cat("== E6 NCAN coverage, size-matched (", test_name, ") ==\n"); print(as.data.frame(st6))
cat("\nGroup means (over", ifelse(unit == "replicate", "preparations", "FOVs"), "):\n"); print(as.data.frame(e6 %>% group_by(size, mix) %>% summarise(mean = mean(y), sd = sd(y), n = n(), .groups = "drop")))
cat("\n== E7 NCAN-ACAN Spearman rho (computed per FOV, averaged per", ifelse(unit == "replicate", "preparation", "FOV"), "),", test_name, "==\nP =", p7, "\n"); print(as.data.frame(e7))
sink()
cat("Driver figures written to", normalizePath(fig_dir), "\n")

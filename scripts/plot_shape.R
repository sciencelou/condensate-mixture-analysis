# =============================================================================
# plot_shape.R
# Called automatically by make_figures.R (or open in RStudio and click "Source").
#
# Size and shape of ACAN condensates: ACAN alone vs ACAN + BCAN vs ACAN + NCAN
# (ACAN channel only; output of shape_analysis.py)
#   Fig D   condensate size and elongation (aspect ratio)
#   Fig S3  number of condensates, % area covered, circularity of size-matched
#           condensates, solidity
#   Fig S4  elongation with different minimum-size cut-offs (robustness)
# Each dot = one field of view (FOV); bar = mean +/- SEM; P values on the graphs.
# Output: vector PDF + SVG + legend (.txt/.docx) + data CSV + stats.txt
# =============================================================================

# ---- 0. settings: EDIT THIS LINE -------------------------------------------
if (!exists("res_dir")) res_dir <- "results/shape"  # folder created by shape_analysis.py (e.g. "D:/LLPS/acan_bcan/results_shape")
cond_order  <- c("ACAN alone", "ACAN+BCAN", "ACAN+NCAN")
cond_labels <- c("ACAN alone" = "ACAN alone", "ACAN+BCAN" = "ACAN + BCAN", "ACAN+NCAN" = "ACAN + NCAN")
size_window <- c(100, 400)    # px; circularity is compared only between condensates of similar size

# ---- 1. packages ------------------------------------------------------------
for (p in c("ggplot2", "dplyr", "readr", "tidyr", "officer", "lme4", "lmerTest")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr); library(tidyr) })
source(file.path(if (exists("scripts")) scripts else "scripts", "style_prism.R"), local = TRUE)  # shared Prism-like style
cond_cols <- fill_cols

fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE)
obj  <- read_csv(file.path(res_dir, "shape_objects.csv"), show_col_types = FALSE)
fov  <- read_csv(file.path(res_dir, "shape_fov_summary.csv"), show_col_types = FALSE)
obj <- obj %>% filter(condition %in% cond_order); fov <- fov %>% filter(condition %in% cond_order)   # triple mix is plotted by plot_triple.R
lev  <- intersect(cond_order, unique(fov$condition))
obj$condition <- factor(obj$condition, levels = lev); fov$condition <- factor(fov$condition, levels = lev)
has_um <- all(!is.na(fov$pixel_um))
area_unit <- if (has_um) "\u00b5m\u00b2" else "px"
# statistical unit: independent preparations if every condition has >= 3 of them, otherwise fields of view
n_rep <- fov %>% group_by(condition) %>% summarise(n = n_distinct(replicate), .groups = "drop")
unit <- if (nrow(n_rep) >= 2 && all(n_rep$n >= 3)) "replicate" else "fov"
to_unit <- function(d) if (unit == "replicate") d %>% group_by(condition, replicate) %>% summarise(y = mean(y), .groups = "drop") else d
dot <- if (unit == "replicate") "Each dot is one independent preparation (mean of its FOVs)" else "Each dot is one field of view (FOV)"

# size-matched circularity per FOV
circ_sm <- obj %>% filter(area_px >= size_window[1], area_px <= size_window[2]) %>%
  group_by(condition, replicate, fov) %>%
  summarise(circ_sizematched = median(circularity), n_sm = n(), .groups = "drop")
fov <- fov %>% left_join(circ_sm, by = c("condition", "replicate", "fov"))
if (has_um) fov$median_area_plot <- fov$median_area_px * fov$pixel_um^2 else fov$median_area_plot <- fov$median_area_px

# ---- 2. helpers (same style as plot_recruitment.R) ---------------------------
compare <- function(d, log = FALSE) {
  d <- droplevels(d); lv <- levels(d$condition); if (length(lv) < 2) return(NULL)
  out <- if (length(lv) == 2) {
    a <- d$y[d$condition == lv[1]]; b <- d$y[d$condition == lv[2]]
    if (log) { a <- log10(a); b <- log10(b) }
    data.frame(group1 = lv[1], group2 = lv[2], p = student_p(a, b), p_adj = student_p(a, b))
  } else anova_tukey(d$y, d$condition, log)
  out
}
fmt_p <- function(p) ifelse(p < 1e-4, "P < 0.0001", paste0("P = ", formatC(signif(p, 2), format = "fg")))
add_p <- function(plt, cmp, ytop, facet = NULL) {
  if (is.null(cmp)) return(plt)
  lv <- levels(droplevels(plt$data$condition))
  cmp$x1 <- match(cmp$group1, lv); cmp$x2 <- match(cmp$group2, lv)
  cmp <- cmp[order(cmp$x2 - cmp$x1), ]                       # short brackets first (lowest)
  cmp$y <- ytop * (1 + 0.11 * (seq_len(nrow(cmp)) - 1)); cmp$lab <- fmt_p(cmp$p_adj)
  if (!is.null(facet)) cmp$metric <- facet
  plt + geom_segment(data = cmp, aes(x = x1 + 0.05, xend = x2 - 0.05, y = y, yend = y),
                     inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = cmp, aes(x = (x1 + x2) / 2, y = y, label = lab), inherit.aes = FALSE,
              vjust = -0.45, size = 2.2) +
    expand_limits(y = max(cmp$y) * 1.1)
}
dotbar <- function(d, ylab) {
  ggplot(d, aes(condition, y, fill = condition)) + prism_layers() +
    scale_fill_manual(values = fill_cols) + scale_x_discrete(labels = cond_labels) +
    labs(x = NULL, y = ylab) + theme_pub
}
save_both <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".svg")), p, width = w, height = h, units = "mm", device = svg)
}
metric_panel <- function(d, ylab, log = FALSE, ymin0 = TRUE) {
  d <- to_unit(d)
  cmp <- compare(d, log)
  p <- dotbar(d, ylab)
  if (ymin0) p <- p + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  list(plot = add_p(p, cmp, ytop = max(d$y) * 1.08), cmp = cmp)
}

# ---- 3. Fig D: size and elongation --------------------------------------------
dSize <- fov %>% transmute(condition, replicate, fov, y = median_area_plot)
dAR   <- fov %>% transmute(condition, replicate, fov, y = median_aspect_ratio)
pD1 <- metric_panel(dSize, paste0("Condensate area (", area_unit, ")"), log = FALSE)   # untransformed, as in the paper
pD2 <- metric_panel(dAR, "Aspect ratio")
pD2$plot <- pD2$plot + geom_hline(yintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey40")
write_csv(fov %>% select(condition, replicate, fov, n_condensates, n_condensates_for_shape,
                         median_area = median_area_plot, median_aspect_ratio),
          file.path(fig_dir, "FigD_data_per_fov.csv"))
save_both(pD1$plot, "FigD1_condensate_size", 50, 60)
save_both(pD2$plot, "FigD2_condensate_elongation", 50, 60)

# ---- 4. Fig S3: other descriptors --------------------------------------------
dN    <- fov %>% transmute(condition, replicate, fov, y = n_condensates)
dCov  <- fov %>% transmute(condition, replicate, fov, y = pct_area_covered)
dCirc <- fov %>% transmute(condition, replicate, fov, y = circ_sizematched) %>% filter(!is.na(y))
dSol  <- fov %>% transmute(condition, replicate, fov, y = median_solidity)
pS3a <- metric_panel(dN,   "Condensates per field")
pS3b <- metric_panel(dCov, "Area fraction (%)")
pS3c <- metric_panel(dCirc, "Circularity")
pS3d <- metric_panel(dSol, "Solidity")
write_csv(fov %>% select(condition, replicate, fov, n_condensates, pct_area_covered,
                         circ_sizematched, n_sm, median_solidity),
          file.path(fig_dir, "FigS3_data_per_fov.csv"))
save_both(pS3a$plot, "FigS3a_number_per_fov", 45, 58)
save_both(pS3b$plot, "FigS3b_area_covered", 45, 58)
save_both(pS3c$plot, "FigS3c_circularity_size_matched", 45, 58)
save_both(pS3d$plot, "FigS3d_solidity", 45, 58)

# ---- 5. Fig S4: elongation vs minimum-size cut-off -------------------------------
cut_cols <- grep("^median_aspect_ratio_min", names(fov), value = TRUE)
dS4 <- fov %>% select(condition, replicate, fov, all_of(cut_cols)) %>%
  pivot_longer(all_of(cut_cols), names_to = "cut", values_to = "y") %>%
  mutate(cut = factor(paste0("\u2265 ", vapply(as.numeric(sub("median_aspect_ratio_min", "", cut)), area_txt, "")),
                      levels = paste0("\u2265 ", vapply(sort(as.numeric(sub("median_aspect_ratio_min", "", cut_cols))), area_txt, ""))))
write_csv(dS4, file.path(fig_dir, "FigS4_data_per_fov.csv"))
if (unit == "replicate") dS4 <- dS4 %>% group_by(condition, replicate, cut) %>% summarise(y = mean(y), .groups = "drop")
pS4 <- ggplot(dS4, aes(cut, y, colour = condition, group = condition)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.4) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.15, linewidth = 0.3) +
  scale_colour_manual(values = line_cols, labels = function(x) cond_labels[x]) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey40") +
  labs(x = "Minimum condensate area included", y = "Aspect ratio", colour = NULL) +
  theme_pub + theme_flat_x + theme(legend.position = "top", legend.key.size = unit(3, "mm"))
save_both(pS4, "FigS4_elongation_cutoff_robustness", 85, 60)

# ---- 6. legends --------------------------------------------------------------
cl <- function(x) gsub("\n", " ", cond_labels[as.character(x)])
n_tab <- fov %>% group_by(condition) %>% summarise(n_fov = n(), n_obj = sum(n_condensates),
                                                   n_shape = sum(n_condensates_for_shape), n_sm = sum(n_sm, na.rm = TRUE), .groups = "drop")
n_fov_txt <- if (length(unique(n_tab$n_fov)) == 1) sprintf("%d FOVs per condition", n_tab$n_fov[1]) else
  paste(sprintf("%s, %d FOVs", cl(n_tab$condition), n_tab$n_fov), collapse = "; ")
n_txt <- if (unit == "replicate") { if (length(unique(n_rep$n)) == 1) sprintf("%d independent preparations per condition", n_rep$n[1]) else
  paste(sprintf("%s, %d independent preparations", cl(n_rep$condition), n_rep$n), collapse = "; ") } else
  paste0(n_fov_txt, " [from a single preparation per condition; FOVs are technical replicates]")
nobj  <- paste(sprintf("%s, %d", cl(n_tab$condition), n_tab$n_obj), collapse = "; ")
nshp  <- paste(sprintf("%s, %d", cl(n_tab$condition), n_tab$n_shape), collapse = "; ")
nsm   <- paste(sprintf("%s, %d", cl(n_tab$condition), n_tab$n_sm), collapse = "; ")
stat  <- paste0("Statistics: ", TEST_ANOVA, " on ", ifelse(unit == "replicate", "preparation", "FOV"), " values.")
common <- paste0("Only the ACAN channel was analysed. All conditions were prepared and imaged identically ",
  "([ACAN concentration, 6% PEG, buffer, time after mixing]; ACAN, [fluorophore]). Image analysis is described in Methods.")
legD1 <- paste0("Area of ACAN condensates", ifelse(is.na(PIXEL_UM), "", paste0(" (pixel size ", PIXEL_UM, " \u00b5m)")), ". ", dot, ifelse(unit == "fov", " (value = median area of all ACAN condensates in that image)", ""), "; ",
  "bars, mean \u00b1 SEM. n = ", n_txt, " (condensates: ", nobj, "). ", stat, " ", common)
legD2 <- paste0("Elongation of ACAN condensates, measured as the aspect ratio (long axis / short axis of the best-fitting ellipse; 1 = round, dashed line). ",
  "Only condensates of at least ", area_txt(50), " (50 px) were included, because the outline of smaller objects is dominated by the pixel grid. ",
  dot, ifelse(unit == "fov", " (median of its condensates)", ""), "; bars, mean \u00b1 SEM. n = ", n_txt, " (condensates analysed: ", nshp, "). ", stat, " ", common)
legS3 <- paste0("Additional descriptors of ACAN condensates. (a) Number of condensates per field of view. (b) Area fraction: percentage of the image area covered by condensates. ",
  "(c) Circularity (4\u03c0 \u00d7 area / perimeter\u00b2; 1 = perfect circle) of condensates of similar size only (", area_range(size_window[1], size_window[2]), "; ",
  "condensates analysed: ", nsm, "), because circularity computed from pixel outlines depends strongly on object size. ",
  "(d) Solidity (area / convex-hull area; 1 = no indentations). ", dot, "; bars, mean \u00b1 SEM. n = ", n_txt, ". ", stat, " ", common)
legS4 <- paste0("Elongation of ACAN condensates (Fig. D2) recalculated with different minimum condensate sizes. Points, mean \u00b1 SEM across ", ifelse(unit == "replicate", "independent preparations", "FOVs"), ". ",
  "n = ", n_txt, ". ", common)

panels <- list(
  list(file = "FigD1_condensate_size",            title = "Figure D1", text = legD1, plot = pD1$plot, h = 60),
  list(file = "FigD2_condensate_elongation",      title = "Figure D2", text = legD2, plot = pD2$plot, h = 60),
  list(file = "FigS3a-d (four files FigS3a..d)",  title = "Figure S3", text = legS3, plot = NULL,     h = 58),
  list(file = "FigS4_elongation_cutoff_robustness", title = "Figure S4", text = legS4, plot = pS4,    h = 60))
doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "Figure legends - ACAN condensate size and shape (auto-generated by plot_shape.R)", style = "heading 1")
doc <- officer::body_add_par(doc, "Items in [brackets] must be completed by hand.", style = "Normal")
all_txt <- character()
for (L in panels) {
  stem <- sub(" .*", "", L$file)
  writeLines(c(paste0(L$title, " (", L$file, ")"), "", L$text), file.path(fig_dir, paste0(stem, "_legend.txt")), useBytes = TRUE)
  all_txt <- c(all_txt, paste0(L$title, " (", L$file, ")"), L$text, "")
  doc <- officer::body_add_par(doc, paste0(L$title, " - ", L$file), style = "heading 2")
  doc <- officer::body_add_par(doc, L$text, style = "Normal")
  if (!is.null(L$plot)) {
    cap <- paste(strwrap(L$text, width = 95), collapse = "\n")
    pc <- L$plot + labs(caption = cap) + theme(plot.caption = element_text(hjust = 0, size = 6.5, lineheight = 1.1),
                                               plot.caption.position = "plot")
    if (grepl("S4", L$file)) pc <- pc + theme(legend.position = "top")
    ggsave(file.path(fig_dir, paste0(stem, "_with_legend.pdf")), pc, width = 180,
           height = L$h + 12 + 3.2 * length(strsplit(cap, "\n")[[1]]), units = "mm", device = cairo_pdf)
  }
}
writeLines(all_txt, file.path(fig_dir, "Figure_legends_shape.txt"), useBytes = TRUE)
print(doc, target = file.path(fig_dir, "Figure_legends_shape.docx"))

# ---- 7. statistics -------------------------------------------------------------
sink(file.path(fig_dir, "stats_shape.txt"))
cat("Unit of analysis:", ifelse(unit == "replicate", "independent preparation", "field of view (technical replicates)"), "-", TEST_ANOVA, "(p_adj = Tukey-adjusted P; ANOVA P in brackets per panel).\n")
for (nm in c("D1 size (median area)", "D2 aspect ratio", "S3a number", "S3b % area",
             "S3c circularity size-matched", "S3d solidity")) {
  cmp <- switch(nm, "D1 size (median area)" = pD1$cmp, "D2 aspect ratio" = pD2$cmp, "S3a number" = pS3a$cmp,
                "S3b % area" = pS3b$cmp, "S3c circularity size-matched" = pS3c$cmp, "S3d solidity" = pS3d$cmp)
  cat("\n==", nm, "==   one-way ANOVA P =", signif(attr(cmp, "anova_p"), 3), "\n"); print(cmp)
}
cat("\nSupporting linear mixed models on individual condensates (random intercept = FOV):\n")
ob <- obj %>% mutate(fov_id = interaction(condition, replicate, fov, drop = TRUE))
cat("\nlog10(area):\n"); print(summary(lmerTest::lmer(log10(area_px) ~ condition + (1 | fov_id), data = ob))$coefficients)
cat("\naspect ratio (>= 50 px):\n")
print(summary(lmerTest::lmer(aspect_ratio ~ condition + (1 | fov_id), data = filter(ob, area_px >= 50)))$coefficients)
sm <- filter(ob, area_px >= size_window[1], area_px <= size_window[2])
cat(sprintf("\nSize-matched condensates (%d-%d px) - circularity:\n", size_window[1], size_window[2]))
print(summary(lmerTest::lmer(circularity ~ condition + (1 | fov_id), data = sm))$coefficients)
cat(sprintf("\nSize-matched condensates (%d-%d px) - aspect ratio:\n", size_window[1], size_window[2]))
print(summary(lmerTest::lmer(aspect_ratio ~ condition + (1 | fov_id), data = sm))$coefficients)
cat("\nPer-condition means of FOV values:\n")
print(as.data.frame(fov %>% group_by(condition) %>% summarise(across(c(n_condensates, pct_area_covered, median_area_plot,
  median_aspect_ratio, circ_sizematched, median_solidity), ~ round(mean(.x, na.rm = TRUE), 3)))))
sink()
cat("Shape figures, legends, data and stats written to", normalizePath(fig_dir), "\n")

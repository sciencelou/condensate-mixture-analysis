# =============================================================================
# plot_triple.R
# Called automatically by make_figures.R (or open in RStudio and click "Source").
#
# Three-protein mix (ACAN + BCAN + NCAN) compared with the two-protein mixes.
#   Fig E1  % client puncta inside ACAN condensates: two-protein vs three-protein mix
#   Fig E2  % ACAN condensates containing the client: two-protein vs three-protein mix
#   Fig E3  where client puncta are (inside / touching / separate)
#   Fig E4  within the three-protein mix: BCAN vs NCAN in the same images (paired)
#   Fig E5  ACAN condensate size and elongation, three-protein vs two-protein mixes
# Each dot = one field of view (FOV); bar = mean +/- SEM; P values on the graphs.
# Output: vector PDF + SVG + legends (.txt/.docx) + data CSV + stats_triple.txt
# =============================================================================

# ---- 0. settings: EDIT THESE LINES ------------------------------------------
if (!exists("res_dir"))   res_dir   <- "results/three_protein"   # output of recruitment_quant.py run on triple_sheet.csv
if (!exists("shape_dir")) shape_dir <- "results/shape"    # output of shape_analysis.py run on shape_sheet.csv
ctx_cols  <- c("ACAN + BCAN" = "#E0413A", "ACAN + NCAN" = "#B8358F", "ACAN + BCAN + NCAN" = "#D98C1F",
               "ACAN alone" = "#2E9E5B")
client_cols <- c(BCAN = "#E0413A", NCAN = "#B8358F")

# ---- 1. packages ------------------------------------------------------------
for (p in c("ggplot2", "dplyr", "readr", "tidyr", "officer")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr); library(tidyr) })
fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE)

summ <- read_csv(file.path(res_dir, "image_summary.csv"), show_col_types = FALSE)
cobj <- read_csv(file.path(res_dir, "objects_client_puncta.csv"), show_col_types = FALSE)
annot <- function(d) d %>% mutate(
  client  = ifelse(grepl("BCAN", condition) & !grepl("NCAN$", condition), "BCAN", "NCAN"),
  client  = case_when(condition == "ACAN+BCAN" ~ "BCAN", condition == "ACAN+NCAN" ~ "NCAN",
                      condition == "Triple|BCAN" ~ "BCAN", condition == "Triple|NCAN" ~ "NCAN"),
  context = case_when(condition == "ACAN+BCAN" ~ "ACAN + BCAN", condition == "ACAN+NCAN" ~ "ACAN + NCAN",
                      TRUE ~ "ACAN + BCAN + NCAN"),
  mix     = factor(ifelse(grepl("^Triple", condition), "Three proteins", "Two proteins"),
                   levels = c("Two proteins", "Three proteins")))
summ <- annot(summ); cobj <- annot(cobj)
# statistical unit: independent preparations if every group has >= 3 of them, otherwise fields of view
n_rep <- summ %>% group_by(condition) %>% summarise(n = n_distinct(replicate), .groups = "drop")
unit <- if (all(n_rep$n >= 3)) "replicate" else "fov"
dot <- if (unit == "replicate") "Each dot is one independent preparation (mean of its FOVs)" else "Each dot is one field of view (FOV)"
xlab_fun <- function(x) ifelse(x == "Two proteins", "Without\nthe other\nclient", "With\nthe other\nclient")

# ---- 2. helpers ----------------------------------------------------------------
fmt_p <- function(p) ifelse(p < 1e-4, "P < 0.0001", ifelse(p > 0.99, "P > 0.99", paste0("P = ", formatC(signif(p, 2), format = "fg"))))
source(file.path(if (exists("scripts")) scripts else "scripts", "style_prism.R"), local = TRUE)  # shared Prism-like style
ctx_cols <- fill_cols; client_cols <- fill_cols
save_both <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".svg")), p, width = w, height = h, units = "mm", device = svg)
}
mwu <- function(a, b) if (length(unique(c(a, b))) == 1) 1 else if (unit == "replicate") student_p(a, b) else
  suppressWarnings(wilcox.test(a, b, exact = TRUE))$p.value
test_name <- if (unit == "replicate") TEST_T else "two-sided Mann-Whitney U test"
# two-protein vs three-protein, one panel per client
two_vs_three <- function(var, ylab, ytop = 108, ylim = c(0, 128), breaks = seq(0, 100, 25)) {
  d <- summ %>% transmute(client, mix, context, replicate, fov, y = .data[[var]])
  if (unit == "replicate") d <- d %>% group_by(client, mix, context, replicate) %>% summarise(y = mean(y), .groups = "drop") %>% mutate(fov = replicate)
  st <- d %>% group_by(client) %>% summarise(
    p = mwu(y[mix == "Two proteins"], y[mix == "Three proteins"]),
    n2 = sum(mix == "Two proteins"), n3 = sum(mix == "Three proteins"), .groups = "drop") %>%
    mutate(p_adj = p, lab = fmt_p(p_adj), y = ytop)
  g <- ggplot(d, aes(mix, y, fill = context)) + prism_layers() +
    geom_segment(data = st, aes(x = 1.05, xend = 1.95, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = st, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE, vjust = -0.45, size = 2.3) +
    facet_wrap(~client, labeller = labeller(client = c(BCAN = "BCAN", NCAN = "NCAN"))) +
    scale_fill_manual(values = fill_cols) +
    scale_x_discrete(labels = c("Two proteins" = "2 proteins", "Three proteins" = "3 proteins")) +
    scale_y_continuous(limits = ylim, breaks = breaks, expand = c(0, 0)) +
    labs(x = NULL, y = ylab) + theme_pub
  list(plot = g, stats = st, data = d)
}

# ---- 3. Fig E1, E2 ---------------------------------------------------------------
E1 <- two_vs_three("pct_client_inside_scaffold_B3", "Puncta within\nACAN condensates (%)")
E2 <- two_vs_three("pct_scaffold_clientpos_B1",     "BCAN/NCAN-positive\nACAN condensates (%)")
write_csv(E1$data, file.path(fig_dir, "FigE1_data_per_fov.csv"))
write_csv(E2$data, file.path(fig_dir, "FigE2_data_per_fov.csv"))
save_both(E1$plot, "FigE1_client_inside_ACAN_2vs3", 62, 58)
# single-client panels (for figure assembly): NCAN (FigE1n) and BCAN (FigE1b)
E1single <- list()
for (cl in c("NCAN", "BCAN")) local({
  d <- filter(E1$data, client == cl); st <- filter(E1$stats, client == cl)
  g <- ggplot(d, aes(mix, y, fill = context)) + prism_layers() +
    annotate("segment", x = 1, xend = 2, y = st$y, yend = st$y, linewidth = 0.3) +
    annotate("text", x = 1.5, y = st$y, label = st$lab, vjust = -0.45) +
    scale_fill_manual(values = fill_cols) +
    scale_x_discrete(labels = c("Two proteins" = paste("ACAN +", cl), "Three proteins" = "ACAN + BCAN + NCAN")) +
    scale_y_continuous(limits = c(0, 128), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    labs(x = NULL, y = paste0(cl, " puncta within\nACAN condensates (%)")) + theme_pub
  save_both(g, paste0("FigE1", tolower(substr(cl, 1, 1)), "_", cl, "_within_condensates_2vs3"), 40, 64)
  E1single[[cl]] <<- g
})
save_both(E2$plot, "FigE2_ACAN_containing_client_2vs3", 62, 58)

# ---- 4. Fig E3: location categories ------------------------------------------
locE <- cobj %>% count(client, mix, context, fov, location) %>%
  complete(nesting(client, mix, context, fov), location, fill = list(n = 0)) %>%
  group_by(client, mix, fov) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
locE <- locE %>% left_join(distinct(cobj, fov, replicate), by = "fov")
unitE <- if (unit == "replicate") locE %>% group_by(client, mix, replicate, location) %>% summarise(pct = mean(pct), .groups = "drop") else locE
meanE <- unitE %>% group_by(client, mix, location) %>% summarise(sd = sd(pct), n = n(), pct = mean(pct), .groups = "drop") %>%
  mutate(location = factor(recode(location, inside = "Inside", contacting = "At edge", separate = "Beyond edge"),
                           levels = c("Beyond edge", "At edge", "Inside")))
write_csv(locE,  file.path(fig_dir, "FigE3_data_per_fov.csv"))
write_csv(meanE, file.path(fig_dir, "FigE3_data_means.csv"))
pE3 <- ggplot(meanE, aes(mix, pct, fill = location)) +
  geom_col(width = 0.6, colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = ifelse(pct >= 6, sprintf("%.0f%%", pct), ""), colour = location == "Inside"),
            position = position_stack(vjust = 0.5), size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_manual(values = loc_cols,
                    breaks = c("Inside", "At edge", "Beyond edge")) +
  facet_wrap(~client) +
  scale_x_discrete(labels = c("Two proteins" = "2 proteins", "Three proteins" = "3 proteins")) +
  scale_y_continuous(limits = c(0, 100.5), expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of BCAN/NCAN\npuncta (%)", fill = "Position relative to\nACAN condensate") +
  theme_pub + theme(legend.position = "right", legend.key.size = unit(3, "mm"),
                    legend.title = element_text(size = 7.5), legend.spacing.y = unit(0.8, "mm")) +
  guides(fill = guide_legend(byrow = TRUE))   # byrow + spacing.y = small gap between keys (ggplot2 >= 3.3)
save_both(pE3, "FigE3_client_location_2vs3", 78, 58)

# ---- 5. Fig E4: BCAN vs NCAN inside the same three-protein images (paired) ------
tri <- summ %>% filter(mix == "Three proteins") %>%
  select(replicate, fov, client, acan_level = B4_median_rel_scaffold_at_client, inside = pct_client_inside_scaffold_B3)
if (unit == "replicate") tri <- tri %>% group_by(replicate, client) %>% summarise(across(c(acan_level, inside), mean), .groups = "drop") %>% mutate(fov = replicate)
tri <- tri %>% select(fov, client, acan_level, inside)
pairs <- tri %>% pivot_longer(c(acan_level, inside), names_to = "metric", values_to = "y")
st4 <- do.call(rbind, lapply(c("acan_level", "inside"), function(m) {
  w <- tri %>% select(fov, client, y = all_of(m)) %>% pivot_wider(names_from = client, values_from = y)
  data.frame(metric = m, n = nrow(w), p = t.test(w$BCAN, w$NCAN, paired = TRUE)$p.value,
             ytop = max(c(w$BCAN, w$NCAN)) * 1.08)
}))
st4$p_adj <- st4$p; st4$lab <- fmt_p(st4$p_adj)
write_csv(tri, file.path(fig_dir, "FigE4_data_per_fov_paired.csv"))
paired_panel <- function(m, ylab, ylim) {
  d <- filter(pairs, metric == m); s <- filter(st4, metric == m)
  ggplot(d, aes(client, y, fill = client)) +
    stat_summary(fun = mean, geom = "col", width = 0.7, colour = NA) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.28, linewidth = 0.4, colour = "black") +
    geom_line(aes(group = fov), colour = "grey35", linewidth = 0.3) +
    geom_point(shape = 21, fill = "white", colour = "black", size = 1.5, stroke = 0.35) +
    annotate("segment", x = 1.05, xend = 1.95, y = s$ytop, yend = s$ytop, linewidth = 0.3) +
    annotate("text", x = 1.5, y = s$ytop, label = s$lab, vjust = -0.45, size = 2.3) +
    scale_fill_manual(values = fill_cols) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    labs(x = "Three-protein mix", y = ylab) + theme_pub + theme_flat_x
}
pE4a <- paired_panel("acan_level", "Normalized ACAN\nintensity at puncta", c(0, 1.45)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey40")
pE4b <- paired_panel("inside", "Puncta within\nACAN condensates (%)", c(0, 115))
save_both(pE4a, "FigE4a_ACAN_level_BCAN_vs_NCAN_paired", 48, 58)
save_both(pE4b, "FigE4b_inside_BCAN_vs_NCAN_paired", 48, 58)

# ---- 6. Fig E5: condensate size and elongation ------------------------------------
sfov <- read_csv(file.path(shape_dir, "shape_fov_summary.csv"), show_col_types = FALSE) %>%
  mutate(context = recode(condition, "ACAN alone" = "ACAN alone", "ACAN+BCAN" = "ACAN + BCAN",
                          "ACAN+NCAN" = "ACAN + NCAN", "ACAN+BCAN+NCAN" = "ACAN + BCAN + NCAN"),
         context = factor(context, levels = c("ACAN alone", "ACAN + BCAN", "ACAN + NCAN", "ACAN + BCAN + NCAN")))
planned <- list(c("ACAN + BCAN", "ACAN + BCAN + NCAN"), c("ACAN + NCAN", "ACAN + BCAN + NCAN"))
shape_panel <- function(var, ylab, log = FALSE) {
  d <- sfov %>% transmute(context, replicate, fov, y = .data[[var]])
  if (unit == "replicate") d <- d %>% group_by(context, replicate) %>% summarise(y = mean(y), .groups = "drop")
  tk <- anova_tukey(d$y, d$context, log)                       # all pairs; the two comparisons with the triple mix are drawn
  st <- do.call(rbind, lapply(planned, function(g) {
    r <- tk[(tk$group1 == g[1] & tk$group2 == g[2]) | (tk$group1 == g[2] & tk$group2 == g[1]), ]
    data.frame(g1 = g[1], g2 = g[2], p = r$p_adj)
  }))
  st$p_adj <- st$p; st$lab <- fmt_p(st$p_adj); attr(st, "all_pairs") <- tk
  lv <- levels(d$context); st$x1 <- match(st$g1, lv); st$x2 <- match(st$g2, lv)
  st <- st[order(st$x2 - st$x1), ]; st$y <- max(d$y) * c(1.08, 1.2)[seq_len(nrow(st))]
  g <- ggplot(d, aes(context, y, fill = context)) + prism_layers() +
    geom_segment(data = st, aes(x = x1 + 0.05, xend = x2 - 0.05, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = st, aes(x = (x1 + x2) / 2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.45, size = 2.2) +
    scale_fill_manual(values = fill_cols) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = ylab) + theme_pub
  list(plot = g, stats = st, data = d)
}
sfov$median_area_plot <- if (is.na(PIXEL_UM)) sfov$median_area_px else sfov$median_area_px * PIXEL_UM^2
E5a <- shape_panel("median_area_plot", paste0("Condensate area (", area_unit_lab, ")"), log = FALSE)   # untransformed, as in the paper
E5b <- shape_panel("median_aspect_ratio", "Aspect ratio")
E5b$plot <- E5b$plot + geom_hline(yintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey40")
write_csv(sfov %>% select(context, fov, n_condensates, median_area_px, median_area_plot, median_aspect_ratio),
          file.path(fig_dir, "FigE5_data_per_fov.csv"))
save_both(E5a$plot, "FigE5a_size_with_triple", 66, 64)
save_both(E5b$plot, "FigE5b_elongation_with_triple", 66, 64)

# ---- 7. legends ----------------------------------------------------------------
n2 <- sum(summ$mix == "Two proteins" & summ$client == "NCAN"); n3 <- sum(summ$mix == "Three proteins" & summ$client == "NCAN")
nfov <- if (unit == "replicate") { if (length(unique(n_rep$n)) == 1) sprintf("n = %d independent preparations per group", n_rep$n[1]) else
  paste0("n = ", paste(sprintf("%s, %d", n_rep$condition, n_rep$n), collapse = "; "), " independent preparations") } else
  sprintf("n = %d FOVs (two-protein mixes) and %d FOVs (three-protein mix) [single preparation per condition; FOVs are technical replicates]", n2, n3)
nfov <- gsub("Triple\\|", "three-protein mix, ", nfov)
npunct <- cobj %>% count(client, mix)
npt <- paste(sprintf("%s, %s: %d puncta", npunct$client, tolower(npunct$mix), npunct$n), collapse = "; ")
common <- paste0("All mixes contained the same protein concentrations and were imaged with identical settings ",
  "([concentrations, 6% PEG, buffer, time after mixing]; ACAN, [fluorophore]; BCAN, [fluorophore]; NCAN, [fluorophore]). ",
  "Image analysis is described in Methods.")
legE1 <- paste0("Recruitment of each client into ACAN condensates with or without the other client. Percentage of BCAN (left) or NCAN (right) puncta ",
  "located inside ACAN condensates (at least half of the punctum area overlapping a condensate), in the two-protein mix ",
  "(ACAN + BCAN or ACAN + NCAN) and in the three-protein mix (ACAN + BCAN + NCAN). ", dot, "; bars, mean \u00b1 SEM. ",
  nfov, " (", npt, "). Statistics: ", test_name, ". ", common)
legE2 <- paste0("BCAN/NCAN-positive ACAN condensates: percentage of ACAN condensates in which BCAN (left) or NCAN (right) is detected (mean client signal inside the condensate above background mean + 3 SD), ",
  "in the two-protein and three-protein mixes. ", dot, "; bars, mean \u00b1 SEM. ", nfov,
  ". Statistics: ", test_name, ". ", common)
legE3 <- paste0("Position of BCAN and NCAN puncta relative to the edge of the ACAN condensate body (pre-specified Otsu threshold) in the two-protein and three-protein mixes. Inside, at least half of the ",
  "punctum overlaps a condensate; at edge, less than half overlaps but the punctum contacts the condensate edge; beyond edge, no contact with the bright condensate body ",
  "(such puncta still lie in dim ACAN signal around condensates; see Fig. S5). ",
  "Bars, mean percentage across ", ifelse(unit == "replicate", "independent preparations", "FOVs"), "; numbers inside bars give the percentages. ", nfov, ". ", common)
legE4 <- paste0("BCAN and NCAN compared within the same three-protein images. (a) ACAN level at client puncta: ACAN intensity at each punctum divided by the ",
  "median ACAN intensity of the ACAN condensates in the same image (1, dashed line = as strong as the median condensate of that image). (b) Percentage of client puncta inside ",
  "ACAN condensates. Each pair of dots joined by a line is one ", ifelse(unit == "replicate", "independent preparation", "FOV"), "; bars, mean \u00b1 SEM. n = ", nrow(tri) / 2, ifelse(unit == "replicate", " preparations. ", " FOVs. "),
  "Statistics: ", TEST_PAIR, ". ", if (TRI_SAT) "ACAN signal was partly saturated in these images, which compresses the measure in (a) towards 1; (a) is therefore used only for this within-image comparison. " else "",
  common)
legE5 <- paste0("Area (a; median per FOV, ", area_unit_lab, ") and aspect ratio (b; condensates of at least ", area_txt(50), ") of ACAN condensates in ACAN alone, the two-protein ",
  "mixes and the three-protein mix. ", dot, "; bars, mean \u00b1 SEM. Two planned comparisons are shown (three-protein mix vs ACAN + BCAN and vs ACAN + NCAN); ",
  "the comparisons among the first three groups are in Fig. D. Statistics: ", TEST_ANOVA, " across all four groups; all pairwise P values are in stats_triple.txt. ",
  "n = ", paste(table(E5a$data$context)[table(E5a$data$context) > 0], collapse = ", "), ifelse(unit == "replicate", " independent preparations", " FOVs"), " per group. ", common)

legE1one <- function(cl) {
  np <- npunct[npunct$client == cl, ]
  st <- E1$stats[E1$stats$client == cl, ]
  paste0("Recruitment of ", cl, " into ACAN condensates without and with ", setdiff(c("BCAN", "NCAN"), cl), ". Percentage of ", cl, " puncta ",
    "located inside ACAN condensates (at least half of the punctum area overlapping a condensate) in the two-protein mix (ACAN + ", cl,
    ") and in the three-protein mix (ACAN + BCAN + NCAN). ", dot, "; bars, mean \u00b1 SEM. ", nfov, " (",
    paste(sprintf("%s: %d puncta", tolower(np$mix), np$n), collapse = "; "), "). Statistics: ", test_name, ", ", st$lab, ". ", common)
}
panels <- list(
  list(file = "FigE1_client_inside_ACAN_2vs3",       title = "Figure E1", text = legE1, plot = E1$plot, h = 58),
  list(file = "FigE1b_BCAN_within_condensates_2vs3", title = "Figure E1 (BCAN panel alone)", text = legE1one("BCAN"), plot = E1single$BCAN, h = 64),
  list(file = "FigE1n_NCAN_within_condensates_2vs3", title = "Figure E1 (NCAN panel alone)", text = legE1one("NCAN"), plot = E1single$NCAN, h = 64),
  list(file = "FigE2_ACAN_containing_client_2vs3",   title = "Figure E2", text = legE2, plot = E2$plot, h = 58),
  list(file = "FigE3_client_location_2vs3",          title = "Figure E3", text = legE3, plot = pE3,     h = 58),
  list(file = "FigE4a_ACAN_level_BCAN_vs_NCAN_paired", title = "Figure E4 (a, b)", text = legE4, plot = pE4a, h = 58),
  list(file = "FigE5a_size_with_triple",             title = "Figure E5 (a, b)", text = legE5, plot = E5a$plot, h = 62))
doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "Figure legends - three-protein mix (auto-generated by plot_triple.R)", style = "heading 1")
doc <- officer::body_add_par(doc, "Items in [brackets] must be completed by hand.", style = "Normal")
all_txt <- character()
for (L in panels) {
  writeLines(c(paste0(L$title, " (", L$file, ")"), "", L$text), file.path(fig_dir, paste0(L$file, "_legend.txt")), useBytes = TRUE)
  all_txt <- c(all_txt, paste0(L$title, " (", L$file, ")"), L$text, "")
  doc <- officer::body_add_par(doc, paste0(L$title, " - ", L$file), style = "heading 2")
  doc <- officer::body_add_par(doc, L$text, style = "Normal")
  cap <- paste(strwrap(L$text, width = 95), collapse = "\n")
  pc <- L$plot + labs(caption = cap) + theme(plot.caption = element_text(hjust = 0, size = 6.5, lineheight = 1.1),
                                             plot.caption.position = "plot")
  if (grepl("E3", L$file)) pc <- pc + theme(legend.position = "right")
  ggsave(file.path(fig_dir, paste0(L$file, "_with_legend.pdf")), pc, width = 180,
         height = L$h + 12 + 3.2 * length(strsplit(cap, "\n")[[1]]), units = "mm", device = cairo_pdf)
}
writeLines(all_txt, file.path(fig_dir, "Figure_legends_triple.txt"), useBytes = TRUE)
print(doc, target = file.path(fig_dir, "Figure_legends_triple.docx"))

# ---- 8. statistics ---------------------------------------------------------------
sink(file.path(fig_dir, "stats_triple.txt"))
cat("Unit of analysis:", ifelse(unit == "replicate", "independent preparation", "field of view (technical replicates)"), "| two-group test:", test_name, "\n")
cat("\n== E1 % client puncta inside ACAN (two vs three proteins;", test_name, ") ==\n"); print(as.data.frame(E1$stats))
cat("\n== E2 % ACAN condensates containing client ==\n"); print(as.data.frame(E2$stats))
cat("\n== E4 BCAN vs NCAN within three-protein images (", TEST_PAIR, ") ==\n"); print(as.data.frame(st4))
cat("\n== E5a size:", TEST_ANOVA, "; ANOVA P =", signif(attr(attr(E5a$stats, "all_pairs"), "anova_p"), 3), "==\n"); print(attr(E5a$stats, "all_pairs"))
cat("\n== E5b aspect ratio:", TEST_ANOVA, "; ANOVA P =", signif(attr(attr(E5b$stats, "all_pairs"), "anova_p"), 3), "==\n"); print(attr(E5b$stats, "all_pairs"))
cat("\nMeans per group:\n")
print(as.data.frame(summ %>% group_by(client, mix) %>% summarise(n_fov = n(),
  inside = mean(pct_client_inside_scaffold_B3), contains = mean(pct_scaffold_clientpos_B1),
  acan_level = mean(B4_median_rel_scaffold_at_client), .groups = "drop")))
sink()
cat("Triple-mix figures, legends, data and stats written to", normalizePath(fig_dir), "\n")

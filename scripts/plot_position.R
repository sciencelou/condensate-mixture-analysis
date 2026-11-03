# =============================================================================
# plot_position.R
# Called automatically by make_figures.R (or open in RStudio and click "Source").
#
# WHERE in the ACAN condensate do BCAN and NCAN sit - core or rim?
#   Fig F1  average intensity profile from condensate centre (0) to edge (1) and just outside
#   Fig F2  client position relative to ACAN (radial ratio), two-protein mixes, BCAN vs NCAN
#   Fig F3  the same in the three-protein mix: BCAN vs NCAN in the SAME condensates (paired)
#   Fig S5  % client puncta inside ACAN: bright body only (Otsu) vs body + dim rim (Li)
#   Fig S6  ACAN level at each client punctum relative to its own condensate (B5), Otsu and Li
#   Fig S7  robustness of the radial ratio (other threshold, no channel-offset correction,
#           simulated client distributed exactly like ACAN)
# Each dot = one independent preparation (mean of its FOVs); bar = mean +/- SEM.
# Output: vector PDF + SVG + legends (.txt/.docx) + data CSV + stats_position.txt
# =============================================================================
if (!exists("pos_dir")) pos_dir <- "results/position"
if (!exists("rec_dir")) rec_dir <- "results/recruitment"
if (!exists("tri_dir")) tri_dir <- "results/three_protein"
for (p in c("ggplot2", "dplyr", "readr", "tidyr", "officer")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr); library(tidyr) })
fig_dir <- file.path(pos_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(if (exists("scripts")) scripts else "scripts", "style_prism.R"), local = TRUE)  # shared Prism-like style
prot_cols <- line_cols
fmt_p <- function(p) ifelse(p < 1e-4, "P < 0.0001", ifelse(p > 0.99, "P > 0.99", paste0("P = ", formatC(signif(p, 2), format = "fg"))))
save_both <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".svg")), p, width = w, height = h, units = "mm", device = svg)
}
bar_layers <- prism_layers()
shift_lab <- "Peripheral shift\nrelative to ACAN (%)"   # (radial ratio - 1) x 100; 0 = distributed like ACAN
paired_layers <- list(
  stat_summary(fun = mean, geom = "col", width = 0.7, colour = NA),
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.28, linewidth = 0.4, colour = "black"),
  geom_line(aes(group = replicate), colour = "grey35", linewidth = 0.3),
  geom_point(shape = 21, fill = "white", colour = "black", size = 1.5, stroke = 0.35))
rd <- function(...) read_csv(file.path(...), show_col_types = FALSE)
client_of <- function(cond) ifelse(grepl("NCAN$", cond), "NCAN", "BCAN")
mix_of <- function(cond) ifelse(grepl("^Triple", cond), "ACAN + BCAN + NCAN", ifelse(cond == "ACAN+BCAN", "ACAN + BCAN", "ACAN + NCAN"))
welch <- function(a, b) student_p(a, b)   # name kept for brevity: two-tailed unpaired Student t-test

# ---- 1. data ---------------------------------------------------------------------
two  <- rd(pos_dir, "two_protein", "radial_fov_summary.csv")
tri  <- rd(pos_dir, "three_protein", "radial_fov_summary.csv") %>% filter(grepl("^Triple", condition))
fovR <- bind_rows(two, tri) %>% mutate(client = client_of(condition), mix = mix_of(condition))
n_rep <- fovR %>% group_by(condition) %>% summarise(n = n_distinct(replicate), .groups = "drop")
unit <- if (all(n_rep$n >= 3)) "replicate" else "fov"
unit_word <- if (unit == "replicate") "independent preparation" else "field of view (FOV)"
dot <- if (unit == "replicate") "Each dot is one independent preparation (mean of its FOVs)" else "Each dot is one field of view (FOV)"
nprep <- if (unit == "replicate" && length(unique(n_rep$n)) == 1) sprintf("n = %d independent preparations per group", n_rep$n[1]) else
  paste0("n = ", paste(sprintf("%s, %d", n_rep$condition, n_rep$n), collapse = "; "), ifelse(unit == "replicate", " independent preparations", " FOVs"))
to_unit <- function(d, ...) if (unit == "replicate") d %>% group_by(..., replicate) %>% summarise(across(where(is.numeric), mean), .groups = "drop") else d
ncond <- fovR %>% group_by(condition) %>% summarise(n = sum(n_clientpos), .groups = "drop")
ncond_txt <- paste(sprintf("%s, %d", recode(ncond$condition, "Triple|BCAN" = "three-protein mix", "Triple|NCAN" = "three-protein mix (same condensates)"), ncond$n), collapse = "; ")

# ---- 2. Fig F1: radial profiles -----------------------------------------------------
prof_raw <- bind_rows(rd(pos_dir, "two_protein", "radial_profiles_binned.csv"),
                      rd(pos_dir, "three_protein", "radial_profiles_binned.csv") %>% filter(grepl("^Triple", condition))) %>%
  filter(client_pos_B1)
prof_img <- prof_raw %>% group_by(condition, replicate, fov, rho) %>%
  summarise(ACAN = mean(scaffold_norm, na.rm = TRUE), client = mean(client_norm, na.rm = TRUE), .groups = "drop")
prof_u <- if (unit == "replicate") prof_img %>% group_by(condition, replicate, rho) %>% summarise(ACAN = mean(ACAN, na.rm = TRUE), client = mean(client, na.rm = TRUE), .groups = "drop") else prof_img
prof_long <- prof_u %>% mutate(client_name = client_of(condition), mix = mix_of(condition)) %>%
  pivot_longer(c(ACAN, client), names_to = "chan", values_to = "y") %>%
  mutate(protein = ifelse(chan == "ACAN", "ACAN", client_name)) %>%
  filter(!(mix == "ACAN + BCAN + NCAN" & chan == "ACAN" & client_name == "NCAN"))   # ACAN drawn once in the triple panel
prof_mean <- prof_long %>% group_by(mix, protein, rho) %>%
  summarise(m = mean(y), se = sd(y) / sqrt(n()), n = n(), .groups = "drop") %>%
  mutate(mix = factor(mix, levels = c("ACAN + BCAN", "ACAN + NCAN", "ACAN + BCAN + NCAN")),
         protein = factor(protein, levels = c("ACAN", "BCAN", "NCAN")))
write_csv(prof_u, file.path(fig_dir, paste0("FigF1_data_per_", unit, ".csv")))
write_csv(prof_mean, file.path(fig_dir, "FigF1_data_means.csv"))
pF1 <- ggplot(prof_mean, aes(rho, m, colour = protein, fill = protein)) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey40") +
  geom_ribbon(aes(ymin = m - se, ymax = m + se), colour = NA, alpha = 0.25) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~mix, nrow = 1) +
  scale_colour_manual(values = prot_cols) + scale_fill_manual(values = prot_cols) +
  annotate("text", x = 1.25, y = Inf, label = "outside", vjust = 1.5, size = 2.3, colour = "grey30") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("centre", "0.5", "edge"), limits = c(0, 1.5)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Normalized radial position", y = "Normalized intensity", colour = NULL, fill = NULL) +
  theme_pub + theme_flat_x + theme(legend.position = "top", legend.key.size = unit(3, "mm"))
save_both(pF1, "FigF1_radial_profiles", 120, 62)

# ---- 3. Fig F2 / F3: radial ratio ---------------------------------------------------
rr_two <- to_unit(fovR %>% filter(mix != "ACAN + BCAN + NCAN") %>% select(condition, client, mix, replicate, fov, y = radial_ratio_median), condition, client, mix)
pF2_p <- welch(rr_two$y[rr_two$client == "BCAN"], rr_two$y[rr_two$client == "NCAN"])
rr_two$shift_pct <- (rr_two$y - 1) * 100
write_csv(rr_two %>% rename(radial_ratio = y), file.path(fig_dir, paste0("FigF2_data_per_", unit, ".csv")))
ytop2 <- max(rr_two$shift_pct) + 2
pF2 <- ggplot(rr_two, aes(mix, shift_pct, fill = client)) + bar_layers + zero_axis() +
  annotate("segment", x = 1, xend = 2, y = ytop2, yend = ytop2, linewidth = 0.3) +
  annotate("text", x = 1.5, y = ytop2, label = fmt_p(pF2_p), vjust = -0.45) +
  scale_fill_manual(values = fill_cols) +
  scale_x_discrete(labels = c("ACAN + BCAN" = "ACAN + BCAN", "ACAN + NCAN" = "ACAN + NCAN")) +
  scale_y_continuous(limits = c(min(0, floor(min(rr_two$shift_pct)) - 1), ceiling(ytop2) + 2), breaks = scales::breaks_pretty(6), expand = c(0, 0)) +
  labs(x = NULL, y = shift_lab) + theme_pub
save_both(pF2, "FigF2_radial_ratio_two_protein", 46, 62)

rr_tri <- to_unit(fovR %>% filter(mix == "ACAN + BCAN + NCAN") %>% select(client, replicate, fov, y = radial_ratio_median), client)
if (unit == "fov") rr_tri <- rr_tri %>% select(client, replicate = fov, y)
w3 <- rr_tri %>% select(replicate, client, y) %>% pivot_wider(names_from = client, values_from = y)
pF3_p <- t.test(w3$NCAN, w3$BCAN, paired = TRUE)$p.value
write_csv(w3 %>% mutate(BCAN_shift_pct = (BCAN - 1) * 100, NCAN_shift_pct = (NCAN - 1) * 100) %>%
            rename(BCAN_radial_ratio = BCAN, NCAN_radial_ratio = NCAN),
          file.path(fig_dir, paste0("FigF3_data_per_", unit, "_paired.csv")))
rr_tri$shift_pct <- (rr_tri$y - 1) * 100
ytop3 <- max(rr_tri$shift_pct) + 2
pF3 <- ggplot(rr_tri, aes(client, shift_pct, fill = client)) + paired_layers + zero_axis() +
  annotate("segment", x = 1, xend = 2, y = ytop3, yend = ytop3, linewidth = 0.3) +
  annotate("text", x = 1.5, y = ytop3, label = fmt_p(pF3_p), vjust = -0.45) +
  scale_fill_manual(values = fill_cols) +
  scale_y_continuous(limits = c(min(0, floor(min(rr_tri$shift_pct)) - 1), ceiling(ytop3) + 2), breaks = scales::breaks_pretty(6), expand = c(0, 0)) +
  labs(x = "ACAN + BCAN + NCAN\n(same condensates)", y = shift_lab) + theme_pub
save_both(pF3, "FigF3_radial_ratio_triple_paired", 52, 66)

# ---- 4. Fig S5: inside ACAN, bright body (Otsu) vs body + dim rim (Li) ---------------
bnd <- bind_rows(
  rd(rec_dir, "image_summary.csv") %>% mutate(boundary = "Bright condensate body\n(main threshold, Otsu)"),
  rd(rec_dir, "check_Li_threshold", "image_summary.csv") %>% mutate(boundary = "Body + dim ACAN rim\n(inclusive threshold, Li)")) %>%
  mutate(boundary = factor(boundary, levels = unique(boundary)), client = client_of(condition))
s5 <- to_unit(bnd %>% select(boundary, client, replicate, fov, y = pct_client_inside_scaffold_B3), boundary, client)
st5 <- s5 %>% group_by(boundary) %>% summarise(p = welch(y[client == "BCAN"], y[client == "NCAN"]),
  BCAN = mean(y[client == "BCAN"]), NCAN = mean(y[client == "NCAN"]), .groups = "drop") %>%
  mutate(p_adj = p, lab = fmt_p(p_adj))
write_csv(s5, file.path(fig_dir, paste0("FigS5_data_per_", unit, ".csv")))
pS5 <- ggplot(s5, aes(client, y, fill = client)) + bar_layers +
  geom_segment(data = st5, aes(x = 1.05, xend = 1.95, y = 108, yend = 108), inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = st5, aes(x = 1.5, y = 108, label = lab), inherit.aes = FALSE, vjust = -0.45, size = 2.3) +
  facet_wrap(~boundary) + scale_fill_manual(values = fill_cols) +
  scale_y_continuous(limits = c(0, 125), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = "Two-protein mixes", y = "Puncta within ACAN (%)") + theme_pub
save_both(pS5, "FigS5_inside_ACAN_body_vs_rim", 80, 62)

# ---- 5. Fig S6: B5, ACAN at punctum / ACAN of its own condensate ---------------------
b5 <- to_unit(bnd %>% filter(!is.na(B5_median_rel_scaffold_to_host)) %>%
                select(boundary, client, replicate, fov, y = B5_median_rel_scaffold_to_host), boundary, client)
st6 <- b5 %>% group_by(boundary) %>% summarise(p = welch(y[client == "BCAN"], y[client == "NCAN"]),
  BCAN = mean(y[client == "BCAN"]), NCAN = mean(y[client == "NCAN"]), .groups = "drop") %>%
  mutate(p_adj = p, lab = fmt_p(p_adj), yy = max(b5$y) * 1.08)
tri_b5 <- rd(tri_dir, "image_summary.csv") %>% filter(grepl("^Triple", condition)) %>%
  mutate(client = client_of(condition)) %>% select(client, replicate, fov, y = B5_median_rel_scaffold_to_host)
tri_b5 <- to_unit(tri_b5, client)
write_csv(b5, file.path(fig_dir, paste0("FigS6_data_per_", unit, ".csv")))
write_csv(tri_b5, file.path(fig_dir, paste0("FigS6_three_protein_values_per_", unit, ".csv")))
pS6 <- ggplot(b5, aes(client, y, fill = client)) + bar_layers +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey40") +
  geom_segment(data = st6, aes(x = 1.05, xend = 1.95, y = yy, yend = yy), inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = st6, aes(x = 1.5, y = yy, label = lab), inherit.aes = FALSE, vjust = -0.45, size = 2.3) +
  facet_wrap(~boundary) + scale_fill_manual(values = fill_cols) +
  scale_y_continuous(limits = c(0, max(b5$y) * 1.2), expand = c(0, 0)) +
  labs(x = "Two-protein mixes", y = "ACAN intensity at puncta\n(norm. to host condensate)") + theme_pub
save_both(pS6, "FigS6_ACAN_at_punctum_vs_own_condensate", 80, 62)

# ---- 6. Fig S7: robustness of the radial ratio ---------------------------------------
variant <- function(sub, label, col = "radial_ratio_median") {
  bind_rows(rd(pos_dir, "two_protein", sub, "radial_fov_summary.csv"),
            rd(pos_dir, "three_protein", sub, "radial_fov_summary.csv") %>% filter(grepl("^Triple", condition))) %>%
    transmute(condition, replicate, fov, y = .data[[col]], analysis = label)
}
s7 <- bind_rows(
  variant(".", "Main\nanalysis"),
  variant("check_Li_threshold", "Inclusive\nthreshold (Li)"),
  variant("check_no_shift_correction", "No offset\ncorrection"),
  variant(".", paste0("Size-matched\n(", area_range(100, 400), ")"), "radial_ratio_median_sizematched"),
  variant(".", "Simulated client\n(same as ACAN)", "radial_ratio_null_median")) %>%
  mutate(analysis = factor(analysis, levels = unique(analysis)), client = client_of(condition),
         group = factor(ifelse(grepl("^Triple", condition), paste0(client, ", 3 prot."), paste0(client, ", 2 prot.")),
                        levels = c("BCAN, 2 prot.", "NCAN, 2 prot.", "BCAN, 3 prot.", "NCAN, 3 prot.")))
s7u <- to_unit(s7, analysis, client, group)
s7u$shift_pct <- (s7u$y - 1) * 100
write_csv(s7u %>% rename(radial_ratio = y), file.path(fig_dir, paste0("FigS7_data_per_", unit, ".csv")))
st7 <- do.call(rbind, lapply(levels(s7u$analysis)[1:4], function(an) {
  d <- s7u[s7u$analysis == an, ]
  pick_g <- function(g) { x <- d[d$group == g, ]; x$y[order(x$replicate)] }
  data.frame(analysis = gsub("\n", " ", an),
             p_two_protein_student = welch(pick_g("BCAN, 2 prot."), pick_g("NCAN, 2 prot.")),
             p_three_protein_paired = t.test(pick_g("NCAN, 3 prot."), pick_g("BCAN, 3 prot."), paired = TRUE)$p.value)
}))
pS7 <- ggplot(s7u, aes(group, shift_pct, fill = client)) + prism_layers(point_size = 1.2, dx = 0.12) + zero_axis() +
  facet_wrap(~analysis, nrow = 1) + scale_fill_manual(values = fill_cols) +
  scale_y_continuous(limits = c(min(0, floor(min(s7u$shift_pct)) - 1), ceiling(max(s7u$shift_pct)) + 2), breaks = scales::breaks_pretty(6), expand = c(0, 0)) +
  labs(x = NULL, y = shift_lab) + theme_pub
save_both(pS7, "FigS7_radial_ratio_robustness", 180, 78)

# ---- 7. legends ----------------------------------------------------------------------
common <- paste0("Same protein concentrations and imaging settings in all mixes ([concentrations, 6% PEG, buffer, time after mixing]; ACAN, [fluorophore]; ",
  "BCAN, [fluorophore]; NCAN, [fluorophore]). Image analysis is described in Methods.")
method <- paste0("For each ACAN condensate the edge was defined from its own ACAN signal (the contour at half-maximum between the surrounding background and ",
  "the condensate's ACAN peak), so the measurement does not depend on a global intensity threshold. Each pixel was assigned a normalized radial position from the centre (0) ",
  "to this edge (1) and up to ", len_txt(3), " beyond it. Only condensates containing the client (client signal above background mean + 3 SD) were analysed.")
legF1 <- paste0("Distribution of BCAN and NCAN within ACAN condensates. ", method, " Curves show the intensity of each protein at each position, ",
  "normalised to that protein's mean inside the condensate, so that the curves compare shape, not brightness. Lines, mean; shading, \u00b1 SEM across ",
  ifelse(unit == "replicate", "independent preparations", "FOVs"), ". Grey area, outside the condensate edge. In ACAN + BCAN the BCAN profile follows ACAN; ",
  "in ACAN + NCAN, NCAN is relatively depleted from the centre and relatively enriched at and beyond the edge. ",
  if (TRI_SAT) "ACAN was partly saturated in the three-protein images, which flattens the ACAN curve in that panel. " else "", nprep, " (condensates analysed: ", ncond_txt, "). ", common)
legF2 <- paste0("Position of the client relative to ACAN in the two-protein mixes. For every condensate the signal-weighted mean position of the client was divided by that of ACAN ",
  "(positions as in Fig. F1; radial ratio) and expressed as the peripheral shift, (radial ratio - 1) x 100%: 0, the client is distributed like ACAN; ",
  "> 0, shifted towards the rim; < 0, towards the centre. ",
  ifelse(unit == "replicate", "Each dot is one independent preparation (median of the condensates in each FOV, averaged over the FOVs of that preparation)", "Each dot is one FOV (median of its condensates)"),
  "; bars, mean \u00b1 SEM. ", nprep, ". Statistics: ", TEST_T, " on the radial ratio. ", common)
legF3 <- paste0("Position of BCAN and NCAN relative to ACAN in the three-protein mix, measured in the same condensates (peripheral shift as in Fig. F2). ",
  "Lines join BCAN and NCAN from the same ", ifelse(unit == "replicate", "preparation", "FOV"), "; bars, mean \u00b1 SEM. ", nprep,
  ". Statistics: ", TEST_PAIR, ". ", if (TRI_SAT) "Because both clients are compared within the same condensates and against the same ACAN signal, partial ACAN saturation affects both equally. " else "",
  common)
legS5 <- paste0("Percentage of BCAN or NCAN puncta lying inside ACAN (at least half of the punctum area), with the ACAN boundary drawn either around the bright condensate body ",
  "(left; the pre-specified threshold used in Figs. B, C) or with a more inclusive threshold that also captures the dim ACAN rim around condensates (right). ",
  "With the inclusive boundary nearly all NCAN puncta lie within ACAN signal, showing that NCAN puncta counted as 'at edge' or 'beyond edge' in Fig. C sit in the dim ACAN periphery of condensates ",
  "rather than free in solution. ", dot, "; bars, mean \u00b1 SEM. ", nprep, ". Statistics: ", TEST_T, " for each boundary. ", common)
legS6 <- paste0("ACAN intensity at each client punctum divided by the mean ACAN intensity of the condensate that the punctum overlaps (1, dashed line = as strong as the average of its own ",
  "condensate; < 1 = in the dimmer, outer part). Left, condensates outlined with the pre-specified threshold (Otsu); right, with the inclusive threshold (Li). ",
  "With the inclusive threshold the condensate average includes the dim rim and is lowered, so the ratio rises above 1 and no longer separates BCAN from NCAN; ",
  "this measure therefore depends on how the condensate boundary is drawn and is shown only as supporting information (the threshold-independent measure is Fig. F). ",
  dot, "; bars, mean \u00b1 SEM. ", nprep, ". Statistics: ", TEST_T, " for each threshold. ", common)
legS7 <- paste0("Robustness of the peripheral shift (Figs. F2, F3). Main analysis; condensates identified with the inclusive threshold (Li) instead of Otsu; ",
  "without correcting the sub-pixel offset between channels; only condensates of ", area_range(100, 400), " (to compare condensates of the same size, because ACAN + NCAN condensates are smaller); and a simulated client distributed exactly like ACAN but with the measured client brightness and ",
  "background noise of each image (expected value if there were no real difference; it stays below 1%, far smaller than the shift measured for NCAN, so low signal-to-noise does not create the shift). ",
  "2 prot., two-protein mix; 3 prot., three-protein mix (same condensates). ", dot, "; bars, mean \u00b1 SEM. ", nprep, ". P values for each variant are in stats_position.txt. ", common)

panels <- list(
  list(file = "FigF1_radial_profiles",                    title = "Figure F1", text = legF1, plot = pF1, h = 62),
  list(file = "FigF2_radial_ratio_two_protein",           title = "Figure F2", text = legF2, plot = pF2, h = 62),
  list(file = "FigF3_radial_ratio_triple_paired",         title = "Figure F3", text = legF3, plot = pF3, h = 62),
  list(file = "FigS5_inside_ACAN_body_vs_rim",            title = "Figure S5", text = legS5, plot = pS5, h = 62),
  list(file = "FigS6_ACAN_at_punctum_vs_own_condensate",  title = "Figure S6", text = legS6, plot = pS6, h = 62),
  list(file = "FigS7_radial_ratio_robustness",            title = "Figure S7", text = legS7, plot = pS7, h = 78))
doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "Figure legends - position of BCAN and NCAN within ACAN condensates (auto-generated by plot_position.R)", style = "heading 1")
doc <- officer::body_add_par(doc, "Items in [brackets] must be completed by hand.", style = "Normal")
all_txt <- character()
for (L in panels) {
  writeLines(c(paste0(L$title, " (", L$file, ")"), "", L$text), file.path(fig_dir, paste0(L$file, "_legend.txt")), useBytes = TRUE)
  all_txt <- c(all_txt, paste0(L$title, " (", L$file, ")"), L$text, "")
  doc <- officer::body_add_par(doc, paste0(L$title, " - ", L$file), style = "heading 2")
  doc <- officer::body_add_par(doc, L$text, style = "Normal")
  cap <- paste(strwrap(L$text, width = 95), collapse = "\n")
  pc <- L$plot + labs(caption = cap) + theme(plot.caption = element_text(hjust = 0, size = 6.5, lineheight = 1.1), plot.caption.position = "plot")
  ggsave(file.path(fig_dir, paste0(L$file, "_with_legend.pdf")), pc, width = 180,
         height = L$h + 12 + 3.2 * length(strsplit(cap, "\n")[[1]]), units = "mm", device = cairo_pdf)
}
writeLines(all_txt, file.path(fig_dir, "Figure_legends_position.txt"), useBytes = TRUE)
print(doc, target = file.path(fig_dir, "Figure_legends_position.docx"))

# ---- 8. statistics -------------------------------------------------------------------
sink(file.path(fig_dir, "stats_position.txt"))
cat("Unit of analysis:", unit_word, "\n\n")
cat("== F2 radial ratio, two-protein mixes, BCAN vs NCAN (", TEST_T, ") ==\nP =", signif(pF2_p, 3), "\n")
print(as.data.frame(rr_two %>% group_by(client) %>% summarise(mean = mean(y), sd = sd(y), n = n())))
cat("\n== F3 radial ratio, three-protein mix, NCAN vs BCAN in the same condensates (", TEST_PAIR, ") ==\nP =", signif(pF3_p, 3), "\n")
print(as.data.frame(w3)); cat("mean difference NCAN - BCAN =", signif(mean(w3$NCAN - w3$BCAN), 3), "\n")
cat("\n== S5 % client puncta inside ACAN by boundary (Student t) ==\n"); print(as.data.frame(st5))
cat("\n== S6 B5 (ACAN at punctum / own condensate) by boundary (Student t) ==\n"); print(as.data.frame(st6))
cat("three-protein mix (Otsu), per", unit_word, ":\n"); print(as.data.frame(tri_b5))
cat("\n== S7 radial ratio robustness: P (two-protein Student t; three-protein paired t) ==\n"); print(as.data.frame(st7))
cat("\nmeans per variant:\n"); print(as.data.frame(s7u %>% group_by(analysis, group) %>% summarise(mean = mean(y), sd = sd(y), n = n(), .groups = "drop")))
sink()
cat("Position figures, legends, data and stats written to", normalizePath(fig_dir), "\n")

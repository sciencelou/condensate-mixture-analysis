# =============================================================================
# plot_recruitment.R
# Called automatically by make_figures.R (or open in RStudio and click "Source").
#
# Reads the output of recruitment_quant.py and makes simple, publication-style plots:
#   Fig A  ACAN level where each client punctum sits (or partition coefficient, if valid)
#   Fig B  % ACAN condensates containing client; % client puncta inside ACAN condensates
#   Fig C  where client puncta are: inside / touching / separate from ACAN condensates
#   Fig S1 the % in Fig B under stricter / looser criteria (robustness)
#   Fig S2 all individual client puncta behind Fig A
# Each dot = one field of view (or one independent preparation, once >= 3 exist);
# bar = mean +/- SEM; P value printed on the graph.
# Output: vector PDF + SVG + legend + the exact data behind every panel (CSV) + stats.txt
# =============================================================================

# ---- 0. settings: EDIT THIS LINE -------------------------------------------
if (!exists("res_dir")) res_dir <- "results/recruitment"  # folder created by recruitment_quant.py (full path is safest,
                              # e.g. "D:/LLPS/acan_bcan/results"; use / not \ )
cond_order  <- c("ACAN+BCAN", "ACAN+NCAN", "ACAN alone")   # x-axis order; absent ones are dropped
cond_labels <- c("ACAN+BCAN" = "ACAN + BCAN", "ACAN+NCAN" = "ACAN + NCAN", "ACAN alone" = "ACAN alone")
stat_unit   <- "auto"  # unit of the statistical test:
                       #  "replicate" = independent preparations (publication standard; needs >= 3 per condition)
                       #  "fov"       = fields of view (technical replicates; stated in the legend)
                       #  "auto"      = "replicate" if every condition has >= 3 preparations, otherwise "fov"

# ---- 1. packages ------------------------------------------------------------
for (p in c("ggplot2", "dplyr", "readr", "tidyr", "officer", "lme4", "lmerTest")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
})
source(file.path(if (exists("scripts")) scripts else "scripts", "style_prism.R"), local = TRUE)  # shared Prism-like style
cond_cols <- fill_cols

fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE)
obj  <- read_csv(file.path(res_dir, "objects_scaffold_condensates.csv"), show_col_types = FALSE)
cobj <- read_csv(file.path(res_dir, "objects_client_puncta.csv"), show_col_types = FALSE)
summ <- read_csv(file.path(res_dir, "image_summary.csv"), show_col_types = FALSE)

lev <- intersect(cond_order, unique(summ$condition))
lev <- c(lev, setdiff(unique(summ$condition), lev))
fix <- function(d) { d$condition <- factor(d$condition, levels = lev)
                     d$replicate <- factor(d$replicate); d }
obj <- fix(obj); cobj <- fix(cobj); summ <- fix(summ)

n_rep_cond <- summ %>% group_by(condition) %>% summarise(n = n_distinct(replicate), .groups = "drop")
unit <- if (stat_unit == "auto") {
  if (nrow(n_rep_cond) >= 2 && all(n_rep_cond$n >= 3)) "replicate" else "fov"
} else stat_unit
message("Statistical unit: ", unit)
# reduce FOV-level values to the plotted/tested unit
to_unit <- function(d, f = mean) {
  if (unit == "replicate") d %>% group_by(condition, replicate, across(any_of("metric"))) %>%
    summarise(y = f(y), .groups = "drop") else d
}

# two conditions: two-tailed unpaired Student t-test; > 2 conditions: one-way ANOVA + Tukey
# ("welch" is only the internal name of the parametric branch; FOV mode uses Mann-Whitney + Holm)
compare <- function(d, method = c("welch", "mwu")) {
  method <- match.arg(method)
  d <- droplevels(d); lv <- levels(d$condition)
  if (length(lv) < 2) return(NULL)
  if (length(lv) > 2 && method == "welch") {
    out <- anova_tukey(d$y, d$condition); out$test <- TEST_ANOVA; return(out)
  }
  pr <- combn(lv, 2, simplify = FALSE)
  out <- do.call(rbind, lapply(pr, function(g) {
    a <- d$y[d$condition == g[1]]; b <- d$y[d$condition == g[2]]
    pv <- if (method == "welch") student_p(a, b) else
      suppressWarnings(wilcox.test(a, b, exact = TRUE))$p.value
    data.frame(group1 = g[1], group2 = g[2], n1 = length(a), n2 = length(b),
               mean1 = mean(a), mean2 = mean(b), p = pv,
               test = ifelse(method == "welch", TEST_T, "two-sided Mann-Whitney U test"),
               p_welch_sensitivity = if (method == "welch") welch_p(a, b) else NA)
  }))
  out$p_adj <- if (nrow(out) > 1) p.adjust(out$p, "holm") else out$p
  out
}
fmt_p <- function(p) ifelse(p < 1e-4, "P < 0.0001", paste0("P = ", formatC(signif(p, 2), format = "fg")))
add_p <- function(plt, cmp, ytop, facet = NULL) {
  if (is.null(cmp) || nrow(cmp) == 0) return(plt)
  lv <- levels(droplevels(plt$data$condition))
  cmp$x1 <- match(cmp$group1, lv); cmp$x2 <- match(cmp$group2, lv)
  cmp$y <- ytop * (1 + 0.09 * (seq_len(nrow(cmp)) - 1)); cmp$lab <- fmt_p(cmp$p_adj)
  if (!is.null(facet)) cmp$metric <- facet
  plt + geom_segment(data = cmp, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE,
                     linewidth = 0.3) +
    geom_text(data = cmp, aes(x = (x1 + x2) / 2, y = y, label = lab), inherit.aes = FALSE,
              vjust = -0.5, size = 2.4) +
    expand_limits(y = max(cmp$y) * 1.1)
}


save_both <- function(p, name, w = 45, h = 55) {   # size in mm
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h,
         units = "mm", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".svg")), p, width = w, height = h,
         units = "mm", device = svg)
}

# simple dot plot: one dot per unit, bar = mean +/- SEM
dotbar <- function(d, ylab) {
  ggplot(d, aes(condition, y, fill = condition)) + prism_layers() +
    scale_fill_manual(values = fill_cols) +
    scale_x_discrete(labels = cond_labels) +
    labs(x = NULL, y = ylab) + theme_pub
}

# ---- 2. Fig A ---------------------------------------------------------------
if (all(summ$PC_valid)) {
  smallA <- obj %>% transmute(condition, replicate, fov, y = PC)
  ylabA <- "Partition coefficient"; nameA <- "FigA_partition_coefficient"
} else {
  message("PC invalid (dilute phase clipped) -> Fig A shows ACAN level at client puncta.")
  smallA <- cobj %>% transmute(condition, replicate, fov, y = rel_scaffold_at_client)
  ylabA <- "Normalized ACAN\nintensity at puncta"
  nameA <- "FigA_ACAN_level_at_client"
}
smallA <- smallA %>% filter(is.finite(y))
fovA  <- smallA %>% group_by(condition, replicate, fov) %>% summarise(y = median(y), .groups = "drop")
unitA <- to_unit(fovA)
write_csv(smallA, file.path(fig_dir, "FigA_FigS2_data_individual_puncta.csv"))
write_csv(unitA,  file.path(fig_dir, paste0("FigA_data_per_", unit, ".csv")))
cmpA <- compare(unitA, "welch")
lmmA <- if (unit == "fov" && nlevels(droplevels(smallA$condition)) >= 2) {
  smallA$fov_id <- interaction(smallA$condition, smallA$replicate, smallA$fov, drop = TRUE)
  lmerTest::lmer(y ~ condition + (1 | fov_id), data = smallA)
} else NULL
pA <- dotbar(unitA, ylabA)
if (!all(summ$PC_valid)) pA <- pA + geom_hline(yintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey40") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
pA <- add_p(pA, cmpA, ytop = max(unitA$y) * 1.12)
save_both(pA, nameA)

# Fig S2: all individual puncta behind Fig A
pS2 <- ggplot(smallA, aes(condition, y, fill = condition)) +
  stat_summary(fun = median, geom = "col", width = 0.7, colour = NA) +
  geom_point(shape = 21, fill = "white", colour = "black", size = 0.7, stroke = 0.2, position = position_swarm(width = 0.6, binfrac = 0.02, dx = 0.03)) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "errorbar", width = 0.5, colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = fill_cols) + scale_x_discrete(labels = cond_labels) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = ylabA) + theme_pub
if (!all(summ$PC_valid)) pS2 <- pS2 + geom_hline(yintercept = 1, linetype = 2, linewidth = 0.25, colour = "grey40")
save_both(pS2, "FigS2_individual_puncta")

# ---- 3. Fig B ---------------------------------------------------------------
labB <- c(pct_scaffold_clientpos_B1     = "BCAN/NCAN-positive\nACAN condensates",
          pct_client_inside_scaffold_B3 = "BCAN/NCAN puncta within\nACAN condensates")
fovB <- summ %>% select(condition, replicate, fov, all_of(names(labB))) %>%
  pivot_longer(all_of(names(labB)), names_to = "metric", values_to = "y") %>%
  mutate(metric = factor(labB[metric], levels = labB))
unitB <- to_unit(fovB)
write_csv(fovB,  file.path(fig_dir, "FigB_data_per_fov.csv"))
if (unit == "replicate") write_csv(unitB, file.path(fig_dir, "FigB_data_per_replicate.csv"))
cmpB <- do.call(rbind, lapply(levels(unitB$metric), function(m) {
  r <- compare(filter(unitB, metric == m), if (unit == "replicate") "welch" else "mwu")
  if (!is.null(r)) r$metric <- m; r }))
pB <- dotbar(unitB, "%") + facet_wrap(~metric) +
  scale_y_continuous(limits = c(0, 128), breaks = seq(0, 100, 25), expand = c(0, 0))
if (!is.null(cmpB)) for (m in unique(cmpB$metric))
  pB <- add_p(pB, cmpB[cmpB$metric == m, setdiff(names(cmpB), "metric")], ytop = 108,
              facet = factor(m, levels = labB))
save_both(pB, "FigB_percent_colocalized", w = 75)
# single panels (for figure assembly): left = counted per condensate; right = counted per punctum
for (k in seq_along(labB)) {
  m <- labB[[k]]; dk <- droplevels(filter(unitB, metric == m))
  ylab_k <- c("BCAN/NCAN-positive\nACAN condensates (%)", "BCAN/NCAN puncta within\nACAN condensates (%)")[k]
  pk <- dotbar(dk, ylab_k) + ggtitle(c("Per condensate", "Per punctum")[k]) +
    theme(plot.title = element_text(size = 8, hjust = 0.5, face = "plain")) +
    scale_y_continuous(limits = c(0, 128), breaks = seq(0, 100, 25), expand = c(0, 0))
  if (!is.null(cmpB)) pk <- add_p(pk, cmpB[cmpB$metric == m, setdiff(names(cmpB), "metric")], ytop = 108)
  save_both(pk, c("FigB1_ACAN_condensates_containing_client", "FigB2_puncta_within_condensates")[k], w = 47, h = 60)
}

# ---- 4. Fig C: where are the client puncta? ---------------------------------
locC <- cobj %>% count(condition, replicate, fov, location) %>%
  complete(nesting(condition, replicate, fov), location, fill = list(n = 0)) %>%
  group_by(condition, replicate, fov) %>% mutate(y = 100 * n / sum(n)) %>% ungroup()
unitC <- if (unit == "replicate") locC %>% group_by(condition, replicate, location) %>% summarise(y = mean(y), .groups = "drop") else locC
meanC <- unitC %>% group_by(condition, location) %>% summarise(pct = mean(y), sd = sd(y), n = n(), .groups = "drop") %>%
  mutate(location = factor(recode(location, inside = "Inside", contacting = "At edge", separate = "Beyond edge"),
                           levels = c("Beyond edge", "At edge", "Inside")))
write_csv(locC,  file.path(fig_dir, "FigC_data_per_fov.csv"))
if (unit == "replicate") write_csv(unitC, file.path(fig_dir, "FigC_data_per_replicate.csv"))
write_csv(meanC, file.path(fig_dir, "FigC_data_means.csv"))
pC <- ggplot(meanC, aes(condition, pct, fill = location)) +
  geom_col(width = 0.6, colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = ifelse(pct >= 6, sprintf("%.0f%%", pct), "")), position = position_stack(vjust = 0.5),
            size = 2.2, colour = ifelse(meanC$location == "Inside", "white", "black")) +
  scale_fill_manual(values = loc_cols,
                    breaks = c("Inside", "At edge", "Beyond edge")) +
  scale_x_discrete(labels = cond_labels) +
  scale_y_continuous(limits = c(0, 100.5), expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of BCAN/NCAN\npuncta (%)", fill = "Position relative to\nACAN condensate") +
  theme_pub + theme(legend.position = "right", legend.key.size = unit(3, "mm"),
                    legend.title = element_text(size = 7.5), legend.spacing.y = unit(0.8, "mm")) +
  guides(fill = guide_legend(byrow = TRUE))   # byrow + spacing.y = small gap between keys (ggplot2 >= 3.3)
save_both(pC, "FigC_client_location", w = 66)

# ---- 5. Fig S1: robustness of Fig B (left panel) -----------------------------
sens_cols <- grep("^pctB[12]_", names(summ), value = TRUE)
fovS <- summ %>% select(condition, replicate, fov, all_of(sens_cols)) %>%
  pivot_longer(all_of(sens_cols), names_to = "rule", values_to = "y") %>%
  mutate(rule_type = factor(ifelse(grepl("B1", rule), "Client signal above\nbackground + k SD",
                                             "Client area covering\n\u2265 x% of condensate"),
                            levels = c("Client signal above\nbackground + k SD", "Client area covering\n\u2265 x% of condensate")),
         level = recode(rule, pctB1_k2 = "2 SD", pctB1_k3 = "3 SD", pctB1_k5 = "5 SD",
                        `pctB2_f0.1` = "10%", `pctB2_f0.25` = "25%", `pctB2_f0.5` = "50%"))
write_csv(fovS, file.path(fig_dir, "FigS1_data_per_fov.csv"))
unitS <- if (unit == "replicate") fovS %>% group_by(condition, replicate, rule_type, level) %>% summarise(y = mean(y), .groups = "drop") else fovS
if (unit == "replicate") write_csv(unitS, file.path(fig_dir, "FigS1_data_per_replicate.csv"))
pS1 <- ggplot(unitS, aes(level, y, colour = condition)) +
  stat_summary(aes(group = condition), fun = mean, geom = "line", linewidth = 0.4) +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.15, linewidth = 0.3) +
  facet_wrap(~rule_type, scales = "free_x") +
  scale_colour_manual(values = line_cols, labels = function(x) cond_labels[x]) +
  scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(x = "Criterion", y = "BCAN/NCAN-positive\nACAN condensates (%)", colour = NULL) +
  theme_pub + theme_flat_x + theme(legend.position = "top", legend.key.size = unit(3, "mm"))
save_both(pS1, "FigS1_criterion_robustness", w = 80, h = 60)

# ---- 5b. apparent partition coefficient (raw single-channel images only) -----------
has_pc <- "rgb_export" %in% names(summ) && !any(as.logical(summ$rgb_export))
if (has_pc) {
  off <- if (is.na(DETECTOR_OFFSET)) 0 else DETECTOR_OFFSET
  pc_obj <- obj %>% left_join(summ %>% select(image_id, dil = dilute_client_median), by = "image_id") %>%
    mutate(y = (client_mean - off) / (dil - off)) %>% filter(is.finite(y), dil - off > 0)
  fovPC <- pc_obj %>% group_by(condition, replicate, fov) %>% summarise(y = median(y), .groups = "drop")
  unitPC <- to_unit(fovPC)
  write_csv(fovPC, file.path(fig_dir, "FigPC_data_per_fov.csv"))
  write_csv(unitPC, file.path(fig_dir, paste0("FigPC_data_per_", unit, ".csv")))
  cmpPC <- compare(unitPC, if (unit == "replicate") "welch" else "mwu")
  pPC <- dotbar(unitPC, if (is.na(DETECTOR_OFFSET)) "Apparent partition\ncoefficient" else "Partition coefficient") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  pPC <- add_p(pPC, cmpPC, ytop = max(unitPC$y) * 1.12)
  save_both(pPC, "FigPC_apparent_partition_coefficient")
}

# ---- 6. figure legends (generated from the data, so n is always current) ----
# Items in [brackets] must be filled in by hand before submission.
n_tab <- summ %>% group_by(condition) %>%
  summarise(n_fov = n(), n_rep = n_distinct(replicate),
            n_cond = sum(n_scaffold_condensates), n_punct = sum(n_client_puncta), .groups = "drop")
cl <- function(x) gsub("\n", " ", cond_labels[as.character(x)])
n_join <- function(f) paste(sprintf("%s, %s", cl(n_tab$condition), f(n_tab)), collapse = "; ")
n_fov_txt <- if (length(unique(n_tab$n_fov)) == 1) sprintf("%d FOVs per condition", n_tab$n_fov[1]) else
  n_join(function(t) sprintf("%d FOVs", t$n_fov))
n_rep_txt <- if (length(unique(n_tab$n_rep)) == 1) sprintf("%d independent preparations per condition", n_tab$n_rep[1]) else
  n_join(function(t) sprintf("%d independent preparation%s", t$n_rep, ifelse(t$n_rep == 1, "", "s")))
n_punct   <- n_join(function(t) sprintf("%d client puncta", t$n_punct))
n_cond    <- n_join(function(t) sprintf("%d ACAN condensates", t$n_cond))
dot_txt <- if (unit == "replicate") "Each dot is one independent preparation" else "Each dot is one field of view (FOV)"
n_txt   <- if (unit == "replicate") n_rep_txt else
  paste0(n_fov_txt, " [from a single preparation per condition; FOVs are technical replicates]")
statA <- if (!is.null(cmpA)) paste0("Statistics: ", cmpA$test[1], if (nrow(cmpA) > 1 && cmpA$test[1] != TEST_ANOVA) " with Holm correction" else "",
  if (!is.null(lmmA)) "; the result was confirmed by a linear mixed model on all individual puncta with FOV as a random effect (stats.txt)" else "",
  ".") else ""
statB <- if (!is.null(cmpB)) paste0("Statistics: ", cmpB$test[1], if (length(unique(paste(cmpB$group1, cmpB$group2))) > 1) " with Holm correction" else "", ".") else ""
common <- paste0("ACAN, [fluorophore]; BCAN, [fluorophore]; NCAN, [fluorophore]; [protein concentrations, 6% PEG, buffer, time after mixing]. ",
  "Image analysis is described in Methods.")

legA <- if (all(summ$PC_valid)) {
  paste0("Partition coefficient of the client into ACAN condensates (client intensity inside condensates / client intensity in the surrounding solution). ",
    dot_txt, " (median of all condensates in it); bars, mean \u00b1 SEM. n = ", n_txt, ". ", statA, " ", common)
} else {
  paste0("ACAN level at client puncta. For every BCAN or NCAN punctum, the ACAN intensity at that position was divided by the ",
    "median ACAN intensity of the ACAN condensates in the same image, so that 1 (dashed line) means the punctum sits in ACAN signal as strong as ",
    "the median condensate of that image, values above 1 mean stronger ACAN signal (e.g. larger, brighter condensates) and values below 1 weaker ACAN signal. ", sub(")$", "", dot_txt), ifelse(unit == "replicate", "", "; value = median of all client puncta in that image"), ifelse(grepl("FOV", dot_txt), ")", ""), "; bars, mean \u00b1 SEM. ",
    "n = ", n_txt, " (", n_punct, "). ", statA, " ", common)
}
legB <- paste0("Co-localization of ACAN condensates and client puncta, counted in two directions. Left (per condensate), BCAN/NCAN-positive ACAN condensates: ",
  "percentage of ACAN condensates in which client signal is detected (mean client intensity inside the condensate above the background mean + 3 SD). ",
  "Right (per punctum), percentage of BCAN or NCAN puncta that lie within ACAN condensates (at least half of the punctum area overlapping a condensate). ",
  "The two measures have different denominators: almost every ACAN condensate contains some NCAN (left), yet many NCAN puncta lie at or beyond the edge of the bright condensate body (right; see Fig. C). ", dot_txt, "; bars, mean \u00b1 SEM. n = ", n_txt,
  " (", n_cond, "; ", n_punct, "). ", statB, " ", common)
legC <- paste0("Position of client puncta relative to the edge of the ACAN condensate body (edge defined by the pre-specified Otsu threshold). ",
  "Inside, at least half of the punctum overlaps a condensate; at edge, less than half overlaps but the punctum contacts the condensate edge; ",
  "beyond edge, no contact with the bright condensate body. With a more inclusive threshold (Li) that also captures the dim ACAN rim around condensates, ",
  "\u2265 99% of BCAN and NCAN puncta lie within ACAN signal (Fig. S5), i.e. 'at edge' and 'beyond edge' puncta sit in dim ACAN material at the condensate periphery rather than free in solution. ",
  "Bars show the mean percentage across ", ifelse(unit == "replicate", "independent preparations", "FOVs"),
  "; numbers inside bars give the percentages. n = ", n_txt, " (", n_punct, "). ", common)
legS1 <- paste0("The percentage of BCAN/NCAN-positive ACAN condensates (Fig. B, left) recalculated with stricter or looser criteria. ",
  "Left, client signal above background mean + 2, 3 (used in Fig. B) or 5 SD. Right, client area covering at least 10%, 25% or 50% ",
  "of the condensate. Points, mean \u00b1 SEM across ", ifelse(unit == "replicate", "independent preparations", "FOVs"), ". n = ", n_txt, ". ", common)
legS2 <- paste0("All individual client puncta underlying Fig. A. Each dot is one punctum (", n_punct, "); bar and black line, median. Dashed line, ACAN level of the median ACAN condensate in the same image. ", common)

legPC <- if (has_pc) paste0(ifelse(is.na(DETECTOR_OFFSET), "Apparent partition coefficient", "Partition coefficient"),
  " of BCAN and NCAN into ACAN condensates: mean client intensity inside each ACAN condensate divided by the median client intensity of the surrounding dilute phase ",
  "(pixels more than 3 px from any condensate or punctum), calculated on raw single-channel images",
  ifelse(is.na(DETECTOR_OFFSET), "; values are not corrected for the detector offset and are therefore lower bounds (any positive offset increases the ratio)",
         paste0("; a detector offset of ", DETECTOR_OFFSET, " grey levels (buffer-only image) was subtracted")), ". ",
  dot_txt, " (median of the condensates in each field of view, averaged per preparation); bars, mean \u00b1 SEM. n = ", n_txt, ". ",
  if (!is.null(cmpPC)) paste0("Statistics: ", cmpPC$test[1], ". ") else "", common) else NULL
legends <- list(
  list(file = nameA,                        title = "Figure A",  text = legA,  plot = pA,  w = 45, h = 55),
  list(file = "FigB_percent_colocalized",   title = "Figure B",  text = legB,  plot = pB,  w = 75, h = 55),
  list(file = "FigC_client_location",       title = "Figure C",  text = legC,  plot = pC,  w = 60, h = 55),
  list(file = "FigS1_criterion_robustness", title = "Figure S1", text = legS1, plot = pS1, w = 80, h = 60),
  list(file = "FigS2_individual_puncta",    title = "Figure S2", text = legS2, plot = pS2, w = 45, h = 55))
if (has_pc) legends <- c(legends, list(list(file = "FigPC_apparent_partition_coefficient", title = "Figure PC", text = legPC, plot = pPC, w = 45, h = 55)))

doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "Figure legends (auto-generated by plot_recruitment.R)", style = "heading 1")
doc <- officer::body_add_par(doc, paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"),
                                        "| Items in [brackets] must be completed by hand."), style = "Normal")
all_txt <- character()
for (L in legends) {
  # (1) one plain-text legend next to each figure
  writeLines(c(paste0(L$title, " (", L$file, ".pdf)"), "", L$text),
             file.path(fig_dir, paste0(L$file, "_legend.txt")), useBytes = TRUE)
  all_txt <- c(all_txt, paste0(L$title, " (", L$file, ".pdf)"), L$text, "")
  # (2) Word document with all legends
  doc <- officer::body_add_par(doc, paste0(L$title, " - ", L$file, ".pdf"), style = "heading 2")
  doc <- officer::body_add_par(doc, L$text, style = "Normal")
  # (3) review copy: figure + legend on one page (vector PDF)
  cap <- paste(strwrap(L$text, width = 95), collapse = "\n")
  p_cap <- L$plot + labs(caption = cap) +
    theme(plot.caption = element_text(hjust = 0, size = 6.5, lineheight = 1.1),
          plot.caption.position = "plot")
  if (L$file == "FigC_client_location") p_cap <- p_cap + theme(legend.position = "right")
  ggsave(file.path(fig_dir, paste0(L$file, "_with_legend.pdf")), p_cap,
         width = 180, height = L$h + 12 + 3.2 * length(strsplit(cap, "\n")[[1]]),
         units = "mm", device = cairo_pdf)
}
writeLines(all_txt, file.path(fig_dir, "Figure_legends_all.txt"), useBytes = TRUE)
print(doc, target = file.path(fig_dir, "Figure_legends_all.docx"))

# ---- 7. statistics -----------------------------------------------------------
sink(file.path(fig_dir, "stats.txt"))
cat("Statistical unit used:", unit, "\n")
cat("  replicate = independent preparations (publication standard, needs >= 3 per condition)\n")
cat("  fov       = fields of view = technical replicates; must be stated as such in the legend\n\n")
print(as.data.frame(n_rep_cond %>% rename(n_independent_preparations = n)))
cat("\n== Fig A ==\n"); print(cmpA)
if (!is.null(lmmA)) { cat("\nFig A supporting LMM (puncta, random intercept = FOV):\n")
  print(summary(lmmA)$coefficients) }
cat("\n== Fig B ==\n"); print(cmpB)
if (has_pc) { cat("\n== Apparent partition coefficient ==\n"); print(cmpPC); print(as.data.frame(unitPC)) }
cat("\nValues entering the tests:\n"); print(as.data.frame(unitA)); print(as.data.frame(unitB))
cat("\nClient location (mean % across", ifelse(unit == "replicate", "independent preparations", "FOVs"), "):\n"); print(as.data.frame(meanC))
sink()
cat("Figures, panel data and stats written to", normalizePath(fig_dir), "\n")

# =============================================================================
# style_prism.R - shared figure style (GraphPad-Prism-like), used by every plot_*.R
#   filled bars without outline, open circles (white fill, black edge) spread sideways,
#   black mean +/- SEM error bars, thick left/bottom axes only, outward ticks,
#   Arial, x labels rotated 45 degrees.
# Change colours / sizes HERE and every figure follows.
# =============================================================================
suppressPackageStartupMessages(library(ggplot2))
prism_font <- "Arial"     # Windows/macOS: Arial. Linux: fontconfig uses the metric-identical Liberation Sans.

# bar fill colours (pastel, as in the template) and darker line colours for line plots
fill_cols <- c("ACAN alone" = "#C8C8C8", "ACAN+BCAN" = "#F8998B", "ACAN+NCAN" = "#BFA2DB", "ACAN+BCAN+NCAN" = "#F5C47E",
               "ACAN + BCAN" = "#F8998B", "ACAN + NCAN" = "#BFA2DB", "ACAN + BCAN + NCAN" = "#F5C47E",
               "Two proteins" = "#F8998B", "Three proteins" = "#F5C47E",
               ACAN = "#C8C8C8", BCAN = "#F8998B", NCAN = "#BFA2DB")
line_cols <- c("ACAN alone" = "#7F7F7F", "ACAN+BCAN" = "#E5533F", "ACAN+NCAN" = "#8A5BB8", "ACAN+BCAN+NCAN" = "#DB8E1F",
               "ACAN + BCAN" = "#E5533F", "ACAN + NCAN" = "#8A5BB8", "ACAN + BCAN + NCAN" = "#DB8E1F",
               ACAN = "#2E9E5B", BCAN = "#E5533F", NCAN = "#8A5BB8")

# position of client puncta relative to the ACAN condensate (stacked bars, Fig C / E3): warm sequential ramp,
# dark = inside the condensate -> light = beyond its edge
loc_cols <- c(Inside = "#A63D2F", `At edge` = "#E59866", `Beyond edge` = "#FBE5CF")

ax_lwd <- 0.45   # axis line width (mm) ~ 1.3 pt
theme_pub <- theme_classic(base_size = 9, base_family = prism_font) +
  theme(axis.line = element_line(colour = "black", linewidth = ax_lwd, lineend = "square"),
        axis.ticks = element_line(colour = "black", linewidth = ax_lwd),
        axis.ticks.length = unit(1.5, "mm"),
        axis.text = element_text(colour = "black", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(colour = "black", size = 9),
        strip.background = element_blank(), strip.text = element_text(colour = "black", size = 8),
        legend.position = "none", legend.text = element_text(size = 8),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
theme_flat_x <- theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))   # for numeric / short x axes
update_geom_defaults("text", list(family = prism_font, size = 2.6))

# ---- open circles spread sideways like a Prism scatter-dot plot (deterministic) ----
PositionSwarm <- ggproto("PositionSwarm", Position,
  required_aes = c("x", "y"),
  setup_params = function(self, data) list(width = self$width, binfrac = self$binfrac, dx = self$dx),
  compute_layer = function(self, data, params, layout) {
    ok <- is.finite(data$y)
    yr <- diff(range(data$y[ok])); if (!is.finite(yr) || yr == 0) yr <- 1
    bw <- yr * params$binfrac
    grp <- interaction(data$PANEL, round(data$x, 6), drop = TRUE)
    off <- numeric(nrow(data))
    for (g in levels(grp)) {
      idx <- which(grp == g & ok)
      if (!length(idx)) next
      bins <- floor(data$y[idx] / bw)
      for (b in unique(bins)) {
        j <- idx[bins == b]; j <- j[order(data$y[j])]; k <- length(j)
        dx <- min(params$dx, params$width / max(k - 1, 1))
        off[j] <- (seq_len(k) - (k + 1) / 2) * dx
      }
    }
    data$x <- data$x + off
    data
  })
position_swarm <- function(width = 0.55, binfrac = 0.035, dx = 0.09)
  ggproto(NULL, PositionSwarm, width = width, binfrac = binfrac, dx = dx)

# ---- bar + SEM + open circles; map fill in the main aes() of the plot ----------------
prism_layers <- function(point_size = 1.5, bar_width = 0.7, dx = 0.09) list(
  stat_summary(fun = mean, geom = "col", width = bar_width, colour = NA),
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.28, linewidth = 0.4, colour = "black"),
  geom_point(shape = 21, fill = "white", colour = "black", size = point_size, stroke = 0.35,
             position = position_swarm(dx = dx)))
# x axis drawn at y = 0 (for data with negative values, as in the template)
zero_axis <- function() list(geom_hline(yintercept = 0, linewidth = ax_lwd, colour = "black"),
                             theme(axis.line.x = element_blank(), axis.ticks.x = element_blank()))
one_line <- function(x) gsub("\\s*\\n\\s*", " ", x)
prism_style_loaded <- TRUE

# ---- statistics used by all figures (lab standard) ------------------------------------
# two groups: two-tailed unpaired Student's t-test; >= 3 groups: one-way ANOVA + Tukey; paired: paired t-test
student_p <- function(a, b) if (length(unique(c(a, b))) == 1) 1 else t.test(a, b, var.equal = TRUE)$p.value
welch_p   <- function(a, b) if (length(unique(c(a, b))) == 1) 1 else t.test(a, b)$p.value     # sensitivity only
anova_tukey <- function(y, g, log = FALSE) {
  g <- droplevels(factor(g)); if (log) y <- log10(y)
  fit <- aov(y ~ g); tk <- TukeyHSD(fit)$g; lv <- levels(g)
  pr <- combn(seq_along(lv), 2, simplify = FALSE)
  out <- do.call(rbind, lapply(pr, function(k) {
    nm <- paste0(lv[k[2]], "-", lv[k[1]])
    data.frame(group1 = lv[k[1]], group2 = lv[k[2]], p_adj = unname(tk[nm, "p adj"]))
  }))
  out$p <- out$p_adj
  attr(out, "anova_p") <- summary(fit)[[1]][["Pr(>F)"]][1]
  out
}
TEST_T     <- "two-tailed unpaired Student's t-test"
TEST_ANOVA <- "one-way ANOVA with Tukey's multiple comparisons test"
TEST_PAIR  <- "two-tailed paired t-test"

# ---- pixel size (from settings.txt) and area labels ------------------------------------
read_px_um <- function(f = "settings.txt") {
  if (!file.exists(f)) return(NA_real_)
  l <- sub("#.*", "", readLines(f, warn = FALSE, encoding = "UTF-8"))
  m <- grep("^\\s*pixel_um\\s*=", l, value = TRUE)
  if (!length(m)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(sub(".*=", "", m[1]))))
}
PIXEL_UM <- read_px_um()
UM2 <- "\u00b5m\u00b2"
fmt_um2 <- function(px) formatC(px * PIXEL_UM^2, digits = 2, format = "fg")
area_txt   <- function(px) if (is.na(PIXEL_UM)) paste0(px, " px") else paste0(fmt_um2(px), " ", UM2)
area_range <- function(a, b) if (is.na(PIXEL_UM)) sprintf("%d-%d px", a, b) else paste0(fmt_um2(a), "-", fmt_um2(b), " ", UM2)
len_txt    <- function(px) if (is.na(PIXEL_UM)) paste0(px, " px") else paste0(px, " px (", formatC(px * PIXEL_UM, digits = 2, format = "fg"), " \u00b5m)")
area_unit_lab <- if (is.na(PIXEL_UM)) "px" else UM2

# ---- was ACAN saturated in the three-protein images? (legend caveats only if > 1 % of pixels) ----
TRI_SAT <- local({
  f <- "results/three_protein/image_summary.csv"
  if (!file.exists(f)) FALSE else {
    d <- utils::read.csv(f); d <- d[grepl("^Triple", d$condition), ]
    nrow(d) > 0 && mean(d$saturated_frac_scaffold, na.rm = TRUE) > 0.01
  }
})
DETECTOR_OFFSET <- local({
  f <- "settings.txt"; v <- NA_real_
  if (file.exists(f)) {
    l <- sub("#.*", "", readLines(f, warn = FALSE, encoding = "UTF-8"))
    m <- grep("^\\s*detector_offset\\s*=", l, value = TRUE)
    if (length(m)) v <- suppressWarnings(as.numeric(trimws(sub(".*=", "", m[1]))))
  }
  v
})

# =============================================================================
# make_figures.R  -  makes every figure, legend and statistics file.
# Normally started automatically by run_analysis.bat.
# By hand: open this file in RStudio and click "Source" (no editing needed).
# =============================================================================
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
need <- c("ggplot2", "dplyr", "readr", "tidyr", "officer", "lme4", "lmerTest", "patchwork", "svglite")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)

# --- find the analysis folder (the folder that contains 'scripts' and 'results') ---
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
if (is.null(this_file)) {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) this_file <- normalizePath(sub("^--file=", "", a[1]))
}
if (is.null(this_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  this_file <- rstudioapi::getSourceEditorContext()$path
base <- dirname(dirname(this_file))
setwd(base)
scripts <- file.path(base, "scripts")
message("Analysis folder: ", base)

run_part <- function(script, needs, label) {
  if (!all(file.exists(needs))) { message("  (skipped ", label, ": no data)"); return(invisible(FALSE)) }
  message(">>> ", label)
  e <- new.env(parent = globalenv())
  ok <- tryCatch({ sys.source(file.path(scripts, script), envir = e); TRUE },
                 error = function(err) { message("  !! ERROR in ", script, ": ", conditionMessage(err)); FALSE })
  invisible(ok)
}
run_part("plot_recruitment.R", "results/recruitment/image_summary.csv", "Fig A-C, S1-S2: recruitment (two-protein mixes)")
run_part("plot_shape.R",       "results/shape/shape_fov_summary.csv",   "Fig D, S3-S4: condensate size and shape")
run_part("plot_triple.R",      c("results/three_protein/image_summary.csv", "results/shape/shape_fov_summary.csv"),
         "Fig E1-E5: three-protein mix")
run_part("plot_driver.R",      "results/three_protein/driver_analysis_objects.csv", "Fig E6-E7: which protein NCAN follows")
run_part("plot_position.R",    c("results/position/two_protein/radial_fov_summary.csv", "results/position/three_protein/radial_fov_summary.csv"),
         "Fig F1-F3, S5-S7: where in the condensate BCAN / NCAN sit")

# --- one Word file with all legends + an overview of all figures ---------------------
docs <- c("results/recruitment/figures/Figure_legends_all.docx", "results/shape/figures/Figure_legends_shape.docx",
          "results/three_protein/figures/Figure_legends_triple.docx", "results/three_protein/figures/Figure_legends_driver.docx",
          "results/position/figures/Figure_legends_position.docx")
docs <- docs[file.exists(docs)]
if (length(docs)) {
  d <- officer::read_docx()
  for (x in docs) d <- officer::body_add_docx(d, src = x)
  print(d, target = "results/ALL_figure_legends.docx")
}
figs <- list.files("results", pattern = "^Fig.*\\.pdf$", recursive = TRUE)
figs <- figs[!grepl("_with_legend|check_", figs)]
writeLines(c("All figures (vector PDF; an SVG, a *_legend.txt and the data CSV sit next to each):", "",
             paste0("  results/", figs), "",
             "Figure + legend on one page: the *_with_legend.pdf files.",
             "All legends in one Word file: results/ALL_figure_legends.docx",
             "Statistics: stats*.txt in each figures folder."),
           "results/figures_overview.txt")
message("Done. See results/figures_overview.txt")

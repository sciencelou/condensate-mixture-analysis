#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_source_data.py - one Excel workbook with the exact numbers behind every figure.

  results/Figure_source_data.xlsx
    sheet "Index"  : figure -> sheet -> which columns to copy (yellow) -> what each value is -> test, P
    one sheet per figure panel, laid out like a GraphPad Prism table:
      * one column per group (yellow = copy these), one row per independent preparation
      * paired figures: the same row = the same preparation (Prism "paired" column table)
      * line figures: X column (blue) + Y columns (yellow)
      * below each table: mean, SD, SEM, n as live Excel formulas, then the statistical test and P value

Runs automatically at the end of run_analysis.bat (after the R figures).
By hand (cmd, in the LLPS_analysis folder):   python scripts\\make_source_data.py
Author: Bin Luo, 2026. Written with AI assistance (Claude, Anthropic).
"""
import os
import sys

import warnings
import numpy as np
import pandas as pd
from scipy import stats
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter as L
warnings.filterwarnings("ignore", category=RuntimeWarning)

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(BASE, "results")
OUT = os.path.join(RES, "Figure_source_data.xlsx")


def _px_um():
    try:
        for line in open(os.path.join(BASE, "settings.txt"), encoding="utf-8-sig"):
            line = line.split("#")[0]
            if line.strip().lower().startswith("pixel_um") and "=" in line:
                return float(line.split("=", 1)[1].strip())
    except (OSError, ValueError):
        pass
    return None


PX_UM = _px_um()
AREA_LAB = "Condensate area (\u00b5m\u00b2; median per field)" if PX_UM else "Condensate area (px; median per field)"

FONT = "Arial"
YELLOW = PatternFill("solid", fgColor="FFF2A8")
BLUE = PatternFill("solid", fgColor="CFE2F3")
GREY = PatternFill("solid", fgColor="E7E6E6")
THIN = Side(style="thin", color="999999")
BOX = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
F = lambda **k: Font(**{"name": FONT, "size": 10, **k})


def rd(*p):
    path = os.path.join(RES, *p)
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, encoding="utf-8")


def to_rep(d, keys, col="y"):
    """FOV values -> mean per independent preparation (same as the R scripts)."""
    return d.groupby(keys + ["replicate"], as_index=False)[col].mean()


def holm(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    adj = np.empty(m); run = 0
    for i, k in enumerate(o):
        run = max(run, (m - i) * p[k]); adj[k] = min(1, run)
    return adj


def welch(a, b, log=False):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if log:
        a, b = np.log10(a), np.log10(b)
    if np.allclose(np.r_[a, b], a[0]):
        return 1.0
    return float(stats.ttest_ind(a, b, equal_var=False).pvalue)


def student(a, b, log=False):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if log:
        a, b = np.log10(a), np.log10(b)
    if np.allclose(np.r_[a, b], a[0]):
        return 1.0
    return float(stats.ttest_ind(a, b, equal_var=True).pvalue)


def anova_tukey(groups, log=False):
    g = [np.log10(np.asarray(x, float)) if log else np.asarray(x, float) for x in groups]
    return float(stats.f_oneway(*g).pvalue), stats.tukey_hsd(*g).pvalue


TEST_T = "two-tailed unpaired Student's t-test"
TEST_ANOVA = "one-way ANOVA with Tukey's multiple comparisons test"
TEST_PAIR = "two-tailed paired t-test"


def paired(a, b):
    return float(stats.ttest_rel(np.asarray(a, float), np.asarray(b, float)).pvalue)


def fmt_p(p):
    return "P < 0.0001" if p < 1e-4 else ("P > 0.99" if p > 0.99 else f"P = {p:.2g}")


# --------------------------------------------------------------------------------------------
wb = Workbook()
index_rows = []


def new_sheet(name, title, yaxis, unit_txt, source):
    ws = wb.create_sheet(name[:31])
    ws["A1"] = title; ws["A1"].font = F(bold=True, size=12)
    ws["A2"] = "Y axis: " + yaxis; ws["A2"].font = F()
    ws["A3"] = unit_txt; ws["A3"].font = F()
    ws["A4"] = "Copy the YELLOW columns (with their header) into Prism / Excel to redraw this panel. Source: " + source
    ws["A4"].font = F(italic=True, color="555555")
    ws.column_dimensions["A"].width = 18
    return ws


def put_table(ws, top, left, rowlab_name, rowlabs, cols, fill=YELLOW, xfill=None, summary=True):
    """cols: list of (header, list_of_values). Returns (first_col_letter, last_col_letter, last_row)."""
    c0 = left
    h = ws.cell(row=top, column=c0, value=rowlab_name); h.font = F(bold=True); h.fill = GREY; h.border = BOX
    for i, r in enumerate(rowlabs):
        c = ws.cell(row=top + 1 + i, column=c0, value=r); c.font = F(); c.fill = GREY; c.border = BOX
    nmax = max(len(v) for _, v in cols) if cols else 0
    for j, (hdr, vals) in enumerate(cols):
        col = c0 + 1 + j
        f = xfill if (xfill is not None and j == 0) else fill
        c = ws.cell(row=top, column=col, value=hdr); c.font = F(bold=True); c.fill = f; c.border = BOX
        c.alignment = Alignment(wrap_text=True, horizontal="center")
        for i, v in enumerate(vals):
            cc = ws.cell(row=top + 1 + i, column=col, value=(None if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)))
            cc.font = F(); cc.fill = f; cc.border = BOX; cc.number_format = "0.000"
        ws.column_dimensions[L(col)].width = max(14, min(26, len(str(hdr)) + 2))
    last = top + max(nmax, len(rowlabs))
    if summary and cols:
        r0 = last + 1
        for k, (lab, fn) in enumerate([("Mean", "AVERAGE"), ("SD", "STDEV"), ("SEM", None), ("n", "COUNT")]):
            rr = r0 + k
            c = ws.cell(row=rr, column=c0, value=lab); c.font = F(bold=True)
            for j in range(len(cols)):
                col = c0 + 1 + j
                rng = f"{L(col)}{top + 1}:{L(col)}{last}"
                if xfill is not None and j == 0:
                    continue
                if fn:
                    ws.cell(row=rr, column=col, value=f"={fn}({rng})")
                else:
                    ws.cell(row=rr, column=col, value=f"=STDEV({rng})/SQRT(COUNT({rng}))")
                ws.cell(row=rr, column=col).font = F(); ws.cell(row=rr, column=col).number_format = "0.000"
        last = r0 + 3
    return L(c0 + 1), L(c0 + len(cols)), last


def put_stats(ws, row, lines):
    for i, t in enumerate(lines):
        c = ws.cell(row=row + i, column=1, value=t); c.font = F(bold=(i == 0))
    return row + len(lines)


def add_index(fig, sheet, cols, what, test):
    index_rows.append((fig, sheet, cols, what, test))


UNIT = "Each value = one independent preparation (mean of its fields of view). Bars = mean +/- SEM."
reps = ["exp1", "exp2", "exp3"]


def wide(d, group_col, groups, ycol="y"):
    out = []
    for g in groups:
        s = d[d[group_col] == g].set_index("replicate")[ycol]
        out.append((g, [s.get(r, np.nan) for r in reps]))
    return out


# ============================ recruitment (Fig A, B, C, S1, S2) ============================
A = rd("recruitment", "figures", "FigA_data_per_replicate.csv")
if A is not None:
    ws = new_sheet("FigA", "Fig A - ACAN level at client puncta", "Normalized ACAN intensity at BCAN/NCAN puncta (1 = median ACAN condensate of the same image)",
                   UNIT, "results/recruitment/figures/FigA_data_per_replicate.csv")
    cols = wide(A, "condition", ["ACAN+BCAN", "ACAN+NCAN"])
    a, b, last = put_table(ws, 6, 1, "Preparation", reps, cols)
    p = student(cols[0][1], cols[1][1])
    put_stats(ws, last + 2, ["Statistics", f"{TEST_T}, ACAN+BCAN vs ACAN+NCAN: {fmt_p(p)} (p = {p:.4g})",
                             f"Sensitivity, Welch's t-test: {fmt_p(welch(cols[0][1], cols[1][1]))}"])
    add_index("Fig A", "FigA", f"{a}:{b}", "ACAN level at client puncta", f"Student's t-test, {fmt_p(p)}")

    S2 = rd("recruitment", "figures", "FigA_FigS2_data_individual_puncta.csv")
    ws = new_sheet("FigS2", "Fig S2 - all individual puncta behind Fig A", "Normalized ACAN intensity at BCAN/NCAN puncta",
                   "Each value = one client punctum. Bar and line = median.", "results/recruitment/figures/FigA_FigS2_data_individual_puncta.csv")
    cols = [(g, S2.loc[S2.condition == g, "y"].tolist()) for g in ["ACAN+BCAN", "ACAN+NCAN"]]
    n = max(len(v) for _, v in cols)
    a, b, last = put_table(ws, 6, 1, "Punctum #", list(range(1, n + 1)), cols, summary=False)
    put_stats(ws, 6, [])
    for j, (g, v) in enumerate(cols):
        ws.cell(row=5, column=2 + j, value=f"median = {np.median(v):.3f} (n = {len(v)})").font = F(italic=True)
    add_index("Fig S2", "FigS2", f"{a}:{b}", "individual puncta (descriptive)", "none (descriptive)")

PCd = rd("recruitment", "figures", "FigPC_data_per_replicate.csv")
if PCd is not None:
    ws = new_sheet("FigPC", "Fig PC - apparent partition coefficient", "Apparent partition coefficient (client inside / dilute phase; lower bound unless a detector offset is set)",
                   UNIT, "results/recruitment/figures/FigPC_data_per_replicate.csv")
    cols = wide(PCd, "condition", ["ACAN+BCAN", "ACAN+NCAN"])
    a_, b_, last = put_table(ws, 6, 1, "Preparation", reps, cols)
    p = student(cols[0][1], cols[1][1])
    put_stats(ws, last + 2, ["Statistics", f"{TEST_T}: {fmt_p(p)} (p = {p:.4g})", f"Sensitivity, Welch's t-test: {fmt_p(welch(cols[0][1], cols[1][1]))}"])
    add_index("Fig PC", "FigPC", f"{a_}:{b_}", "apparent partition coefficient", f"Student's t-test, {fmt_p(p)}")

B = rd("recruitment", "figures", "FigB_data_per_replicate.csv")
if B is not None:
    B["metric"] = B["metric"].str.replace("\n", " ")
    ws = new_sheet("FigB", "Fig B - co-localization (two panels)", "Left (per condensate): BCAN/NCAN-positive ACAN condensates (%); right (per punctum): BCAN/NCAN puncta within ACAN condensates (%)", UNIT, "results/recruitment/figures/FigB_data_per_replicate.csv")
    stats_lines = [f"Statistics ({TEST_T}, ACAN+BCAN vs ACAN+NCAN; Welch's t-test in brackets as sensitivity check)"]
    rngs = []
    left = 1
    for m in B["metric"].unique():
        ws.cell(row=5, column=left + 1, value=m).font = F(bold=True)
        cols = wide(B[B.metric == m], "condition", ["ACAN+BCAN", "ACAN+NCAN"])
        a, b, last = put_table(ws, 6, left, "Preparation", reps, cols)
        p = student(cols[0][1], cols[1][1])
        stats_lines.append(f"{m}: {fmt_p(p)} (p = {p:.4g}) [Welch {fmt_p(welch(cols[0][1], cols[1][1]))}]")
        rngs.append(f"{a}:{b} ({m})")
        left += 4
    put_stats(ws, last + 2, stats_lines)
    add_index("Fig B", "FigB", "; ".join(rngs), "BCAN/NCAN-positive ACAN condensates (%) / BCAN/NCAN puncta within ACAN condensates (%)", "Student's t-test: " + "; ".join(l.split(" (p =")[0] for l in stats_lines[1:]))

C = rd("recruitment", "figures", "FigC_data_per_replicate.csv")
if C is not None:
    C["location"] = C["location"].map({"inside": "Inside", "contacting": "At edge", "separate": "Beyond edge"})
    ws = new_sheet("FigC", "Fig C - position of BCAN/NCAN puncta (stacked bars)", "Fraction of BCAN/NCAN puncta (%) at each position relative to the ACAN condensate",
                   "Top table = values plotted (mean of the 3 preparations). Lower table = each preparation.", "results/recruitment/figures/FigC_data_per_replicate.csv")
    locs = ["Inside", "At edge", "Beyond edge"]
    cols = [(g, [C[(C.condition == g) & (C.location == l)]["y"].mean() for l in locs]) for g in ["ACAN+BCAN", "ACAN+NCAN"]]
    a, b, last = put_table(ws, 6, 1, "Location", locs, cols, summary=False)
    ws.cell(row=last + 2, column=1, value="Per preparation").font = F(bold=True)
    cols2 = []
    for g in ["ACAN+BCAN", "ACAN+NCAN"]:
        for r in reps:
            cols2.append((f"{g} {r}", [C[(C.condition == g) & (C.location == l) & (C.replicate == r)]["y"].mean() for l in locs]))
    put_table(ws, last + 3, 1, "Location", locs, cols2, fill=GREY, summary=False)
    add_index("Fig C", "FigC", f"{a}:{b} (rows 7-9)", "BCAN/NCAN puncta inside / at edge / beyond edge (%)", "descriptive (Inside % is tested in Fig B)")

S1 = rd("recruitment", "figures", "FigS1_data_per_replicate.csv")
if S1 is not None:
    S1["rule_type"] = S1["rule_type"].str.replace("\n", " ")
    ws = new_sheet("FigS1", "Fig S1 - criterion robustness (grouped table: rows = criterion, sub-columns = preparations)",
                   "BCAN/NCAN-positive ACAN condensates (%)", "Prism: Grouped table, 3 replicate values per cell.", "results/recruitment/figures/FigS1_data_per_replicate.csv")
    levels = ["2 SD", "3 SD", "5 SD", "10%", "25%", "50%"]
    cols = []
    for g in ["ACAN+BCAN", "ACAN+NCAN"]:
        for r in reps:
            cols.append((f"{g} {r}", [S1[(S1.condition == g) & (S1.level == lv) & (S1.replicate == r)]["y"].mean() for lv in levels]))
    a, b, last = put_table(ws, 6, 1, "Criterion", levels, cols, summary=False)
    add_index("Fig S1", "FigS1", f"{a}:{b}", "% BCAN/NCAN-positive condensates under each criterion", "descriptive")

# ============================ shape (Fig D, S3, S4) ============================
D = rd("shape", "figures", "FigD_data_per_fov.csv")
S3 = rd("shape", "figures", "FigS3_data_per_fov.csv")
SG = ["ACAN alone", "ACAN+BCAN", "ACAN+NCAN"]


def three_group_sheet(name, title, yaxis, d, col, log=False, source=""):
    r = to_rep(d.rename(columns={col: "y"}).dropna(subset=["y"]), ["condition"])
    ws = new_sheet(name, title, yaxis, UNIT, source)
    cols = wide(r, "condition", SG)
    a, b, last = put_table(ws, 6, 1, "Preparation", reps, cols)
    pairs = [(0, 1), (0, 2), (1, 2)]
    pan, tk = anova_tukey([c[1] for c in cols], log)
    sens = holm([welch(cols[i][1], cols[j][1], log) for i, j in pairs])
    lines = [f"Statistics: {TEST_ANOVA}" + (" on log10 values" if log else "") + f"; one-way ANOVA {fmt_p(pan)}"]
    lines += [f"{cols[i][0]} vs {cols[j][0]}: Tukey {fmt_p(tk[i, j])}" for i, j in pairs]
    lines += ["Sensitivity (pairwise Welch's t-test, Holm): " + "; ".join(f"{cols[i][0]} vs {cols[j][0]} {fmt_p(q)}" for (i, j), q in zip(pairs, sens))]
    put_stats(ws, last + 2, lines)
    add_index(title.split(" - ")[0], name, f"{a}:{b}", yaxis, "ANOVA + Tukey" + (" (log10)" if log else "") + f"; ANOVA {fmt_p(pan)}; " + "; ".join(l.replace("Tukey ", "") for l in lines[1:4]))


if D is not None:
    three_group_sheet("FigD1", "Fig D1 - condensate area", AREA_LAB, D, "median_area", False, "results/shape/figures/FigD_data_per_fov.csv")
    three_group_sheet("FigD2", "Fig D2 - aspect ratio", "Aspect ratio", D, "median_aspect_ratio", False, "results/shape/figures/FigD_data_per_fov.csv")
if S3 is not None:
    src = "results/shape/figures/FigS3_data_per_fov.csv"
    three_group_sheet("FigS3a", "Fig S3a - condensates per field", "Condensates per field", S3, "n_condensates", False, src)
    three_group_sheet("FigS3b", "Fig S3b - area fraction", "Area fraction (%)", S3, "pct_area_covered", False, src)
    three_group_sheet("FigS3c", "Fig S3c - circularity (size-matched condensates)", "Circularity", S3, "circ_sizematched", False, src)
    three_group_sheet("FigS3d", "Fig S3d - solidity", "Solidity", S3, "median_solidity", False, src)
S4 = rd("shape", "figures", "FigS4_data_per_fov.csv")
if S4 is not None:
    r4 = S4.groupby(["condition", "cut", "replicate"], as_index=False)["y"].mean()
    import re
    cuts = sorted(r4["cut"].unique(), key=lambda s: float(re.search(r"[0-9.]+", s).group()))
    ws = new_sheet("FigS4", "Fig S4 - aspect ratio vs minimum condensate area (grouped table)", "Aspect ratio",
                   "Prism: Grouped table, 3 replicate values per cell.", "results/shape/figures/FigS4_data_per_fov.csv")
    cols = [(f"{g} {r}", [r4[(r4.condition == g) & (r4.cut == c) & (r4.replicate == r)]["y"].mean() for c in cuts]) for g in SG for r in reps]
    a, b, last = put_table(ws, 6, 1, "Minimum size", cuts, cols, summary=False)
    add_index("Fig S4", "FigS4", f"{a}:{b}", "aspect ratio per size cut-off", "descriptive")

# ============================ three-protein (Fig E1-E7) ============================
def two_vs_three(name, title, yaxis, file):
    d = rd("three_protein", "figures", file)
    if d is None:
        return
    ws = new_sheet(name, title, yaxis, UNIT, "results/three_protein/figures/" + file)
    left, rngs, ps, labs, sens = 1, [], [], [], []
    for cl in ["BCAN", "NCAN"]:
        ws.cell(row=5, column=left + 1, value=cl + " panel").font = F(bold=True)
        sub = d[d.client == cl]
        cols = wide(sub, "mix", ["Two proteins", "Three proteins"])
        cols = [("2 proteins (ACAN + " + cl + ")", cols[0][1]), ("3 proteins (ACAN + BCAN + NCAN)", cols[1][1])]
        a, b, last = put_table(ws, 6, left, "Preparation", reps, cols)
        ps.append(student(cols[0][1], cols[1][1])); rngs.append(f"{a}:{b} ({cl})"); labs.append(cl); sens.append(welch(cols[0][1], cols[1][1]))
        left += 4
    pa = ps
    put_stats(ws, last + 2, [f"Statistics: {TEST_T}, 2 vs 3 proteins, for each client"] +
              [f"{l}: {fmt_p(q)} (p = {p:.4g})" for l, p, q in zip(labs, ps, pa)] +
              ["Sensitivity (Welch's t-test, Holm over the 2 clients): " + "; ".join(f"{l} {fmt_p(q)}" for l, q in zip(labs, holm(sens)))])
    add_index(title.split(" - ")[0], name, "; ".join(rngs), yaxis, "Student's t-test: " + "; ".join(f"{l} {fmt_p(q)}" for l, q in zip(labs, pa)))


two_vs_three("FigE1", "Fig E1 - puncta within ACAN condensates, 2 vs 3 proteins", "Puncta within ACAN condensates (%)", "FigE1_data_per_fov.csv")
two_vs_three("FigE2", "Fig E2 - BCAN/NCAN-positive ACAN condensates, 2 vs 3 proteins", "BCAN/NCAN-positive ACAN condensates (%)", "FigE2_data_per_fov.csv")

E3 = rd("three_protein", "figures", "FigE3_data_means.csv")
if E3 is not None:
    ws = new_sheet("FigE3", "Fig E3 - position of BCAN/NCAN puncta, 2 vs 3 proteins (stacked bars)", "Fraction of BCAN/NCAN puncta (%) at each position relative to the ACAN condensate",
                   "Values plotted = mean of the 3 preparations.", "results/three_protein/figures/FigE3_data_means.csv")
    locs = ["Inside", "At edge", "Beyond edge"]
    cols = []
    for cl in ["BCAN", "NCAN"]:
        for mx, lab in [("Two proteins", "2 proteins"), ("Three proteins", "3 proteins")]:
            s = E3[(E3.client == cl) & (E3["mix"] == mx)].set_index("location")["pct"]
            cols.append((f"{cl}, {lab}", [s.get(l, np.nan) for l in locs]))
    a, b, last = put_table(ws, 6, 1, "Location", locs, cols, summary=False)
    add_index("Fig E3", "FigE3", f"{a}:{b} (rows 7-9)", "% inside / at edge / beyond edge", "descriptive (Inside % tested in Fig E1)")

E4 = rd("three_protein", "figures", "FigE4_data_per_fov_paired.csv")
if E4 is not None:
    ws = new_sheet("FigE4", "Fig E4a, b - BCAN vs NCAN in the same three-protein images (paired)", "(a) Relative ACAN level at client puncta; (b) client puncta inside ACAN (%)",
                   "Same row = same preparation (Prism: Column table, paired / repeated measures).", "results/three_protein/figures/FigE4_data_per_fov_paired.csv")
    left, rngs, ps, labs, sens = 1, [], [], [], []
    for m, lab in [("acan_level", "E4a ACAN level"), ("inside", "E4b % inside")]:
        ws.cell(row=5, column=left + 1, value=lab).font = F(bold=True)
        w = E4.pivot(index="fov", columns="client", values=m).reindex(reps)
        cols = [("BCAN", w["BCAN"].tolist()), ("NCAN", w["NCAN"].tolist())]
        a, b, last = put_table(ws, 6, left, "Preparation", reps, cols)
        ps.append(paired(cols[0][1], cols[1][1])); rngs.append(f"{a}:{b} ({lab})"); labs.append(lab)
        left += 4
    pa = ps
    put_stats(ws, last + 2, [f"Statistics: {TEST_PAIR}"] + [f"{l}: {fmt_p(q)} (p = {p:.4g})" for l, p, q in zip(labs, ps, pa)] +
              ["Sensitivity (paired t, Holm over the 2 measures): " + "; ".join(f"{l} {fmt_p(q)}" for l, q in zip(labs, holm(ps)))])
    add_index("Fig E4a, b", "FigE4", "; ".join(rngs), "BCAN vs NCAN in the same images", "paired t: " + "; ".join(f"{l} {fmt_p(q)}" for l, q in zip(labs, pa)))

E5 = rd("three_protein", "figures", "FigE5_data_per_fov.csv")
if E5 is not None:
    E5["replicate"] = E5["fov"].str.split("_").str[0]
    G4 = ["ACAN alone", "ACAN + BCAN", "ACAN + NCAN", "ACAN + BCAN + NCAN"]
    for nm, col, yl, lg in [("FigE5a", "median_area_plot" if "median_area_plot" in E5.columns else "median_area_px", AREA_LAB, False),
                            ("FigE5b", "median_aspect_ratio", "Aspect ratio", False)]:
        r = to_rep(E5.rename(columns={col: "y"}), ["context"])
        ws = new_sheet(nm, f"Fig {nm[3:]} - {yl.split(' (')[0].lower()} incl. three-protein mix", yl, UNIT, "results/three_protein/figures/FigE5_data_per_fov.csv")
        cols = wide(r, "context", G4)
        a, b, last = put_table(ws, 6, 1, "Preparation", reps, cols)
        pan, tk = anova_tukey([c[1] for c in cols], lg)
        allp = [f"{G4[i]} vs {G4[j]}: Tukey {fmt_p(tk[i, j])}" for i in range(4) for j in range(i + 1, 4)]
        put_stats(ws, last + 2, [f"Statistics: {TEST_ANOVA}" + (" on log10 values" if lg else "") + f" across the 4 groups; one-way ANOVA {fmt_p(pan)}",
                                 "Shown on the figure: ACAN + BCAN vs ACAN + BCAN + NCAN and ACAN + NCAN vs ACAN + BCAN + NCAN", "All pairs:"] + allp +
                                ["Sensitivity (Welch's t-test, 2 planned comparisons, Holm): " + "; ".join(f"{G4[i]} vs triple {fmt_p(q)}" for i, q in zip([1, 2], holm([welch(cols[1][1], cols[3][1], lg), welch(cols[2][1], cols[3][1], lg)])))])
        add_index(f"Fig {nm[3:]}", nm, f"{a}:{b}", yl, f"ANOVA + Tukey; ANOVA {fmt_p(pan)}; vs ACAN + BCAN {fmt_p(tk[1, 3])}; vs ACAN + NCAN {fmt_p(tk[2, 3])}")

E6 = rd("three_protein", "figures", "FigE6_data_per_replicate.csv")
if E6 is not None:
    ws = new_sheet("FigE6", "Fig E6 - NCAN in ACAN condensates of the same size", "NCAN-positive area of ACAN condensates (%)", UNIT,
                   "results/three_protein/figures/FigE6_data_per_replicate.csv")
    left, rngs, ps, labs, sens = 1, [], [], [], []
    for sz in list(dict.fromkeys(E6["size"])):
        ws.cell(row=5, column=left + 1, value="Condensates " + sz).font = F(bold=True)
        cols = wide(E6[E6["size"] == sz], "mix", ["ACAN + NCAN", "ACAN + BCAN + NCAN"])
        a, b, last = put_table(ws, 6, left, "Preparation", reps, cols)
        ps.append(student(cols[0][1], cols[1][1])); rngs.append(f"{a}:{b} ({sz})"); labs.append(sz); sens.append(welch(cols[0][1], cols[1][1]))
        left += 4
    pa = ps
    put_stats(ws, last + 2, [f"Statistics: {TEST_T} for each size class"] + [f"{l}: {fmt_p(q)} (p = {p:.4g})" for l, p, q in zip(labs, ps, pa)] +
              ["Sensitivity (Welch's t-test, Holm over the 2 sizes): " + "; ".join(f"{l} {fmt_p(q)}" for l, q in zip(labs, holm(sens)))])
    add_index("Fig E6", "FigE6", "; ".join(rngs), "NCAN-positive area, size-matched", "Student's t-test: " + "; ".join(f"{l} {fmt_p(q)}" for l, q in zip(labs, pa)))

E7 = rd("three_protein", "figures", "FigE7_data_per_replicate.csv")
if E7 is not None:
    ws = new_sheet("FigE7", "Fig E7 - NCAN-ACAN correlation across condensates", "Spearman's rho (NCAN vs ACAN; computed per FOV, averaged per preparation)", UNIT,
                   "results/three_protein/figures/FigE7_data_per_replicate.csv")
    cols = wide(E7, "mix", ["ACAN + NCAN", "ACAN + BCAN + NCAN"])
    a, b, last = put_table(ws, 6, 1, "Preparation", reps, cols)
    p = student(cols[0][1], cols[1][1])
    put_stats(ws, last + 2, ["Statistics", f"{TEST_T}: {fmt_p(p)} (p = {p:.4g})", f"Sensitivity, Welch's t-test: {fmt_p(welch(cols[0][1], cols[1][1]))}"])
    add_index("Fig E7", "FigE7", f"{a}:{b}", "Spearman rho", f"Student's t-test, {fmt_p(p)}")

# ============================ position (Fig F1-F3, S5-S7) ============================
F1 = rd("position", "figures", "FigF1_data_means.csv")
if F1 is not None:
    ws = new_sheet("FigF1", "Fig F1 - centre-to-edge intensity profiles", "Normalized intensity (1 = mean inside condensate)",
                   "X (blue) = position: 0 centre, 1 edge, > 1 outside. Y = mean of 3 preparations, with SEM and n (Prism: XY table, 'Mean, SEM, N').",
                   "results/position/figures/FigF1_data_means.csv")
    left, rngs = 1, []
    for mx in ["ACAN + BCAN", "ACAN + NCAN", "ACAN + BCAN + NCAN"]:
        sub = F1[F1["mix"] == mx]
        rho = sorted(sub["rho"].unique())
        ws.cell(row=5, column=left + 1, value=mx + " panel").font = F(bold=True)
        cols = [("Position (X)", [round(x, 2) for x in rho])]
        for pr in ["ACAN", "BCAN", "NCAN"]:
            s = sub[sub.protein == pr].set_index("rho")
            if not len(s):
                continue
            cols += [(f"{pr} mean", [s["m"].get(x, np.nan) for x in rho]), (f"{pr} SEM", [s["se"].get(x, np.nan) for x in rho]),
                     (f"{pr} N", [s["n"].get(x, np.nan) for x in rho])]
        a, b, last = put_table(ws, 6, left, "Bin", list(range(1, len(rho) + 1)), cols, xfill=BLUE, summary=False)
        rngs.append(f"{a}:{b} ({mx}; first column = X)")
        left += len(cols) + 2
    add_index("Fig F1", "FigF1", "; ".join(rngs), "profile curves (mean +/- SEM)", "descriptive")

F2 = rd("position", "figures", "FigF2_data_per_replicate.csv")
if F2 is not None:
    ws = new_sheet("FigF2", "Fig F2 - peripheral shift, two-protein mixes", "Peripheral shift relative to ACAN (%) = (radial ratio - 1) x 100",
                   UNIT + " Grey columns = the same values as radial ratio (used for the test).", "results/position/figures/FigF2_data_per_replicate.csv")
    cols = wide(F2, "mix", ["ACAN + BCAN", "ACAN + NCAN"], "shift_pct")
    a, b, last = put_table(ws, 6, 1, "Preparation", reps, cols)
    rr = wide(F2, "mix", ["ACAN + BCAN", "ACAN + NCAN"], "radial_ratio")
    put_table(ws, 6, 5, "Preparation", reps, [(g + " radial ratio", v) for g, v in rr], fill=GREY, summary=False)
    p = student(rr[0][1], rr[1][1])
    put_stats(ws, last + 2, ["Statistics", f"{TEST_T} on the radial ratio: {fmt_p(p)} (p = {p:.4g})", f"Sensitivity, Welch's t-test: {fmt_p(welch(rr[0][1], rr[1][1]))}"])
    add_index("Fig F2", "FigF2", f"{a}:{b}", "shift towards rim (%)", f"Student's t-test, {fmt_p(p)}")

F3 = rd("position", "figures", "FigF3_data_per_replicate_paired.csv")
if F3 is not None:
    ws = new_sheet("FigF3", "Fig F3 - BCAN vs NCAN in the same three-protein condensates (paired)", "Peripheral shift relative to ACAN (%)",
                   "Same row = same preparation (Prism: Column table, paired).", "results/position/figures/FigF3_data_per_replicate_paired.csv")
    F3 = F3.set_index("replicate").reindex(reps)
    cols = [("BCAN", F3["BCAN_shift_pct"].tolist()), ("NCAN", F3["NCAN_shift_pct"].tolist())]
    a, b, last = put_table(ws, 6, 1, "Preparation", reps, cols)
    p = paired(F3["NCAN_radial_ratio"], F3["BCAN_radial_ratio"])
    put_stats(ws, last + 2, ["Statistics", f"{TEST_PAIR}: {fmt_p(p)} (p = {p:.4g})"])
    add_index("Fig F3", "FigF3", f"{a}:{b}", "shift towards rim (%), paired", f"paired t, {fmt_p(p)}")

for nm, file, yl, title in [("FigS5", "FigS5_data_per_replicate.csv", "Puncta within ACAN (%)", "Fig S5 - inside ACAN: bright body vs body + dim rim"),
                            ("FigS6", "FigS6_data_per_replicate.csv", "ACAN intensity at puncta (normalized to host condensate)", "Fig S6 - B5 under two thresholds")]:
    d = rd("position", "figures", file)
    if d is None:
        continue
    d["boundary"] = d["boundary"].str.replace("\n", " ")
    ws = new_sheet(nm, title, yl, UNIT, "results/position/figures/" + file)
    left, rngs, ps, labs, sens = 1, [], [], [], []
    for bd in d["boundary"].unique():
        ws.cell(row=5, column=left + 1, value=bd).font = F(bold=True)
        cols = wide(d[d.boundary == bd], "client", ["BCAN", "NCAN"])
        a, b, last = put_table(ws, 6, left, "Preparation", reps, cols)
        ps.append(student(cols[0][1], cols[1][1])); rngs.append(f"{a}:{b} ({bd.split(' (')[0]})"); labs.append(bd); sens.append(welch(cols[0][1], cols[1][1]))
        left += 4
    pa = ps
    put_stats(ws, last + 2, [f"Statistics: {TEST_T} for each threshold"] + [f"{l}: {fmt_p(q)} (p = {p:.4g})" for l, p, q in zip(labs, ps, pa)] +
              ["Sensitivity (Welch's t-test, Holm over the 2 thresholds): " + "; ".join(f"{l.split(' (')[0]} {fmt_p(q)}" for l, q in zip(labs, holm(sens)))])
    add_index(title.split(" - ")[0], nm, "; ".join(rngs), yl, "Student's t-test: " + "; ".join(f"{l.split(' (')[0]} {fmt_p(q)}" for l, q in zip(labs, pa)))

S7 = rd("position", "figures", "FigS7_data_per_replicate.csv")
if S7 is not None:
    S7["analysis"] = S7["analysis"].str.replace("\n", " ")
    ws = new_sheet("FigS7", "Fig S7 - robustness of the peripheral shift", "Peripheral shift relative to ACAN (%)", UNIT,
                   "results/position/figures/FigS7_data_per_replicate.csv")
    left, rngs, lines = 1, [], ["Statistics: two-protein Student's t-test (BCAN vs NCAN); three-protein paired t-test (NCAN vs BCAN)"]
    groups = ["BCAN, 2 prot.", "NCAN, 2 prot.", "BCAN, 3 prot.", "NCAN, 3 prot."]
    for an in S7["analysis"].unique():
        sub = S7[S7.analysis == an]
        ws.cell(row=5, column=left + 1, value=an).font = F(bold=True)
        cols = wide(sub, "group", groups, "shift_pct")
        a, b, last = put_table(ws, 6, left, "Preparation", reps, cols)
        rr = dict(wide(sub, "group", groups, "radial_ratio"))
        if "Simulated" not in an:
            lines.append(f"{an}: two-protein {fmt_p(student(rr[groups[0]], rr[groups[1]]))}; three-protein {fmt_p(paired(rr[groups[3]], rr[groups[2]]))}")
        rngs.append(f"{a}:{b}")
        left += 6
    put_stats(ws, last + 2, lines)
    add_index("Fig S7", "FigS7", "; ".join(rngs) + " (one block per analysis variant)", "shift towards rim (%)", " | ".join(lines[1:]))

# ============================ index sheet ============================
ix = wb["Sheet"]; ix.title = "Index"
ix["A1"] = "Source data for every figure - which columns to use"; ix["A1"].font = F(bold=True, size=13)
ix["A2"] = ("Yellow columns = the values plotted (one column per group; one row per independent preparation unless noted). "
            "Blue column = X values (line plots). Grey = labels / extra information. Mean, SD, SEM and n under each table are live Excel formulas.")
ix["A2"].font = F(); ix["A2"].alignment = Alignment(wrap_text=True)
ix.merge_cells("A2:F2"); ix.row_dimensions[2].height = 45
hdr = ["Figure", "Sheet", "Columns to copy (yellow)", "What the values are", "Statistics on the figure", "Figure file"]
for j, h in enumerate(hdr, 1):
    c = ix.cell(row=4, column=j, value=h); c.font = F(bold=True); c.fill = GREY; c.border = BOX
figfiles = {}
for root, _, files in os.walk(RES):
    for f in files:
        if f.startswith("Fig") and f.endswith(".pdf") and "with_legend" not in f:
            figfiles.setdefault(f.split("_")[0], []).append(os.path.relpath(os.path.join(root, f), BASE).replace("\\", "/"))
for i, (fig, sheet, cols, what, test) in enumerate(index_rows, 5):
    key = sheet.replace("FigE4", "FigE4a")
    ff = "; ".join(sorted(figfiles.get(key, figfiles.get(sheet, []))))
    vals = [fig, sheet, cols, what, test, ff]
    for j, v in enumerate(vals, 1):
        c = ix.cell(row=i, column=j, value=v); c.font = F(); c.border = BOX; c.alignment = Alignment(wrap_text=True, vertical="top")
        if j == 3:
            c.fill = YELLOW
    ix.cell(row=i, column=2).hyperlink = f"#'{sheet}'!A1"; ix.cell(row=i, column=2).font = F(color="0563C1", underline="single")
for col, w in zip("ABCDEF", [12, 10, 38, 42, 30, 60]):
    ix.column_dimensions[col].width = w
ix.freeze_panes = "A5"
for ws in wb.worksheets:
    ws.sheet_view.showGridLines = ws.title != "Index"
wb.save(OUT)
print("Source-data workbook written:", OUT)

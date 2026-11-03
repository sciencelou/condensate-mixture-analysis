#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
recruitment_quant.py
====================
Quantify client recruitment into scaffold condensates (e.g. BCAN or NCAN into ACAN)
from two-channel fluorescence images.

Metric A  Partition coefficient (PC)
          PC_i = (mean client intensity inside scaffold condensate i - bg_client)
                 / (median client intensity in the dilute phase - bg_client)
          Only valid when the dilute phase is measurable (raw, non-clipped data).

Metric B  Object-based co-positivity
          B1  % scaffold condensates that are client-positive
              (intensity rule: mean client inside object > dilute mean + k * dilute SD)
          B2  % scaffold condensates that are client-positive
              (overlap rule: >= f of object area overlaps the client mask)
          B3  % client puncta that lie inside scaffold condensates
              (>= 50 % of client-object area inside the scaffold mask)
          B4  relative scaffold (ACAN) intensity at each client punctum, normalised to the median
              intensity of the image's scaffold condensates (1 = sits in a typical ACAN condensate;
              < 1 = sits in dimmer / peripheral ACAN material). Does not depend on the scaffold mask.
          plus M_client: fraction of total client signal (above bg) inside the scaffold mask
          (thresholded Manders coefficient).

Run in: cmd (Windows) / Terminal (macOS/Linux)
    python recruitment_quant.py --sheet sample_sheet.csv --out results

All parameters are written to <out>/params_used.json so that the run is reproducible.
Author: Bin Luo, 2026. Written with AI assistance (Claude, Anthropic).
"""
import argparse
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage import filters, measure, morphology, segmentation
from skimage.registration import phase_cross_correlation

warnings.filterwarnings("ignore", category=FutureWarning)

# ----------------------------------------------------------------------------------
# Default parameters (fixed BEFORE analysis; change only with justification and re-run
# every image with the same values)
# ----------------------------------------------------------------------------------
DEFAULTS = {
    "smooth_sigma_px": 1.0,          # Gaussian sigma used ONLY for segmentation; intensities are measured on unsmoothed data
    "scaffold_threshold": "otsu",    # otsu | li | triangle | <number>
    "client_threshold": "otsu",      # threshold for client puncta (used for B2, B3, dilute exclusion)
    "min_area_px": 9,                # objects smaller than this are discarded (sub-resolution)
    "exclude_border": True,          # drop objects touching the image edge (truncated)
    "split_touching": False,         # watershed split of touching condensates
    "dilute_exclusion_px": 3,        # dilate all condensate masks by this before defining the dilute phase (removes PSF tails)
    "erode_inside_px": 0,            # erode each object before measuring 'inside' intensity (0 = use full mask)
    "positivity_k_sd": 3.0,          # B1: client+ if mean_inside > dilute_mean + k * dilute_sd
    "positivity_min_abs": 0.0,       # B1: also require mean_inside - bg >= this (guards against SD ~ 0)
    "overlap_frac": 0.25,            # B2: client+ if >= this fraction of object area is in the client mask
    "client_inside_frac": 0.5,       # B3: client object counts as 'inside' if >= this fraction lies in scaffold mask
    "max_zero_frac_for_pc": 0.05,    # if > this fraction of dilute pixels are exactly 0 -> PC flagged invalid
    "correct_chromatic_shift": True,   # client channel is shifted by the sub-pixel offset measured by phase correlation
                                       # (ACAN+NCAN data showed a consistent ~1.2 px offset across all FOVs)
    "contact_dilate_px": 1,            # B3 categories: a client punctum <50 % inside ACAN is 'contacting' if, dilated by
                                       # this many px, it touches the ACAN mask; otherwise 'separate'
    "sensitivity_k_sd": [2, 3, 5],
    "sensitivity_overlap_frac": [0.10, 0.25, 0.50],
}


def load_channel(path):
    """Load a single-channel image. RGB exports are collapsed by taking the max over colour."""
    a = tifffile.imread(path)
    info = {"dtype": str(a.dtype), "raw_shape": list(a.shape)}
    a = np.squeeze(a)
    if a.ndim == 3 and a.shape[-1] in (3, 4):          # RGB(A) export
        a = a[..., :3].max(axis=-1)
        info["rgb_export"] = True
    elif a.ndim == 3:                                    # z-stack / multi-page: take max projection
        a = a.max(axis=0)
        info["projected"] = True
    else:
        info["rgb_export"] = False
    info["max_possible"] = float(np.iinfo(a.dtype).max) if np.issubdtype(a.dtype, np.integer) else np.nan
    return a.astype(np.float64), info


def get_threshold(img, rule):
    if isinstance(rule, (int, float)):
        return float(rule)
    try:
        return float(rule)
    except (TypeError, ValueError):
        pass
    rule = str(rule).lower()
    if rule == "otsu":
        return float(filters.threshold_otsu(img))
    if rule == "li":
        return float(filters.threshold_li(img))
    if rule == "triangle":
        return float(filters.threshold_triangle(img))
    raise ValueError(f"unknown threshold rule: {rule}")


def segment(img, p, rule):
    sm = filters.gaussian(img, p["smooth_sigma_px"], preserve_range=True)
    thr = get_threshold(sm, rule)
    mask = sm > thr
    mask = ndi.binary_fill_holes(mask)
    lab = measure.label(mask)
    if p["split_touching"]:
        dist = ndi.distance_transform_edt(mask)
        peaks = morphology.h_maxima(filters.gaussian(dist, 1), 1)
        lab = segmentation.watershed(-dist, measure.label(peaks), mask=mask)
    # size filter
    sizes = np.bincount(lab.ravel())
    small = sizes < p["min_area_px"]
    small[0] = False
    lab[small[lab]] = 0
    full_mask = lab > 0                      # mask kept for dilute exclusion (includes border objects)
    if p["exclude_border"]:
        lab = segmentation.clear_border(lab)
    lab, _, _ = segmentation.relabel_sequential(lab)
    return lab, full_mask, thr


def analyse_image(row, p, qc_dir):
    S, s_info = load_channel(row["scaffold_path"])
    C, c_info = load_channel(row["client_path"])
    if S.shape != C.shape:
        raise ValueError(f"{row['image_id']}: channel shapes differ {S.shape} vs {C.shape}")

    bg_s = float(row.get("bg_scaffold", 0) or 0)
    bg_c = float(row.get("bg_client", 0) or 0)
    px = row.get("pixel_um", np.nan)
    px = float(px) if px not in (None, "") and not pd.isna(px) else np.nan

    # chromatic shift between channels (sub-pixel, phase correlation)
    shift, _, _ = phase_cross_correlation(S, C, upsample_factor=10, normalization=None)
    if p["correct_chromatic_shift"]:
        C = ndi.shift(C, shift, order=1, mode="nearest")

    lab_s, full_s, thr_s = segment(S, p, p["scaffold_threshold"])
    lab_c, full_c, thr_c = segment(C, p, p["client_threshold"])

    # ---------------- dilute phase ----------------
    dense = morphology.dilation(full_s | full_c, morphology.disk(p["dilute_exclusion_px"]))
    dil = ~dense
    dC = C[dil]
    dS = S[dil]
    dil_c_mean, dil_c_med, dil_c_sd = float(dC.mean()), float(np.median(dC)), float(dC.std())
    dil_s_med = float(np.median(dS))
    zero_frac_c = float((dC == 0).mean())
    zero_frac_s = float((dS == 0).mean())
    pc_valid = (zero_frac_c <= p["max_zero_frac_for_pc"]) and (dil_c_med - bg_c > 0)

    sat_c = float((C >= c_info["max_possible"]).mean()) if not np.isnan(c_info["max_possible"]) else np.nan
    sat_s = float((S >= s_info["max_possible"]).mean()) if not np.isnan(s_info["max_possible"]) else np.nan

    # ---------------- per scaffold-condensate measurements ----------------
    objs = []
    for r in measure.regionprops(lab_s):
        rr, cc = r.coords[:, 0], r.coords[:, 1]
        if p["erode_inside_px"] > 0:
            sub = np.zeros(S.shape, bool)
            sub[rr, cc] = True
            sub = morphology.erosion(sub, morphology.disk(p["erode_inside_px"]))
            if sub.sum() >= 1:
                rr, cc = np.nonzero(sub)
        c_in = C[rr, cc]
        s_in = S[rr, cc]
        overlap = float(full_c[r.coords[:, 0], r.coords[:, 1]].mean())
        mean_c = float(c_in.mean())
        thr_pos = dil_c_mean + p["positivity_k_sd"] * dil_c_sd
        pos_B1 = (mean_c > thr_pos) and (mean_c - bg_c >= p["positivity_min_abs"])
        pc = (mean_c - bg_c) / (dil_c_med - bg_c) if (dil_c_med - bg_c) > 0 else np.inf
        perim = r.perimeter if r.perimeter > 0 else np.nan
        objs.append({
            "image_id": row["image_id"], "condition": row["condition"],
            "replicate": row["replicate"], "fov": row["fov"],
            "object_id": r.label,
            "area_px": int(r.area),
            "area_um2": r.area * px * px if not np.isnan(px) else np.nan,
            "eq_diam_um": r.equivalent_diameter_area * px if not np.isnan(px) else np.nan,
            "centroid_y": r.centroid[0], "centroid_x": r.centroid[1],
            "circularity": 4 * np.pi * r.area / perim ** 2 if perim == perim else np.nan,
            "aspect_ratio": r.major_axis_length / r.minor_axis_length if r.minor_axis_length > 0 else np.nan,
            "solidity": r.solidity,
            "scaffold_mean": float(s_in.mean()),
            "client_mean": mean_c,
            "client_max": float(c_in.max()),
            # client/scaffold intensity ratio inside the condensate ("how yellow is it");
            # comparable between conditions ONLY if acquisition + export settings are identical
            "client_to_scaffold_ratio": (mean_c - bg_c) / (float(s_in.mean()) - bg_s)
            if (float(s_in.mean()) - bg_s) > 0 else np.nan,
            "client_overlap_frac": overlap,
            "PC": pc if pc_valid else np.nan,
            "PC_raw_unvalidated": pc,
            "client_pos_B1": bool(pos_B1),
            "client_pos_B2": bool(overlap >= p["overlap_frac"]),
        })
    obj_df = pd.DataFrame(objs)

    # ---------------- per client-puncta measurements (B3) ----------------
    cobjs = []
    s_touch = morphology.dilation(full_s, morphology.disk(p["contact_dilate_px"]))
    # B4 reference: median scaffold intensity of the image's scaffold condensates
    s_ref = float(obj_df["scaffold_mean"].median()) - bg_s if len(obj_df) else np.nan
    # host condensate of each client punctum (all scaffold objects, incl. border-touching ones)
    host_lab = measure.label(full_s)
    host_mean = ndi.mean(S, host_lab, index=np.arange(1, host_lab.max() + 1)) if host_lab.max() else np.array([])
    for r in measure.regionprops(lab_c):
        inside = float(full_s[r.coords[:, 0], r.coords[:, 1]].mean())
        touching = bool(s_touch[r.coords[:, 0], r.coords[:, 1]].any())
        loc = "inside" if inside >= p["client_inside_frac"] else ("contacting" if touching else "separate")
        hosts = host_lab[r.coords[:, 0], r.coords[:, 1]]
        hosts = hosts[hosts > 0]
        s_at = float(S[r.coords[:, 0], r.coords[:, 1]].mean())
        if len(hosts):
            h = np.bincount(hosts).argmax()
            hm = float(host_mean[h - 1])
            rel_host = (s_at - bg_s) / (hm - bg_s) if hm - bg_s > 0 else np.nan
            host_rel = (hm - bg_s) / s_ref if s_ref and s_ref > 0 else np.nan
        else:
            rel_host, host_rel = np.nan, np.nan
        cobjs.append({
            "image_id": row["image_id"], "condition": row["condition"],
            "replicate": row["replicate"], "fov": row["fov"],
            "client_object_id": r.label, "area_px": int(r.area),
            "client_mean": float(C[r.coords[:, 0], r.coords[:, 1]].mean()),
            "scaffold_mean_at_client": float(S[r.coords[:, 0], r.coords[:, 1]].mean()),
            # B4: ACAN density where the client punctum sits, relative to a typical ACAN condensate (threshold-light)
            "rel_scaffold_at_client": (float(S[r.coords[:, 0], r.coords[:, 1]].mean()) - bg_s) / s_ref
            if s_ref and s_ref > 0 else np.nan,
            # B5: ACAN intensity at the punctum relative to the mean ACAN intensity of the condensate it
            # overlaps (1 = like the condensate average; < 1 = in the dimmer, outer part of its own condensate).
            # Independent of the image-wide reference used for B4. NaN for puncta that touch no condensate.
            "rel_scaffold_to_host": rel_host,
            "host_rel_scaffold": host_rel,
            "frac_inside_scaffold": inside,
            "inside_scaffold": bool(inside >= p["client_inside_frac"]),
            "location": loc,
        })
    cobj_df = pd.DataFrame(cobjs)

    # ---------------- Manders-type fraction ----------------
    Cb = np.clip(C - bg_c, 0, None)
    M_client = float(Cb[full_s].sum() / Cb.sum()) if Cb.sum() > 0 else np.nan

    # ---------------- sensitivity analysis ----------------
    sens = {}
    if len(obj_df):
        for k in p["sensitivity_k_sd"]:
            t = dil_c_mean + k * dil_c_sd
            sens[f"pctB1_k{k:g}"] = 100 * float((obj_df["client_mean"] > t).mean())
        for f in p["sensitivity_overlap_frac"]:
            sens[f"pctB2_f{f:g}"] = 100 * float((obj_df["client_overlap_frac"] >= f).mean())

    n = len(obj_df)
    summ = {
        "image_id": row["image_id"], "condition": row["condition"],
        "replicate": row["replicate"], "fov": row["fov"],
        "n_scaffold_condensates": n,
        "n_client_puncta": len(cobj_df),
        "pct_scaffold_clientpos_B1": 100 * obj_df["client_pos_B1"].mean() if n else np.nan,
        "pct_scaffold_clientpos_B2": 100 * obj_df["client_pos_B2"].mean() if n else np.nan,
        "pct_client_inside_scaffold_B3": 100 * cobj_df["inside_scaffold"].mean() if len(cobj_df) else np.nan,
        "pct_client_contacting_scaffold": 100 * (cobj_df["location"] == "contacting").mean() if len(cobj_df) else np.nan,
        "pct_client_separate": 100 * (cobj_df["location"] == "separate").mean() if len(cobj_df) else np.nan,
        "B4_median_rel_scaffold_at_client": float(cobj_df["rel_scaffold_at_client"].median()) if len(cobj_df) else np.nan,
        "B5_median_rel_scaffold_to_host": float(cobj_df["rel_scaffold_to_host"].median()) if len(cobj_df) and cobj_df["rel_scaffold_to_host"].notna().any() else np.nan,
        "M_client_frac_signal_in_scaffold": M_client,
        "PC_median": float(obj_df["PC"].median()) if n and pc_valid else np.nan,
        "PC_mean": float(obj_df["PC"].mean()) if n and pc_valid else np.nan,
        "PC_valid": pc_valid,
        "dilute_client_mean": dil_c_mean, "dilute_client_median": dil_c_med,
        "dilute_client_sd": dil_c_sd, "dilute_client_zero_frac": zero_frac_c,
        "dilute_scaffold_median": dil_s_med, "dilute_scaffold_zero_frac": zero_frac_s,
        "B1_threshold": dil_c_mean + p["positivity_k_sd"] * dil_c_sd,
        "scaffold_seg_threshold": thr_s, "client_seg_threshold": thr_c,
        "chromatic_shift_y_px": float(shift[0]), "chromatic_shift_x_px": float(shift[1]),
        "saturated_frac_scaffold": sat_s, "saturated_frac_client": sat_c,
        "dtype": s_info["dtype"], "rgb_export": bool(s_info.get("rgb_export", False)),
        "pixel_um": px,
    }
    summ.update(sens)

    # ---------------- warnings ----------------
    warn = []
    if s_info.get("rgb_export") or c_info.get("rgb_export"):
        warn.append("RGB display export: intensities may be contrast-stretched; use raw single-channel data for PC.")
    if s_info["dtype"] == "uint8":
        warn.append("8-bit data: limited dynamic range.")
    if not pc_valid:
        warn.append(f"PC INVALID: {zero_frac_c:.0%} of dilute client pixels are exactly 0 "
                    f"(black level clipped) -> denominator not measurable.")
    if (sat_c == sat_c and sat_c > 0.001) or (sat_s == sat_s and sat_s > 0.001):
        warn.append(f"Saturation: scaffold {sat_s:.3%}, client {sat_c:.3%} of pixels at max value.")
    if np.hypot(*shift) > 0.5:
        warn.append(f"Chromatic shift {shift[0]:.2f}, {shift[1]:.2f} px (y, x) between channels.")
    if np.isnan(px):
        warn.append("No pixel size given: areas reported in px only.")
    summ["warnings"] = " | ".join(warn)

    # ---------------- QC overlay ----------------
    save_qc(S, C, lab_s, full_c, obj_df, row["image_id"], qc_dir)
    return obj_df, cobj_df, summ


def save_qc(S, C, lab_s, mask_c, obj_df, image_id, qc_dir):
    """QC image: left = merge, right = merge + outlines (cyan = client+ by B1, white = client-, magenta = client mask)."""
    from PIL import Image
    def norm(x):
        lo, hi = np.percentile(x, [0.5, 99.8])
        return np.clip((x - lo) / (hi - lo + 1e-9), 0, 1)
    rgb = np.zeros(S.shape + (3,))
    rgb[..., 0] = norm(C)          # client in red (magenta-safe alternative possible)
    rgb[..., 1] = norm(S)          # scaffold in green
    base = (rgb * 255).astype(np.uint8)
    out = base.copy()
    pos_ids = set(obj_df.loc[obj_df["client_pos_B1"], "object_id"]) if len(obj_df) else set()
    bnd_c = segmentation.find_boundaries(mask_c, mode="outer")
    out[bnd_c] = [255, 0, 255]
    for lab_id in np.unique(lab_s)[1:]:
        b = segmentation.find_boundaries(lab_s == lab_id, mode="outer")
        out[b] = [0, 255, 255] if lab_id in pos_ids else [255, 255, 255]
    img = np.concatenate([base, np.full((S.shape[0], 6, 3), 128, np.uint8), out], axis=1)
    img = np.kron(img, np.ones((2, 2, 1), np.uint8))  # 2x upscale for viewing
    Image.fromarray(img).save(os.path.join(qc_dir, f"QC_{image_id}.png"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True, help="sample sheet CSV")
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--params", default=None, help="optional JSON file overriding DEFAULTS")
    args = ap.parse_args()

    p = dict(DEFAULTS)
    if args.params:
        with open(args.params, encoding="utf-8") as fh:
            p.update(json.load(fh))

    sheet = pd.read_csv(args.sheet)
    need = {"image_id", "condition", "replicate", "fov", "scaffold_path", "client_path"}
    miss = need - set(sheet.columns)
    if miss:
        sys.exit(f"sample sheet missing columns: {miss}")
    base = os.path.dirname(os.path.abspath(args.sheet))
    for col in ("scaffold_path", "client_path"):
        sheet[col] = [x if os.path.isabs(x) else os.path.join(base, x) for x in sheet[col]]

    os.makedirs(args.out, exist_ok=True)
    qc_dir = os.path.join(args.out, "qc")
    os.makedirs(qc_dir, exist_ok=True)

    all_obj, all_cobj, all_summ = [], [], []
    for _, row in sheet.iterrows():
        print(f"[{row['image_id']}] ...", flush=True)
        o, c, s = analyse_image(row, p, qc_dir)
        all_obj.append(o); all_cobj.append(c); all_summ.append(s)
        if s["warnings"]:
            print("   WARNING:", s["warnings"])

    obj = pd.concat(all_obj, ignore_index=True)
    cobj = pd.concat(all_cobj, ignore_index=True)
    summ = pd.DataFrame(all_summ)
    obj.to_csv(os.path.join(args.out, "objects_scaffold_condensates.csv"), index=False)
    cobj.to_csv(os.path.join(args.out, "objects_client_puncta.csv"), index=False)
    summ.to_csv(os.path.join(args.out, "image_summary.csv"), index=False)

    # replicate-level means (the statistical unit) -> computed from FOV-level values
    num = [c for c in summ.columns if c.startswith(("pct", "M_client", "PC_median", "n_", "B4_", "B5_"))]
    rep = summ.groupby(["condition", "replicate"])[num].mean().reset_index()
    rep["n_fov"] = summ.groupby(["condition", "replicate"]).size().values
    rep.to_csv(os.path.join(args.out, "replicate_summary.csv"), index=False)

    with open(os.path.join(args.out, "params_used.json"), "w", encoding="utf-8") as fh:
        json.dump(p, fh, indent=2)
    print("\nDone. Results in", os.path.abspath(args.out))
    cols = ["image_id", "n_scaffold_condensates", "pct_scaffold_clientpos_B1",
            "pct_scaffold_clientpos_B2", "pct_client_inside_scaffold_B3", "B4_median_rel_scaffold_at_client",
            "M_client_frac_signal_in_scaffold", "PC_median", "PC_valid"]
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(summ[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()

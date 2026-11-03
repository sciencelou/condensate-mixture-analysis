#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
radial_profile.py  -  WHERE inside an ACAN condensate does the client (BCAN / NCAN) sit: core or rim?

Idea (threshold-light; the global Otsu/Li threshold is used ONLY to find the condensates)
  For every ACAN condensate its own edge is defined from its own ACAN signal:
      edge = the contour where ACAN falls to 50 % of the way between the dilute background and the
             condensate's ACAN peak (95th percentile)  -> "half-maximum edge" (per condensate).
  Every pixel then gets a relative position
      rho = 0 at the deepest point (centre), rho = 1 on the half-maximum edge, rho > 1 just outside
      (outside pixels: rho = 1 + distance-outside / max-depth; up to 3 px or 0.5 x max-depth outside).
  Pixels closer to a different condensate are excluded (each pixel belongs to its nearest condensate).

Outputs per condensate
  mean_rho_scaffold, mean_rho_client : signal-weighted mean position (weights = intensity - dilute bg, >= 0)
  radial_ratio = mean_rho_client / mean_rho_scaffold
       = 1   client distributed like ACAN
       > 1   client shifted to the rim / outside   (peripheral)
       < 1   client concentrated in the core
  frac_client_outside_core: fraction of client signal at rho > 0.7 minus the same fraction for ACAN
  binned profiles (rho bins) of both channels, each normalised to its own mean inside the edge.

Run in: cmd (Windows) / Terminal (macOS/Linux)
    python radial_profile.py --sheet sample_table.csv --out results/radial [--params params.json]
Author: Bin Luo, 2026. Written with AI assistance (Claude, Anthropic).
"""
import argparse
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage import filters, measure, morphology
from skimage.registration import phase_cross_correlation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from recruitment_quant import DEFAULTS as RQ_DEFAULTS, load_channel, segment  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)

DEFAULTS = dict(RQ_DEFAULTS)
DEFAULTS.update({
    "edge_frac": 0.5,            # per-condensate edge = bg + edge_frac * (peak - bg)
    "peak_percentile": 95,       # condensate 'peak' ACAN = this percentile of smoothed ACAN inside the object
    "outside_px": 3,             # how far outside the edge pixels are still counted (px; also capped at 0.5 x max depth)
    "min_depth_px": 2.0,         # condensates whose half-max region is thinner than this are skipped (too small)
    "core_rho": 0.7,             # rho > core_rho = 'outer shell + rim'
    "rho_bins": [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5],
    "size_window_px": [100, 400],
    "null_reps": 20,             # repeats of the SNR-matched null per condensate  # size-matched subset (same window as the shape analysis)
})


def analyse_image(row, p):
    S, _ = load_channel(row["scaffold_path"])
    C, _ = load_channel(row["client_path"])
    shift, _, _ = phase_cross_correlation(S, C, upsample_factor=10, normalization=None)
    if p["correct_chromatic_shift"]:
        C = ndi.shift(C, shift, order=1, mode="nearest")

    lab_s, full_s, _ = segment(S, p, p["scaffold_threshold"])
    lab_c, full_c, _ = segment(C, p, p["client_threshold"])
    dense = morphology.dilation(full_s | full_c, morphology.disk(p["dilute_exclusion_px"]))
    bgS, bgC = float(np.median(S[~dense])), float(np.median(C[~dense]))
    dC = C[~dense]
    noise_pool = dC - bgC                      # dilute-phase client residuals (for the SNR-matched null)
    rng = np.random.default_rng(2026)
    pos_thr = float(dC.mean()) + p["positivity_k_sd"] * float(dC.std())

    # nearest-condensate map (all condensates incl. border ones) -> each pixel belongs to one condensate
    all_lab = measure.label(full_s)
    _, (iy, ix) = ndi.distance_transform_edt(all_lab == 0, return_indices=True)
    owner = all_lab[iy, ix]

    Ssm = filters.gaussian(S, p["smooth_sigma_px"], preserve_range=True)
    edges = np.asarray(p["rho_bins"], float)
    mids = (edges[:-1] + edges[1:]) / 2
    objs, prof = [], []
    for r in measure.regionprops(lab_s):
        rr, cc = r.coords[:, 0], r.coords[:, 1]
        own_id = np.bincount(all_lab[rr, cc]).argmax()
        mean_c = float(C[rr, cc].mean())
        pad = int(max(8, r.equivalent_diameter_area))
        y0, x0 = max(r.bbox[0] - pad, 0), max(r.bbox[1] - pad, 0)
        y1, x1 = min(r.bbox[2] + pad, S.shape[0]), min(r.bbox[3] + pad, S.shape[1])
        s_w, c_w, sm_w = S[y0:y1, x0:x1], C[y0:y1, x0:x1], Ssm[y0:y1, x0:x1]
        mine = owner[y0:y1, x0:x1] == own_id
        obj_mask = np.zeros_like(mine)
        obj_mask[rr - y0, cc - x0] = True
        peak = float(np.percentile(sm_w[obj_mask], p["peak_percentile"]))
        if peak - bgS <= 0:
            continue
        hm = (sm_w > bgS + p["edge_frac"] * (peak - bgS)) & mine
        hm = ndi.binary_fill_holes(hm)
        hl = measure.label(hm)
        ids = hl[obj_mask]
        ids = ids[ids > 0]
        if not len(ids):
            continue
        hm = hl == np.bincount(ids).argmax()
        depth_in = ndi.distance_transform_edt(hm)
        dmax = float(depth_in.max())
        if dmax < p["min_depth_px"]:
            continue
        d_out = ndi.distance_transform_edt(~hm)
        out_lim = min(p["outside_px"], 0.5 * dmax)
        rho = np.where(hm, 1 - (depth_in - 1) / max(dmax - 1, 1e-9), 1 + d_out / dmax)
        use = mine & (hm | (d_out <= out_lim))
        rho_u = rho[use]
        ws = np.clip(s_w[use] - bgS, 0, None)
        wc = np.clip(c_w[use] - bgC, 0, None)
        if ws.sum() <= 0 or wc.sum() <= 0:
            mrs = mrc = np.nan
        else:
            mrs = float((rho_u * ws).sum() / ws.sum())
            mrc = float((rho_u * wc).sum() / wc.sum())
        # SNR-matched null: a client that is distributed EXACTLY like ACAN, scaled to this client's total
        # signal, plus noise drawn from this image's dilute-phase client pixels. Any ratio > 1 in the null
        # is bias from low signal-to-noise (clipped noise spreads over the larger outer area).
        null_ratio = np.nan
        if ws.sum() > 0 and wc.sum() > 0 and mrs == mrs:
            vals = []
            for _ in range(p["null_reps"]):
                sim = (s_w[use] - bgS) * (np.clip(c_w[use] - bgC, 0, None).sum() / np.clip(s_w[use] - bgS, 0, None).sum())
                sim = np.clip(sim + rng.choice(noise_pool, size=sim.size), 0, None)
                if sim.sum() > 0:
                    vals.append(float((rho_u * sim).sum() / sim.sum()) / mrs)
            null_ratio = float(np.mean(vals)) if vals else np.nan
        shell = rho_u > p["core_rho"]
        fs = float(ws[shell].sum() / ws.sum()) if ws.sum() > 0 else np.nan
        fc = float(wc[shell].sum() / wc.sum()) if wc.sum() > 0 else np.nan
        objs.append({
            "image_id": row["image_id"], "condition": row["condition"], "replicate": row["replicate"],
            "fov": row["fov"], "object_id": r.label, "area_px": int(r.area),
            "halfmax_area_px": int(hm.sum()), "max_depth_px": dmax,
            "client_pos_B1": bool(mean_c > pos_thr),
            "mean_rho_scaffold": mrs, "mean_rho_client": mrc,
            "radial_ratio": mrc / mrs if mrs and mrs == mrs and mrs > 0 else np.nan,
            "radial_ratio_null": null_ratio,
            "shell_frac_scaffold": fs, "shell_frac_client": fc,
            "shell_frac_diff": fc - fs if fc == fc and fs == fs else np.nan,
        })
        # binned profile, each channel normalised to its own mean inside the edge
        inside = rho_u <= 1
        ns = ws[inside].mean() if inside.any() else np.nan
        nc = wc[inside].mean() if inside.any() else np.nan
        b = np.digitize(rho_u, edges) - 1
        for k, m in enumerate(mids):
            sel = b == k
            prof.append({"image_id": row["image_id"], "condition": row["condition"], "replicate": row["replicate"],
                         "fov": row["fov"], "object_id": r.label, "area_px": int(r.area),
                         "client_pos_B1": bool(mean_c > pos_thr), "rho": m,
                         "scaffold_norm": ws[sel].mean() / ns if sel.any() and ns > 0 else np.nan,
                         "client_norm": wc[sel].mean() / nc if sel.any() and nc > 0 else np.nan})
    return pd.DataFrame(objs), pd.DataFrame(prof), shift


def summarise(obj, p, label_cols=("condition", "replicate", "fov", "image_id")):
    lo, hi = p["size_window_px"]
    rows = []
    for key, g in obj.groupby(list(label_cols)):
        d = dict(zip(label_cols, key))
        gp = g[g["client_pos_B1"] & g["radial_ratio"].notna()]
        gs = gp[(gp["area_px"] >= lo) & (gp["area_px"] <= hi)]
        d.update({
            "n_condensates_measured": len(g), "n_clientpos": len(gp), "n_clientpos_sizematched": len(gs),
            "radial_ratio_median": gp["radial_ratio"].median() if len(gp) else np.nan,
            "radial_ratio_median_sizematched": gs["radial_ratio"].median() if len(gs) else np.nan,
            "radial_ratio_null_median": gp["radial_ratio_null"].median() if len(gp) else np.nan,
            "radial_excess_median": (gp["radial_ratio"] - gp["radial_ratio_null"]).median() if len(gp) else np.nan,
            "shell_frac_diff_median": gp["shell_frac_diff"].median() if len(gp) else np.nan,
        })
        rows.append(d)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--params", default=None)
    a = ap.parse_args()
    p = dict(DEFAULTS)
    if a.params:
        with open(a.params, encoding="utf-8") as fh:
            p.update(json.load(fh))
    sheet = pd.read_csv(a.sheet)
    os.makedirs(a.out, exist_ok=True)
    objs, profs = [], []
    for _, row in sheet.iterrows():
        print(f"[radial] {row['image_id']}", flush=True)
        o, pr, _ = analyse_image(row, p)
        objs.append(o); profs.append(pr)
    obj = pd.concat(objs, ignore_index=True)
    prof = pd.concat(profs, ignore_index=True)
    obj.to_csv(os.path.join(a.out, "radial_objects.csv"), index=False)
    prof.to_csv(os.path.join(a.out, "radial_profiles_binned.csv"), index=False)
    fov = summarise(obj, p)
    fov.to_csv(os.path.join(a.out, "radial_fov_summary.csv"), index=False)
    num = [c for c in fov.columns if c.startswith(("n_", "radial_", "shell_"))]
    rep = fov.groupby(["condition", "replicate"])[num].mean().reset_index()
    rep["n_fov"] = fov.groupby(["condition", "replicate"]).size().values
    rep.to_csv(os.path.join(a.out, "radial_replicate_summary.csv"), index=False)
    with open(os.path.join(a.out, "params_used.json"), "w", encoding="utf-8") as fh:
        json.dump(p, fh, indent=2)
    with pd.option_context("display.width", 200):
        print(rep.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

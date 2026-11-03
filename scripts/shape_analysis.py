#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
shape_analysis.py
=================
Size and shape of ACAN condensates (ACAN channel only) in ACAN alone vs ACAN+BCAN vs ACAN+NCAN.
Uses exactly the same segmentation as recruitment_quant.py (imported from it).

Called automatically by run_all.py.
    python scripts\\shape_analysis.py --sheet shape_sheet.csv --out results_shape

Per condensate : area, equivalent diameter, circularity (4*pi*area/perimeter^2), aspect ratio
                 (major/minor axis), solidity (area/convex-hull area)
Per FOV        : number of condensates, % image area covered, medians of the per-condensate values
Shape metrics are computed only for condensates >= shape_min_area_px, because the outline of very
small objects is dominated by the pixel grid. Sensitivity runs use other cut-offs.
Also writes a gallery of randomly chosen condensates (fixed random seed) for the high-mag panel.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from skimage import measure

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recruitment_quant as rq  # noqa: E402

SHAPE_DEFAULTS = {
    "shape_min_area_px": 50,            # primary cut-off for shape metrics (~8 px diameter)
    "sensitivity_min_area_px": [25, 50, 100],
    "gallery_n_per_condition": 8,
    "gallery_crop_px": 40,
    "gallery_seed": 2026,
}


def analyse(row, p):
    S, info = rq.load_channel(row["scaffold_path"])
    lab, full, thr = rq.segment(S, p, p["scaffold_threshold"])
    px = row.get("pixel_um", np.nan)
    px = float(px) if not pd.isna(px) and px != "" else np.nan
    objs = []
    for r in measure.regionprops(lab):
        per = r.perimeter
        objs.append({
            "image_id": row["image_id"], "condition": row["condition"],
            "replicate": row["replicate"], "fov": row["fov"], "object_id": r.label,
            "area_px": int(r.area),
            "area_um2": r.area * px * px if px == px else np.nan,
            "eq_diam_px": r.equivalent_diameter_area,
            "circularity": min(4 * np.pi * r.area / per ** 2, 1.0) if per > 0 else np.nan,
            "aspect_ratio": r.major_axis_length / r.minor_axis_length if r.minor_axis_length > 0 else np.nan,
            "solidity": r.solidity,
            "centroid_y": r.centroid[0], "centroid_x": r.centroid[1],
            "scaffold_mean": float(S[r.coords[:, 0], r.coords[:, 1]].mean()),
            "img_h": S.shape[0], "img_w": S.shape[1],
        })
    o = pd.DataFrame(objs)
    big = o[o.area_px >= p["shape_min_area_px"]]
    summ = {
        "image_id": row["image_id"], "condition": row["condition"],
        "replicate": row["replicate"], "fov": row["fov"],
        "n_condensates": len(o),
        "n_condensates_for_shape": len(big),
        "pct_area_covered": 100 * float(full.mean()),
        "median_area_px": float(o.area_px.median()) if len(o) else np.nan,
        "median_eq_diam_px": float(o.eq_diam_px.median()) if len(o) else np.nan,
        "median_circularity": float(big.circularity.median()) if len(big) else np.nan,
        "median_aspect_ratio": float(big.aspect_ratio.median()) if len(big) else np.nan,
        "median_solidity": float(big.solidity.median()) if len(big) else np.nan,
        "seg_threshold": thr, "pixel_um": px,
    }
    for c in p["sensitivity_min_area_px"]:
        b = o[o.area_px >= c]
        summ[f"median_circularity_min{c}"] = float(b.circularity.median()) if len(b) else np.nan
        summ[f"median_aspect_ratio_min{c}"] = float(b.aspect_ratio.median()) if len(b) else np.nan
    return o, summ


COLOUR_RGB = {"green": (0, 1, 0), "red": (1, 0, 0), "magenta": (1, 0, 1), "blue": (0, 0, 1),
              "cyan": (0, 1, 1), "yellow": (1, 1, 0), "grey": (1, 1, 1), "gray": (1, 1, 1)}


def as_rgb(path, fallback="grey"):
    """RGB display image for the gallery. Display exports are used as they are; raw single-channel
    images (from Fiji / Bio-Formats) are coloured by the colour word in their file name
    (..._ch00_red.tif) and scaled 0.1-99.9 percentile -> 0-255 for DISPLAY ONLY."""
    import re
    a = np.squeeze(rq.tifffile.imread(path))
    if a.ndim == 3 and a.shape[-1] in (3, 4):
        return a[..., :3].astype(np.uint8)
    if a.ndim == 3:
        a = a.max(axis=0)
    m = re.search(r"_(red|green|blue|magenta|cyan|yellow|grey|gray)\.tif", os.path.basename(path), flags=re.I)
    col = COLOUR_RGB[m.group(1).lower()] if m else COLOUR_RGB[fallback]
    lo, hi = np.percentile(a, [0.1, 99.9])
    g = np.clip((a.astype(float) - lo) / (hi - lo + 1e-9), 0, 1)
    return (np.stack([g * c for c in col], axis=-1) * 255).astype(np.uint8)


def gallery(sheet, objs, p, out_dir):
    """Random condensates (>= shape cut-off) per condition, same crop size, ACAN channel
    (merged with client channel when present). Saved as PNG + TIFF (raster images)."""
    from PIL import Image, ImageDraw, ImageFont
    import tifffile
    rng = np.random.default_rng(p["gallery_seed"])
    half = p["gallery_crop_px"] // 2
    rows_img, picked = [], []
    for cond in sheet.condition.unique():
        cand = objs[(objs.condition == cond) & (objs.area_px >= p["shape_min_area_px"])]
        # keep only objects whose crop fits inside the image
        cand = cand[(cand.centroid_y >= half) & (cand.centroid_y < cand.img_h - half) &
                    (cand.centroid_x >= half) & (cand.centroid_x < cand.img_w - half)]
        sel = cand.sample(n=min(p["gallery_n_per_condition"], len(cand)),
                          random_state=int(rng.integers(1e9)))
        tiles = []
        for _, o in sel.iterrows():
            r = sheet[sheet.image_id == o.image_id].iloc[0]
            rgb = as_rgb(r.scaffold_path, "green")                         # ACAN
            if "client_path" in r and isinstance(r.client_path, str) and r.client_path:
                rgb = np.maximum(rgb, as_rgb(r.client_path, "red" if "BCAN" in cond else "magenta"))
            if "client2_path" in r and isinstance(r.client2_path, str) and r.client2_path:
                rgb = np.maximum(rgb, as_rgb(r.client2_path, "magenta"))     # NCAN in the three-protein mix
            y, x = int(round(o.centroid_y)), int(round(o.centroid_x))
            crop = rgb[y - half:y + half, x - half:x + half, :3]
            tiles.append(np.kron(crop, np.ones((4, 4, 1), np.uint8)))   # 4x nearest-neighbour enlargement
            picked.append({"condition": cond, "image_id": o.image_id, "object_id": o.object_id,
                           "area_px": o.area_px, "circularity": o.circularity, "aspect_ratio": o.aspect_ratio})
        gap = np.full((tiles[0].shape[0], 6, 3), 255, np.uint8)
        row_img = np.concatenate(sum([[t, gap] for t in tiles], [])[:-1], axis=1)
        rows_img.append(row_img)
    # rows in the same order as all other figures (sampling order above unchanged -> same tiles)
    canon = ["ACAN alone", "ACAN+BCAN", "ACAN+NCAN", "ACAN+BCAN+NCAN"]
    conds = list(sheet.condition.unique())
    order = sorted(range(len(conds)), key=lambda k: canon.index(conds[k]) if conds[k] in canon else 99)
    rows_img = [rows_img[k] for k in order]
    conds = [conds[k] for k in order]
    w = max(r.shape[1] for r in rows_img)
    rows_img = [np.pad(r, ((0, 0), (0, w - r.shape[1]), (0, 0)), constant_values=255) for r in rows_img]
    vgap = np.full((10, w, 3), 255, np.uint8)
    img = np.concatenate(sum([[r, vgap] for r in rows_img], [])[:-1], axis=0)
    # white label column on the left + scale bar under the first tile
    font = None
    for f in ("arial.ttf", "Arial.ttf", "DejaVuSans.ttf", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        try:
            font = ImageFont.truetype(f, 22); break
        except OSError:
            continue
    font = font or ImageFont.load_default()
    lab_w, bottom = 215, 40
    canvas = Image.new("RGB", (img.shape[1] + lab_w, img.shape[0] + bottom), "white")
    canvas.paste(Image.fromarray(img), (lab_w, 0))
    d = ImageDraw.Draw(canvas)
    tile_h = rows_img[0].shape[0]
    for i, cond in enumerate(conds):
        label = cond.replace("+", " + ")
        if len(label) > 12:                       # e.g. "ACAN + BCAN + NCAN" -> two lines
            parts = label.split(" + ")
            label = parts[0] + "\n+ " + " + ".join(parts[1:])
        n_lines = label.count("\n") + 1
        y = i * (tile_h + 10) + tile_h // 2 - 13 * n_lines
        d.multiline_text((8, y), label, fill=(0, 0, 0), font=font, spacing=4)
    px_um = sheet["pixel_um"].dropna().astype(float)
    px_um = float(px_um.iloc[0]) if len(px_um) else np.nan
    if px_um == px_um and px_um > 0:               # scale bar of 1 um when the pixel size is known
        bar, bar_lab = int(round(1.0 / px_um * 4)), "1 \u00b5m"
    else:
        bar, bar_lab = 10 * 4, "10 px"             # 10 image px, displayed 4x enlarged
    d.rectangle([lab_w, img.shape[0] + 12, lab_w + bar, img.shape[0] + 17], fill=(0, 0, 0))
    d.text((lab_w + bar + 8, img.shape[0] + 4), bar_lab, fill=(0, 0, 0), font=font)
    fig_dir = os.path.join(out_dir, "figures"); os.makedirs(fig_dir, exist_ok=True)
    canvas.save(os.path.join(fig_dir, "FigD3_gallery_random_condensates.png"))
    tifffile.imwrite(os.path.join(fig_dir, "FigD3_gallery_random_condensates.tif"), np.array(canvas))
    sel_df = pd.DataFrame(picked)
    sel_df.to_csv(os.path.join(fig_dir, "FigD3_data_selected_condensates.csv"), index=False)
    n = p["gallery_n_per_condition"]
    legend = ("Figure D3 (FigD3_gallery_random_condensates.png/.tif)\n\n"
              f"Representative ACAN condensates at high magnification. In each row, {n} condensates of at least "
              f"{p['shape_min_area_px']} px were chosen at random by the analysis script (fixed random seed {p['gallery_seed']}), "
              "so the examples are not hand-picked. Green, ACAN; red, BCAN; magenta, NCAN; yellow/white, overlap. "
              f"All tiles are {p['gallery_crop_px']} x {p['gallery_crop_px']} px crops shown at the same scale. "
              + (f"Each tile is {p['gallery_crop_px'] * px_um:.2f} x {p['gallery_crop_px'] * px_um:.2f} \u00b5m; scale bar, 1 \u00b5m. " if px_um == px_um else "Scale bar, 10 px. ") +
              ("Raw intensities, scaled for display only (0.1-99.9th percentile of each channel in each image); "
                 if np.squeeze(rq.tifffile.imread(sheet.iloc[0].scaffold_path)).ndim == 2 else
                 "Display settings as exported from the microscope software; ") +
              "channels are shown without the ~1 px inter-channel offset correction used for quantification. "
              "The identity of each condensate is listed in FigD3_data_selected_condensates.csv.")
    with open(os.path.join(fig_dir, "FigD3_gallery_legend.txt"), "w", encoding="utf-8") as fh:
        fh.write(legend)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--params", default=None)
    a = ap.parse_args()
    p = dict(rq.DEFAULTS); p.update(SHAPE_DEFAULTS)
    if a.params:
        with open(a.params, encoding="utf-8") as fh:
            p.update(json.load(fh))
    sheet = pd.read_csv(a.sheet)
    base = os.path.dirname(os.path.abspath(a.sheet))
    for col in ("scaffold_path", "client_path", "client2_path"):
        if col in sheet:
            sheet[col] = [x if (not isinstance(x, str)) or os.path.isabs(x) else os.path.join(base, x)
                          for x in sheet[col]]
    os.makedirs(a.out, exist_ok=True)
    objs, summs = [], []
    for _, row in sheet.iterrows():
        o, s = analyse(row, p)
        objs.append(o); summs.append(s)
        print(f"[{row['image_id']}] n={s['n_condensates']}  median area={s['median_area_px']:.0f} px  "
              f"circ={s['median_circularity']:.3f}  AR={s['median_aspect_ratio']:.2f}")
    objs = pd.concat(objs, ignore_index=True)
    summ = pd.DataFrame(summs)
    objs.to_csv(os.path.join(a.out, "shape_objects.csv"), index=False)
    summ.to_csv(os.path.join(a.out, "shape_fov_summary.csv"), index=False)
    gallery(sheet, objs, p, a.out)
    with open(os.path.join(a.out, "params_used.json"), "w", encoding="utf-8") as fh:
        json.dump(p, fh, indent=2)
    print("\nPer-condition means of FOV medians:")
    cols = ["n_condensates", "pct_area_covered", "median_area_px", "median_circularity",
            "median_aspect_ratio", "median_solidity"]
    print(summ.groupby("condition")[cols].mean().round(3).T.to_string())


if __name__ == "__main__":
    main()

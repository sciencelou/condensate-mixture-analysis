#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_all.py  -  one-command analysis of ACAN / BCAN / NCAN condensate images.
Run:  python scripts/run_all.py        (Windows: double-click 2_run_analysis.bat)

What it does
  1. Scans every single-channel .tif in images/ (one sub-folder per independent preparation,
     e.g. images/exp1, images/exp2, images/exp3). Channels of the same field of view are grouped by
     file name (Fiji/Bio-Formats export "<name>.lif - SeriesNNN ... - C=<n>.tif").
  2. Assigns each channel to a protein from the 'channels <lif name> = ...' lines in settings.txt and
     from that the sample type (ACAN alone / ACAN+BCAN / ACAN+NCAN / ACAN+BCAN+NCAN).
  3. Runs the image analyses (recruitment_quant.py, radial_profile.py, shape_analysis.py,
     driver_analysis.py) and writes the sample tables it built to results/sample_tables/.
  4. Makes figures, legends and statistics with R (make_figures.R, plot_F3_prism_style.R) and the
     source-data workbook (make_source_data.py).
Author: Bin Luo, 2026. Written with AI assistance (Claude, Anthropic).
"""
import glob
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

COLOURS = {(0, 1, 0): "green", (1, 0, 0): "red", (0, 0, 1): "blue", (1, 0, 1): "magenta",
           (0, 1, 1): "cyan", (1, 1, 0): "yellow", (1, 1, 1): "grey"}


def say(msg):
    print(msg, flush=True)


def norm_name(t):
    return re.sub(r"\s+", " ", t).strip().lower()


def read_settings():
    s = {"green": "ACAN", "red": "BCAN", "magenta": "NCAN", "pixel_um": "", "detector_offset": "", "channels": {}}
    path = os.path.join(BASE, "settings.txt")
    if os.path.exists(path):
        for line in open(path, encoding="utf-8-sig"):
            line = line.split("#")[0].strip()
            if "=" in line:
                k, v = [x.strip() for x in line.split("=", 1)]
                if k.lower().startswith("channels "):        # channels <lif name> = protein of C=0, C=1, ...
                    s["channels"][norm_name(k[9:])] = [x.strip().upper() for x in v.split(",") if x.strip()]
                else:
                    s[k.lower()] = v
    return s


def channel_colour(path):
    import numpy as np
    import tifffile
    a = tifffile.imread(path)
    a = np.squeeze(a)
    if a.ndim == 3 and a.shape[-1] in (3, 4):
        on = tuple(int(x > 0) for x in a[..., :3].reshape(-1, 3).max(0))
        return COLOURS.get(on, "unknown")
    # raw single-channel image: the colour may be written in the file name,
    # e.g. "..._Series014_ch00_red.tif"
    m = re.search(r"_(red|green|blue|magenta|cyan|yellow|grey|gray)\.tif", os.path.basename(path), flags=re.I)
    return m.group(1).lower().replace("gray", "grey") if m else "grey"


def scan_images(settings):
    """Return a list of fields of view: dict(experiment, fov, files={protein: path})."""
    img_dir = os.path.join(BASE, "images")
    tifs = [f for f in glob.glob(os.path.join(img_dir, "**", "*.tif*"), recursive=True)
            if "overlay" not in os.path.basename(f).lower()]
    if not tifs:
        sys.exit("No .tif files found in the 'images' folder.")
    groups = defaultdict(dict)
    problems = []
    for f in sorted(tifs):
        name = os.path.basename(f)
        rel = os.path.relpath(os.path.dirname(f), img_dir)
        experiment = "exp1" if rel in (".", "") else rel.split(os.sep)[0]
        fiji = re.search(r"^(.*?)\s*-\s*C=(\d+)\.tiff?$", name, flags=re.I)     # Fiji default: "x.lif - Series014 4X - C=0.tif"
        m = fiji or re.search(r"^(.*?)[_ ]?ch(\d+)", name, flags=re.I)
        if not m:
            problems.append(f"{name}: no '_chNN' or ' - C=N' in the file name, skipped"); continue
        key_name, ch = m.group(1).rstrip("_ "), int(m.group(2))
        colour = channel_colour(f)
        protein = None
        if colour in ("grey", "unknown") and settings["channels"]:
            lif = norm_name(name.split(".lif")[0]) if ".lif" in name.lower() else norm_name(key_name)
            hits = [k for k in settings["channels"] if lif == k] or \
                   sorted([k for k in settings["channels"] if lif.startswith(k)], key=len, reverse=True)[:1]
            if hits and ch < len(settings["channels"][hits[0]]):
                protein = settings["channels"][hits[0]][ch]
        if protein is None and colour not in ("grey", "unknown"):
            protein = settings.get(colour)
        protein = protein or settings.get(f"ch{ch:02d}") or settings.get(f"ch{ch}")
        if not protein:
            problems.append(f"{name}: cannot tell which protein this channel is (colour '{colour}'); "
                            f"add a 'channels {name.split('.lif')[0]} = ...' line to settings.txt"); continue
        groups[(experiment, key_name)][protein.upper()] = f
    fovs = []
    for (experiment, key_name), files in sorted(groups.items()):
        if "ACAN" not in files:
            problems.append(f"{key_name}: no ACAN channel found, skipped"); continue
        fovs.append(dict(experiment=experiment, fov=key_name, files=files))
    return fovs, problems


def sample_type(files):
    clients = [p for p in ("BCAN", "NCAN") if p in files]
    return "ACAN alone" if not clients else "ACAN+" + "+".join(clients)


def short_id(fov, used):
    m = re.search(r"series\s*0*(\d+)", fov, flags=re.I)
    base = ("S" + m.group(1).zfill(3)) if m else re.sub(r"\W+", "", fov)[-12:]
    sid, i = base, 2
    while sid in used:
        sid = f"{base}_{i}"; i += 1
    used.add(sid)
    return sid


def build_tables(fovs, settings, out_dir):
    import pandas as pd
    px = settings.get("pixel_um", "")
    px = float(px) if px not in ("", None) else None
    off = settings.get("detector_offset", "")
    off = float(off) if off not in ("", None) else 0.0      # dark-level of raw images (subtracted before ratios / PC)
    two, triple, shape = [], [], []
    used = defaultdict(set)
    for f in fovs:
        typ = sample_type(f["files"])
        sid = short_id(f["fov"], used[(f["experiment"], typ)])
        fov_id = f"{f['experiment']}_{sid}"
        A = f["files"]["ACAN"]
        base = dict(replicate=f["experiment"], fov=fov_id, scaffold_path=A, pixel_um=px, bg_scaffold=off, bg_client=off)
        if typ in ("ACAN+BCAN", "ACAN+NCAN"):
            client = f["files"]["BCAN" if typ == "ACAN+BCAN" else "NCAN"]
            row = dict(image_id=f"{typ}_{fov_id}", condition=typ, client_path=client, **base)
            two.append(row); triple.append(row)
            shape.append(dict(image_id=row["image_id"], condition=typ, client_path=client, client2_path="", **base))
        elif typ == "ACAN+BCAN+NCAN":
            B, N = f["files"]["BCAN"], f["files"]["NCAN"]
            triple.append(dict(image_id=f"TRI_BCAN_{fov_id}", condition="Triple|BCAN", client_path=B, other_path="", **base))
            triple.append(dict(image_id=f"TRI_NCAN_{fov_id}", condition="Triple|NCAN", client_path=N, other_path=B, **base))
            shape.append(dict(image_id=f"TRI_{fov_id}", condition=typ, client_path=B, client2_path=N, **base))
        else:
            shape.append(dict(image_id=f"ALONE_{fov_id}", condition="ACAN alone", client_path="", client2_path="", **base))
    os.makedirs(out_dir, exist_ok=True)
    paths = {}
    for name, rows in (("two_protein", two), ("three_protein", triple), ("shape", shape)):
        if rows:
            p = os.path.join(out_dir, f"sample_table_{name}.csv")
            pd.DataFrame(rows).to_csv(p, index=False, encoding="utf-8-sig")
            paths[name] = p
    return paths, two, triple, shape


def run(cmd, label):
    say(f"\n>>> {label}")
    r = subprocess.run(cmd, cwd=BASE)
    if r.returncode != 0:
        sys.exit(f"Step failed: {label}")


def find_rscript():
    exe = shutil.which("Rscript")
    if exe:
        return exe
    cands = sorted(glob.glob(r"C:\Program Files\R\R-*\bin\Rscript.exe")) + \
            sorted(glob.glob(os.path.expanduser(r"~\AppData\Local\Programs\R\R-*\bin\Rscript.exe")))
    return cands[-1] if cands else None


def main():
    say("=" * 70 + "\nACAN / BCAN / NCAN condensate analysis\n" + "=" * 70)
    settings = read_settings()
    fovs, problems = scan_images(settings)
    for p in problems:
        say("  NOTE: " + p)
    say(f"\nFound {len(fovs)} fields of view:")
    count = defaultdict(lambda: defaultdict(int))
    for f in fovs:
        count[f["experiment"]][sample_type(f["files"])] += 1
    for exp, d in count.items():
        say(f"  experiment '{exp}': " + ", ".join(f"{k} = {v}" for k, v in sorted(d.items())))

    res = os.path.join(BASE, "results")
    tables, two, triple, shape = build_tables(fovs, settings, os.path.join(res, "sample_tables"))
    py = sys.executable
    sc = lambda n: os.path.join(HERE, n)
    li = os.path.join(res, "sample_tables", "params_li.json")
    nosh = os.path.join(res, "sample_tables", "params_no_shift_correction.json")
    open(li, "w").write('{"scaffold_threshold": "li"}')
    open(nosh, "w").write('{"correct_chromatic_shift": false}')

    if two:
        out = os.path.join(res, "recruitment")
        run([py, sc("recruitment_quant.py"), "--sheet", tables["two_protein"], "--out", out], "1/4 Recruitment of BCAN / NCAN into ACAN condensates")
        run([py, sc("recruitment_quant.py"), "--sheet", tables["two_protein"], "--out", os.path.join(out, "check_Li_threshold"), "--params", li], "    robustness check: other threshold")
        run([py, sc("recruitment_quant.py"), "--sheet", tables["two_protein"], "--out", os.path.join(out, "check_no_shift_correction"), "--params", nosh], "    robustness check: no channel-offset correction")
        pos = os.path.join(res, "position", "two_protein")
        run([py, sc("radial_profile.py"), "--sheet", tables["two_protein"], "--out", pos], "    where in the condensate: centre-to-edge profiles")
        run([py, sc("radial_profile.py"), "--sheet", tables["two_protein"], "--out", os.path.join(pos, "check_Li_threshold"), "--params", li], "    robustness check: other threshold")
        run([py, sc("radial_profile.py"), "--sheet", tables["two_protein"], "--out", os.path.join(pos, "check_no_shift_correction"), "--params", nosh], "    robustness check: no channel-offset correction")
    if shape:
        run([py, sc("shape_analysis.py"), "--sheet", tables["shape"], "--out", os.path.join(res, "shape")], "2/4 Size and shape of ACAN condensates")
    if any(r["condition"].startswith("Triple") for r in triple):
        out = os.path.join(res, "three_protein")
        run([py, sc("recruitment_quant.py"), "--sheet", tables["three_protein"], "--out", out], "3/4 Three-protein mix")
        run([py, sc("driver_analysis.py"), tables["three_protein"], out], "    which protein NCAN follows")
        tri_only = os.path.join(res, "sample_tables", "sample_table_three_protein_only.csv")
        import pandas as pd
        t = pd.read_csv(tables["three_protein"]); t[t["condition"].str.startswith("Triple")].to_csv(tri_only, index=False, encoding="utf-8-sig")
        pos = os.path.join(res, "position", "three_protein")
        run([py, sc("radial_profile.py"), "--sheet", tri_only, "--out", pos], "    where in the condensate (three-protein mix)")
        run([py, sc("radial_profile.py"), "--sheet", tri_only, "--out", os.path.join(pos, "check_Li_threshold"), "--params", li], "    robustness check: other threshold")
        run([py, sc("radial_profile.py"), "--sheet", tri_only, "--out", os.path.join(pos, "check_no_shift_correction"), "--params", nosh], "    robustness check: no channel-offset correction")

    rs = find_rscript()
    if rs:
        run([rs, sc("make_figures.R")], "4/4 Figures, legends and statistics (R)")
        say("\n>>>     Excel workbook with the numbers behind every figure")
        if subprocess.run([py, sc("make_source_data.py")], cwd=BASE).returncode != 0:
            say("  (Excel workbook not made - see the message above. If it says 'No module named openpyxl', double-click 1_install_packages.bat once, then run again.)")
        run([rs, sc("plot_F3_prism_style.R")], "    Prism-style version of the paired peripheral-shift panel")
        say("\nDONE. All results are in the 'results' folder:")
        say("  results/figures_overview.txt lists every figure.")
        say("  results/Figure_source_data.xlsx = the numbers behind every figure (yellow columns).")
    else:
        say("\nR (Rscript) was not found automatically.")
        say("Last step by hand: open scripts/make_figures.R in RStudio and click 'Source'.")


if __name__ == "__main__":
    main()

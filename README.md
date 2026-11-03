# LLPS condensate analysis (ACAN, BCAN, NCAN)

Image-analysis code for the in vitro condensate experiments in [Author list, Title, Journal, Year]
(Figure 8 and Figure S11). The pipeline segments ACAN condensates and BCAN/NCAN puncta in two- and
three-protein mixtures and measures:

- whether each ACAN condensate contains the client protein,
- where each BCAN/NCAN punctum lies relative to the ACAN condensate (inside, at the edge, beyond the edge),
- the radial position of the client relative to ACAN within each condensate (peripheral shift),
- the size of ACAN condensates.

Raw images are not included (see `images/README.md`).

## Repository layout

```
settings.txt                pixel size and channel order of each .lif file
1_install_packages.bat      installs the Python packages (Windows)
2_run_analysis.bat          runs the whole pipeline (Windows)
requirements.txt            Python packages
images/                     put raw single-channel TIFFs here (not tracked)
scripts/
  run_all.py                entry point: finds images, builds sample tables, runs every step
  recruitment_quant.py      channel registration, segmentation, client positivity and punctum location
  radial_profile.py         per-condensate half-maximum edge, radial profiles, peripheral shift
  shape_analysis.py         condensate area and shape
  driver_analysis.py        NCAN-ACAN correlation in the three-protein mixture (supplementary analysis)
  make_figures.R            runs all plot_*.R scripts (figures, legends, statistics)
  plot_recruitment.R        client positivity and location (two-protein mixtures)
  plot_triple.R             two- vs three-protein mixtures, condensate area with all four sample types
  plot_position.R           radial profiles and peripheral shift
  plot_shape.R              condensate size and shape (ACAN alone and two-protein mixtures)
  plot_driver.R             NCAN-ACAN correlation
  plot_F3_prism_style.R     paired peripheral-shift panel in the style of the other figure panels
  style_prism.R             shared plot style and statistical tests
  make_source_data.py       Excel workbook with the values behind every plot
```

## Requirements

- Python >= 3.10 with the packages in `requirements.txt` (tested with Python 3.13.16, NumPy 2.5.3,
  SciPy 1.18.1, scikit-image 0.26.0, pandas 3.0.5).
- R >= 4.2 with ggplot2, dplyr, readr, tidyr, officer, lme4, lmerTest, patchwork and svglite
  (installed automatically on the first run; tested with R 4.3.3, ggplot2 3.4.4).
- Fiji with Bio-Formats, only to export the raw images from the Leica .lif files.

## Usage

1. Export the raw images and place them in `images/exp1`, `images/exp2`, ... (one folder per independent
   preparation; see `images/README.md`).
2. Check `settings.txt` (pixel size and the channel order of each .lif file).
3. Install the Python packages: `pip install -r requirements.txt` (Windows: double-click `1_install_packages.bat`).
4. Run the pipeline: `python scripts/run_all.py` (Windows: double-click `2_run_analysis.bat`).

All output is written to `results/`: per-image and per-object tables (CSV), figures (PDF/SVG), legends,
statistics (`stats*.txt`) and `Figure_source_data.xlsx`.

## Correspondence between paper panels and outputs

| Paper panel | Measurement | Output (in `results/`) |
| --- | --- | --- |
| Figure 8H | ACAN condensates positive for BCAN/NCAN (%) | `recruitment/figures/FigB1_*` |
| Figure 8I | BCAN/NCAN puncta within ACAN condensates (%) | `recruitment/figures/FigB2_*` |
| Figure 8J | Puncta inside / at edge / beyond edge (%) | `recruitment/figures/FigC_*` |
| Figure 8L, 8M | Puncta within condensates, two vs three proteins | `three_protein/figures/FigE1b_*`, `FigE1n_*` |
| Figure 8N | Condensate area, four sample types | `three_protein/figures/FigE5a_*` |
| Figure S11E | Centre-to-edge intensity profiles | `position/figures/FigF1_*` |
| Figure S11F | Peripheral shift, two-protein mixtures | `position/figures/FigF2_*` |
| Figure S11G | Peripheral shift, BCAN vs NCAN in the same condensates | `position/figures/FigF3_*` |

The pipeline also produces robustness checks (Li threshold, no channel registration, alternative
positivity criteria, size-matched condensates) and additional analyses that are not shown in the paper.

## Method summary

1. **Channel registration.** The offset between ACAN and the client channel is measured in each image by
   sub-pixel phase cross-correlation (upsampling factor 10; Guizar-Sicairos et al., 2008) and corrected
   before analysis.
2. **Segmentation.** Each channel is smoothed with a Gaussian filter (sigma = 1 pixel) and thresholded with
   Otsu's method (Otsu, 1979); holes are filled and objects < 9 pixels are removed. Condensates touching the
   image border are excluded from per-condensate values. Intensities are measured on unsmoothed images.
   The dilute phase is all pixels > 3 pixels from any object.
3. **Client-positive condensate.** Mean client intensity inside the ACAN condensate > dilute-phase mean + 3 SD.
4. **Punctum location.** Inside: >= 50% of the punctum area overlaps the ACAN mask; at the edge: < 50% but
   touches the mask after a 1-pixel dilation; beyond the edge: otherwise.
5. **Radial distribution.** For each client-positive condensate, the edge is the half-maximum contour of its
   own ACAN signal, background + 0.5 x (95th percentile - background). Each pixel gets a normalized radial
   position rho from its Euclidean distance to this edge (0 = centre, 1 = edge, > 1 outside, up to
   min(3 pixels, 0.5 x maximum depth)). Peripheral shift = (intensity-weighted mean rho of the client /
   that of ACAN - 1) x 100%.
6. **Condensate area.** Median area of all ACAN condensates in each field.
7. **Statistics.** Field values are averaged per independent preparation (n = 3), the statistical unit.
   Two groups: two-tailed unpaired Student's t-test; same condensates: paired t-test; four groups: one-way
   ANOVA with Tukey's multiple comparisons test. Mean +/- SEM.

## References

- Guizar-Sicairos M, Thurman ST, Fienup JR (2008). Efficient subpixel image registration algorithms.
  Opt Lett 33, 156-158. https://doi.org/10.1364/OL.33.000156
- Otsu N (1979). A threshold selection method from gray-level histograms. IEEE Trans Syst Man Cybern 9,
  62-66. https://doi.org/10.1109/TSMC.1979.4310076
- van der Walt S et al. (2014). scikit-image: image processing in Python. PeerJ 2, e453.
  https://doi.org/10.7717/peerj.453
- Virtanen P et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python.
  Nat Methods 17, 261-272. https://doi.org/10.1038/s41592-019-0686-2

## Citation

If you use this code, please cite: [Author list, Title, Journal, Year, DOI].

## License

MIT (see `LICENSE`). Code written with AI assistance (Claude, Anthropic).

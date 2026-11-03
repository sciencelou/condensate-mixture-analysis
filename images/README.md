# Input images (not included)

Place the raw single-channel TIFFs here, one sub-folder per independent protein preparation:

```
images/
  exp1/
    ACAN FL 6%.lif - Series005 4x - C=0.tif
    BCAN-ACAN FL 6% PEG.lif - Series012 4X - C=0.tif
    BCAN-ACAN FL 6% PEG.lif - Series012 4X - C=1.tif
    ...
  exp2/
  exp3/
```

Export: open each Leica .lif file in Fiji with Bio-Formats (one series at a time, "Colorized" off),
then File > Save As > Tiff for each channel (Image > Stacks > Stack to Images). Pixel values must not be
changed (no contrast adjustment, no RGB conversion). The channel order of every .lif file is declared in
`settings.txt`.

The raw images are available from the corresponding author upon reasonable request [or: at <data repository DOI>].

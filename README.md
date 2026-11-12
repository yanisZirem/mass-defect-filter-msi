# Mass defect filter for lipid ion lists in mass spectrometry imaging

R script that flags which ions of a MALDI-MSI peak list are **lipid-like** and which are **non-lipid** (matrix clusters, salts, peptides, small metabolites, ...), based on their mass defect.

For each observed *m/z*, the script checks whether its mass defect falls inside the envelope of theoretical lipid ions at the same nominal mass. Ions outside this window are reported as non-lipid and can be listed, plotted and discussed separately.

## Method

1. **Theoretical lipid ions** are generated for the following classes, over a range of total chain length, number of double bonds and (for some classes) additional oxygens:
   FA, LPC, LPE, LPS, LPI, LPA, LPG, PC, PE, PS, PI, PA, PG, PC-O, PE-O, SM, Cer, HexCer, DG, TG, CE, sterols (ST) and sulfated sterols.
2. Each species is converted into ions with the adducts `[M+H]+`, `[M+Na]+`, `[M+K]+`, `[M+NH4]+`, `[M+H-H2O]+` and `[M+2Na-H]+`.
3. The **mass defect** is computed as `m/z - floor(m/z)`. For every nominal mass, the **lipid envelope** is the min-max mass defect of the theoretical ions, widened by a tolerance (default ±0.010 Da).
4. Observed ions are classified as:
   - `Lipid-like`: defect inside the window
   - `Non-lipid (defect below lipid window)`
   - `Non-lipid (defect above lipid window)`
5. Additional outputs: relative mass defect (ppm), Kendrick mass defect (CH2 base), and a hint when an ion matches a DHB matrix cluster `[nDHB - mH2O + H/Na/K]+` (within 5 mDa; indicative only).
6. **Optional (imzML)**: the same filter is applied to all peaks detected in the mean spectrum of an imzML file, and the share of the summed intensity attributed to lipid-like vs non-lipid peaks is reported.

## Requirements

- R >= 4.0
- CRAN: `readxl`, `ggplot2`, `writexl`
- Bioconductor, only for the imzML part: [`Cardinal`](https://bioconductor.org/packages/Cardinal) >= 3.0

```r
install.packages(c("readxl", "ggplot2", "writexl"))
# only if you use an imzML file
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Cardinal")
```

## Input

An Excel file (`.xlsx`) with one ion per row. Defaults (editable in the `CONFIG` block at the top of the script):

| Setting | Default | Meaning |
|---|---|---|
| `col_mz` | 1 | column with the *m/z* values |
| `col_annotation` | 7 | column with the annotation (optional content) |
| `col_msms` | 8 | column with the MS/MS confirmation (optional content) |
| `first_row` | 3 | first data row (rows above are title/headers) |

Rows whose *m/z* is not numeric (e.g. a `Total` row) are ignored. Annotation columns are only used to compare annotated ions with the filter result; empty cells are treated as unannotated.

## Usage

```bash
# ion list only
Rscript mass_defect_filter.R ion_list.xlsx

# ion list + imzML (mean spectrum, all detected peaks)
Rscript mass_defect_filter.R ion_list.xlsx sample.imzML
```

Or from R / RStudio: set `xlsx_file` (and `imzml_file`) in the `CONFIG` block and run `source("mass_defect_filter.R")`.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `md_tol` | `0.010` | tolerance (Da) added on each side of the lipid envelope; increase it to be more permissive |
| `mz_range` | `c(400, 1000)` | *m/z* range for the envelope and plots |
| `ppm_tol` | `10` | tolerance (ppm) to find an ion of the list in the imzML mean spectrum |
| `snr_min` | `3` | S/N threshold for peak detection in the mean spectrum |
| `out_dir` | `mass_defect_output` | output folder |

The theoretical lipid set (classes, chain ranges, adducts) is defined in the `classes` and `adducts` objects and can be extended.

## Outputs

Written in `out_dir`:

| File | Content |
|---|---|
| `mass_defect_results.xlsx` | sheets: all ions with class, defect, window, KMD and check; non-lipid ions only; summary; annotation vs filter; lipid envelope |
| `mass_defect_plot.png` / `.pdf` | mass defect vs *m/z*, theoretical lipid window in green, annotated vs unannotated ions |
| `kendrick_plot.png` | Kendrick mass defect (CH2) plot |
| `mass_defect_imzML.xlsx` | (imzML only) detected peaks with class, intensity share per class, intensity of the listed ions |
| `mean_spectrum_mdf.png` | (imzML only) mean spectrum coloured by class |

The `check` column cross-references annotation and filter:

- `Annotated, passes MDF`
- `Annotated but fails MDF (to re-check)`: the annotation is not compatible with a lipid mass defect and should be verified
- `Unannotated, lipid-like (candidate)`
- `Unannotated, non-lipid`

## Limitations

- The filter is based on mass defect only: it indicates that an ion is **compatible** with a lipid, not that it **is** a lipid. Non-lipid compounds with a lipid-like defect pass the filter, and lipid classes absent from the theoretical set (e.g. cardiolipins, gangliosides) may be flagged as non-lipid. Add them to `classes` if needed.
- The window width depends on the lipid classes included; a broader set gives a wider, less selective window.
- Defects are computed with the `floor` convention, which is appropriate for lipids in the *m/z* 400-1000 range. Species with a negative defect relative to their nominal mass (heavily oxidised or phosphorylated lipids) would be misclassified.
- DHB cluster matches are indicative and may be coincidental at the 5 mDa tolerance.
- Peak detection in the imzML part is a simple local-maximum search on the mean spectrum; for publication-grade peak picking, provide your own peak list through the Excel input.

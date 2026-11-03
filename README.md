# Diagnostic bile acid MS/MS spectrum filter

R script that filters MS/MS spectra of **bile acid–alkylamine conjugates** (made by reductive amination at the 3-position of the bile acid) using diagnostic fragment ions of the bile acid core.

A spectrum is kept only if it contains the diagnostic core ion **and** that ion is more intense than at least one of its water-loss fragments.

## Diagnostic ions

| Bile acid core (3-alkylated) | Diagnostic ion (m/z) | –H₂O | –2H₂O | Library keyword |
|---|---|---|---|---|
| Lithocholic acid (from dehydrolithocholic acid) | 359.2945 (C₂₄H₃₉O₂⁺) | 341.2839 (C₂₄H₃₇O⁺) | 323.2733 (C₂₄H₃₅⁺) | `Dehydrolithocholicacid` |
| Deoxycholic acid (from 3-oxo-deoxycholic acid) | 357.2788 (C₂₄H₃₇O₂⁺) | 339.2682 (C₂₄H₃₅O⁺) | 321.2577 (C₂₄H₃₃⁺) | `3-oxodeoxycholicacid` |
| Cholic acid (from 3-oxo-cholic acid) | 355.2632 (C₂₄H₃₅O₂⁺) | 337.2526 (C₂₄H₃₃O⁺) | 319.2420 (C₂₄H₃₁⁺) | `3-oxocholicacid` |

Peaks are matched within **±20 ppm**.

## Filtering rule

For each spectrum whose `NAME` contains the library keyword, the spectrum **passes** if:

1. the diagnostic ion (359 / 357 / 355) is present with intensity ≥ 1, **and**
2. at least one water-loss ion is present with intensity ≥ 1, **and**
3. diagnostic ion / water-loss ion intensity ratio is **> 1** for that pair.

Each spectrum gets a `cutoff_reason` (e.g. `PASS via: 359/341` or `FAIL: intensity_359 absent or < 1`) so the decision can be audited.

## Input

An MGF spectral library, e.g. `GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf` from GNPS, with standard `BEGIN IONS … END IONS` blocks containing `NAME=`, `SCANS=`, `SPECTRUMID=` and peak lines (`m/z intensity`).

## Usage

1. Install R (≥ 4.0). No extra packages are needed (base R only).
2. Put the MGF file in a folder and change the `setwd(...)` line in each section of the script to that folder.
3. Run the script:

```r
source("MASSQL_bile_acid_alkamines filter.R")
```

## Output

For each bile acid core, two CSV files are written:

| File | Content |
|---|---|
| `<bileacid>_all_results.csv` | every matching spectrum with intensities, ratios, pass/fail and reason |
| `<bileacid>_filtered_ratio_above1.csv` | only spectra that pass |

Columns: `NAME, SCANS, SPECTRUMID, intensity_<ion>, intensity_<ion-H2O>, intensity_<ion-2H2O>, ratio_*, passes_cutoff, cutoff_reason`.

A summary (spectra extracted / passing / failing) is printed to the console.

## Related resources

- GNPS spectral library: GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY
- MassQL: https://mwang87.github.io/MassQueryLanguage_Documentation/

## Author

Julius Agongo (Dorrestein Lab)

## Citation

If you use this code, please cite this repository. <!-- add paper citation / DOI when available -->

## License

<!-- e.g. MIT — add a LICENSE file via "Add file → Create new file → LICENSE" and choose a template -->

# Diagnostic bile acid MS/MS spectrum filter

R script that filters MS/MS spectra of **bile acid–alkylamine conjugates** (made by reductive amination at the 3-position of the bile acid) using diagnostic fragment ions of the bile acid core.

A spectrum is kept only if it contains the diagnostic core ion **and** that ion is more intense than at least one of its water-loss fragments.

## Diagnostic ions

| Bile acid core (3-alkylated) | Diagnostic ion (m/z) | –H₂O | –2H₂O |
|---|---|---|---|
| Lithocholic acid (from dehydrolithocholic acid) | 359.2930 | 341.2831 | 323.2722 |
| Deoxycholic acid (from 3-keto-deoxycholic acid) | 357.2773 | 339.2667 | 321.2561 |
| Cholic acid (from 3-keto-cholic acid) | 355.2616 | 337.2510 | 319.2405 |

Peaks are matched within **±20 ppm**.

## Filtering rule

For each spectrum whose `NAME` contains the bile acid keyword, the spectrum **passes** if:

1. the diagnostic ion (359 / 357 / 355) is present with intensity ≥ 1, **and**
2. at least one water-loss ion is present with intensity ≥ 1, **and**
3. the diagnostic ion / water-loss ion intensity ratio is **> 1** for that pair.

Each spectrum gets a `cutoff_reason` (e.g. `PASS via: 359/341`) so every decision can be audited.

## Input

An MGF spectral library, e.g. `GNPS-ALKYLAMINES-BILE-ACIDS-LIBRARY.mgf` from GNPS.

## Usage

1. Install R (≥ 4.0). Only base R is needed.
2. Put the MGF file in a folder and set `setwd(...)` in the script to that folder.
3. Run:

```r
source("MASSQL_bile_acid_alkamines filter.R")
```

## Output

For each bile acid core, two CSV files:

| File | Content |
|---|---|
| `<bileacid>_all_results.csv` | every matching spectrum with intensities, ratios, pass/fail and reason |
| `<bileacid>_filtered_ratio_above1.csv` | only spectra that pass |

A summary (extracted / passing / failing) is printed to the console.

## Author

Julius Agongo (Dorrestein Lab)

## Citation

If you use this code, please cite this repository.

## License

Add a license via **Add file → Create new file → LICENSE** and choose a template.

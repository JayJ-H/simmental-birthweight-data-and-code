# Simmental Calf Birth Weight Data and Reproduction Code

Data and R code supporting **Low and high birth weight phenotypes in singleton
Simmental calves: distributional characteristics and associated factors**.

This release corresponds to the revised analyses, not the initial submission.
The study uses routine calving records from one commercial herd, covering
January 1, 2023 through June 22, 2026. The raw export contains 9,480 records;
9,478 have birth weight, and the primary analytical cohort contains 9,384
singleton records from 4,624 dams. It includes calves that survived the day
of birth. Sex-specific P10/P90 classification yields 948 lower-tail and 955
upper-tail records. These are distributional categories, not clinical risk
cut-offs.

## Contents

| Path | Content |
| --- | --- |
| `original_analysis/data/raw/calving_records.xlsx` | Original export with document properties removed |
| `original_analysis/data/processed/` | Cleaned records, singleton and twin datasets, derived variables |
| `original_analysis/scripts/` | Data cleaning, original model and exploratory-figure code |
| `revision_analysis/scripts/` | Final threshold, restricted-sample, ancestry and bootstrap analyses |
| `expected_results/revision/` | Numerical results used in the revised manuscript |
| `original_analysis/outputs/figures/` | CSV source data for figure panels |
| `figure_source/` | Updated Figure 3E and 3F source data |
| `figures/` | Four final main figures |
| `tables/` | Main tables 1-4 and seven supplementary worksheets |
| `environment/` | Recorded software versions |
| `manifest.json` | SHA256 file checksums |

## Run

Clone the GitHub repository or extract the release archive, then run from its
root directory:

```sh
python verify_release.py
Rscript run_reproduction.R --clean-only
Rscript run_reproduction.R
```

The first command checks archived files. The second rebuilds the analytical
dataset from the raw workbook and compares all columns and rows. The third
reruns model fits and the 999-replicate dam-level bootstrap; it may take
substantial time. See `REPRODUCIBILITY.md` for dependencies and output locations.

Animal IDs and recorded sire/ancestry labels are deliberately retained under
the data owner's permission for public release. Computer usernames, local
paths and document author/editor properties are not part of the release.
Manuscripts, reviewer letters, signed forms and internal revision records are
excluded.

## Availability and Licensing

GitHub repository:
https://github.com/JayJ-H/simmental-birthweight-data-and-code

Data and derived research outputs are licensed under CC BY 4.0.
R and Python code is licensed under MIT. See `LICENSE.md` and `RIGHTS.md`
for scope and attribution, and `LICENSE_CODE.txt` for the MIT terms.
The Zenodo archival record is being finalized; no registered DOI is claimed
until the public record has been verified.

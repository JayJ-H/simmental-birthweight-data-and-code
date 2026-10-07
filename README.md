# Simmental Calf Birth Weight: Data and R Code

Data and reproduction code supporting *Low and high birth weight phenotypes in
singleton Simmental calves: distributional characteristics and associated
factors*.

## Download and Run

Download `GitHub_release.zip` from this repository and extract it. The archive
contains 94 files, including raw and processed records, original and revised
R analyses, numerical reference outputs, figure source data, four main figures,
main tables, and seven supplementary worksheets.

From the extracted archive directory:

```sh
python verify_release.py
Rscript run_reproduction.R --clean-only
Rscript run_reproduction.R
```

See the archived `README.md`, `REPRODUCIBILITY.md`, and `DATA_DICTIONARY.md`
for the detailed directory structure, dependencies, and data definitions.
The archived documentation describes the preparation state of this fixed
snapshot; this GitHub repository now provides public access to that snapshot.

## Study and Verification

Routine calving records cover January 1, 2023 through June 22, 2026.
The primary analysis includes 9,384 singleton records from 4,624 dams.
The cohort includes calves that survived the day of birth.
Sex-specific P10/P90 groups are distributional categories, not validated
clinical risk cut-offs.

The raw-to-analysis reconstruction was run and matched the existing cohort.
The complete model and 999-replicate bootstrap code is included; all models
were not freshly refitted during release packaging.

The uploaded archive was downloaded from GitHub and verified against the local
file. Its SHA-256 is:

```text
b36b89d80d52d3cc2cbb65f14bd844bc9a01d4b6900c81e17ea7d1766489b986
```

## Rights and Privacy

Animal identifiers and sire/ancestry labels are retained with the data owner's
permission for public release. Computer paths and document author/editor
properties were removed. Manuscripts, reviewer correspondence, signed forms,
and internal revision records are excluded.

Public availability does not grant a general reuse license. Data and code
reuse licenses remain to be confirmed; see the archived `RIGHTS.md`.
A Zenodo archival record is being prepared. No published DOI is claimed here.

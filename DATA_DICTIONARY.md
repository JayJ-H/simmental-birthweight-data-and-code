# Data Dictionary

## Raw Fields

| Original field | Analysis field | Meaning |
| --- | --- | --- |
| 牛号 | `dam_id` | Recorded dam identifier |
| 月龄 | `dam_age_month` | Recorded dam age in months |
| 胎次 | `parity` | Recorded parity |
| 产犊日期 | `calving_date` | Calving date |
| 受孕公牛 | `service_sire` | Recorded semen-source bull code |
| 初生重(犊牛) | `birth_weight` | Birth weight in kg |
| 犊牛牛号(犊牛) | `calf_id` | Recorded calf identifier |
| 外祖父 | `maternal_grandsire` | Historical maternal-grandsire grouping label |
| 外祖母 | `maternal_granddam` | Recorded maternal-granddam label |

IDs should be read as text, not converted to numeric values. They are original
animal identifiers, not human participant identifiers.

## Derived Variables

- `row_id`: one-based row identifier in the imported raw export.
- `calf_sex`: `male` for calf IDs beginning with G; `female` for other nonmissing
  IDs under the farm's existing numbering convention.
- `calf_sex_code`: the corresponding G-coded/non-G-coded field.
- `calf_birth_year`, `calf_birth_month`: derived from calving date.
- `calf_birth_season`: Spring March-May, Summer June-August, Autumn
  September-November, Winter December-February.
- `parity_group`: 1 through 6 and >=7.
- `birth_type`: identified from dam ID and calving date after duplicate checks.
- `lower_tail`, `upper_tail`: <= within-sex P10 and >= within-sex P90.
- `lower_tail_p5`, `upper_tail_p95`, `lower_tail_p15`, `upper_tail_p85`:
  alternative within-sex percentile definitions.
- Other absolute/overall tail flags are legacy comparison fields, not the
  revised study's primary clinical definitions.
- `dam_birth_date_estimated` and related birth-year fields are calculated
  dates, not independently recorded dam birth dates.
- Export-date-based age fields and impossible-age flags are QC comparisons.
  The primary models use parity, not reconstructed age.
- Duplicate and extreme-value flags are documented in
  `original_analysis/scripts/01_data_cleaning.R`; the primary models retain
  extreme weights passing record checks. The separate 25-65 kg analysis uses
  inclusive limits.

`data_columns.csv` lists every column in the primary processed dataset, with
missingness counts and an example. Empty cells and `NA` indicate missing
values. Spreadsheet date formatting is retained in the raw workbook; processed
CSV dates use year-month-day.

## Population and Ancestry

The record extract covers calves surviving the day of birth; it is not a dataset
of all births including same-day deaths. Twins are described separately.
Ancestry fields are used as recorded groups, not as a newly validated pedigree.
Missing ancestry labels are not fabricated. Each ancestry comparison uses the
same available-record subset for its baseline and adjusted model.

## Numerical Outputs

`OR`, `CI_low`, `CI_high`, `p` denote odds ratios, 95% confidence limits and
P values. `log_OR` is the logarithm of the odds ratio. `OR_ratio` compares the
named analysis with its reference model, not a risk ratio.
Quantile-regression `estimate` values are in kg, with empirical bootstrap
intervals. `model_diagnostics.csv` and `numerical_checks.csv` document numerical
messages, boundary variance estimates and optimizer comparisons.

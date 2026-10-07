source("scripts/00_config.R")

write_log("03_descriptive_stats.R started")

all_valid_records <- read_csv_safe(file.path(PROCESSED_DIR, "all_valid_records.csv")) |> prepare_analysis_types()
singleton_analysis <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis.csv")) |> prepare_analysis_types()
twin_records <- read_csv_safe(file.path(PROCESSED_DIR, "twin_records.csv")) |> prepare_analysis_types()

sex_var <- sex_variable_name()

dataset_summary <- dplyr::bind_rows(
  summary_stats(all_valid_records) |> dplyr::mutate(dataset = "all_valid_records", .before = 1),
  summary_stats(singleton_analysis) |> dplyr::mutate(dataset = "singleton_analysis", .before = 1),
  summary_stats(twin_records) |> dplyr::mutate(dataset = "twin_records", .before = 1)
)

grouped_summary <- function(data, group_var) {
  data |>
    dplyr::group_by(.data[[group_var]], .drop = FALSE) |>
    dplyr::group_modify(~ summary_stats(.x)) |>
    dplyr::ungroup() |>
    dplyr::rename(group = 1) |>
    dplyr::mutate(grouping_variable = group_var, .before = 1)
}

year_season_summary <- singleton_analysis |>
  dplyr::group_by(calf_birth_year, calf_birth_season, .drop = FALSE) |>
  dplyr::group_modify(~ summary_stats(.x)) |>
  dplyr::ungroup()

service_sire_summary <- singleton_analysis |>
  dplyr::group_by(service_sire) |>
  dplyr::summarise(
    n = dplyr::n(),
    missing_label = all(is.na(service_sire)),
    mean_birth_weight = mean(birth_weight, na.rm = TRUE),
    SD = stats::sd(birth_weight, na.rm = TRUE),
    lower_extreme_n = sum(extreme_low_flag, na.rm = TRUE),
    upper_extreme_n = sum(extreme_high_flag, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n), service_sire)

maternal_grandsire_summary <- singleton_analysis |>
  dplyr::group_by(maternal_grandsire) |>
  dplyr::summarise(
    n = dplyr::n(),
    missing_label = all(is.na(maternal_grandsire)),
    mean_birth_weight = mean(birth_weight, na.rm = TRUE),
    SD = stats::sd(birth_weight, na.rm = TRUE),
    lower_extreme_n = sum(extreme_low_flag, na.rm = TRUE),
    upper_extreme_n = sum(extreme_high_flag, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(n), maternal_grandsire)

safe_write_xlsx(
  list(
    dataset_summary = dataset_summary,
    singleton_by_sex = grouped_summary(singleton_analysis, sex_var),
    singleton_by_sex_code = grouped_summary(singleton_analysis, "calf_sex_code"),
    singleton_by_parity = grouped_summary(singleton_analysis, "parity_group"),
    singleton_by_dam_age = grouped_summary(singleton_analysis, "dam_age_group"),
    singleton_by_year = grouped_summary(singleton_analysis, "calf_birth_year"),
    singleton_by_season = grouped_summary(singleton_analysis, "calf_birth_season"),
    singleton_year_season = year_season_summary,
    supplementary_dam_birth_year = grouped_summary(singleton_analysis, "dam_birth_year"),
    supplementary_dam_birth_season = grouped_summary(singleton_analysis, "dam_birth_season"),
    service_sire_counts = service_sire_summary,
    maternal_grandsire_counts = maternal_grandsire_summary
  ),
  file.path(TABLE_DIR, "descriptive_birth_weight.xlsx")
)

write_log("03_descriptive_stats.R completed")

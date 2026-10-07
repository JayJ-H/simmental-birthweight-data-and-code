source("scripts/00_config.R")

write_log("01_data_cleaning.R started")

required_original_fields <- c(
  "牛号", "月龄", "胎次", "产犊日期", "受孕公牛",
  "初生重(犊牛)", "犊牛牛号(犊牛)", "外祖父", "外祖母"
)

raw <- safe_read_excel(RAW_FILE)
missing_fields <- setdiff(required_original_fields, names(raw))
if (length(missing_fields) > 0) {
  stop("Missing required original field(s): ", paste(missing_fields, collapse = ", "), call. = FALSE)
}

field_overview <- tibble::tibble(
  field = names(raw),
  class = vapply(raw, function(x) paste(class(x), collapse = ";"), character(1)),
  non_missing = vapply(raw, function(x) sum(!is.na(x)), integer(1)),
  missing = vapply(raw, function(x) sum(is.na(x)), integer(1))
)

preview <- utils::head(raw, 10)
safe_write_xlsx(
  list(field_overview = field_overview, raw_preview_first10 = preview),
  file.path(TABLE_DIR, "raw_field_overview.xlsx")
)

dat <- raw |>
  dplyr::transmute(
    dam_id = clean_id(.data[["牛号"]]),
    dam_age_month = to_numeric_safe(.data[["月龄"]]),
    parity = to_numeric_safe(.data[["胎次"]]),
    calving_date = parse_excel_date(.data[["产犊日期"]]),
    service_sire = clean_id(.data[["受孕公牛"]]),
    birth_weight = to_numeric_safe(.data[["初生重(犊牛)"]]),
    calf_id = clean_id(.data[["犊牛牛号(犊牛)"]]),
    maternal_grandsire = clean_id(.data[["外祖父"]]),
    maternal_granddam = clean_id(.data[["外祖母"]])
  )

conversion_issues <- dplyr::bind_rows(
  raw |>
    dplyr::mutate(row_id = dplyr::row_number(), parsed = parse_excel_date(.data[["产犊日期"]])) |>
    dplyr::filter(is.na(parsed) & !is.na(.data[["产犊日期"]])) |>
    dplyr::transmute(row_id, field = "产犊日期", raw_value = as.character(.data[["产犊日期"]]), issue = "date_parse_failed"),
  raw |>
    dplyr::mutate(row_id = dplyr::row_number(), parsed = to_numeric_safe(.data[["月龄"]])) |>
    dplyr::filter(is.na(parsed) & !is.na(.data[["月龄"]])) |>
    dplyr::transmute(row_id, field = "月龄", raw_value = as.character(.data[["月龄"]]), issue = "numeric_parse_failed"),
  raw |>
    dplyr::mutate(row_id = dplyr::row_number(), parsed = to_numeric_safe(.data[["胎次"]])) |>
    dplyr::filter(is.na(parsed) & !is.na(.data[["胎次"]])) |>
    dplyr::transmute(row_id, field = "胎次", raw_value = as.character(.data[["胎次"]]), issue = "numeric_parse_failed"),
  raw |>
    dplyr::mutate(row_id = dplyr::row_number(), parsed = to_numeric_safe(.data[["初生重(犊牛)"]])) |>
    dplyr::filter(is.na(parsed) & !is.na(.data[["初生重(犊牛)"]])) |>
    dplyr::transmute(row_id, field = "初生重(犊牛)", raw_value = as.character(.data[["初生重(犊牛)"]]), issue = "numeric_parse_failed")
)

if (nrow(conversion_issues) > 0) {
  safe_write_xlsx(conversion_issues, file.path(TABLE_DIR, "conversion_issues.xlsx"))
}

dat <- dat |>
  dplyr::mutate(
    row_id = dplyr::row_number(),
    calf_birth_year = lubridate::year(calving_date),
    calf_birth_month = lubridate::month(calving_date),
    calf_birth_season = season_from_month(calf_birth_month),
    dam_birth_date_estimated = subtract_months(calving_date, dam_age_month),
    dam_birth_date_export_based = subtract_months(DATA_EXPORT_DATE, dam_age_month),
    dam_age_month_export_based_at_calving = lubridate::time_length(
      lubridate::interval(dam_birth_date_export_based, calving_date),
      unit = "month"
    ),
    export_based_impossible_age_flag = !is.na(dam_age_month_export_based_at_calving) &
      dam_age_month_export_based_at_calving < 15,
    dam_birth_year = lubridate::year(dam_birth_date_estimated),
    dam_birth_month = lubridate::month(dam_birth_date_estimated),
    dam_birth_season = season_from_month(dam_birth_month),
    calf_sex_code = dplyr::case_when(
      is.na(calf_id) ~ NA_character_,
      stringr::str_detect(calf_id, "^G") ~ "G-coded",
      TRUE ~ "non-G-coded"
    ),
    calf_sex = dplyr::case_when(
      is.na(calf_id) ~ NA_character_,
      SEX_RULE_CONFIRMED & stringr::str_detect(calf_id, "^G") ~ "male",
      SEX_RULE_CONFIRMED ~ "female",
      TRUE ~ NA_character_
    ),
    parity_group = dplyr::case_when(
      is.na(parity) ~ NA_character_,
      parity >= 7 ~ ">=7",
      TRUE ~ as.character(as.integer(parity))
    ),
    parity_group = factor(parity_group, levels = c("1", "2", "3", "4", "5", "6", ">=7"), ordered = FALSE),
    dam_age_month_cont = dam_age_month,
    dam_age_group = dplyr::case_when(
      is.na(dam_age_month) ~ NA_character_,
      dam_age_month <= 30 ~ "<=30",
      dam_age_month <= 42 ~ "31-42",
      dam_age_month <= 60 ~ "43-60",
      TRUE ~ ">60"
    ),
    dam_age_group = factor(dam_age_group, levels = c("<=30", "31-42", "43-60", ">60"), ordered = FALSE),
    extreme_low_flag = !is.na(birth_weight) & birth_weight < 25,
    extreme_high_flag = !is.na(birth_weight) & birth_weight > 65,
    implausible_flag = !is.na(birth_weight) & (birth_weight < 20 | birth_weight > 80)
  )

core_fields <- c(
  "dam_id", "dam_age_month", "parity", "calving_date", "service_sire",
  "birth_weight", "calf_id", "maternal_grandsire", "maternal_granddam"
)

dat <- dat |>
  dplyr::mutate(
    exact_duplicate_record_flag = duplicated(dplyr::pick(dplyr::all_of(core_fields))) |
      duplicated(dplyr::pick(dplyr::all_of(core_fields)), fromLast = TRUE),
    exact_duplicate_exclude_flag = duplicated(dplyr::pick(dplyr::all_of(core_fields))),
    duplicate_dam_date_calf_id_flag = !is.na(dam_id) & !is.na(calving_date) & !is.na(calf_id) &
      (duplicated(dplyr::pick(dam_id, calving_date, calf_id)) |
         duplicated(dplyr::pick(dam_id, calving_date, calf_id), fromLast = TRUE)),
    duplicate_dam_date_calf_id_exclude_flag = !is.na(dam_id) & !is.na(calving_date) & !is.na(calf_id) &
      duplicated(dplyr::pick(dam_id, calving_date, calf_id)),
    duplicate_calf_id_flag = !is.na(calf_id) &
      (duplicated(calf_id) | duplicated(calf_id, fromLast = TRUE))
  )

same_birth_weight_rows <- dat |>
  dplyr::filter(!is.na(dam_id), !is.na(calving_date), !is.na(birth_weight)) |>
  dplyr::group_by(dam_id, calving_date) |>
  dplyr::mutate(dam_date_n = dplyr::n()) |>
  dplyr::group_by(dam_id, calving_date, birth_weight, .add = FALSE) |>
  dplyr::mutate(same_weight_n = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::filter(dam_date_n > 1, same_weight_n > 1) |>
  dplyr::pull(row_id)

dat <- dat |>
  dplyr::mutate(
    same_dam_date_same_birth_weight_flag = row_id %in% same_birth_weight_rows,
    obvious_duplicate_exclude_flag = exact_duplicate_exclude_flag | duplicate_dam_date_calf_id_exclude_flag
  )

dedup <- dat |>
  dplyr::filter(!obvious_duplicate_exclude_flag) |>
  dplyr::group_by(dam_id, calving_date) |>
  dplyr::mutate(
    dam_date_record_n = dplyr::n(),
    dam_date_unique_calf_n = dplyr::n_distinct(calf_id, na.rm = TRUE),
    birth_type = dplyr::case_when(
      dam_date_record_n == 1 ~ "singleton",
      dam_date_record_n == 2 & dam_date_unique_calf_n == 2 ~ "twin",
      dam_date_record_n >= 3 ~ "multiple_or_duplicate_check",
      TRUE ~ "multiple_or_duplicate_check"
    )
  ) |>
  dplyr::ungroup()

dat <- dat |>
  dplyr::left_join(
    dedup |> dplyr::select(row_id, dam_date_record_n, dam_date_unique_calf_n, birth_type),
    by = "row_id"
  ) |>
  dplyr::mutate(
    birth_type = dplyr::if_else(obvious_duplicate_exclude_flag, "duplicate_excluded", birth_type),
    birth_type = factor(birth_type, levels = c("singleton", "twin", "multiple_or_duplicate_check", "duplicate_excluded"), ordered = FALSE),
    calf_birth_season = factor(calf_birth_season, levels = c("Spring", "Summer", "Autumn", "Winter"), ordered = FALSE),
    dam_birth_season = factor(dam_birth_season, levels = c("Spring", "Summer", "Autumn", "Winter"), ordered = FALSE),
    calf_sex_code = factor(calf_sex_code, levels = c("G-coded", "non-G-coded"), ordered = FALSE),
    calf_sex = factor(calf_sex, levels = c("male", "female"), ordered = FALSE)
  )

all_valid_records <- dat |>
  dplyr::filter(
    !is.na(birth_weight), !is.na(calving_date), !is.na(dam_id),
    !is.na(calf_id), !is.na(calf_sex_code),
    !obvious_duplicate_exclude_flag,
    birth_type != "multiple_or_duplicate_check"
  )

singleton_analysis <- all_valid_records |>
  dplyr::filter(birth_type == "singleton")

twin_records <- all_valid_records |>
  dplyr::filter(birth_type == "twin")

duplicate_or_multiple_check <- dat |>
  dplyr::filter(
    obvious_duplicate_exclude_flag |
      duplicate_calf_id_flag |
      same_dam_date_same_birth_weight_flag |
      birth_type %in% c("multiple_or_duplicate_check", "duplicate_excluded")
  )

write_csv_safe(dat, file.path(PROCESSED_DIR, "calving_clean.csv"))
write_csv_safe(all_valid_records, file.path(PROCESSED_DIR, "all_valid_records.csv"))
write_csv_safe(singleton_analysis, file.path(PROCESSED_DIR, "singleton_analysis.csv"))
write_csv_safe(twin_records, file.path(PROCESSED_DIR, "twin_records.csv"))
write_csv_safe(duplicate_or_multiple_check, file.path(PROCESSED_DIR, "duplicate_or_multiple_check.csv"))

potential_duplicate_records <- duplicate_or_multiple_check |>
  dplyr::arrange(dam_id, calving_date, calf_id, row_id)
safe_write_xlsx(potential_duplicate_records, file.path(TABLE_DIR, "potential_duplicate_records.xlsx"))

birth_type_distribution <- dat |>
  dplyr::count(birth_type, name = "n") |>
  dplyr::mutate(percent = n / sum(n))
safe_write_xlsx(birth_type_distribution, file.path(TABLE_DIR, "birth_type_distribution.xlsx"))

analysis_dataset_flow <- tibble::tibble(
  step = c(
    "Raw records",
    "Obvious duplicate copies excluded",
    "Valid records for all-record analyses",
    "Singleton main-analysis records",
    "Twin records for supplementary description",
    "Multiple or duplicate-check records"
  ),
  n = c(
    nrow(dat),
    sum(dat$obvious_duplicate_exclude_flag, na.rm = TRUE),
    nrow(all_valid_records),
    nrow(singleton_analysis),
    nrow(twin_records),
    sum(dat$birth_type == "multiple_or_duplicate_check", na.rm = TRUE)
  )
)
safe_write_xlsx(analysis_dataset_flow, file.path(TABLE_DIR, "analysis_dataset_flow.xlsx"))

age_interpretation_check <- tibble::tibble(
  metric = c(
    "Records with non-missing original dam_age_month",
    "Export-date based estimated age at calving < 0 months",
    "Export-date based estimated age at calving < 15 months",
    "Primary dam birth-date rule used in analysis"
  ),
  value = c(
    sum(!is.na(dat$dam_age_month)),
    sum(dat$dam_age_month_export_based_at_calving < 0, na.rm = TRUE),
    sum(dat$dam_age_month_export_based_at_calving < 15, na.rm = TRUE),
    "calving_date minus original dam_age_month"
  )
)
safe_write_xlsx(age_interpretation_check, file.path(TABLE_DIR, "dam_age_interpretation_check.xlsx"))

write_log(paste("01_data_cleaning.R completed. Raw n =", nrow(dat),
                "singleton n =", nrow(singleton_analysis)))

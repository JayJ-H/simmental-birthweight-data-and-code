source("scripts/00_config.R")

write_log("06_tail_phenotypes.R started")

singleton_analysis <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis.csv")) |> prepare_analysis_types()
sex_var <- sex_variable_name()

threshold_probs <- c(P5 = 0.05, P10 = 0.10, P15 = 0.15, P85 = 0.85, P90 = 0.90, P95 = 0.95)

sex_thresholds <- singleton_analysis |>
  dplyr::group_by(.data[[sex_var]]) |>
  dplyr::summarise(
    n = dplyr::n(),
    P5 = stats::quantile(birth_weight, 0.05, na.rm = TRUE, names = FALSE),
    P10 = stats::quantile(birth_weight, 0.10, na.rm = TRUE, names = FALSE),
    P15 = stats::quantile(birth_weight, 0.15, na.rm = TRUE, names = FALSE),
    P85 = stats::quantile(birth_weight, 0.85, na.rm = TRUE, names = FALSE),
    P90 = stats::quantile(birth_weight, 0.90, na.rm = TRUE, names = FALSE),
    P95 = stats::quantile(birth_weight, 0.95, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )
names(sex_thresholds)[1] <- sex_var

overall_thresholds <- singleton_analysis |>
  dplyr::summarise(
    n = dplyr::n(),
    P5 = stats::quantile(birth_weight, 0.05, na.rm = TRUE, names = FALSE),
    P10 = stats::quantile(birth_weight, 0.10, na.rm = TRUE, names = FALSE),
    P15 = stats::quantile(birth_weight, 0.15, na.rm = TRUE, names = FALSE),
    P85 = stats::quantile(birth_weight, 0.85, na.rm = TRUE, names = FALSE),
    P90 = stats::quantile(birth_weight, 0.90, na.rm = TRUE, names = FALSE),
    P95 = stats::quantile(birth_weight, 0.95, na.rm = TRUE, names = FALSE)
  )

threshold_lookup <- sex_thresholds |>
  dplyr::select(dplyr::all_of(sex_var), P5, P10, P15, P85, P90, P95)

singleton_tail <- singleton_analysis |>
  dplyr::left_join(threshold_lookup, by = sex_var) |>
  dplyr::mutate(
    lower_tail = as.integer(birth_weight <= P10),
    upper_tail = as.integer(birth_weight >= P90),
    lower_tail_p5 = as.integer(birth_weight <= P5),
    upper_tail_p95 = as.integer(birth_weight >= P95),
    lower_tail_p15 = as.integer(birth_weight <= P15),
    upper_tail_p85 = as.integer(birth_weight >= P85),
    lower_tail_overall = as.integer(birth_weight <= overall_thresholds$P10[[1]]),
    upper_tail_overall = as.integer(birth_weight >= overall_thresholds$P90[[1]]),
    low_absolute_30 = as.integer(birth_weight < 30),
    high_absolute_55 = as.integer(birth_weight >= 55),
    high_absolute_60 = as.integer(birth_weight > 60)
  ) |>
  dplyr::select(-P5, -P10, -P15, -P85, -P90, -P95)

event_summary <- tibble::tibble(
  phenotype = c("lower_tail", "upper_tail", "lower_tail_p5", "upper_tail_p95", "lower_tail_p15", "upper_tail_p85", "lower_tail_overall", "upper_tail_overall", "low_absolute_30", "high_absolute_55", "high_absolute_60"),
  events = vapply(phenotype, function(v) sum(singleton_tail[[v]], na.rm = TRUE), numeric(1)),
  n = nrow(singleton_tail),
  event_rate = events / n
)

year_season_event_check <- singleton_tail |>
  dplyr::group_by(calf_birth_year, calf_birth_season, .drop = FALSE) |>
  dplyr::summarise(
    n = dplyr::n(),
    lower_events = sum(lower_tail, na.rm = TRUE),
    lower_nonevents = n - lower_events,
    upper_events = sum(upper_tail, na.rm = TRUE),
    upper_nonevents = n - upper_events,
    lower_sparse_flag = lower_events < 5 | lower_nonevents < 5 | n < 30,
    upper_sparse_flag = upper_events < 5 | upper_nonevents < 5 | n < 30,
    .groups = "drop"
  )

write_csv_safe(singleton_tail, file.path(PROCESSED_DIR, "singleton_analysis_with_tail_phenotypes.csv"))
safe_write_xlsx(
  list(
    sex_specific_thresholds = sex_thresholds,
    overall_thresholds = overall_thresholds,
    event_summary = event_summary,
    notes = tibble::tibble(
      item = c("primary_lower_tail", "primary_upper_tail", "absolute_cutoffs"),
      definition = c(
        "birth_weight <= sex-specific P10",
        "birth_weight >= sex-specific P90",
        "Supplementary descriptive thresholds only; not primary phenotype definitions"
      )
    )
  ),
  file.path(TABLE_DIR, "sex_specific_thresholds.xlsx")
)
safe_write_xlsx(year_season_event_check, file.path(TABLE_DIR, "year_season_event_check.xlsx"))

write_log("06_tail_phenotypes.R completed")

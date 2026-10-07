source("scripts/00_config.R")

write_log("02_qc_summary.R started", file.path(LOG_DIR, "qc_log.txt"))

calving_clean <- read_csv_safe(file.path(PROCESSED_DIR, "calving_clean.csv")) |> prepare_analysis_types()
all_valid_records <- read_csv_safe(file.path(PROCESSED_DIR, "all_valid_records.csv")) |> prepare_analysis_types()
singleton_analysis <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis.csv")) |> prepare_analysis_types()

qc_fields <- c(
  "dam_id", "dam_age_month", "parity", "calving_date", "service_sire",
  "birth_weight", "calf_id", "maternal_grandsire", "maternal_granddam"
)

missingness <- purrr::map_dfr(qc_fields, function(field) {
  tibble::tibble(
    field = field,
    missing_n = sum(is.na(calving_clean[[field]])),
    total_n = nrow(calving_clean),
    missing_percent = missing_n / total_n
  )
})
safe_write_xlsx(missingness, file.path(TABLE_DIR, "missingness.xlsx"))

extreme_birth_weight_records <- calving_clean |>
  dplyr::filter(extreme_low_flag | extreme_high_flag | implausible_flag) |>
  dplyr::arrange(birth_weight)
safe_write_xlsx(extreme_birth_weight_records, file.path(TABLE_DIR, "extreme_birth_weight_records.xlsx"))

date_summary <- tibble::tibble(
  metric = c("earliest_calving_date", "latest_calving_date", "raw_records", "all_valid_records", "singleton_records", "twin_records"),
  value = c(
    as.character(min(calving_clean$calving_date, na.rm = TRUE)),
    as.character(max(calving_clean$calving_date, na.rm = TRUE)),
    nrow(calving_clean),
    nrow(all_valid_records),
    nrow(singleton_analysis),
    sum(all_valid_records$birth_type == "twin", na.rm = TRUE)
  )
)

year_counts <- calving_clean |>
  dplyr::count(calf_birth_year, name = "n") |>
  dplyr::arrange(calf_birth_year)

month_counts <- calving_clean |>
  dplyr::count(calf_birth_month, name = "n") |>
  dplyr::arrange(calf_birth_month)

unique_counts <- purrr::map_dfr(c("dam_id", "calf_id", "service_sire", "maternal_grandsire", "maternal_granddam"), function(field) {
  tibble::tibble(
    field = field,
    unique_non_missing = dplyr::n_distinct(calving_clean[[field]], na.rm = TRUE),
    missing_n = sum(is.na(calving_clean[[field]]))
  )
})

duplicate_summary <- tibble::tibble(
  metric = c(
    "exact_duplicate_record_members",
    "exact_duplicate_excluded_copies",
    "duplicate_dam_date_calf_id_members",
    "duplicate_dam_date_calf_id_excluded_copies",
    "duplicate_calf_id_members",
    "same_dam_date_same_birth_weight_members",
    "export_based_estimated_age_at_calving_lt0",
    "export_based_estimated_age_at_calving_lt15"
  ),
  n = c(
    sum(calving_clean$exact_duplicate_record_flag, na.rm = TRUE),
    sum(calving_clean$exact_duplicate_exclude_flag, na.rm = TRUE),
    sum(calving_clean$duplicate_dam_date_calf_id_flag, na.rm = TRUE),
    sum(calving_clean$duplicate_dam_date_calf_id_exclude_flag, na.rm = TRUE),
    sum(calving_clean$duplicate_calf_id_flag, na.rm = TRUE),
    sum(calving_clean$same_dam_date_same_birth_weight_flag, na.rm = TRUE),
    sum(calving_clean$dam_age_month_export_based_at_calving < 0, na.rm = TRUE),
    sum(calving_clean$dam_age_month_export_based_at_calving < 15, na.rm = TRUE)
  )
)

safe_write_xlsx(
  list(
    date_summary = date_summary,
    year_counts = year_counts,
    month_counts = month_counts,
    unique_counts = unique_counts,
    duplicate_summary = duplicate_summary
  ),
  file.path(TABLE_DIR, "qc_summary.xlsx")
)

sex_var <- sex_variable_name()
sex_distribution <- all_valid_records |>
  dplyr::count(.data[[sex_var]], name = "n") |>
  dplyr::mutate(percent = n / sum(n), sex_rule_confirmed = SEX_RULE_CONFIRMED)
names(sex_distribution)[1] <- sex_var
safe_write_xlsx(sex_distribution, file.path(TABLE_DIR, "sex_distribution.xlsx"))

dam_age_quantiles <- singleton_analysis |>
  dplyr::summarise(
    min = min(dam_age_month, na.rm = TRUE),
    P5 = stats::quantile(dam_age_month, 0.05, na.rm = TRUE, names = FALSE),
    P10 = stats::quantile(dam_age_month, 0.10, na.rm = TRUE, names = FALSE),
    P25 = stats::quantile(dam_age_month, 0.25, na.rm = TRUE, names = FALSE),
    median = stats::median(dam_age_month, na.rm = TRUE),
    P75 = stats::quantile(dam_age_month, 0.75, na.rm = TRUE, names = FALSE),
    P90 = stats::quantile(dam_age_month, 0.90, na.rm = TRUE, names = FALSE),
    P95 = stats::quantile(dam_age_month, 0.95, na.rm = TRUE, names = FALSE),
    max = max(dam_age_month, na.rm = TRUE)
  )

parity_counts <- singleton_analysis |>
  dplyr::count(parity_group, name = "n") |>
  dplyr::mutate(percent = n / sum(n))

dam_age_group_counts <- singleton_analysis |>
  dplyr::count(dam_age_group, name = "n") |>
  dplyr::mutate(percent = n / sum(n))

parity_age_cross <- singleton_analysis |>
  dplyr::count(parity_group, dam_age_group, name = "n") |>
  tidyr::pivot_wider(names_from = dam_age_group, values_from = n, values_fill = 0)

spearman_data <- singleton_analysis |>
  dplyr::filter(!is.na(parity), !is.na(dam_age_month))
spearman_test <- suppressWarnings(stats::cor.test(spearman_data$parity, spearman_data$dam_age_month, method = "spearman"))
spearman_summary <- tibble::tibble(
  comparison = "parity vs dam_age_month",
  spearman_rho = unname(spearman_test$estimate),
  p_value = spearman_test$p.value,
  interpretation = dplyr::case_when(
    abs(spearman_rho) >= 0.80 ~ "strong collinearity; prefer parity_group in main model and move dam-age term to sensitivity analysis",
    abs(spearman_rho) >= 0.60 ~ "moderate collinearity; compare model variants",
    TRUE ~ "limited collinearity"
  )
)

fit_lm_variant <- function(terms, name) {
  terms <- c(sex_var, terms, "calf_birth_year", "calf_birth_season")
  formula <- safe_formula("birth_weight", terms)
  data <- singleton_analysis |>
    dplyr::select(dplyr::all_of(all.vars(formula))) |>
    stats::na.omit()
  if (nrow(data) < 100) {
    return(list(summary = tibble::tibble(model = name, formula = deparse(formula), n = nrow(data), AIC = NA_real_, BIC = NA_real_, max_vif = NA_real_, status = "insufficient complete cases"), vif = tibble::tibble()))
  }
  fit_model <- stats::lm(formula, data = data)
  model_aic <- stats::AIC(fit_model)
  model_bic <- stats::BIC(fit_model)
  vif_raw <- tryCatch(car::vif(fit_model), error = function(e) e)
  if (inherits(vif_raw, "error")) {
    vif_table <- tibble::tibble(term = NA_character_, vif_metric = NA_real_, note = vif_raw$message)
  } else if (is.matrix(vif_raw)) {
    vif_table <- as.data.frame(vif_raw) |>
      tibble::rownames_to_column("term") |>
      dplyr::mutate(vif_metric = GVIF^(1 / (2 * Df)), note = "GVIF adjusted") |>
      dplyr::select(term, vif_metric, note)
  } else {
    vif_table <- tibble::tibble(term = names(vif_raw), vif_metric = as.numeric(vif_raw), note = "VIF")
  }
  list(
    summary = tibble::tibble(
      model = name,
      formula = paste(deparse(formula), collapse = " "),
      n = nrow(data),
      AIC = model_aic,
      BIC = model_bic,
      max_vif = suppressWarnings(max(vif_table$vif_metric, na.rm = TRUE)),
      status = "fit"
    ),
    vif = dplyr::mutate(vif_table, model = name, .before = 1)
  )
}

variants <- list(
  A_parity_plus_age_group = c("parity_group", "dam_age_group"),
  B_parity_plus_age_cont = c("parity_group", "dam_age_month_cont"),
  C_parity_only = c("parity_group"),
  D_age_group_only = c("dam_age_group")
)
variant_fits <- purrr::imap(variants, fit_lm_variant)
variant_summary <- purrr::map_dfr(variant_fits, "summary") |>
  dplyr::mutate(
    recommendation = dplyr::if_else(
      model == "C_parity_only" & abs(spearman_summary$spearman_rho[[1]]) >= 0.80,
      "Recommended for main model because parity and dam age are strongly collinear",
      "Sensitivity or alternative model"
    )
  )
variant_vif <- purrr::map_dfr(variant_fits, "vif")

safe_write_xlsx(
  list(
    dam_age_quantiles = dam_age_quantiles,
    parity_counts = parity_counts,
    dam_age_group_counts = dam_age_group_counts,
    parity_age_cross = parity_age_cross,
    spearman = spearman_summary
  ),
  file.path(TABLE_DIR, "dam_age_parity_summary.xlsx")
)

safe_write_xlsx(
  list(
    model_variant_summary = variant_summary,
    vif_by_model = variant_vif,
    recommendation = spearman_summary
  ),
  file.path(TABLE_DIR, "collinearity_report.xlsx")
)

write_log(paste("Raw records:", nrow(calving_clean)), file.path(LOG_DIR, "qc_log.txt"))
write_log(paste("Singleton analysis records:", nrow(singleton_analysis)), file.path(LOG_DIR, "qc_log.txt"))
write_log(paste("Missing dam_age_month:", sum(is.na(calving_clean$dam_age_month))), file.path(LOG_DIR, "qc_log.txt"))
write_log(paste("Parity-age Spearman rho:", round(spearman_summary$spearman_rho, 3)), file.path(LOG_DIR, "qc_log.txt"))
write_log("02_qc_summary.R completed", file.path(LOG_DIR, "qc_log.txt"))

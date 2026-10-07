source("scripts/00_config.R")

log_file <- file.path(LOG_DIR, "lmm_model_log.txt")
write_log("05_lmm_birth_weight.R started", log_file)

singleton_analysis <- read_csv_safe(file.path(PROCESSED_DIR, "singleton_analysis.csv")) |> prepare_analysis_types()
sex_var <- sex_variable_name()

lmm_terms <- list(
  main_parity_only = c(sex_var, "parity_group", "calf_birth_year", "calf_birth_season"),
  parity_plus_age_group = c(sex_var, "parity_group", "dam_age_group", "calf_birth_year", "calf_birth_season"),
  parity_plus_age_cont = c(sex_var, "parity_group", "dam_age_month_cont", "calf_birth_year", "calf_birth_season"),
  age_group_only = c(sex_var, "dam_age_group", "calf_birth_year", "calf_birth_season")
)

fit_lmer_safe <- function(name, fixed_terms, REML = FALSE) {
  formula <- safe_formula("birth_weight", fixed_terms, random_terms = "(1 | dam_id)")
  model_data <- singleton_analysis |>
    dplyr::select(dplyr::all_of(all.vars(formula))) |>
    stats::na.omit()
  result <- tryCatch(
    lmerTest::lmer(
      formula,
      data = model_data,
      REML = REML,
      control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
    ),
    error = function(e) e,
    warning = function(w) {
      write_log(paste(name, "warning:", conditionMessage(w)), log_file)
      invokeRestart("muffleWarning")
    }
  )
  list(name = name, formula = formula, data = model_data, model = result)
}

ml_fits <- purrr::imap(lmm_terms, ~ fit_lmer_safe(.y, .x, REML = FALSE))

comparison <- purrr::map_dfr(ml_fits, function(fit) {
  if (inherits(fit$model, "error")) {
    return(tibble::tibble(model = fit$name, formula = paste(deparse(fit$formula), collapse = " "), n = nrow(fit$data), AIC = NA_real_, BIC = NA_real_, logLik = NA_real_, singular = NA, status = fit$model$message))
  }
  tibble::tibble(
    model = fit$name,
    formula = paste(deparse(fit$formula), collapse = " "),
    n = nrow(fit$data),
    AIC = stats::AIC(fit$model),
    BIC = stats::BIC(fit$model),
    logLik = as.numeric(stats::logLik(fit$model)),
    singular = lme4::isSingular(fit$model, tol = 1e-4),
    status = "fit"
  )
}) |>
  dplyr::arrange(AIC)

selected_model_name <- if ("main_parity_only" %in% comparison$model && any(comparison$model == "main_parity_only" & comparison$status == "fit")) {
  "main_parity_only"
} else {
  comparison$model[which(comparison$status == "fit")][[1]]
}

selected_terms <- lmm_terms[[selected_model_name]]
final_fit <- fit_lmer_safe(paste0(selected_model_name, "_REML"), selected_terms, REML = TRUE)
if (inherits(final_fit$model, "error")) stop("Final LMM failed: ", final_fit$model$message, call. = FALSE)

fixed_effects <- broom.mixed::tidy(final_fit$model, effects = "fixed", conf.int = TRUE) |>
  dplyr::mutate(model = selected_model_name, .before = 1)

random_effects <- broom.mixed::tidy(final_fit$model, effects = "ran_pars", conf.int = FALSE) |>
  dplyr::mutate(model = selected_model_name, .before = 1)

anova_table <- as.data.frame(stats::anova(final_fit$model, type = 3)) |>
  tibble::rownames_to_column("term") |>
  dplyr::mutate(model = selected_model_name, .before = 1)

emm_terms <- selected_terms[!selected_terms %in% c("dam_age_month_cont")]
emm_table <- purrr::map_dfr(emm_terms, function(term) {
  if (is.numeric(final_fit$data[[term]])) return(tibble::tibble())
  emm <- tryCatch(emmeans::emmeans(final_fit$model, stats::as.formula(paste("~", term))), error = function(e) e)
  if (inherits(emm, "error")) return(tibble::tibble(term = term, status = emm$message))
  as.data.frame(emm) |>
    dplyr::mutate(term = term, status = "ok", .before = 1)
})

diagnostics <- tibble::tibble(
  metric = c("selected_model", "n_complete_cases", "is_singular", "sigma", "REML_logLik"),
  value = c(
    selected_model_name,
    nrow(final_fit$data),
    as.character(lme4::isSingular(final_fit$model, tol = 1e-4)),
    as.character(stats::sigma(final_fit$model)),
    as.character(as.numeric(stats::logLik(final_fit$model)))
  )
)

safe_write_xlsx(comparison, file.path(MODEL_DIR, "lmm_model_comparison.xlsx"))
safe_write_xlsx(
  list(
    diagnostics = diagnostics,
    fixed_effects = fixed_effects,
    type_III_ANOVA = anova_table,
    random_effects = random_effects,
    estimated_marginal_means = emm_table
  ),
  file.path(MODEL_DIR, "lmm_birth_weight_results.xlsx")
)

write_csv_safe(emm_table, file.path(PROCESSED_DIR, "lmm_emmeans_for_figures.csv"))
saveRDS(final_fit$model, file.path(MODEL_DIR, "lmm_birth_weight_final.rds"))

write_log(paste("Selected LMM:", selected_model_name), log_file)
write_log("05_lmm_birth_weight.R completed", log_file)

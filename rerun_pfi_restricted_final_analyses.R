#!/usr/bin/env Rscript

# Final PFI-restricted analyses using the 20-gene V4 subject-to-XCI proxy score.
#
# Analyses:
#   1. TP53 mutation logistic models (unadjusted, adjusted, quartile trend, LRTs)
#   2. TP53-by-XCI genomic-instability analyses for four outcomes
#   3. PFI Cox models using PFI_time/PFI_event only (never OS_time/OS_event)
#
# Patients are classified as TP53 wild-type only if they are represented in
# the GDC MAF. All analysis cohorts require a valid PFI time and event.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(broom)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(MASS)
})

options(stringsAsFactors = FALSE)
set.seed(2026)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else getwd()
project_root <- normalizePath(project_root, mustWork = TRUE)
output_dir <- if (length(args) >= 2) args[[2]] else file.path(
  project_root,
  "outputs",
  "pfi_restricted_final_analyses"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- function(...) file.path(project_root, ...)

required_files <- c(
  xci = input_path("outputs", "RQ1_V4_XCI_linked", "rq1_v4_plotting_analysis_dataset.csv"),
  pfi = input_path("data", "tcga_brca_cdr_survival_endpoints.csv"),
  maf = input_path("data", "tcga_brca_mutect2_maf.rds"),
  rq3 = input_path("outputs", "Aim3_RQ3_genomic_instability", "rq3_complete_analysis_dataset.csv"),
  clinical = input_path("data", "tcga_brca_clinical_raw.csv"),
  pam50 = input_path("data", "xist_pam50_boxplot_data.csv")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

clean_id <- function(x) str_to_upper(substr(str_trim(as.character(x)), 1, 12))

clean_stage <- function(x) {
  y <- str_to_lower(str_squish(str_replace_all(as.character(x), "[_-]", " ")))
  case_when(
    is.na(y) | y == "" | y %in% c("not reported", "not available", "unknown") ~ NA_character_,
    str_detect(y, "stage iv") | y %in% c("iv", "4") ~ "Stage IV",
    str_detect(y, "stage iii") | y %in% c("iii", "3") ~ "Stage III",
    str_detect(y, "stage ii") | y %in% c("ii", "2") ~ "Stage II",
    str_detect(y, "stage i") | y %in% c("i", "1") ~ "Stage I",
    TRUE ~ NA_character_
  )
}

clean_pam50 <- function(x) {
  y <- str_to_lower(str_squish(str_replace_all(as.character(x), "_", " ")))
  case_when(
    y %in% c("luma", "luminal a", "luminal-a") ~ "LumA",
    y %in% c("lumb", "luminal b", "luminal-b") ~ "LumB",
    y %in% c("her2", "her2 enriched", "her2-enriched") ~ "HER2",
    y %in% c("basal", "basal like", "basal-like") ~ "Basal",
    y %in% c("normal", "normal like", "normal-like") ~ "Normal",
    TRUE ~ NA_character_
  )
}

scale_numeric <- function(x) as.numeric(scale(as.numeric(x)))

wald_or_table <- function(model, model_name) {
  model_n <- stats::nobs(model)
  model_events <- sum(model$y == 1, na.rm = TRUE)
  broom::tidy(model) %>%
    mutate(
      model = model_name,
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error),
      n = model_n,
      events = model_events
    ) %>%
    dplyr::select(model, term, n, events, estimate, std.error, statistic, p.value,
                  OR, CI_low, CI_high)
}

robust_table <- function(model, model_id, outcome, xci_parameterisation,
                         adjustment_set, exponentiate = FALSE) {
  vc <- sandwich::vcovHC(model, type = "HC3")
  ct <- lmtest::coeftest(model, vcov. = vc)
  ci <- lmtest::coefci(model, vcov. = vc)
  out <- tibble(
    model_id = model_id,
    outcome = outcome,
    xci_parameterisation = xci_parameterisation,
    adjustment_set = adjustment_set,
    term = rownames(ct),
    estimate = unname(ct[, 1]),
    robust_se = unname(ct[, 2]),
    statistic = unname(ct[, 3]),
    p.value = unname(ct[, 4]),
    CI_low = unname(ci[, 1]),
    CI_high = unname(ci[, 2]),
    n = stats::nobs(model)
  )
  if (exponentiate) {
    out <- out %>%
      mutate(effect = exp(estimate), effect_CI_low = exp(CI_low),
             effect_CI_high = exp(CI_high), effect_scale = "IRR")
  } else {
    out <- out %>%
      mutate(effect = estimate, effect_CI_low = CI_low,
             effect_CI_high = CI_high, effect_scale = "difference")
  }
  out
}

extract_interaction <- function(tab) {
  tab %>% filter(str_detect(term, "TP53_mut:|:TP53_mut"))
}

pairwise_wilcox_tidy <- function(data, outcome) {
  f <- reformulate("combined_group", response = outcome)
  pw <- pairwise.wilcox.test(
    x = data[[outcome]],
    g = data$combined_group,
    p.adjust.method = "BH",
    exact = FALSE
  )
  mat <- pw$p.value
  if (is.null(mat)) return(tibble())
  idx <- which(!is.na(mat), arr.ind = TRUE)
  tibble(
    outcome = outcome,
    group1 = colnames(mat)[idx[, "col"]],
    group2 = rownames(mat)[idx[, "row"]],
    p_adjusted_BH = mat[idx]
  )
}

write_model_lrt <- function(model0, model1, filename, comparison) {
  tab <- as.data.frame(anova(model0, model1, test = "Chisq")) %>%
    rownames_to_column("model_row") %>%
    mutate(comparison = comparison)
  write_csv(tab, file.path(output_dir, filename))
  tab
}

# -----------------------------------------------------------------------------
# Load and harmonise core inputs
# -----------------------------------------------------------------------------

xci_raw <- read_csv(required_files[["xci"]], show_col_types = FALSE)
pfi_raw <- read_csv(required_files[["pfi"]], show_col_types = FALSE)
maf_raw <- readRDS(required_files[["maf"]])
rq3_raw <- read_csv(required_files[["rq3"]], show_col_types = FALSE)
clinical_raw <- read_csv(required_files[["clinical"]], show_col_types = FALSE)
pam50_raw <- read_csv(required_files[["pam50"]], show_col_types = FALSE)

non_silent_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
  "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site",
  "Translation_Start_Site", "Nonstop_Mutation"
)

maf_patients <- maf_raw %>%
  transmute(patient_barcode = clean_id(Tumor_Sample_Barcode)) %>%
  filter(!is.na(patient_barcode), patient_barcode != "") %>%
  distinct()

tp53_patients <- maf_raw %>%
  filter(Hugo_Symbol == "TP53", Variant_Classification %in% non_silent_classes) %>%
  transmute(patient_barcode = clean_id(Tumor_Sample_Barcode), TP53_mut = 1L) %>%
  distinct()

pfi <- pfi_raw %>%
  transmute(
    patient_barcode = clean_id(patient_barcode),
    age_years = as.numeric(age_years),
    PFI_time = as.numeric(PFI_time),
    PFI_event = as.integer(PFI_event)
  ) %>%
  filter(!is.na(PFI_time), PFI_time > 0, PFI_event %in% c(0L, 1L)) %>%
  distinct(patient_barcode, .keep_all = TRUE)

clinical_covariates <- clinical_raw %>%
  transmute(
    patient_barcode = clean_id(submitter_id),
    stage_clinical_raw = clean_stage(ajcc_pathologic_stage)
  ) %>%
  distinct(patient_barcode, .keep_all = TRUE)

pam50_covariates <- pam50_raw %>%
  transmute(
    patient_barcode = clean_id(patient_barcode),
    pam50_independent = clean_pam50(PAM50)
  ) %>%
  filter(!is.na(patient_barcode), patient_barcode != "") %>%
  distinct(patient_barcode, .keep_all = TRUE)

xci <- xci_raw %>%
  transmute(
    patient_barcode = clean_id(patient_barcode),
    xci_score = as.numeric(xci_erosion_score_v4_xci_linked_gene_only),
    xist_log_cpm = as.numeric(xist_log_cpm),
    stage_embedded = clean_stage(stage_clean),
    pam50_embedded = clean_pam50(pam50)
  ) %>%
  distinct(patient_barcode, .keep_all = TRUE) %>%
  left_join(clinical_covariates, by = "patient_barcode") %>%
  left_join(pam50_covariates, by = "patient_barcode") %>%
  mutate(
    # Preserve the previously used values where present, but recover patients
    # lost by the old combined stage/PAM50 merge from the independent sources.
    stage = coalesce(stage_embedded, stage_clinical_raw),
    pam50 = coalesce(pam50_embedded, pam50_independent)
  ) %>%
  dplyr::select(-stage_embedded, -pam50_embedded,
                -stage_clinical_raw, -pam50_independent)

input_audit <- tibble(
  item = c(
    "V4 score rows", "Valid PFI rows", "MAF-covered patients",
    "TP53-mutant patients in MAF", "RQ3 input rows",
    "Independent clinical-stage patients", "Independent PAM50 patients"
  ),
  n = c(nrow(xci), nrow(pfi), nrow(maf_patients), nrow(tp53_patients), nrow(rq3_raw),
        sum(!is.na(clinical_covariates$stage_clinical_raw)),
        sum(!is.na(pam50_covariates$pam50_independent)))
)
write_csv(input_audit, file.path(output_dir, "00_input_audit.csv"))

# -----------------------------------------------------------------------------
# 1. PFI-restricted TP53 logistic analyses
# -----------------------------------------------------------------------------

tp53_pfi_full <- xci %>%
  filter(!is.na(xci_score)) %>%
  inner_join(maf_patients, by = "patient_barcode") %>%
  inner_join(pfi, by = "patient_barcode") %>%
  left_join(tp53_patients, by = "patient_barcode") %>%
  mutate(TP53_mut = coalesce(TP53_mut, 0L))

tp53_cc <- tp53_pfi_full %>%
  filter(!is.na(age_years), !is.na(stage), !is.na(pam50)) %>%
  mutate(
    xci_z = scale_numeric(xci_score),
    xci_quartile_numeric = ntile(xci_score, 4),
    xci_quartile = factor(xci_quartile_numeric, levels = 1:4),
    stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    pam50 = factor(pam50, levels = c("LumA", "LumB", "HER2", "Basal", "Normal"))
  ) %>%
  droplevels()

tp53_pfi_full <- tp53_pfi_full %>% mutate(xci_z = scale_numeric(xci_score))

tp53_cohort_summary <- bind_rows(
  tibble(
    cohort = "PFI-valid V4 + MAF cohort",
    n = nrow(tp53_pfi_full),
    TP53_mutant_n = sum(tp53_pfi_full$TP53_mut == 1),
    PFI_events = sum(tp53_pfi_full$PFI_event == 1)
  ),
  tibble(
    cohort = "Age/stage/PAM50 complete-case cohort",
    n = nrow(tp53_cc),
    TP53_mutant_n = sum(tp53_cc$TP53_mut == 1),
    PFI_events = sum(tp53_cc$PFI_event == 1)
  )
)
write_csv(tp53_cohort_summary, file.path(output_dir, "01_tp53_cohort_summary.csv"))
write_csv(tp53_cc, file.path(output_dir, "02_tp53_complete_case_dataset.csv"))

tp53_full_unadj <- glm(TP53_mut ~ xci_z, data = tp53_pfi_full, family = binomial())
tp53_cc_unadj <- glm(TP53_mut ~ xci_z, data = tp53_cc, family = binomial())
tp53_base <- glm(TP53_mut ~ age_years + stage + pam50,
                 data = tp53_cc, family = binomial())
tp53_adj <- glm(TP53_mut ~ xci_z + age_years + stage + pam50,
                data = tp53_cc, family = binomial())
tp53_quartile <- glm(TP53_mut ~ xci_quartile_numeric + age_years + stage + pam50,
                     data = tp53_cc, family = binomial())

tp53_model_results <- bind_rows(
  wald_or_table(tp53_full_unadj, "Unadjusted full PFI/MAF cohort"),
  wald_or_table(tp53_cc_unadj, "Unadjusted complete-case cohort"),
  wald_or_table(tp53_base, "Clinical base model"),
  wald_or_table(tp53_adj, "Adjusted continuous V4"),
  wald_or_table(tp53_quartile, "Adjusted ordinal V4 quartile trend")
)
write_csv(tp53_model_results, file.path(output_dir, "03_tp53_logistic_all_coefficients.csv"))
write_csv(
  tp53_model_results %>% filter(term %in% c("xci_z", "xci_quartile_numeric")),
  file.path(output_dir, "04_tp53_logistic_key_results.csv")
)

tp53_lrt_cont <- write_model_lrt(
  tp53_base, tp53_adj, "05_tp53_lrt_continuous_v4.csv",
  "Clinical base vs base + continuous V4"
)
tp53_lrt_quartile <- write_model_lrt(
  tp53_base, tp53_quartile, "06_tp53_lrt_quartile_trend.csv",
  "Clinical base vs base + ordinal V4 quartile"
)

# -----------------------------------------------------------------------------
# 2. PFI-restricted TP53-by-XCI genomic-instability analyses
# -----------------------------------------------------------------------------

rq3 <- rq3_raw %>%
  transmute(
    patient_barcode = clean_id(patient_barcode),
    xci_score = as.numeric(xci_score),
    stage_embedded = clean_stage(stage_covariate),
    pam50_embedded = clean_pam50(pam50),
    xist_log_cpm = as.numeric(xist_log_cpm),
    TP53_mut = as.integer(TP53_mut),
    fga_autosomal = as.numeric(fga_autosomal),
    weighted_abs_segment_mean = as.numeric(weighted_abs_segment_mean),
    cnv_segment_count = as.numeric(cnv_segment_count),
    tmb_per_mb = as.numeric(tmb_per_mb)
  ) %>%
  left_join(clinical_covariates, by = "patient_barcode") %>%
  left_join(pam50_covariates, by = "patient_barcode") %>%
  mutate(
    stage = coalesce(stage_embedded, stage_clinical_raw),
    pam50 = coalesce(pam50_embedded, pam50_independent)
  ) %>%
  dplyr::select(-stage_embedded, -pam50_embedded,
                -stage_clinical_raw, -pam50_independent) %>%
  inner_join(pfi, by = "patient_barcode") %>%
  filter(
    !is.na(xci_score), TP53_mut %in% c(0L, 1L),
    is.finite(fga_autosomal), is.finite(weighted_abs_segment_mean),
    is.finite(cnv_segment_count), is.finite(tmb_per_mb)
  ) %>%
  mutate(
    xci_z = scale_numeric(xci_score),
    xci_quartile_numeric = ntile(xci_score, 4),
    high_xci = as.integer(xci_quartile_numeric == 4),
    TP53_status = factor(TP53_mut, levels = c(0, 1), labels = c("WT", "Mutant")),
    combined_group = case_when(
      TP53_mut == 0 & high_xci == 0 ~ "TP53 WT / lower XCI",
      TP53_mut == 0 & high_xci == 1 ~ "TP53 WT / high XCI",
      TP53_mut == 1 & high_xci == 0 ~ "TP53 mutant / lower XCI",
      TP53_mut == 1 & high_xci == 1 ~ "TP53 mutant / high XCI"
    ),
    combined_group = factor(
      combined_group,
      levels = c(
        "TP53 WT / lower XCI", "TP53 WT / high XCI",
        "TP53 mutant / lower XCI", "TP53 mutant / high XCI"
      )
    ),
    fga_percent = 100 * fga_autosomal,
    log1p_tmb = log1p(tmb_per_mb)
  )

# Confirm that the RQ3 score is the same final V4 score.
score_audit <- rq3 %>%
  dplyr::select(patient_barcode, rq3_score = xci_score) %>%
  inner_join(xci %>% dplyr::select(patient_barcode, v4_score = xci_score), by = "patient_barcode")
score_correlation <- cor(score_audit$rq3_score, score_audit$v4_score,
                         use = "complete.obs", method = "pearson")
if (!is.finite(score_correlation) || score_correlation < 0.999999) {
  stop("RQ3 xci_score does not match the final V4 score; correlation = ", score_correlation)
}

write_csv(rq3, file.path(output_dir, "07_rq3_pfi_analysis_dataset.csv"))
write_csv(
  rq3 %>% count(combined_group, name = "n") %>%
    mutate(total_n = nrow(rq3), PFI_events = sum(rq3$PFI_event == 1)),
  file.path(output_dir, "08_rq3_four_group_counts.csv")
)

outcome_map <- tribble(
  ~outcome, ~label,
  "fga_autosomal", "Autosomal FGA",
  "weighted_abs_segment_mean", "Weighted absolute CNV magnitude",
  "cnv_segment_count", "Autosomal CNV segment count",
  "tmb_per_mb", "Tumour mutation burden"
)

rq3_long <- rq3 %>%
  pivot_longer(cols = all_of(outcome_map$outcome), names_to = "outcome", values_to = "value") %>%
  left_join(outcome_map, by = "outcome")

rq3_descriptive <- rq3_long %>%
  group_by(outcome, label, combined_group) %>%
  summarise(
    n = n(), median = median(value), q1 = quantile(value, 0.25),
    q3 = quantile(value, 0.75), mean = mean(value), .groups = "drop"
  )
write_csv(rq3_descriptive, file.path(output_dir, "09_rq3_four_group_descriptive_summary.csv"))

rq3_kw <- map_dfr(outcome_map$outcome, function(outcome) {
  test <- kruskal.test(reformulate("combined_group", response = outcome), data = rq3)
  tibble(
    outcome = outcome,
    n = nrow(rq3),
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p.value = test$p.value,
    epsilon_squared = max(0, (unname(test$statistic) - 4 + 1) / (nrow(rq3) - 4))
  )
}) %>% mutate(FDR_q = p.adjust(p.value, method = "BH"))
write_csv(rq3_kw, file.path(output_dir, "10_rq3_kruskal_wallis_results.csv"))

rq3_pairwise <- map_dfr(outcome_map$outcome, ~pairwise_wilcox_tidy(rq3, .x))
write_csv(rq3_pairwise, file.path(output_dir, "11_rq3_pairwise_wilcoxon_BH.csv"))

# Model-specific complete-case cohorts prevent XIST missingness from removing
# patients from models that do not include XIST.
factor_rq3_covariates <- function(data) data %>%
  mutate(
    stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    pam50 = factor(pam50, levels = c("LumA", "LumB", "HER2", "Basal", "Normal"))
  ) %>%
  droplevels()

rq3_unadjusted <- factor_rq3_covariates(rq3)
rq3_clinical <- rq3 %>%
  filter(!is.na(age_years), !is.na(stage), !is.na(pam50)) %>%
  factor_rq3_covariates()
rq3_clinical_xist <- rq3_clinical %>%
  filter(!is.na(xist_log_cpm)) %>%
  droplevels()

write_csv(rq3_unadjusted, file.path(output_dir, "12a_rq3_interaction_unadjusted_dataset.csv"))
write_csv(rq3_clinical, file.path(output_dir, "12b_rq3_interaction_clinical_dataset.csv"))
write_csv(rq3_clinical_xist, file.path(output_dir, "12c_rq3_interaction_clinical_xist_dataset.csv"))
write_csv(
  tibble(
    adjustment_set = c("unadjusted", "clinical_adjusted", "clinical_plus_xist"),
    n = c(nrow(rq3_unadjusted), nrow(rq3_clinical), nrow(rq3_clinical_xist)),
    PFI_events = c(sum(rq3_unadjusted$PFI_event), sum(rq3_clinical$PFI_event),
                   sum(rq3_clinical_xist$PFI_event))
  ),
  file.path(output_dir, "12_rq3_interaction_cohort_summary.csv")
)

model_specs <- tribble(
  ~outcome, ~model_type, ~exponentiate,
  "fga_percent", "lm", FALSE,
  "weighted_abs_segment_mean", "lm", FALSE,
  "cnv_segment_count", "negative_binomial", TRUE,
  "log1p_tmb", "lm", FALSE
)

adjustment_sets <- list(
  unadjusted = list(data = rq3_unadjusted, covariates = character(0)),
  clinical_adjusted = list(
    data = rq3_clinical, covariates = c("age_years", "stage", "pam50")
  ),
  clinical_plus_xist = list(
    data = rq3_clinical_xist,
    covariates = c("age_years", "stage", "pam50", "xist_log_cpm")
  )
)

fit_rq3_model <- function(data, outcome, model_type, xci_var, covariates) {
  rhs <- c(paste0("TP53_mut * ", xci_var), covariates)
  formula <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  if (model_type == "lm") {
    lm(formula, data = data)
  } else {
    MASS::glm.nb(formula, data = data)
  }
}

rq3_model_tables <- list()
rq3_interactions <- list()
rq3_models <- list()
counter <- 1L

for (i in seq_len(nrow(model_specs))) {
  spec <- model_specs[i, ]
  for (xci_var in c("xci_z", "high_xci")) {
    for (adj_name in names(adjustment_sets)) {
      model_id <- paste(spec$outcome, xci_var, adj_name, sep = "__")
      analysis_set <- adjustment_sets[[adj_name]]
      fit <- fit_rq3_model(
        analysis_set$data, spec$outcome, spec$model_type, xci_var,
        analysis_set$covariates
      )
      tab <- robust_table(
        fit, model_id, spec$outcome,
        ifelse(xci_var == "xci_z", "continuous", "top_quartile"),
        adj_name, spec$exponentiate
      )
      rq3_model_tables[[counter]] <- tab
      rq3_interactions[[counter]] <- extract_interaction(tab)
      rq3_models[[model_id]] <- fit
      counter <- counter + 1L
    }
  }
}

rq3_coefficients <- bind_rows(rq3_model_tables)
rq3_interaction_results <- bind_rows(rq3_interactions)
write_csv(rq3_coefficients, file.path(output_dir, "13_rq3_all_model_coefficients_robust.csv"))
write_csv(rq3_interaction_results, file.path(output_dir, "14_rq3_interaction_terms_robust.csv"))
saveRDS(rq3_models, file.path(output_dir, "15_rq3_fitted_models.rds"))

# Continuous XCI slopes within TP53 strata, using the robust covariance matrix.
rq3_slopes <- map_dfr(names(rq3_models)[str_detect(names(rq3_models), "__xci_z__")], function(model_id) {
  fit <- rq3_models[[model_id]]
  vc <- sandwich::vcovHC(fit, type = "HC3")
  b <- coef(fit)
  x_term <- "xci_z"
  int_term <- names(b)[str_detect(names(b), "TP53_mut:xci_z|xci_z:TP53_mut")]
  if (!(x_term %in% names(b)) || length(int_term) != 1) return(tibble())
  int_term <- int_term[[1]]
  contrast_wt <- setNames(rep(0, length(b)), names(b)); contrast_wt[x_term] <- 1
  contrast_mut <- contrast_wt; contrast_mut[int_term] <- 1
  make_slope <- function(contrast, stratum) {
    est <- sum(contrast * b)
    se <- sqrt(drop(t(contrast) %*% vc %*% contrast))
    is_nb <- inherits(fit, "negbin")
    tibble(
      model_id = model_id, TP53_stratum = stratum, estimate = est, robust_se = se,
      statistic = est / se, p.value = 2 * pnorm(-abs(est / se)),
      CI_low = est - 1.96 * se, CI_high = est + 1.96 * se,
      effect = if (is_nb) exp(est) else est,
      effect_CI_low = if (is_nb) exp(est - 1.96 * se) else est - 1.96 * se,
      effect_CI_high = if (is_nb) exp(est + 1.96 * se) else est + 1.96 * se,
      effect_scale = if (is_nb) "IRR per 1-SD V4" else "difference per 1-SD V4",
      n = nobs(fit)
    )
  }
  bind_rows(make_slope(contrast_wt, "WT"), make_slope(contrast_mut, "Mutant"))
})
write_csv(rq3_slopes, file.path(output_dir, "16_rq3_continuous_v4_slopes_by_tp53.csv"))

# -----------------------------------------------------------------------------
# 3. V4 prognostic analysis using PFI only
# -----------------------------------------------------------------------------

pfi_core <- xci %>%
  filter(!is.na(xci_score)) %>%
  inner_join(pfi, by = "patient_barcode") %>%
  mutate(
    xci_z = scale_numeric(xci_score),
    age_z = scale_numeric(age_years),
    xci_group = factor(
      if_else(xci_score >= median(xci_score), "High V4", "Low V4"),
      levels = c("Low V4", "High V4")
    ),
    stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    pam50 = factor(pam50, levels = c("LumA", "LumB", "HER2", "Basal", "Normal"))
  ) %>%
  droplevels()

pfi_clinical <- pfi_core %>%
  filter(!is.na(age_years), !is.na(stage), !is.na(pam50)) %>%
  mutate(xci_z = scale_numeric(xci_score), age_z = scale_numeric(age_years)) %>%
  droplevels()

pfi_xist_data <- pfi_clinical %>%
  filter(!is.na(xist_log_cpm)) %>%
  mutate(
    xci_z = scale_numeric(xci_score), age_z = scale_numeric(age_years),
    xist_z = scale_numeric(xist_log_cpm)
  ) %>%
  droplevels()

write_csv(pfi_core, file.path(output_dir, "17a_v4_pfi_core_dataset.csv"))
write_csv(pfi_clinical, file.path(output_dir, "17b_v4_pfi_clinical_dataset.csv"))
write_csv(pfi_xist_data, file.path(output_dir, "17c_v4_pfi_clinical_xist_dataset.csv"))
write_csv(
  tibble(
    endpoint = "PFI",
    cohort = c("V4 score + valid PFI", "Age/stage/PAM50 complete case",
               "Age/stage/PAM50/XIST complete case"),
    n = c(nrow(pfi_core), nrow(pfi_clinical), nrow(pfi_xist_data)),
    events = c(sum(pfi_core$PFI_event), sum(pfi_clinical$PFI_event),
               sum(pfi_xist_data$PFI_event)),
    censored = n - events,
    OS_columns_used = FALSE
  ),
  file.path(output_dir, "18_v4_pfi_cohort_summary.csv")
)

pfi_univ <- coxph(Surv(PFI_time, PFI_event) ~ xci_z,
                  data = pfi_core, x = TRUE)
pfi_base <- coxph(Surv(PFI_time, PFI_event) ~ age_z + stage + pam50,
                  data = pfi_clinical, x = TRUE)
pfi_clinical_v4 <- coxph(Surv(PFI_time, PFI_event) ~ age_z + stage + pam50 + xci_z,
                         data = pfi_clinical, x = TRUE)
pfi_clinical_xist <- coxph(Surv(PFI_time, PFI_event) ~ age_z + stage + pam50 + xist_z,
                           data = pfi_xist_data, x = TRUE)
pfi_full <- coxph(Surv(PFI_time, PFI_event) ~ age_z + stage + pam50 + xist_z + xci_z,
                  data = pfi_xist_data, x = TRUE)
pfi_clinical_spline <- coxph(
  Surv(PFI_time, PFI_event) ~ age_z + stage + pam50 + splines::ns(xci_z, df = 3),
  data = pfi_clinical, x = TRUE
)
pfi_full_spline <- coxph(
  Surv(PFI_time, PFI_event) ~ age_z + stage + pam50 + xist_z +
    splines::ns(xci_z, df = 3),
  data = pfi_xist_data, x = TRUE
)

tidy_cox <- function(model, model_name) {
  # nobs.coxph returns the number of events, not the number of fitted subjects.
  model_n <- length(model$linear.predictors)
  model_events <- model$nevent
  broom::tidy(model) %>%
    mutate(
      model = model_name,
      HR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error),
      n = model_n,
      events = model_events
    ) %>%
    dplyr::select(model, term, n, events, estimate, std.error, statistic, p.value,
                  HR, CI_low, CI_high)
}

pfi_cox_results <- bind_rows(
  tidy_cox(pfi_univ, "Univariable V4"),
  tidy_cox(pfi_base, "Age + stage + PAM50"),
  tidy_cox(pfi_clinical_v4, "Age + stage + PAM50 + V4"),
  tidy_cox(pfi_clinical_xist, "Age + stage + PAM50 + XIST"),
  tidy_cox(pfi_full, "Age + stage + PAM50 + XIST + V4")
)
write_csv(pfi_cox_results, file.path(output_dir, "19_v4_pfi_cox_all_coefficients.csv"))
write_csv(
  pfi_cox_results %>% filter(term == "xci_z"),
  file.path(output_dir, "20_v4_pfi_cox_key_v4_results.csv")
)

pfi_lrt_clinical <- write_model_lrt(
  pfi_base, pfi_clinical_v4, "21_v4_pfi_lrt_add_v4_to_clinical.csv",
  "Age/stage/PAM50 vs age/stage/PAM50 + V4"
)
pfi_lrt_full <- write_model_lrt(
  pfi_clinical_xist, pfi_full, "22_v4_pfi_lrt_add_v4_after_xist.csv",
  "Age/stage/PAM50/XIST vs age/stage/PAM50/XIST + V4"
)
write_model_lrt(
  pfi_base, pfi_clinical_spline, "22a_v4_pfi_lrt_add_spline_to_clinical.csv",
  "Age/stage/PAM50 vs age/stage/PAM50 + 3-df V4 spline"
)
write_model_lrt(
  pfi_clinical_v4, pfi_clinical_spline,
  "22b_v4_pfi_lrt_spline_vs_linear_clinical.csv",
  "Age/stage/PAM50 + linear V4 vs age/stage/PAM50 + 3-df V4 spline"
)
write_model_lrt(
  pfi_clinical_xist, pfi_full_spline,
  "22c_v4_pfi_lrt_add_spline_after_xist.csv",
  "Age/stage/PAM50/XIST vs age/stage/PAM50/XIST + 3-df V4 spline"
)
write_model_lrt(
  pfi_full, pfi_full_spline,
  "22d_v4_pfi_lrt_spline_vs_linear_after_xist.csv",
  "Age/stage/PAM50/XIST + linear V4 vs age/stage/PAM50/XIST + 3-df V4 spline"
)

pfi_model_fit <- tibble(
  model = c("Univariable V4", "Age + stage + PAM50", "Age + stage + PAM50 + V4",
            "Age + stage + PAM50 + XIST", "Age + stage + PAM50 + XIST + V4"),
  n = c(length(pfi_univ$linear.predictors), length(pfi_base$linear.predictors),
        length(pfi_clinical_v4$linear.predictors),
        length(pfi_clinical_xist$linear.predictors), length(pfi_full$linear.predictors)),
  events = c(pfi_univ$nevent, pfi_base$nevent, pfi_clinical_v4$nevent,
             pfi_clinical_xist$nevent, pfi_full$nevent),
  AIC = c(AIC(pfi_univ), AIC(pfi_base), AIC(pfi_clinical_v4),
          AIC(pfi_clinical_xist), AIC(pfi_full)),
  concordance = c(summary(pfi_univ)$concordance[[1]], summary(pfi_base)$concordance[[1]],
                  summary(pfi_clinical_v4)$concordance[[1]],
                  summary(pfi_clinical_xist)$concordance[[1]],
                  summary(pfi_full)$concordance[[1]])
)
write_csv(pfi_model_fit, file.path(output_dir, "23_v4_pfi_model_fit.csv"))

km_fit <- survfit(Surv(PFI_time, PFI_event) ~ xci_group, data = pfi_core)
km_diff <- survdiff(Surv(PFI_time, PFI_event) ~ xci_group, data = pfi_core)
km_cox <- coxph(Surv(PFI_time, PFI_event) ~ xci_group,
                data = pfi_core, x = TRUE)
km_cox_clinical <- coxph(
  Surv(PFI_time, PFI_event) ~ age_z + stage + pam50 + xci_group,
  data = pfi_clinical, x = TRUE
)
km_cox_full <- coxph(
  Surv(PFI_time, PFI_event) ~ age_z + stage + pam50 + xist_z + xci_group,
  data = pfi_xist_data, x = TRUE
)
km_p <- pchisq(km_diff$chisq, df = length(km_diff$n) - 1, lower.tail = FALSE)
write_csv(
  tibble(endpoint = "PFI", comparison = "High vs low V4 (median split)",
         n = nrow(pfi_core), events = sum(pfi_core$PFI_event),
         chisq = unname(km_diff$chisq), p.value = km_p),
  file.path(output_dir, "24_v4_pfi_km_logrank.csv")
)
write_csv(
  bind_rows(
    tidy_cox(km_cox, "Univariable V4 median split"),
    tidy_cox(km_cox_clinical, "Age + stage + PAM50 + V4 median split"),
    tidy_cox(km_cox_full, "Age + stage + PAM50 + XIST + V4 median split")
  ) %>% filter(str_detect(term, "xci_group")),
  file.path(output_dir, "24a_v4_pfi_median_split_cox.csv")
)

ph_test <- cox.zph(pfi_full)
ph_df <- as.data.frame(ph_test$table) %>% rownames_to_column("term")
write_csv(ph_df, file.path(output_dir, "25_v4_pfi_full_model_ph_test.csv"))
ph_clinical_test <- cox.zph(pfi_clinical_v4)
ph_clinical_df <- as.data.frame(ph_clinical_test$table) %>% rownames_to_column("term")
write_csv(
  ph_clinical_df,
  file.path(output_dir, "25a_v4_pfi_clinical_model_ph_test.csv")
)

# -----------------------------------------------------------------------------
# Compact machine-readable conclusions and run metadata
# -----------------------------------------------------------------------------

key_tp53 <- tp53_model_results %>%
  filter(
    (model == "Unadjusted complete-case cohort" & term == "xci_z") |
      (model == "Adjusted continuous V4" & term == "xci_z") |
      (model == "Adjusted ordinal V4 quartile trend" & term == "xci_quartile_numeric")
  )

key_interactions <- rq3_interaction_results %>%
  filter(adjustment_set == "clinical_adjusted") %>%
  dplyr::select(outcome, xci_parameterisation, n, effect, effect_CI_low,
                effect_CI_high, effect_scale, p.value)

key_pfi <- pfi_cox_results %>%
  filter(term == "xci_z") %>%
  dplyr::select(model, n, events, HR, CI_low, CI_high, p.value)

write_csv(key_interactions, file.path(output_dir, "26_key_adjusted_interactions.csv"))
write_csv(key_pfi, file.path(output_dir, "27_key_pfi_results.csv"))

capture.output(sessionInfo(), file = file.path(output_dir, "28_sessionInfo.txt"))

summary_lines <- c(
  "PFI-RESTRICTED FINAL ANALYSIS SUMMARY",
  paste0("Generated: ", Sys.time()),
  "OS_time and OS_event were not used.",
  "",
  "TP53 logistic cohort summary:",
  capture.output(print(tp53_cohort_summary)),
  "",
  "Key TP53 logistic results:",
  capture.output(print(key_tp53)),
  "",
  "RQ3 PFI-restricted group counts:",
  capture.output(print(rq3 %>% count(combined_group))),
  "",
  "Key clinically adjusted interaction terms:",
  capture.output(print(key_interactions)),
  "",
  "PFI survival cohorts:",
  paste0("Core n = ", nrow(pfi_core), "; PFI events = ", sum(pfi_core$PFI_event)),
  paste0("Clinical n = ", nrow(pfi_clinical), "; PFI events = ",
         sum(pfi_clinical$PFI_event)),
  paste0("Clinical + XIST n = ", nrow(pfi_xist_data), "; PFI events = ",
         sum(pfi_xist_data$PFI_event)),
  "",
  "Key PFI Cox results:",
  capture.output(print(key_pfi)),
  "",
  paste0("PFI median-split log-rank p = ", signif(km_p, 6)),
  paste0("RQ3/V4 score correlation audit = ", signif(score_correlation, 8))
)
writeLines(summary_lines, file.path(output_dir, "29_run_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")

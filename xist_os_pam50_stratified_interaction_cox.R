############################################################
## xist_os_pam50_stratified_interaction_cox.R
##
## TCGA-BRCA
## XIST and OS:
## 1. PAM50-stratified Cox models
## 2. XIST × PAM50 interaction model
##
## Uses existing local files only:
## - xist_patient_level_expression.csv
## - tcga_brca_clinical_raw.csv
## - xist_stage_pam50_metadata.csv
############################################################

############################
## 0. Load packages
############################

cran_pkgs <- c("dplyr", "stringr", "readr", "survival", "broom", "purrr")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(stringr)
library(readr)
library(survival)
library(broom)
library(purrr)

############################
## 1. Set directory
############################

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, mustWork = TRUE)
outdir <- file.path(project_dir, "data")
if (!dir.exists(outdir)) stop("Cannot find project data/ directory")
setwd(outdir)

############################
## 2. Helper functions
############################

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

pick_col_regex <- function(df, patterns) {
  nm <- colnames(df)
  for (p in patterns) {
    hit <- grep(p, nm, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  return(NA_character_)
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

clean_stage4 <- function(x) {
  x <- as.character(x)
  case_when(
    str_detect(x, regex("^stage iv", ignore_case = TRUE)) ~ "Stage IV",
    str_detect(x, regex("^stage iii", ignore_case = TRUE)) ~ "Stage III",
    str_detect(x, regex("^stage ii", ignore_case = TRUE)) ~ "Stage II",
    str_detect(x, regex("^stage i($|[^v])", ignore_case = TRUE)) ~ "Stage I",
    TRUE ~ NA_character_
  )
}

############################
## 3. Read XIST expression
############################

xist_df <- read_csv(
  file.path(outdir, "xist_patient_level_expression.csv"),
  show_col_types = FALSE
)

############################
## 4. Read clinical survival data
############################

clinical_df <- read_csv(
  file.path(outdir, "tcga_brca_clinical_raw.csv"),
  show_col_types = FALSE
)

patient_col <- pick_col(clinical_df, c("submitter_id", "bcr_patient_barcode", "case_submitter_id"))
if (is.na(patient_col)) {
  stop("Could not find patient barcode column in clinical data.")
}

age_col <- pick_col(clinical_df, c("age_at_diagnosis", "age_at_index"))
if (is.na(age_col)) {
  age_col <- pick_col_regex(clinical_df, c("^age_"))
}

vital_col <- pick_col(clinical_df, c("vital_status"))
if (is.na(vital_col)) {
  vital_col <- pick_col_regex(clinical_df, c("vital_status"))
}

death_col <- pick_col(clinical_df, c("days_to_death"))
if (is.na(death_col)) {
  death_col <- pick_col_regex(clinical_df, c("days_to_death"))
}

lfu_col <- pick_col(clinical_df, c("days_to_last_follow_up", "days_to_last_known_alive"))
if (is.na(lfu_col)) {
  lfu_col <- pick_col_regex(clinical_df, c("days_to_last_follow"))
}
if (is.na(lfu_col)) {
  lfu_col <- pick_col_regex(clinical_df, c("days_to_last_known_alive"))
}

stage_col <- pick_col(clinical_df, c(
  "ajcc_pathologic_stage",
  "pathologic_stage",
  "tumor_stage",
  "ajcc_clinical_stage"
))
if (is.na(stage_col)) {
  stage_col <- pick_col_regex(clinical_df, c("ajcc.*stage", "pathologic_stage", "tumor_stage"))
}

clinical_use <- clinical_df %>%
  mutate(
    patient_barcode = substr(.data[[patient_col]], 1, 12),
    
    age_raw = if (!is.na(age_col)) .data[[age_col]] else NA,
    age_num = safe_num(age_raw),
    
    age_years = case_when(
      is.na(age_num) ~ NA_real_,
      !is.na(age_col) && age_col == "age_at_diagnosis" ~ abs(age_num) / 365.25,
      !is.na(age_num) & age_num > 200 ~ age_num / 365.25,
      TRUE ~ age_num
    ),
    
    vital_status_raw = if (!is.na(vital_col)) as.character(.data[[vital_col]]) else NA_character_,
    vital_status_low = tolower(vital_status_raw),
    
    days_to_death_num = if (!is.na(death_col)) safe_num(.data[[death_col]]) else NA_real_,
    days_to_last_fu_num = if (!is.na(lfu_col)) safe_num(.data[[lfu_col]]) else NA_real_,
    
    OS_time = case_when(
      !is.na(days_to_death_num) ~ days_to_death_num,
      is.na(days_to_death_num) & !is.na(days_to_last_fu_num) ~ days_to_last_fu_num,
      TRUE ~ NA_real_
    ),
    
    OS_event = case_when(
      vital_status_low %in% c("dead", "deceased") ~ 1,
      vital_status_low %in% c("alive") ~ 0,
      TRUE ~ NA_real_
    ),
    
    stage_raw = if (!is.na(stage_col)) as.character(.data[[stage_col]]) else NA_character_,
    stage4 = clean_stage4(stage_raw)
  ) %>%
  select(patient_barcode, age_years, OS_time, OS_event, stage4) %>%
  distinct(patient_barcode, .keep_all = TRUE)

############################
## 5. Read PAM50 metadata
############################

pam50_meta <- read_csv(
  file.path(outdir, "xist_stage_pam50_metadata.csv"),
  show_col_types = FALSE
)

pam50_df <- pam50_meta %>%
  transmute(
    patient_barcode = patient_barcode,
    PAM50_subtype = case_when(
      pam50 %in% c("LumA", "LumB") ~ as.character(pam50),
      pam50 %in% c("HER2-enriched", "Her2", "HER2") ~ "HER2-enriched",
      pam50 %in% c("Basal-like", "Basal") ~ "Basal",
      pam50 %in% c("Normal-like", "Normal") ~ "Normal-like",
      TRUE ~ as.character(pam50)
    )
  ) %>%
  filter(!is.na(patient_barcode), !is.na(PAM50_subtype)) %>%
  distinct(patient_barcode, .keep_all = TRUE)

############################
## 6. Build full analysis dataset
############################

analysis_df <- clinical_use %>%
  inner_join(xist_df, by = "patient_barcode") %>%
  left_join(pam50_df, by = "patient_barcode") %>%
  filter(
    !is.na(OS_time),
    !is.na(OS_event),
    !is.na(age_years),
    !is.na(stage4),
    !is.na(XIST_logCPM),
    !is.na(PAM50_subtype),
    OS_time > 0,
    stage4 %in% c("Stage I", "Stage II", "Stage III", "Stage IV")
  ) %>%
  mutate(
    stage4 = factor(stage4, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    PAM50_subtype = factor(
      PAM50_subtype,
      levels = c("LumA", "LumB", "HER2-enriched", "Basal", "Normal-like")
    )
  )

write_csv(
  analysis_df,
  file.path(outdir, "xist_os_pam50_stratified_dataset.csv")
)

############################
## 7. PAM50-stratified Cox models
############################
## Model within each subtype:
## Surv(OS_time, OS_event) ~ age_years + stage4 + XIST_logCPM

subtype_levels <- levels(droplevels(analysis_df$PAM50_subtype))

fit_one_subtype <- function(st_name) {
  df_sub <- analysis_df %>%
    filter(PAM50_subtype == st_name) %>%
    droplevels()
  
  n_patients <- nrow(df_sub)
  n_events <- sum(df_sub$OS_event, na.rm = TRUE)
  
  # Need enough events to fit a Cox model sensibly
  if (n_patients < 20 || n_events < 10) {
    return(data.frame(
      PAM50_subtype = st_name,
      term = "XIST_logCPM",
      hazard_ratio = NA_real_,
      conf_low_95 = NA_real_,
      conf_high_95 = NA_real_,
      p_value = NA_real_,
      n_patients = n_patients,
      n_events = n_events,
      note = "Too few patients/events for stable Cox model"
    ))
  }
  
  fit <- tryCatch(
    coxph(Surv(OS_time, OS_event) ~ age_years + stage4 + XIST_logCPM, data = df_sub),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.frame(
      PAM50_subtype = st_name,
      term = "XIST_logCPM",
      hazard_ratio = NA_real_,
      conf_low_95 = NA_real_,
      conf_high_95 = NA_real_,
      p_value = NA_real_,
      n_patients = n_patients,
      n_events = n_events,
      note = "Model failed to converge"
    ))
  }
  
  tidy_fit <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  
  xist_row <- tidy_fit %>%
    filter(term == "XIST_logCPM")
  
  if (nrow(xist_row) == 0) {
    return(data.frame(
      PAM50_subtype = st_name,
      term = "XIST_logCPM",
      hazard_ratio = NA_real_,
      conf_low_95 = NA_real_,
      conf_high_95 = NA_real_,
      p_value = NA_real_,
      n_patients = n_patients,
      n_events = n_events,
      note = "XIST row not found"
    ))
  }
  
  data.frame(
    PAM50_subtype = st_name,
    term = "XIST_logCPM",
    hazard_ratio = xist_row$estimate,
    conf_low_95 = xist_row$conf.low,
    conf_high_95 = xist_row$conf.high,
    p_value = xist_row$p.value,
    n_patients = n_patients,
    n_events = n_events,
    note = ""
  )
}

stratified_results <- map_dfr(subtype_levels, fit_one_subtype)

write_csv(
  stratified_results,
  file.path(outdir, "xist_os_pam50_stratified_cox_results.csv")
)

############################
## 8. Interaction model
############################
## Main-effects model:
## Surv(OS_time, OS_event) ~ age_years + stage4 + PAM50_subtype + XIST_logCPM
##
## Interaction model:
## Surv(OS_time, OS_event) ~ age_years + stage4 + PAM50_subtype * XIST_logCPM

cox_main <- coxph(
  Surv(OS_time, OS_event) ~ age_years + stage4 + PAM50_subtype + XIST_logCPM,
  data = analysis_df
)

cox_interaction <- coxph(
  Surv(OS_time, OS_event) ~ age_years + stage4 + PAM50_subtype * XIST_logCPM,
  data = analysis_df
)

saveRDS(cox_main, file.path(outdir, "cox_pam50_main_effects.rds"))
saveRDS(cox_interaction, file.path(outdir, "cox_pam50_interaction.rds"))

interaction_results <- broom::tidy(
  cox_interaction,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  mutate(
    n_patients = nrow(analysis_df),
    n_events = sum(analysis_df$OS_event)
  )

write_csv(
  interaction_results,
  file.path(outdir, "xist_os_pam50_interaction_model_results.csv")
)

############################
## 9. Likelihood ratio test for interaction
############################

lrt_interaction <- anova(cox_main, cox_interaction, test = "LRT")

interaction_lrt_df <- data.frame(
  comparison = "Add XIST × PAM50 interaction terms",
  chisq = if ("Chisq" %in% colnames(lrt_interaction)) lrt_interaction[2, "Chisq"] else NA_real_,
  df = if ("Df" %in% colnames(lrt_interaction)) lrt_interaction[2, "Df"] else NA_real_,
  p_value = if ("Pr(>|Chi|)" %in% colnames(lrt_interaction)) {
    lrt_interaction[2, "Pr(>|Chi|)"]
  } else if ("P(>|Chi|)" %in% colnames(lrt_interaction)) {
    lrt_interaction[2, "P(>|Chi|)"]
  } else {
    NA_real_
  },
  n_patients = nrow(analysis_df),
  n_events = sum(analysis_df$OS_event)
)

write_csv(
  interaction_lrt_df,
  file.path(outdir, "xist_os_pam50_interaction_lrt.csv")
)

############################
## 10. Console summary
############################

cat("\n============================\n")
cat("PAM50-stratified and interaction Cox analysis complete\n")
cat("============================\n\n")

cat("Overall dataset:\n")
cat("N patients:", nrow(analysis_df), "\n")
cat("N events:", sum(analysis_df$OS_event), "\n\n")

cat("Subtype counts:\n")
print(table(analysis_df$PAM50_subtype))

cat("\nStratified Cox results for XIST:\n")
print(stratified_results)

cat("\nInteraction LRT result:\n")
print(interaction_lrt_df)

cat("\nFiles saved to:\n", outdir, "\n")
cat("- xist_os_pam50_stratified_dataset.csv\n")
cat("- xist_os_pam50_stratified_cox_results.csv\n")
cat("- xist_os_pam50_interaction_model_results.csv\n")
cat("- xist_os_pam50_interaction_lrt.csv\n")
cat("- cox_pam50_main_effects.rds\n")
cat("- cox_pam50_interaction.rds\n")

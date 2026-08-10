############################################################
## xist_os_multivariable_cox.R
##
## TCGA-BRCA
## Association between XIST expression and overall survival
##
## Analyses:
## 1. Univariable Cox: OS ~ XIST
## 2. Multivariable Cox: OS ~ age + stage + XIST
## 3. Multivariable Cox: OS ~ age + stage + PAM50 + XIST
##
## Uses existing local files only:
## - xist_patient_level_expression.csv
## - tcga_brca_clinical_raw.csv
## - xist_stage_pam50_metadata.csv
############################################################

############################
## 0. Load packages
############################

cran_pkgs <- c("dplyr", "stringr", "readr", "survival", "broom")

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
## 5. Read PAM50 subtype metadata
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
## 6. Build datasets
############################

analysis_df_basic <- clinical_use %>%
  inner_join(xist_df, by = "patient_barcode") %>%
  filter(
    !is.na(OS_time),
    !is.na(OS_event),
    !is.na(age_years),
    !is.na(stage4),
    !is.na(XIST_logCPM),
    OS_time > 0,
    stage4 %in% c("Stage I", "Stage II", "Stage III", "Stage IV")
  ) %>%
  mutate(
    stage4 = factor(stage4, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
  )

analysis_df_pam50 <- analysis_df_basic %>%
  left_join(pam50_df, by = "patient_barcode") %>%
  filter(!is.na(PAM50_subtype)) %>%
  mutate(
    PAM50_subtype = factor(
      PAM50_subtype,
      levels = c("LumA", "LumB", "HER2-enriched", "Basal", "Normal-like")
    )
  )

write_csv(
  analysis_df_basic,
  file.path(outdir, "xist_os_cox_dataset_basic.csv")
)

write_csv(
  analysis_df_pam50,
  file.path(outdir, "xist_os_cox_dataset_with_pam50.csv")
)

############################
## 7. Fit Cox models
############################

# Model 1: univariable
cox_uni <- coxph(
  Surv(OS_time, OS_event) ~ XIST_logCPM,
  data = analysis_df_basic
)

# Model 2: adjusted for age + stage
cox_adj_basic <- coxph(
  Surv(OS_time, OS_event) ~ age_years + stage4 + XIST_logCPM,
  data = analysis_df_basic
)

# Model 3: adjusted for age + stage + PAM50
cox_adj_pam50 <- coxph(
  Surv(OS_time, OS_event) ~ age_years + stage4 + PAM50_subtype + XIST_logCPM,
  data = analysis_df_pam50
)

############################
## 8. Save model objects
############################

saveRDS(cox_uni, file.path(outdir, "cox_univariable_xist.rds"))
saveRDS(cox_adj_basic, file.path(outdir, "cox_age_stage_xist.rds"))
saveRDS(cox_adj_pam50, file.path(outdir, "cox_age_stage_pam50_xist.rds"))

############################
## 9. Extract tidy results
############################

tidy_uni <- broom::tidy(cox_uni, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    model = "Univariable: XIST only",
    n_patients = nrow(analysis_df_basic),
    n_events = sum(analysis_df_basic$OS_event)
  )

tidy_adj_basic <- broom::tidy(cox_adj_basic, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    model = "Multivariable: age + stage + XIST",
    n_patients = nrow(analysis_df_basic),
    n_events = sum(analysis_df_basic$OS_event)
  )

tidy_adj_pam50 <- broom::tidy(cox_adj_pam50, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    model = "Multivariable: age + stage + PAM50 + XIST",
    n_patients = nrow(analysis_df_pam50),
    n_events = sum(analysis_df_pam50$OS_event)
  )

all_results <- bind_rows(
  tidy_uni,
  tidy_adj_basic,
  tidy_adj_pam50
)

write_csv(
  all_results,
  file.path(outdir, "xist_os_cox_all_results.csv")
)

############################
## 10. Extract XIST rows only
############################

xist_results <- all_results %>%
  filter(term == "XIST_logCPM") %>%
  select(
    model, term, estimate, conf.low, conf.high, p.value, n_patients, n_events
  ) %>%
  rename(
    hazard_ratio = estimate,
    conf_low_95 = conf.low,
    conf_high_95 = conf.high,
    p_value = p.value
  )

write_csv(
  xist_results,
  file.path(outdir, "xist_os_cox_xist_only_results.csv")
)

############################
## 11. Likelihood ratio tests for added XIST
############################

# For age + stage model
cox_adj_basic_no_xist <- coxph(
  Surv(OS_time, OS_event) ~ age_years + stage4,
  data = analysis_df_basic
)

lrt_basic <- anova(cox_adj_basic_no_xist, cox_adj_basic, test = "LRT")

# For age + stage + PAM50 model
cox_adj_pam50_no_xist <- coxph(
  Surv(OS_time, OS_event) ~ age_years + stage4 + PAM50_subtype,
  data = analysis_df_pam50
)

lrt_pam50 <- anova(cox_adj_pam50_no_xist, cox_adj_pam50, test = "LRT")

extract_lrt <- function(tbl, model_name) {
  out <- data.frame(
    model = model_name,
    chisq = NA_real_,
    df = NA_real_,
    p_value = NA_real_
  )
  
  if ("Chisq" %in% colnames(tbl)) {
    out$chisq <- tbl[2, "Chisq"]
  }
  if ("Df" %in% colnames(tbl)) {
    out$df <- tbl[2, "Df"]
  }
  if ("Pr(>|Chi|)" %in% colnames(tbl)) {
    out$p_value <- tbl[2, "Pr(>|Chi|)"]
  } else if ("P(>|Chi|)" %in% colnames(tbl)) {
    out$p_value <- tbl[2, "P(>|Chi|)"]
  }
  
  out
}

lrt_results <- bind_rows(
  extract_lrt(lrt_basic, "Add XIST to age + stage"),
  extract_lrt(lrt_pam50, "Add XIST to age + stage + PAM50")
)

write_csv(
  lrt_results,
  file.path(outdir, "xist_os_cox_lrt_results.csv")
)

############################
## 12. Console summary
############################

cat("\n============================\n")
cat("XIST OS multivariable Cox analysis complete\n")
cat("============================\n\n")

cat("Dataset 1: age + stage + XIST\n")
cat("N patients:", nrow(analysis_df_basic), "\n")
cat("N events:", sum(analysis_df_basic$OS_event), "\n\n")

cat("Dataset 2: age + stage + PAM50 + XIST\n")
cat("N patients:", nrow(analysis_df_pam50), "\n")
cat("N events:", sum(analysis_df_pam50$OS_event), "\n\n")

cat("XIST rows from all models:\n")
print(xist_results)

cat("\nLikelihood ratio tests:\n")
print(lrt_results)

cat("\nFiles saved to:\n", outdir, "\n")
cat("- xist_os_cox_dataset_basic.csv\n")
cat("- xist_os_cox_dataset_with_pam50.csv\n")
cat("- cox_univariable_xist.rds\n")
cat("- cox_age_stage_xist.rds\n")
cat("- cox_age_stage_pam50_xist.rds\n")
cat("- xist_os_cox_all_results.csv\n")
cat("- xist_os_cox_xist_only_results.csv\n")
cat("- xist_os_cox_lrt_results.csv\n")

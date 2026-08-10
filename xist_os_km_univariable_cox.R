############################################################
## xist_os_km_univariable_cox.R
##
## TCGA-BRCA
## Overall association between XIST expression and OS
##
## Analyses:
## 1. Kaplan-Meier curve: XIST high vs low
## 2. Univariable Cox model: continuous XIST
##
## Input files (already existing locally):
## - xist_patient_level_expression.csv
## - tcga_brca_clinical_raw.csv
##
## Output files:
## - xist_os_analysis_dataset.csv
## - XIST_KM_OS_high_vs_low.pdf
## - xist_km_logrank_result.csv
## - xist_univariable_cox_result.csv
############################################################

############################
## 0. Load packages
############################

cran_pkgs <- c("dplyr", "stringr", "readr", "survival", "survminer")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(stringr)
library(readr)
library(survival)
library(survminer)

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

############################
## 3. Read existing XIST data
############################

xist_df <- read_csv(
  file.path(outdir, "xist_patient_level_expression.csv"),
  show_col_types = FALSE
)

############################
## 4. Read existing clinical data
############################

clinical_df <- read_csv(
  file.path(outdir, "tcga_brca_clinical_raw.csv"),
  show_col_types = FALSE
)

############################
## 5. Detect clinical columns
############################

patient_col <- pick_col(clinical_df, c("submitter_id", "bcr_patient_barcode", "case_submitter_id"))
if (is.na(patient_col)) {
  stop("Could not find patient barcode column in clinical data.")
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

############################
## 6. Build clinical survival table
############################

clinical_use <- clinical_df %>%
  mutate(
    patient_barcode = substr(.data[[patient_col]], 1, 12),
    
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
    )
  ) %>%
  select(patient_barcode, OS_time, OS_event) %>%
  distinct(patient_barcode, .keep_all = TRUE)

############################
## 7. Merge XIST with survival data
############################

analysis_df <- clinical_use %>%
  inner_join(xist_df, by = "patient_barcode") %>%
  filter(
    !is.na(OS_time),
    !is.na(OS_event),
    !is.na(XIST_logCPM),
    OS_time > 0
  )

############################
## 8. Define XIST high vs low groups
############################
## Here we use median split

xist_median <- median(analysis_df$XIST_logCPM, na.rm = TRUE)

analysis_df <- analysis_df %>%
  mutate(
    XIST_group = ifelse(XIST_logCPM >= xist_median, "High XIST", "Low XIST"),
    XIST_group = factor(XIST_group, levels = c("Low XIST", "High XIST"))
  )

write_csv(
  analysis_df,
  file.path(outdir, "xist_os_analysis_dataset.csv")
)

############################
## 9. Kaplan-Meier analysis
############################

km_fit <- survfit(Surv(OS_time, OS_event) ~ XIST_group, data = analysis_df)

logrank_test <- survdiff(Surv(OS_time, OS_event) ~ XIST_group, data = analysis_df)

logrank_chisq <- unname(logrank_test$chisq)
logrank_p <- 1 - pchisq(logrank_chisq, df = length(logrank_test$n) - 1)

km_result_df <- data.frame(
  n_patients = nrow(analysis_df),
  n_events = sum(analysis_df$OS_event),
  median_XIST_cutoff = xist_median,
  logrank_chisq = logrank_chisq,
  logrank_df = length(logrank_test$n) - 1,
  logrank_p_value = logrank_p
)

write_csv(
  km_result_df,
  file.path(outdir, "xist_km_logrank_result.csv")
)

############################
## 10. Plot Kaplan-Meier curve
############################

km_plot <- ggsurvplot(
  km_fit,
  data = analysis_df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("#D55E00", "#0072B2"),
  legend.title = "XIST group",
  legend.labs = c("Low XIST", "High XIST"),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  title = "Kaplan-Meier curve of overall survival by XIST expression group",
  break.time.by = 365.25
)

ggsave(
  filename = file.path(outdir, "XIST_KM_OS_high_vs_low.pdf"),
  plot = km_plot$plot,
  width = 8,
  height = 6
)

############################
## 11. Univariable Cox model
############################

cox_uni <- coxph(Surv(OS_time, OS_event) ~ XIST_logCPM, data = analysis_df)

cox_summary <- summary(cox_uni)

cox_result_df <- data.frame(
  term = "XIST_logCPM",
  hazard_ratio = unname(cox_summary$coefficients[, "exp(coef)"]),
  conf_low_95 = unname(cox_summary$conf.int[, "lower .95"]),
  conf_high_95 = unname(cox_summary$conf.int[, "upper .95"]),
  p_value = unname(cox_summary$coefficients[, "Pr(>|z|)"]),
  n_patients = nrow(analysis_df),
  n_events = sum(analysis_df$OS_event)
)

write_csv(
  cox_result_df,
  file.path(outdir, "xist_univariable_cox_result.csv")
)

############################
## 12. Console summary
############################

cat("\n============================\n")
cat("XIST-OS overall association analysis complete\n")
cat("============================\n")
cat("N patients:", nrow(analysis_df), "\n")
cat("N OS events:", sum(analysis_df$OS_event), "\n")
cat("Median XIST cutoff:", round(xist_median, 4), "\n\n")

cat("Kaplan-Meier log-rank test:\n")
cat("Chi-square =", round(logrank_chisq, 4), "\n")
cat("p-value    =", signif(logrank_p, 4), "\n\n")

cat("Univariable Cox model:\n")
cat("HR for XIST_logCPM =", round(cox_result_df$hazard_ratio, 4), "\n")
cat("95% CI =", paste0("(", round(cox_result_df$conf_low_95, 4), ", ",
                       round(cox_result_df$conf_high_95, 4), ")"), "\n")
cat("p-value =", signif(cox_result_df$p_value, 4), "\n\n")

cat("Files saved to:\n", outdir, "\n")
cat("- xist_os_analysis_dataset.csv\n")
cat("- XIST_KM_OS_high_vs_low.pdf\n")
cat("- xist_km_logrank_result.csv\n")
cat("- xist_univariable_cox_result.csv\n")

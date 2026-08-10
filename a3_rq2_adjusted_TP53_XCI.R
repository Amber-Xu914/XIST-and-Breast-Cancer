# ============================================================
# Aim 3 RQ2
# Is the XCI erosion–TP53 mutation association independent of
# PAM50 subtype and clinical covariates?
#
# Main model:
# TP53_mut ~ xci_score_z + age_years + stage + PAM50
# ============================================================

library(tidyverse)
library(broom)

# -------------------------
# 0. Output folder
# -------------------------

out_dir <- "outputs/aim3_rq2_tp53_adjusted_no_xist"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# -------------------------
# 1. Helper functions
# -------------------------

pick_col <- function(df, exact = character(), regex = NULL, required = TRUE, label = "column") {
  nms <- names(df)
  lower_nms <- tolower(nms)
  
  if (length(exact) > 0) {
    idx <- match(tolower(exact), lower_nms)
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) {
      return(nms[idx[1]])
    }
  }
  
  if (!is.null(regex)) {
    idx <- which(stringr::str_detect(lower_nms, regex))
    if (length(idx) > 0) {
      return(nms[idx[1]])
    }
  }
  
  if (required) {
    stop(paste0("Could not find required ", label, ". Available columns are:\n",
                paste(nms, collapse = ", ")))
  } else {
    return(NA_character_)
  }
}


detect_age_col <- function(df) {
  nms <- names(df)
  lower_nms <- tolower(nms)
  
  exact_age_cols <- c(
    "age_years",
    "age_at_diagnosis",
    "age_at_index",
    "diagnosis_age",
    "days_to_birth",
    "diagnoses.age_at_diagnosis"
  )
  
  idx <- match(tolower(exact_age_cols), lower_nms)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) > 0) {
    return(nms[idx[1]])
  }
  
  # Avoid accidentally matching "stage"
  idx2 <- which(
    (stringr::str_detect(lower_nms, "age") |
       stringr::str_detect(lower_nms, "birth")) &
      !stringr::str_detect(lower_nms, "stage")
  )
  
  if (length(idx2) > 0) {
    return(nms[idx2[1]])
  }
  
  return(NA_character_)
}


make_age_years <- function(x, colname) {
  v <- suppressWarnings(as.numeric(x))
  med_abs <- median(abs(v), na.rm = TRUE)
  
  # TCGA age_at_diagnosis or days_to_birth is often in days.
  # If values are much larger than normal human age, convert days to years.
  if (is.finite(med_abs) && med_abs > 150) {
    v <- abs(v) / 365.25
  }
  
  return(v)
}


read_any_table <- function(file) {
  if (grepl("\\.rds$", file, ignore.case = TRUE)) {
    obj <- readRDS(file)
    if (is.data.frame(obj)) return(obj)
    return(NULL)
  }
  
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    return(read.csv(file, stringsAsFactors = FALSE, check.names = FALSE))
  }
  
  if (grepl("\\.tsv$|\\.txt$|\\.tsv\\.gz$|\\.txt\\.gz$", file, ignore.case = TRUE)) {
    return(read.delim(file, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#"))
  }
  
  if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      message("Skipping xlsx file because readxl is not installed: ", file)
      return(NULL)
    }
    return(as.data.frame(readxl::read_excel(file)))
  }
  
  return(NULL)
}


load_age_metadata <- function() {
  metadata_files <- list.files(
    path = c("data", "outputs"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE,
    pattern = "clinical|metadata|cdr|survival|analysis_dataset|proxy|plotting"
  )
  
  metadata_files <- metadata_files[
    grepl("\\.csv$|\\.tsv$|\\.txt$|\\.rds$|\\.xlsx$", metadata_files, ignore.case = TRUE)
  ]
  
  metadata_files <- metadata_files[
    !grepl("maf|mutation|methylation|counts|segments|cnv|star", metadata_files, ignore.case = TRUE)
  ]
  
  cat("\nSearching for age metadata in candidate files:\n")
  print(metadata_files)
  
  for (f in metadata_files) {
    cat("\nChecking file:\n")
    print(f)
    
    tmp <- tryCatch(read_any_table(f), error = function(e) NULL)
    
    if (is.null(tmp) || !is.data.frame(tmp)) {
      next
    }
    
    patient_col <- pick_col(
      tmp,
      exact = c("patient_barcode", "bcr_patient_barcode", "submitter_id", "case_submitter_id"),
      regex = "patient|barcode|submitter",
      required = FALSE,
      label = "patient barcode"
    )
    
    age_col <- detect_age_col(tmp)
    
    if (!is.na(patient_col) && !is.na(age_col)) {
      cat("\nFound age metadata:\n")
      cat("File:", f, "\n")
      cat("Patient column:", patient_col, "\n")
      cat("Age column:", age_col, "\n")
      
      age_df <- tmp %>%
        transmute(
          patient_barcode = substr(as.character(.data[[patient_col]]), 1, 12),
          age_years = make_age_years(.data[[age_col]], age_col)
        ) %>%
        filter(!is.na(patient_barcode), !is.na(age_years)) %>%
        group_by(patient_barcode) %>%
        summarise(
          age_years = first(age_years),
          .groups = "drop"
        )
      
      return(age_df)
    }
  }
  
  stop("No usable age metadata was found. Please provide or create a file with patient_barcode and age_years.")
}


tidy_or <- function(model, model_name) {
  broom::tidy(model) %>%
    mutate(
      model = model_name,
      odds_ratio = exp(estimate),
      conf.low = exp(estimate - 1.96 * std.error),
      conf.high = exp(estimate + 1.96 * std.error)
    ) %>%
    select(
      model,
      term,
      odds_ratio,
      std.error,
      statistic,
      p.value,
      conf.low,
      conf.high
    )
}


# -------------------------
# 2. Load XCI erosion dataset
# -------------------------

xci_file <- "outputs/RQ1_V4_XCI_linked/rq1_v4_plotting_analysis_dataset.csv"

if (!file.exists(xci_file)) {
  stop("XCI file not found: ", xci_file)
}

xci_raw <- read.csv(xci_file, stringsAsFactors = FALSE, check.names = FALSE)

cat("\nXCI dataset loaded:\n")
cat("File:", xci_file, "\n")
cat("Dimension:", dim(xci_raw), "\n")
cat("\nColumns:\n")
print(colnames(xci_raw))


patient_col <- pick_col(
  xci_raw,
  exact = c("patient_barcode", "bcr_patient_barcode", "submitter_id"),
  regex = "patient|barcode|submitter",
  label = "patient barcode"
)

score_col <- pick_col(
  xci_raw,
  exact = c(
    "xci_erosion_score_v4_xci_linked_gene_only",
    "XCI_erosion_score_mean_chrX",
    "xci_erosion_score",
    "xci_score"
  ),
  regex = "xci.*erosion.*score",
  label = "XCI erosion score"
)

pam50_col <- pick_col(
  xci_raw,
  exact = c("pam50", "PAM50", "PAM50_subtype", "subtype"),
  regex = "pam50",
  label = "PAM50 subtype"
)

stage_col <- pick_col(
  xci_raw,
  exact = c("stage_clean", "stage4", "pathologic_stage", "ajcc_pathologic_stage"),
  regex = "stage",
  label = "tumour stage"
)

age_col_xci <- detect_age_col(xci_raw)


cat("\nDetected columns:\n")
cat("Patient:", patient_col, "\n")
cat("XCI score:", score_col, "\n")
cat("PAM50:", pam50_col, "\n")
cat("Stage:", stage_col, "\n")
cat("Age in XCI file:", age_col_xci, "\n")


xci_df <- xci_raw %>%
  transmute(
    patient_barcode = substr(as.character(.data[[patient_col]]), 1, 12),
    xci_score = as.numeric(.data[[score_col]]),
    pam50 = na_if(as.character(.data[[pam50_col]]), ""),
    stage_covariate = na_if(as.character(.data[[stage_col]]), "")
  ) %>%
  distinct(patient_barcode, .keep_all = TRUE)

if (!is.na(age_col_xci)) {
  xci_df$age_years <- make_age_years(xci_raw[[age_col_xci]], age_col_xci)
}


# -------------------------
# 3. Add age if not in XCI dataset
# -------------------------

if (!"age_years" %in% colnames(xci_df) || all(is.na(xci_df$age_years))) {
  cat("\nAge is not available in the XCI file. Searching external metadata files...\n")
  age_df <- load_age_metadata()
  
  xci_df <- xci_df %>%
    left_join(age_df, by = "patient_barcode")
}


# -------------------------
# 4. Load MAF and define TP53 mutation status
# -------------------------

maf_rds <- "data/tcga_brca_mutect2_maf.rds"
maf_tsv <- "data/tcga_brca_mutect2.maf.tsv.gz"

if (file.exists(maf_rds)) {
  maf_df <- readRDS(maf_rds)
  cat("\nLoaded cached MAF RDS:\n")
  cat(maf_rds, "\n")
} else {
  if (!file.exists(maf_tsv)) {
    stop("MAF file not found: ", maf_tsv)
  }
  
  maf_df <- read.delim(
    maf_tsv,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "#"
  )
  
  saveRDS(maf_df, maf_rds)
  
  cat("\nLoaded MAF TSV and saved RDS cache:\n")
  cat(maf_tsv, "\n")
}

cat("\nMAF dimension:\n")
print(dim(maf_df))

non_silent_classes <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Splice_Site",
  "Translation_Start_Site",
  "Nonstop_Mutation"
)

tp53_df <- maf_df %>%
  mutate(
    patient_barcode = substr(as.character(Tumor_Sample_Barcode), 1, 12)
  ) %>%
  filter(
    Hugo_Symbol == "TP53",
    Variant_Classification %in% non_silent_classes
  ) %>%
  distinct(patient_barcode) %>%
  mutate(TP53_mut = 1L)

cat("\nNumber of TP53-mutated patients in MAF:\n")
print(nrow(tp53_df))


# -------------------------
# 5. Merge datasets
# -------------------------

analysis_df <- xci_df %>%
  left_join(tp53_df, by = "patient_barcode") %>%
  mutate(
    TP53_mut = replace_na(TP53_mut, 0L),
    xci_score_z = as.numeric(scale(xci_score)),
    xci_quartile_numeric = ntile(xci_score, 4),
    xci_quartile = factor(
      xci_quartile_numeric,
      levels = c(1, 2, 3, 4),
      labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest")
    )
  )

cat("\nMerged analysis dataset:\n")
print(dim(analysis_df))

cat("\nMissing values before complete-case filtering:\n")
missing_summary <- analysis_df %>%
  summarise(
    n_total = n(),
    missing_xci_score_z = sum(is.na(xci_score_z)),
    missing_age_years = sum(is.na(age_years)),
    missing_stage = sum(is.na(stage_covariate)),
    missing_pam50 = sum(is.na(pam50)),
    missing_tp53 = sum(is.na(TP53_mut))
  )

print(missing_summary)


# -------------------------
# 6. Complete-case dataset for RQ2
# -------------------------

rq2_df <- analysis_df %>%
  filter(
    !is.na(TP53_mut),
    !is.na(xci_score_z),
    !is.na(xci_quartile_numeric),
    !is.na(age_years),
    !is.na(stage_covariate),
    !is.na(pam50)
  ) %>%
  mutate(
    TP53_mut = as.integer(TP53_mut),
    pam50 = factor(pam50),
    stage_covariate = factor(stage_covariate)
  )

# Optional reference levels
if ("LumA" %in% levels(rq2_df$pam50)) {
  rq2_df$pam50 <- relevel(rq2_df$pam50, ref = "LumA")
}

if ("Stage I" %in% levels(rq2_df$stage_covariate)) {
  rq2_df$stage_covariate <- relevel(rq2_df$stage_covariate, ref = "Stage I")
}

cat("\nRQ2 complete-case sample size:\n")
print(nrow(rq2_df))

cat("\nTP53 mutation count in RQ2 complete-case dataset:\n")
print(table(rq2_df$TP53_mut))

cat("\nPAM50 distribution:\n")
print(table(rq2_df$pam50))

cat("\nStage distribution:\n")
print(table(rq2_df$stage_covariate))

write.csv(
  rq2_df,
  file.path(out_dir, "aim3_rq2_complete_case_dataset.csv"),
  row.names = FALSE
)


# -------------------------
# 7. Main RQ2 models
# -------------------------

# Unadjusted model using the same complete-case sample
model_unadj_continuous <- glm(
  TP53_mut ~ xci_score_z,
  data = rq2_df,
  family = binomial
)

# Adjusted main model
# XIST expression is intentionally excluded to avoid potential
# overadjustment because it is biologically closely related to XCI erosion.
model_adj_continuous <- glm(
  TP53_mut ~ xci_score_z + age_years + stage_covariate + pam50,
  data = rq2_df,
  family = binomial
)

# Quartile trend model, adjusted
model_adj_quartile <- glm(
  TP53_mut ~ xci_quartile_numeric + age_years + stage_covariate + pam50,
  data = rq2_df,
  family = binomial
)

# Base model without XCI, useful for likelihood-ratio test
model_base <- glm(
  TP53_mut ~ age_years + stage_covariate + pam50,
  data = rq2_df,
  family = binomial
)


# -------------------------
# 8. Extract odds ratios
# -------------------------

results_all <- bind_rows(
  tidy_or(model_unadj_continuous, "Unadjusted: XCI z-score"),
  tidy_or(model_adj_continuous, "Adjusted: XCI z-score"),
  tidy_or(model_adj_quartile, "Adjusted: XCI quartile trend")
)

results_labelled <- results_all %>%
  mutate(
    variable = case_when(
      term == "xci_score_z" ~ "XCI erosion score (z-score)",
      term == "xci_quartile_numeric" ~ "XCI erosion quartile (ordinal trend)",
      term == "age_years" ~ "Age at diagnosis (years)",
      str_detect(term, "^stage_covariate") ~ str_replace(term, "^stage_covariate", "Stage: "),
      str_detect(term, "^pam50") ~ str_replace(term, "^pam50", "PAM50: "),
      term == "(Intercept)" ~ "Intercept",
      TRUE ~ term
    )
  )

main_xci_results <- results_labelled %>%
  filter(term %in% c("xci_score_z", "xci_quartile_numeric")) %>%
  select(
    model,
    variable,
    odds_ratio,
    std.error,
    statistic,
    p.value,
    conf.low,
    conf.high
  )

cat("\nMain RQ2 XCI results:\n")
print(
  main_xci_results %>%
    mutate(
      odds_ratio = round(odds_ratio, 3),
      std.error = round(std.error, 4),
      statistic = round(statistic, 3),
      p.value = signif(p.value, 3),
      conf.low = round(conf.low, 3),
      conf.high = round(conf.high, 3)
    )
)

write.csv(
  results_labelled,
  file.path(out_dir, "aim3_rq2_all_logistic_results.csv"),
  row.names = FALSE
)

write.csv(
  main_xci_results,
  file.path(out_dir, "aim3_rq2_main_xci_adjusted_results.csv"),
  row.names = FALSE
)


# -------------------------
# 9. Likelihood-ratio test:
# Does XCI improve the adjusted model?
# -------------------------

lrt_continuous <- anova(
  model_base,
  model_adj_continuous,
  test = "Chisq"
)

lrt_quartile <- anova(
  model_base,
  model_adj_quartile,
  test = "Chisq"
)

cat("\nLikelihood-ratio test: base model vs adjusted model with XCI z-score\n")
print(lrt_continuous)

cat("\nLikelihood-ratio test: base model vs adjusted model with XCI quartile trend\n")
print(lrt_quartile)

write.csv(
  as.data.frame(lrt_continuous),
  file.path(out_dir, "aim3_rq2_lrt_xci_z_score.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(lrt_quartile),
  file.path(out_dir, "aim3_rq2_lrt_xci_quartile_trend.csv"),
  row.names = TRUE
)


# -------------------------
# 10. Forest plot for main adjusted XCI effects
# -------------------------

plot_df <- main_xci_results %>%
  filter(str_detect(model, "^Adjusted")) %>%
  mutate(
    label = case_when(
      variable == "XCI erosion score (z-score)" ~ "XCI erosion score\n(z-score)",
      variable == "XCI erosion quartile (ordinal trend)" ~ "XCI erosion quartile\n(ordinal trend)",
      TRUE ~ variable
    )
  )

p <- ggplot(plot_df, aes(x = odds_ratio, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    orientation = "y",
    width = 0.15
  ) +
  scale_x_log10() +
  labs(
    title = "Adjusted association between XCI erosion and TP53 mutation",
    x = "Odds ratio, log scale",
    y = NULL
  ) +
  theme_classic(base_size = 13)

print(p)

ggsave(
  file.path(out_dir, "aim3_rq2_adjusted_xci_tp53_forest_plot.png"),
  p,
  width = 7,
  height = 4,
  dpi = 300
)

ggsave(
  file.path(out_dir, "aim3_rq2_adjusted_xci_tp53_forest_plot.pdf"),
  p,
  width = 7,
  height = 4
)


cat("\nRQ2 analysis complete. Outputs saved to:\n")
cat(out_dir, "\n")
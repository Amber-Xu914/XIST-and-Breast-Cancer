############################################################
# RQ2: V4 XCI erosion score vs genomic instability
#      and major breast cancer driver alterations
############################################################

rm(list = ls())

############################################################
# 0. Load packages
############################################################

packages <- c(
  "dplyr", "tidyr", "readr", "data.table", "ggplot2",
  "stringr", "purrr", "broom", "forcats"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(ggplot2)
library(stringr)
library(purrr)
library(broom)
library(forcats)
library(TCGAbiolinks)
library(readr)

dir.create("data", showWarnings = FALSE)

maf_path <- "data/tcga_brca_mutect2.maf.tsv.gz"

if (!file.exists(maf_path)) {
  query_maf <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Simple Nucleotide Variation",
    access = "open",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  
  GDCdownload(
    query_maf,
    method = "api",
    files.per.chunk = 10
  )
  maf <- GDCprepare(query_maf)
  
  readr::write_tsv(
    maf,
    maf_path
  )
} else {
  message("Using cached MAF file: ", maf_path)
}

############################################################
# 1. User settings: change these paths
############################################################

# Main V4 score table from RQ1
# Must contain patient_barcode and V4 XCI erosion score
score_path <- "./data/rq1_v4_xci_linked_gene_only_integrated_score.csv"

clinical_path <- "./outputs/RQ1_V4_XCI_linked/rq1_v4_plotting_analysis_dataset.csv"

# TCGA MAF file for somatic mutations
# Set to NA if you do not have it yet
maf_path <- "data/tcga_brca_mutect2.maf.tsv.gz"

# TCGA CNV segment file
# Used to calculate fraction genome altered
# Set to NA if unavailable
seg_path <- "./data/tcga_brca_cnv_masked_segments.rds"

gistic_path <- NA


out_dir <- "results/rq2_v4"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# 2. Helper functions
############################################################

clean_colnames <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    tolower()
}

standardise_tcga_patient <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\.", "-")
  substr(x, 1, 12)
}

read_table_auto <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(NULL)
  }
  
  message("Reading: ", path)
  
  if (str_detect(path, "\\.rds$")) {
    readRDS(path)
  } else if (str_detect(path, "\\.csv$")) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    data.table::fread(path, data.table = FALSE)
  }
}

find_first_col <- function(df, candidates) {
  candidates <- clean_colnames(candidates)
  matched <- intersect(candidates, names(df))
  if (length(matched) == 0) return(NA_character_)
  matched[1]
}

safe_rename <- function(df, new_name, candidates) {
  candidates <- clean_colnames(candidates)
  
  # If the target column already exists, use it directly instead of renaming it.
  if (new_name %in% names(df)) {
    return(df)
  }
  
  old_name <- find_first_col(df, candidates)
  
  if (!is.na(old_name) && old_name != new_name) {
    names(df)[names(df) == old_name] <- new_name
  }
  
  df
}

############################################################
# 3. Load V4 score and clinical metadata
############################################################

score_df <- read_table_auto(score_path)
if (is.null(score_df)) stop("Score file not found. Please check score_path.")

names(score_df) <- clean_colnames(names(score_df))

score_df <- safe_rename(
  score_df,
  "patient_barcode",
  c("patient_barcode", "case_submitter_id", "submitter_id", "barcode", "sample")
)

score_df <- safe_rename(
  score_df,
  "xci_score",
  c(
    "xci_erosion_score_v4_xci_linked_gene_only",
    "v4_xci_linked_gene_only_erosion_score",
    "v4_integrated_score",
    "xci_erosion_score_v4",
    "xci_score"
  )
)

if (!"patient_barcode" %in% names(score_df)) {
  stop("Cannot find patient barcode column in score table.")
}

if (!"xci_score" %in% names(score_df)) {
  stop("Cannot find V4 XCI erosion score column in score table.")
}

score_df <- score_df %>%
  mutate(
    patient_barcode = standardise_tcga_patient(patient_barcode),
    xci_score = as.numeric(xci_score),
    xci_score_z = as.numeric(scale(xci_score))
  ) %>%
  filter(!is.na(patient_barcode), !is.na(xci_score))

clinical_df <- read_table_auto(clinical_path)

if (!is.null(clinical_df)) {
  names(clinical_df) <- clean_colnames(names(clinical_df))
  
  clinical_df <- safe_rename(
    clinical_df,
    "patient_barcode",
    c("patient_barcode", "case_submitter_id", "submitter_id", "barcode", "sample")
  )
  
  clinical_df <- safe_rename(
    clinical_df,
    "pam50",
    c("pam50", "pam50_subtype", "subtype", "brca_subtype", "molecular_subtype")
  )
  
  clinical_df <- safe_rename(
    clinical_df,
    "stage",
    c("stage", "tumour_stage", "tumor_stage", "pathologic_stage", "stage4")
  )
  
  clinical_df <- safe_rename(
    clinical_df,
    "grade",
    c("grade", "tumour_grade", "tumor_grade", "neoplasm_histologic_grade")
  )
  
  clinical_df <- safe_rename(
    clinical_df,
    "age",
    c("age", "age_years", "age_at_diagnosis", "diagnosis_age")
  )
  
  clinical_df <- safe_rename(
    clinical_df,
    "purity",
    c("purity", "tumour_purity", "tumor_purity", "abs_purity", "estimate_purity")
  )
  
  clinical_df <- clinical_df %>%
    mutate(patient_barcode = standardise_tcga_patient(patient_barcode)) %>%
    distinct(patient_barcode, .keep_all = TRUE)
  
  df <- score_df %>%
    left_join(clinical_df, by = "patient_barcode")
} else {
  warning("Clinical file not found. Analysis will only use covariates already present in score table.")
  df <- score_df
}

############################################################
# 4. Define four V4 erosion quartile groups
#    using 25%, 50%, and 75% cut-offs of the V4 XCI erosion score
############################################################

quartile_levels <- c(
  "Q1 low erosion",
  "Q2 low-intermediate erosion",
  "Q3 high-intermediate erosion",
  "Q4 high erosion"
)

q25 <- quantile(df$xci_score, 0.25, na.rm = TRUE)
q50 <- quantile(df$xci_score, 0.50, na.rm = TRUE)
q75 <- quantile(df$xci_score, 0.75, na.rm = TRUE)

df <- df %>%
  mutate(
    xci_erosion_group = case_when(
      is.na(xci_score) ~ NA_character_,
      xci_score <= q25 ~ "Q1 low erosion",
      xci_score <= q50 ~ "Q2 low-intermediate erosion",
      xci_score <= q75 ~ "Q3 high-intermediate erosion",
      TRUE ~ "Q4 high erosion"
    ),
    xci_erosion_group = factor(
      xci_erosion_group,
      levels = quartile_levels
    )
  )

message("V4 XCI erosion quartile cut-offs:")
print(c(
  q1_to_q2 = q25,
  q2_to_q3 = q50,
  q3_to_q4 = q75
))

message("V4 XCI erosion quartile group counts:")
print(table(df$xci_erosion_group, useNA = "ifany"))

############################################################
# 5. Load MAF and create mutation burden + driver mutations
############################################################

driver_mut_genes <- c(
  "TP53", "PIK3CA", "GATA3", "MAP3K1", "CDH1",
  "PTEN", "AKT1", "BRCA1", "BRCA2", "RB1",
  "ESR1", "ERBB2", "MAP2K4"
)

non_syn_classes <- c(
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

maf_df <- read_table_auto(maf_path)

if (!is.null(maf_df)) {
  names(maf_df) <- clean_colnames(names(maf_df))
  
  maf_df <- safe_rename(
    maf_df,
    "patient_barcode",
    c("tumor_sample_barcode", "tumour_sample_barcode", "sample_barcode", "sample")
  )
  
  # Important:
  # In GDC MAF, "gene" can be Ensembl ID, e.g. ENSG00000141510.
  # Driver genes are HGNC symbols, e.g. TP53 / PIK3CA.
  # So force gene to use hugo_symbol or symbol.
  if ("hugo_symbol" %in% names(maf_df)) {
    maf_df$gene <- maf_df$hugo_symbol
  } else if ("symbol" %in% names(maf_df)) {
    maf_df$gene <- maf_df$symbol
  } else {
    stop("Cannot find hugo_symbol or symbol column in MAF file.")
  }
  
  maf_df <- safe_rename(
    maf_df,
    "variant_classification",
    c("variant_classification", "classification")
  )
  
  maf_df <- maf_df %>%
    mutate(
      patient_barcode = standardise_tcga_patient(patient_barcode),
      gene = as.character(gene),
      variant_classification = as.character(variant_classification)
    )
  
  maf_nonsyn <- maf_df %>%
    filter(variant_classification %in% non_syn_classes)
  
  mutation_burden <- maf_nonsyn %>%
    count(patient_barcode, name = "n_nonsynonymous_mutations") %>%
    mutate(log1p_n_nonsynonymous_mutations = log1p(n_nonsynonymous_mutations))
  
  driver_mut_matrix <- maf_nonsyn %>%
    filter(gene %in% driver_mut_genes) %>%
    distinct(patient_barcode, gene) %>%
    mutate(value = 1, alteration = paste0(gene, "_mut")) %>%
    select(patient_barcode, alteration, value) %>%
    tidyr::pivot_wider(
      names_from = alteration,
      values_from = value,
      values_fill = 0
    )
  
  df <- df %>%
    left_join(mutation_burden, by = "patient_barcode") %>%
    left_join(driver_mut_matrix, by = "patient_barcode")
  
  df$n_nonsynonymous_mutations[is.na(df$n_nonsynonymous_mutations)] <- 0
  df$log1p_n_nonsynonymous_mutations[is.na(df$log1p_n_nonsynonymous_mutations)] <- 0
  
  mut_cols <- paste0(driver_mut_genes, "_mut")
  mut_cols <- intersect(mut_cols, names(df))
  df[mut_cols] <- lapply(df[mut_cols], function(x) replace_na(x, 0))
  
} else {
  warning("MAF file not found. Mutation burden and mutation-driver analyses will be skipped.")
}

############################################################
# 6. Load CNV segments and calculate genomic instability
#    Fixed version:
#    1) calculate CNV burden at sample level first
#    2) merge overlapping intervals
#    3) select one representative sample per patient
############################################################

standardise_tcga_sample <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\.", "-")
  substr(x, 1, 16)
}

sum_merged_interval_bp <- function(start, end) {
  ok <- !is.na(start) & !is.na(end) & end >= start
  
  if (!any(ok)) {
    return(0)
  }
  
  intervals <- tibble::tibble(
    start = as.numeric(start[ok]),
    end = as.numeric(end[ok])
  ) %>%
    arrange(start, end)
  
  current_start <- intervals$start[1]
  current_end <- intervals$end[1]
  total_bp <- 0
  
  if (nrow(intervals) > 1) {
    for (i in 2:nrow(intervals)) {
      this_start <- intervals$start[i]
      this_end <- intervals$end[i]
      
      if (this_start <= current_end + 1) {
        current_end <- max(current_end, this_end)
      } else {
        total_bp <- total_bp + current_end - current_start + 1
        current_start <- this_start
        current_end <- this_end
      }
    }
  }
  
  total_bp <- total_bp + current_end - current_start + 1
  total_bp
}

seg_df <- read_table_auto(seg_path)

if (!is.null(seg_df)) {
  names(seg_df) <- clean_colnames(names(seg_df))
  
  seg_df <- safe_rename(
    seg_df,
    "sample",
    c("sample", "tumor_sample_barcode", "tumour_sample_barcode", "sample_id")
  )
  
  seg_df <- safe_rename(
    seg_df,
    "chromosome",
    c("chromosome", "chrom", "chr")
  )
  
  seg_df <- safe_rename(
    seg_df,
    "start",
    c("start", "start_position", "loc_start")
  )
  
  seg_df <- safe_rename(
    seg_df,
    "end",
    c("end", "end_position", "loc_end")
  )
  
  seg_df <- safe_rename(
    seg_df,
    "segment_mean",
    c("segment_mean", "seg_mean", "segmean", "log2_copy_ratio")
  )
  
  needed_seg_cols <- c("sample", "chromosome", "start", "end", "segment_mean")
  
  if (all(needed_seg_cols %in% names(seg_df))) {
    
    autosomal_genome_bp <- 2880000000
    
    cnv_segment_clean <- seg_df %>%
      mutate(
        sample_barcode = standardise_tcga_sample(.data$sample),
        patient_barcode = standardise_tcga_patient(.data$sample),
        sample_type = substr(sample_barcode, 14, 15),
        chromosome = str_replace_all(as.character(chromosome), "chr|CHR", ""),
        chromosome = as.character(chromosome),
        start = as.numeric(start),
        end = as.numeric(end),
        segment_mean = as.numeric(segment_mean),
        altered = abs(segment_mean) >= 0.2
      ) %>%
      filter(
        chromosome %in% as.character(1:22),
        !is.na(sample_barcode),
        !is.na(patient_barcode),
        !is.na(start),
        !is.na(end),
        !is.na(segment_mean),
        end >= start
      ) %>%
      distinct(
        sample_barcode,
        chromosome,
        start,
        end,
        segment_mean,
        .keep_all = TRUE
      )
    
    cnv_by_chr <- cnv_segment_clean %>%
      group_by(patient_barcode, sample_barcode, sample_type, chromosome) %>%
      summarise(
        total_chr_bp = sum_merged_interval_bp(start, end),
        altered_chr_bp = sum_merged_interval_bp(start[altered], end[altered]),
        cnv_altered_segment_count = sum(altered, na.rm = TRUE),
        .groups = "drop"
      )
    
    cnv_instability_sample_level <- cnv_by_chr %>%
      group_by(patient_barcode, sample_barcode, sample_type) %>%
      summarise(
        cnv_total_autosomal_bp = sum(total_chr_bp, na.rm = TRUE),
        cnv_altered_bp = sum(altered_chr_bp, na.rm = TRUE),
        cnv_altered_mb = cnv_altered_bp / 1e6,
        cnv_altered_segment_count = sum(cnv_altered_segment_count, na.rm = TRUE),
        fraction_genome_altered = if_else(
          cnv_total_autosomal_bp > 0,
          cnv_altered_bp / cnv_total_autosomal_bp,
          NA_real_
        ),
        cnv_coverage_ratio = cnv_total_autosomal_bp / autosomal_genome_bp,
        .groups = "drop"
      )
    
    readr::write_csv(
      cnv_instability_sample_level,
      file.path(out_dir, "rq2_cnv_instability_sample_level_qc.csv")
    )
    
    cnv_instability <- cnv_instability_sample_level %>%
      mutate(
        sample_priority = case_when(
          sample_type == "01" ~ 1L,   # primary solid tumour
          TRUE ~ 2L
        ),
        coverage_distance = abs(cnv_coverage_ratio - 1)
      ) %>%
      arrange(
        patient_barcode,
        sample_priority,
        coverage_distance,
        sample_barcode
      ) %>%
      group_by(patient_barcode) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      select(
        patient_barcode,
        cnv_sample_barcode = sample_barcode,
        cnv_sample_type = sample_type,
        cnv_total_autosomal_bp,
        cnv_coverage_ratio,
        cnv_altered_bp,
        cnv_altered_mb,
        cnv_altered_segment_count,
        fraction_genome_altered
      )
    
    if (any(cnv_instability$fraction_genome_altered > 1, na.rm = TRUE)) {
      warning("Some patient-level fraction_genome_altered values are still > 1. Please inspect rq2_cnv_instability_sample_level_qc.csv.")
    }
    
    message("CNV sample-level records:")
    print(nrow(cnv_instability_sample_level))
    
    message("CNV patient-level records after selecting one sample per patient:")
    print(nrow(cnv_instability))
    
    df <- df %>%
      left_join(cnv_instability, by = "patient_barcode")
    
  } else {
    warning("CNV segment file found, but required columns are missing. CNV instability skipped.")
  }
  
} else {
  warning("CNV segment file not found. FGA / CNV burden analysis will be skipped.")
}

############################################################
# 7. Optional: load GISTIC gene-level CNV driver alterations
############################################################

cnv_driver_def <- tibble::tribble(
  ~gene,    ~direction, ~alteration,
  "ERBB2",  "amp",      "ERBB2_amp",
  "MYC",    "amp",      "MYC_amp",
  "CCND1",  "amp",      "CCND1_amp",
  "FGFR1",  "amp",      "FGFR1_amp",
  "MDM2",   "amp",      "MDM2_amp",
  "PTEN",   "del",      "PTEN_del",
  "RB1",    "del",      "RB1_del",
  "BRCA1",  "del",      "BRCA1_del",
  "CDKN2A", "del",      "CDKN2A_del"
)

gistic_df <- read_table_auto(gistic_path)

if (!is.null(gistic_df)) {
  names(gistic_df) <- clean_colnames(names(gistic_df))
  
  gene_col <- find_first_col(gistic_df, c("gene_symbol", "hugo_symbol", "gene", "genesymbol"))
  
  if (!is.na(gene_col)) {
    gistic_long <- gistic_df %>%
      rename(gene = all_of(gene_col)) %>%
      filter(gene %in% cnv_driver_def$gene) %>%
      pivot_longer(
        cols = -gene,
        names_to = "sample",
        values_to = "gistic_value"
      ) %>%
      mutate(
        patient_barcode = standardise_tcga_patient(sample),
        gistic_value = as.numeric(gistic_value)
      ) %>%
      inner_join(cnv_driver_def, by = "gene") %>%
      mutate(
        value = case_when(
          direction == "amp" & gistic_value >= 2 ~ 1,
          direction == "del" & gistic_value <= -2 ~ 1,
          TRUE ~ 0
        )
      ) %>%
      group_by(patient_barcode, alteration) %>%
      summarise(value = max(value, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(
        names_from = alteration,
        values_from = value,
        values_fill = 0
      )
    
    df <- df %>%
      left_join(gistic_long, by = "patient_barcode")
    
    cnv_driver_cols <- cnv_driver_def$alteration
    cnv_driver_cols <- intersect(cnv_driver_cols, names(df))
    df[cnv_driver_cols] <- lapply(df[cnv_driver_cols], function(x) replace_na(x, 0))
  }
} else {
  message("No GISTIC file provided. Gene-level CNV driver alteration analysis will be skipped.")
}

############################################################
# 8. Prepare covariates
############################################################

candidate_covariates <- c("pam50", "age", "stage", "grade", "purity")
available_covariates <- intersect(candidate_covariates, names(df))

if ("pam50" %in% names(df)) {
  df$pam50 <- as.factor(df$pam50)
}

if ("stage" %in% names(df)) {
  df$stage <- as.factor(df$stage)
}

if ("grade" %in% names(df)) {
  df$grade <- as.factor(df$grade)
}

if ("age" %in% names(df)) {
  df$age <- as.numeric(df$age)
}

if ("purity" %in% names(df)) {
  df$purity <- as.numeric(df$purity)
}

message("Available covariates used in adjusted models:")
print(available_covariates)

############################################################
# 9. Unadjusted Spearman correlations
############################################################

continuous_outcomes <- c(
  "fraction_genome_altered",
  "cnv_altered_mb",
  "cnv_altered_segment_count",
  "n_nonsynonymous_mutations",
  "log1p_n_nonsynonymous_mutations",
  "hrd_score",
  "loh_score",
  "aneuploidy_score"
)

continuous_outcomes <- intersect(continuous_outcomes, names(df))

spearman_results <- map_dfr(continuous_outcomes, function(y) {
  tmp <- df %>%
    select(xci_score, all_of(y)) %>%
    drop_na()
  
  if (nrow(tmp) < 20 || length(unique(tmp[[y]])) < 5) {
    return(NULL)
  }
  
  ct <- suppressWarnings(cor.test(tmp$xci_score, tmp[[y]], method = "spearman"))
  
  tibble(
    outcome = y,
    n = nrow(tmp),
    spearman_rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}) %>%
  mutate(fdr_q = p.adjust(p_value, method = "BH"))

readr::write_csv(
  spearman_results,
  file.path(out_dir, "rq2_unadjusted_spearman_v4_vs_instability.csv")
)

############################################################
# 10. Adjusted linear models for genomic instability
############################################################

run_linear_model <- function(outcome, data, covariates) {
  model_vars <- c(outcome, "xci_score_z", covariates)
  
  model_df <- data %>%
    select(all_of(model_vars)) %>%
    drop_na()
  
  if (nrow(model_df) < 50 || length(unique(model_df[[outcome]])) < 5) {
    return(NULL)
  }
  
  rhs <- paste(c("xci_score_z", covariates), collapse = " + ")
  formula <- as.formula(paste(outcome, "~", rhs))
  
  fit <- lm(formula, data = model_df)
  
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term == "xci_score_z") %>%
    mutate(
      outcome = outcome,
      n = nrow(model_df),
      model = paste(deparse(formula), collapse = "")
    ) %>%
    select(outcome, n, term, estimate, conf.low, conf.high, p.value, model)
}

linear_results <- map_dfr(
  continuous_outcomes,
  run_linear_model,
  data = df,
  covariates = available_covariates
)

if (nrow(linear_results) > 0 && "p.value" %in% names(linear_results)) {
  linear_results <- linear_results %>%
    mutate(fdr_q = p.adjust(p.value, method = "BH"))
} else {
  warning("No adjusted linear model results were generated. Saving empty linear result table.")
  linear_results <- tibble::tibble(
    outcome = character(),
    n = integer(),
    term = character(),
    estimate = numeric(),
    conf.low = numeric(),
    conf.high = numeric(),
    p.value = numeric(),
    model = character(),
    fdr_q = numeric()
  )
}

readr::write_csv(
  linear_results,
  file.path(out_dir, "rq2_adjusted_linear_models_genomic_instability.csv")
)

############################################################
# 11. Adjusted logistic models for driver alterations
############################################################

driver_cols <- c(
  paste0(driver_mut_genes, "_mut"),
  cnv_driver_def$alteration
)

driver_cols <- intersect(driver_cols, names(df))

run_logistic_model <- function(outcome, data, covariates) {
  model_vars <- c(outcome, "xci_score_z", covariates)
  
  model_df <- data %>%
    select(all_of(model_vars)) %>%
    drop_na()
  
  model_df[[outcome]] <- as.numeric(model_df[[outcome]])
  
  n_event <- sum(model_df[[outcome]] == 1, na.rm = TRUE)
  n_nonevent <- sum(model_df[[outcome]] == 0, na.rm = TRUE)
  
  if (nrow(model_df) < 50 || n_event < 10 || n_nonevent < 10) {
    return(NULL)
  }
  
  rhs <- paste(c("xci_score_z", covariates), collapse = " + ")
  formula <- as.formula(paste(outcome, "~", rhs))
  
  fit <- glm(formula, data = model_df, family = binomial())
  
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term == "xci_score_z") %>%
    mutate(
      outcome = outcome,
      n = nrow(model_df),
      n_event = n_event,
      odds_ratio = exp(estimate),
      or_conf_low = exp(conf.low),
      or_conf_high = exp(conf.high),
      model = paste(deparse(formula), collapse = "")
    ) %>%
    select(
      outcome, n, n_event, term,
      estimate, conf.low, conf.high,
      odds_ratio, or_conf_low, or_conf_high,
      p.value, model
    )
}

logistic_results <- map_dfr(
  driver_cols,
  run_logistic_model,
  data = df,
  covariates = available_covariates
)

if (nrow(logistic_results) > 0 && "p.value" %in% names(logistic_results)) {
  logistic_results <- logistic_results %>%
    mutate(fdr_q = p.adjust(p.value, method = "BH"))
} else {
  warning("No adjusted logistic model results were generated. Saving empty logistic result table.")
  logistic_results <- tibble::tibble(
    outcome = character(),
    n = integer(),
    n_event = integer(),
    term = character(),
    estimate = numeric(),
    conf.low = numeric(),
    conf.high = numeric(),
    odds_ratio = numeric(),
    or_conf_low = numeric(),
    or_conf_high = numeric(),
    p.value = numeric(),
    model = character(),
    fdr_q = numeric()
  )
}

readr::write_csv(
  logistic_results,
  file.path(out_dir, "rq2_adjusted_logistic_models_driver_alterations.csv")
)

############################################################
# 12. Four-level erosion quartile descriptive comparison
############################################################

group_df <- df %>%
  filter(
    xci_erosion_group %in% quartile_levels
  )

group_continuous_results <- map_dfr(continuous_outcomes, function(y) {
  tmp <- group_df %>%
    select(xci_erosion_group, all_of(y)) %>%
    drop_na()

  if (nrow(tmp) < 20 || length(unique(tmp[[y]])) < 5) return(NULL)

  kt <- kruskal.test(tmp[[y]] ~ tmp$xci_erosion_group)

  tmp %>%
    group_by(xci_erosion_group) %>%
    summarise(
      n = n(),
      median = median(.data[[y]], na.rm = TRUE),
      iqr = IQR(.data[[y]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      outcome = y,
      p_value = kt$p.value,
      test = "Kruskal-Wallis"
    )
})

if (nrow(group_continuous_results) > 0 && "outcome" %in% names(group_continuous_results)) {
  continuous_p_table <- group_continuous_results %>%
    distinct(outcome, p_value) %>%
    mutate(fdr_q = p.adjust(p_value, method = "BH"))

  group_continuous_results <- group_continuous_results %>%
    left_join(continuous_p_table, by = c("outcome", "p_value"))
} else {
  warning("No quartile-group continuous outcome results were generated. Saving empty continuous comparison table.")

  group_continuous_results <- tibble::tibble(
    xci_erosion_group = factor(
      levels = quartile_levels
    ),
    n = integer(),
    median = numeric(),
    iqr = numeric(),
    outcome = character(),
    p_value = numeric(),
    test = character(),
    fdr_q = numeric()
  )
}

readr::write_csv(
  group_continuous_results,
  file.path(out_dir, "rq2_quartile_erosion_continuous_outcomes.csv")
)

group_driver_results <- map_dfr(driver_cols, function(y) {
  tmp <- group_df %>%
    select(xci_erosion_group, all_of(y)) %>%
    drop_na()

  if (nrow(tmp) < 20) return(NULL)

  tmp[[y]] <- as.numeric(tmp[[y]])

  tab <- table(tmp$xci_erosion_group, tmp[[y]])

  if (ncol(tab) < 2 || nrow(tab) < 2) return(NULL)

  p_val <- tryCatch(
    fisher.test(tab)$p.value,
    error = function(e) {
      suppressWarnings(chisq.test(tab)$p.value)
    }
  )

  tmp %>%
    group_by(xci_erosion_group) %>%
    summarise(
      n = n(),
      n_altered = sum(.data[[y]] == 1, na.rm = TRUE),
      alteration_frequency = n_altered / n,
      .groups = "drop"
    ) %>%
    mutate(
      alteration = y,
      p_value = p_val,
      test = "Fisher exact test"
    )
})

if (nrow(group_driver_results) > 0 && "alteration" %in% names(group_driver_results)) {
  driver_p_table <- group_driver_results %>%
    distinct(alteration, p_value) %>%
    mutate(fdr_q = p.adjust(p_value, method = "BH"))

  group_driver_results <- group_driver_results %>%
    left_join(driver_p_table, by = c("alteration", "p_value"))
} else {
  warning("No quartile-group driver alteration frequency results were generated. Saving empty driver frequency table.")

  group_driver_results <- tibble::tibble(
    xci_erosion_group = factor(
      levels = quartile_levels
    ),
    n = integer(),
    n_altered = integer(),
    alteration_frequency = numeric(),
    alteration = character(),
    p_value = numeric(),
    test = character(),
    fdr_q = numeric()
  )
}

readr::write_csv(
  group_driver_results,
  file.path(out_dir, "rq2_quartile_erosion_driver_frequencies.csv")
)

############################################################
# 13. Plots
############################################################

plot_scatter <- function(y) {
  tmp <- df %>%
    select(xci_score, all_of(y)) %>%
    drop_na()
  
  if (nrow(tmp) < 20 || length(unique(tmp[[y]])) < 5) return(NULL)
  
  p <- ggplot(tmp, aes(x = xci_score, y = .data[[y]])) +
    geom_point(alpha = 0.45) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Association between V4 XCI erosion score and", y),
      x = "V4 XCI-linked gene-only erosion score",
      y = y
    )
  
  ggsave(
    filename = file.path(out_dir, paste0("scatter_v4_xci_vs_", y, ".png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

walk(continuous_outcomes, plot_scatter)

plot_box <- function(y) {
  tmp <- group_df %>%
    select(xci_erosion_group, all_of(y)) %>%
    drop_na()
  
  if (nrow(tmp) < 20 || length(unique(tmp[[y]])) < 5) return(NULL)
  
  p <- ggplot(tmp, aes(x = xci_erosion_group, y = .data[[y]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.18, alpha = 0.45) +
    theme_bw(base_size = 14) +
    labs(
      title = paste(y, "across V4 XCI erosion quartiles"),
      x = "V4 XCI erosion group",
      y = y
    )
  
  ggsave(
    filename = file.path(out_dir, paste0("boxplot_quartile_v4_", y, ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

walk(continuous_outcomes, plot_box)

if (nrow(group_driver_results) > 0) {
  p_driver <- group_driver_results %>%
    ggplot(aes(x = alteration, y = alteration_frequency, fill = xci_erosion_group)) +
    geom_col(position = "dodge") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Driver alteration frequency across V4 XCI erosion quartiles",
      x = "Driver alteration",
      y = "Alteration frequency",
      fill = "V4 erosion quartile"
    )
  
  ggsave(
    filename = file.path(out_dir, "barplot_driver_frequency_quartile_v4_erosion.png"),
    plot = p_driver,
    width = 10,
    height = 5.5,
    dpi = 300
  )
}

############################################################
# 14. Save final merged analysis table
############################################################

readr::write_csv(
  df,
  file.path(out_dir, "rq2_v4_merged_analysis_table.csv")
)

message("RQ2 V4 analysis completed.")
message("Results saved to: ", out_dir)

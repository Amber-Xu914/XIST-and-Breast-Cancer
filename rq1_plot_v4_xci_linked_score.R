############################################################
# RQ1 V4 plotting and validation
#
# Input:
# data/rq1_v4_xci_linked_gene_only_integrated_score.csv
#
# Main purpose:
# 1. Plot V4 XCI-linked gene-only erosion score
# 2. Check component behaviour
# 3. Check association with XIST expression
# 4. Check association with PAM50, tumour stage, and grade if available
# 5. Compare V4 with V3 if V3 file exists
############################################################

############################################################
# 0. packages
############################################################

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "janitor",
  "stringr",
  "broom",
  "forcats",
  "tibble"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(janitor)
  library(stringr)
  library(broom)
  library(forcats)
  library(tibble)
})

############################################################
# 1. paths
############################################################

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, mustWork = TRUE)

data_dir <- file.path(project_dir, "data")

v4_file <- file.path(
  data_dir,
  "rq1_v4_xci_linked_gene_only_integrated_score.csv"
)

v3_file <- file.path(
  data_dir,
  "rq1_v3_expression_methylation_cnv_integrated_score.csv"
)

expr_v1_file <- file.path(
  data_dir,
  "analysis_dataset_with_XCI_erosion_score.csv"
)

metadata_file <- file.path(
  data_dir,
  "xist_stage_pam50_metadata.csv"
)

output_dir <- file.path(
  project_dir,
  "outputs",
  "RQ1_V4_XCI_linked"
)

figure_dir <- file.path(output_dir, "figures")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
# 2. helper functions
############################################################

safe_spearman <- function(data, x_col, y_col, label) {
  
  test_data <- data %>%
    dplyr::filter(
      !is.na(.data[[x_col]]),
      !is.na(.data[[y_col]])
    )
  
  if (nrow(test_data) < 5) {
    return(
      tibble::tibble(
        comparison = label,
        x = x_col,
        y = y_col,
        n = nrow(test_data),
        spearman_rho = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  ct <- suppressWarnings(
    cor.test(
      test_data[[x_col]],
      test_data[[y_col]],
      method = "spearman",
      exact = FALSE
    )
  )
  
  tibble::tibble(
    comparison = label,
    x = x_col,
    y = y_col,
    n = nrow(test_data),
    spearman_rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

kw_test_with_epsilon <- function(data, score_col, group_col, label) {
  
  test_data <- data %>%
    dplyr::filter(
      !is.na(.data[[score_col]]),
      !is.na(.data[[group_col]])
    )
  
  n <- nrow(test_data)
  k <- dplyr::n_distinct(test_data[[group_col]])
  
  if (n < 5 || k < 2) {
    return(
      tibble::tibble(
        comparison = label,
        score = score_col,
        group = group_col,
        n = n,
        groups = k,
        statistic = NA_real_,
        p_value = NA_real_,
        epsilon_squared = NA_real_
      )
    )
  }
  
  kt <- kruskal.test(
    stats::as.formula(paste(score_col, "~", group_col)),
    data = test_data
  )
  
  h <- unname(kt$statistic)
  
  eps <- max((h - k + 1) / (n - k), 0)
  
  tibble::tibble(
    comparison = label,
    score = score_col,
    group = group_col,
    n = n,
    groups = k,
    statistic = h,
    p_value = kt$p.value,
    epsilon_squared = eps
  )
}

nice_p <- function(p) {
  ifelse(
    is.na(p),
    "NA",
    ifelse(p < 0.001, formatC(p, format = "e", digits = 2), signif(p, 3))
  )
}

############################################################
# 3. load V4 data
############################################################

if (!file.exists(v4_file)) {
  stop("Cannot find V4 file: ", v4_file)
}

v4 <- readr::read_csv(v4_file, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    patient_barcode = stringr::str_sub(patient_barcode, 1, 12)
  )

score_col <- "xci_erosion_score_v4_xci_linked_gene_only"

if (!score_col %in% colnames(v4)) {
  stop("Cannot find V4 score column: ", score_col)
}

cat("\nLoaded V4 data:\n")
print(dim(v4))

cat("\nV4 columns:\n")
print(colnames(v4))

############################################################
# 4. add XIST expression if missing
############################################################

if (!"xist_log_cpm" %in% colnames(v4)) {
  
  xist_sources <- list()
  
  if (file.exists(expr_v1_file)) {
    xist_sources[["v1"]] <- readr::read_csv(expr_v1_file, show_col_types = FALSE) %>%
      janitor::clean_names()
  }
  
  if (file.exists(metadata_file)) {
    xist_sources[["metadata"]] <- readr::read_csv(metadata_file, show_col_types = FALSE) %>%
      janitor::clean_names()
  }
  
  xist_lookup_list <- list()
  
  for (nm in names(xist_sources)) {
    
    temp <- xist_sources[[nm]]
    
    if ("patient_barcode" %in% colnames(temp)) {
      temp <- temp %>%
        dplyr::mutate(patient_barcode = stringr::str_sub(patient_barcode, 1, 12))
    } else if ("patient_id" %in% colnames(temp)) {
      temp <- temp %>%
        dplyr::mutate(patient_barcode = stringr::str_sub(patient_id, 1, 12))
    } else if ("sample_barcode" %in% colnames(temp)) {
      temp <- temp %>%
        dplyr::mutate(patient_barcode = stringr::str_sub(sample_barcode, 1, 12))
    } else {
      next
    }
    
    xist_col <- intersect(
      c("xist_log_cpm", "xist_logcpm"),
      colnames(temp)
    )[1]
    
    if (!is.na(xist_col)) {
      xist_lookup_list[[nm]] <- temp %>%
        dplyr::select(
          patient_barcode,
          xist_log_cpm = dplyr::all_of(xist_col)
        )
    }
  }
  
  if (length(xist_lookup_list) > 0) {
    
    xist_lookup <- dplyr::bind_rows(xist_lookup_list) %>%
      dplyr::filter(!is.na(patient_barcode), !is.na(xist_log_cpm)) %>%
      dplyr::group_by(patient_barcode) %>%
      dplyr::summarise(
        xist_log_cpm = mean(xist_log_cpm, na.rm = TRUE),
        .groups = "drop"
      )
    
    v4 <- v4 %>%
      dplyr::left_join(xist_lookup, by = "patient_barcode")
    
    cat("\nAdded XIST expression from external files.\n")
    cat("Patients with XIST expression:\n")
    print(sum(!is.na(v4$xist_log_cpm)))
    
  } else {
    warning("Could not find XIST expression column in available files.")
  }
}

############################################################
# 5. add V3 score if available
############################################################

if (file.exists(v3_file)) {
  
  v3 <- readr::read_csv(v3_file, show_col_types = FALSE) %>%
    janitor::clean_names() %>%
    dplyr::mutate(
      patient_barcode = stringr::str_sub(patient_barcode, 1, 12)
    )
  
  possible_v3_cols <- c(
    "xci_erosion_score_v3_expr_methylation_cnv",
    "xci_erosion_score_v3_expression_methylation_cnv"
  )
  
  v3_score_col <- possible_v3_cols[
    possible_v3_cols %in% colnames(v3)
  ][1]
  
  if (!is.na(v3_score_col)) {
    v3_keep <- v3 %>%
      dplyr::select(
        patient_barcode,
        v3_integrated_score = dplyr::all_of(v3_score_col)
      ) %>%
      dplyr::distinct(patient_barcode, .keep_all = TRUE)
    
    v4 <- v4 %>%
      dplyr::left_join(v3_keep, by = "patient_barcode")
    
    cat("\nAdded V3 integrated score for comparison.\n")
    cat("Patients with V3 score:\n")
    print(sum(!is.na(v4$v3_integrated_score)))
  }
}

############################################################
# 6. clean factors
############################################################

if ("pam50" %in% colnames(v4)) {
  v4 <- v4 %>%
    dplyr::mutate(
      pam50 = as.factor(pam50)
    )
}

if (!"stage_clean" %in% colnames(v4) && "stage4" %in% colnames(v4)) {
  v4 <- v4 %>%
    dplyr::mutate(stage_clean = stage4)
}

if ("stage_clean" %in% colnames(v4)) {
  v4 <- v4 %>%
    dplyr::mutate(
      stage_clean = factor(
        stage_clean,
        levels = c("Stage I", "Stage II", "Stage III", "Stage IV"),
        ordered = TRUE
      )
    )
}

if ("grade_clean" %in% colnames(v4)) {
  v4 <- v4 %>%
    dplyr::mutate(
      grade_clean = as.factor(grade_clean)
    )
}

############################################################
# 7. save cleaned analysis dataset
############################################################

readr::write_csv(
  v4,
  file.path(output_dir, "rq1_v4_plotting_analysis_dataset.csv")
)

############################################################
# 8. missingness and component summary
############################################################

missingness_summary <- v4 %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_v4_score = sum(!is.na(.data[[score_col]])),
    n_xist = if ("xist_log_cpm" %in% colnames(v4)) sum(!is.na(xist_log_cpm)) else NA_integer_,
    n_pam50 = if ("pam50" %in% colnames(v4)) sum(!is.na(pam50)) else NA_integer_,
    n_stage = if ("stage_clean" %in% colnames(v4)) sum(!is.na(stage_clean)) else NA_integer_,
    n_grade = if ("grade_clean" %in% colnames(v4)) sum(!is.na(grade_clean)) else NA_integer_,
    n_v3 = if ("v3_integrated_score" %in% colnames(v4)) sum(!is.na(v3_integrated_score)) else NA_integer_
  )

readr::write_csv(
  missingness_summary,
  file.path(output_dir, "v4_missingness_summary.csv")
)

component_cols <- c(
  "expression_erosion_score_raw",
  "methylation_erosion_score_raw",
  "cnv_erosion_score_raw",
  "expression_erosion_z",
  "methylation_erosion_z",
  "cnv_erosion_z",
  score_col
)

component_summary <- v4 %>%
  dplyr::select(dplyr::any_of(component_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "component",
    values_to = "value"
  ) %>%
  dplyr::group_by(component) %>%
  dplyr::summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  component_summary,
  file.path(output_dir, "v4_component_summary.csv")
)

############################################################
# 9. Spearman correlations
############################################################

spearman_results <- list()

if ("xist_log_cpm" %in% colnames(v4)) {
  
  spearman_results[["xist_v4"]] <- safe_spearman(
    v4,
    "xist_log_cpm",
    score_col,
    "XIST expression vs V4 integrated score"
  )
  
  spearman_results[["xist_expr"]] <- safe_spearman(
    v4,
    "xist_log_cpm",
    "expression_erosion_z",
    "XIST expression vs expression component"
  )
  
  spearman_results[["xist_meth"]] <- safe_spearman(
    v4,
    "xist_log_cpm",
    "methylation_erosion_z",
    "XIST expression vs methylation component"
  )
  
  spearman_results[["xist_cnv"]] <- safe_spearman(
    v4,
    "xist_log_cpm",
    "cnv_erosion_z",
    "XIST expression vs CNV component"
  )
}

spearman_results[["v4_expr"]] <- safe_spearman(
  v4,
  score_col,
  "expression_erosion_z",
  "V4 integrated score vs expression component"
)

spearman_results[["v4_meth"]] <- safe_spearman(
  v4,
  score_col,
  "methylation_erosion_z",
  "V4 integrated score vs methylation component"
)

spearman_results[["v4_cnv"]] <- safe_spearman(
  v4,
  score_col,
  "cnv_erosion_z",
  "V4 integrated score vs CNV component"
)

if ("v3_integrated_score" %in% colnames(v4)) {
  spearman_results[["v3_v4"]] <- safe_spearman(
    v4,
    "v3_integrated_score",
    score_col,
    "V3 all-chrX integrated score vs V4 XCI-linked score"
  )
}

spearman_results <- dplyr::bind_rows(spearman_results)

readr::write_csv(
  spearman_results,
  file.path(output_dir, "v4_spearman_correlation_results.csv")
)

cat("\nSpearman correlation results:\n")
print(spearman_results)

############################################################
# 10. Kruskal-Wallis tests
############################################################

kw_results <- list()

if ("pam50" %in% colnames(v4)) {
  kw_results[["pam50"]] <- kw_test_with_epsilon(
    v4,
    score_col,
    "pam50",
    "V4 score across PAM50 subtypes"
  )
}

if ("stage_clean" %in% colnames(v4)) {
  kw_results[["stage"]] <- kw_test_with_epsilon(
    v4,
    score_col,
    "stage_clean",
    "V4 score across tumour stages"
  )
}

if ("grade_clean" %in% colnames(v4)) {
  kw_results[["grade"]] <- kw_test_with_epsilon(
    v4,
    score_col,
    "grade_clean",
    "V4 score across tumour grades"
  )
}

kw_results <- dplyr::bind_rows(kw_results)

readr::write_csv(
  kw_results,
  file.path(output_dir, "v4_kruskal_wallis_results.csv")
)

cat("\nKruskal-Wallis results:\n")
print(kw_results)

############################################################
# 11. Pairwise Wilcoxon for PAM50 if available
############################################################

if ("pam50" %in% colnames(v4)) {
  
  pam50_test_data <- v4 %>%
    dplyr::filter(
      !is.na(pam50),
      !is.na(.data[[score_col]])
    )
  
  if (nrow(pam50_test_data) > 5 && dplyr::n_distinct(pam50_test_data$pam50) >= 2) {
    
    pairwise_pam50 <- pairwise.wilcox.test(
      x = pam50_test_data[[score_col]],
      g = pam50_test_data$pam50,
      p.adjust.method = "BH",
      exact = FALSE
    )
    
    pairwise_pam50_df <- as.data.frame(as.table(pairwise_pam50$p.value)) %>%
      dplyr::filter(!is.na(Freq)) %>%
      dplyr::rename(
        group_1 = Var1,
        group_2 = Var2,
        adjusted_p_value = Freq
      ) %>%
      dplyr::mutate(
        comparison = "V4 pairwise PAM50 subtype comparison"
      ) %>%
      dplyr::select(
        comparison,
        group_1,
        group_2,
        adjusted_p_value
      )
    
    readr::write_csv(
      pairwise_pam50_df,
      file.path(output_dir, "v4_pairwise_pam50_wilcoxon_results.csv")
    )
  }
}

############################################################
# 12. Plot 1: V4 score distribution
############################################################

p_score_distribution <- ggplot(
  v4,
  aes(x = .data[[score_col]])
) +
  geom_histogram(bins = 40, alpha = 0.75) +
  geom_vline(
    xintercept = median(v4[[score_col]], na.rm = TRUE),
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    title = "Distribution of V4 XCI-linked gene-only erosion score",
    subtitle = "Dashed line shows median score",
    x = "V4 XCI-linked gene-only erosion score",
    y = "Number of patients"
  )

ggsave(
  file.path(figure_dir, "01_v4_score_distribution.png"),
  p_score_distribution,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
# 13. Plot 2: component z-score distributions
############################################################

component_long <- v4 %>%
  dplyr::select(
    patient_barcode,
    expression_erosion_z,
    methylation_erosion_z,
    cnv_erosion_z,
    dplyr::all_of(score_col)
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      expression_erosion_z,
      methylation_erosion_z,
      cnv_erosion_z,
      dplyr::all_of(score_col)
    ),
    names_to = "component",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    component = dplyr::recode(
      component,
      expression_erosion_z = "Expression z",
      methylation_erosion_z = "Methylation erosion z",
      cnv_erosion_z = "CNV z",
      !!score_col := "V4 integrated score"
    )
  )

p_component_distribution <- ggplot(
  component_long,
  aes(x = component, y = score)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
  theme_bw() +
  labs(
    title = "V4 score and component distributions",
    x = "Component",
    y = "Score"
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave(
  file.path(figure_dir, "02_v4_component_distributions.png"),
  p_component_distribution,
  width = 8,
  height = 5,
  dpi = 300
)

############################################################
# 14. Plot 3: component correlation heatmap
############################################################

component_corr_data <- v4 %>%
  dplyr::select(
    expression_erosion_z,
    methylation_erosion_z,
    cnv_erosion_z,
    dplyr::all_of(score_col)
  )

component_corr <- cor(
  component_corr_data,
  use = "pairwise.complete.obs",
  method = "spearman"
)

component_corr_df <- as.data.frame(as.table(component_corr)) %>%
  dplyr::rename(
    component_1 = Var1,
    component_2 = Var2,
    spearman_rho = Freq
  )

readr::write_csv(
  component_corr_df,
  file.path(output_dir, "v4_component_correlation_matrix.csv")
)

p_component_corr <- ggplot(
  component_corr_df,
  aes(x = component_1, y = component_2, fill = spearman_rho)
) +
  geom_tile() +
  geom_text(
    aes(label = round(spearman_rho, 2)),
    size = 3.5
  ) +
  theme_bw() +
  labs(
    title = "Spearman correlation between V4 components",
    x = NULL,
    y = NULL,
    fill = "rho"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  file.path(figure_dir, "03_v4_component_correlation_heatmap.png"),
  p_component_corr,
  width = 7,
  height = 6,
  dpi = 300
)

############################################################
# 15. Plot 4: XIST vs V4 and components
############################################################

if ("xist_log_cpm" %in% colnames(v4)) {
  
  xist_plot_data <- v4 %>%
    dplyr::select(
      patient_barcode,
      xist_log_cpm,
      expression_erosion_z,
      methylation_erosion_z,
      cnv_erosion_z,
      dplyr::all_of(score_col)
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        expression_erosion_z,
        methylation_erosion_z,
        cnv_erosion_z,
        dplyr::all_of(score_col)
      ),
      names_to = "score_type",
      values_to = "score_value"
    ) %>%
    dplyr::mutate(
      score_type = dplyr::recode(
        score_type,
        expression_erosion_z = "Expression component",
        methylation_erosion_z = "Methylation component",
        cnv_erosion_z = "CNV component",
        !!score_col := "V4 integrated score"
      )
    )
  
  p_xist_components <- ggplot(
    xist_plot_data,
    aes(x = xist_log_cpm, y = score_value)
  ) +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~ score_type, scales = "free_y") +
    theme_bw() +
    labs(
      title = "Association between XIST expression and V4 score components",
      x = "XIST expression (logCPM)",
      y = "Score"
    )
  
  ggsave(
    file.path(figure_dir, "04_xist_vs_v4_and_components.png"),
    p_xist_components,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  xist_v4_result <- spearman_results %>%
    dplyr::filter(comparison == "XIST expression vs V4 integrated score")
  
  xist_subtitle <- if (nrow(xist_v4_result) == 1) {
    paste0(
      "Spearman rho = ",
      round(xist_v4_result$spearman_rho, 3),
      ", p = ",
      nice_p(xist_v4_result$p_value),
      ", n = ",
      xist_v4_result$n
    )
  } else {
    ""
  }
  
  p_xist_v4 <- ggplot(
    v4,
    aes(x = xist_log_cpm, y = .data[[score_col]])
  ) +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_bw() +
    labs(
      title = "Association between XIST expression and V4 XCI-linked score",
      subtitle = xist_subtitle,
      x = "XIST expression (logCPM)",
      y = "V4 XCI-linked gene-only erosion score"
    )
  
  ggsave(
    file.path(figure_dir, "05_xist_vs_v4_integrated_score.png"),
    p_xist_v4,
    width = 7,
    height = 5,
    dpi = 300
  )
}

############################################################
# 16. Plot 5: V4 across PAM50
############################################################

if ("pam50" %in% colnames(v4)) {
  
  pam50_result <- kw_results %>%
    dplyr::filter(group == "pam50")
  
  pam50_subtitle <- if (nrow(pam50_result) == 1) {
    paste0(
      "Kruskal-Wallis p = ",
      nice_p(pam50_result$p_value),
      ", epsilon squared = ",
      round(pam50_result$epsilon_squared, 4),
      ", n = ",
      pam50_result$n
    )
  } else {
    ""
  }
  
  p_pam50 <- ggplot(
    v4 %>% dplyr::filter(!is.na(pam50)),
    aes(
      x = forcats::fct_reorder(pam50, .data[[score_col]], .fun = median),
      y = .data[[score_col]]
    )
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.18, alpha = 0.35, size = 1) +
    theme_bw() +
    labs(
      title = "V4 XCI-linked gene-only erosion score across PAM50 subtypes",
      subtitle = pam50_subtitle,
      x = "PAM50 molecular subtype",
      y = "V4 XCI-linked gene-only erosion score"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
  
  ggsave(
    file.path(figure_dir, "06_v4_score_across_pam50.png"),
    p_pam50,
    width = 8,
    height = 5,
    dpi = 300
  )
}

############################################################
# 17. Plot 6: V4 across tumour stage
############################################################

if ("stage_clean" %in% colnames(v4)) {
  
  stage_result <- kw_results %>%
    dplyr::filter(group == "stage_clean")
  
  stage_subtitle <- if (nrow(stage_result) == 1) {
    paste0(
      "Kruskal-Wallis p = ",
      nice_p(stage_result$p_value),
      ", epsilon squared = ",
      round(stage_result$epsilon_squared, 4),
      ", n = ",
      stage_result$n
    )
  } else {
    ""
  }
  
  p_stage <- ggplot(
    v4 %>% dplyr::filter(!is.na(stage_clean)),
    aes(x = stage_clean, y = .data[[score_col]])
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.18, alpha = 0.35, size = 1) +
    theme_bw() +
    labs(
      title = "V4 XCI-linked gene-only erosion score across tumour stages",
      subtitle = stage_subtitle,
      x = "Tumour stage",
      y = "V4 XCI-linked gene-only erosion score"
    )
  
  ggsave(
    file.path(figure_dir, "07_v4_score_across_stage.png"),
    p_stage,
    width = 7,
    height = 5,
    dpi = 300
  )
}

############################################################
# 18. Plot 7: V4 across grade if available
############################################################

if ("grade_clean" %in% colnames(v4)) {
  
  grade_result <- kw_results %>%
    dplyr::filter(group == "grade_clean")
  
  grade_subtitle <- if (nrow(grade_result) == 1) {
    paste0(
      "Kruskal-Wallis p = ",
      nice_p(grade_result$p_value),
      ", epsilon squared = ",
      round(grade_result$epsilon_squared, 4),
      ", n = ",
      grade_result$n
    )
  } else {
    ""
  }
  
  p_grade <- ggplot(
    v4 %>% dplyr::filter(!is.na(grade_clean)),
    aes(x = grade_clean, y = .data[[score_col]])
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.18, alpha = 0.35, size = 1) +
    theme_bw() +
    labs(
      title = "V4 XCI-linked gene-only erosion score across tumour grades",
      subtitle = grade_subtitle,
      x = "Tumour grade",
      y = "V4 XCI-linked gene-only erosion score"
    )
  
  ggsave(
    file.path(figure_dir, "08_v4_score_across_grade.png"),
    p_grade,
    width = 7,
    height = 5,
    dpi = 300
  )
}

############################################################
# 19. Plot 8: V4 vs V3 if available
############################################################

if ("v3_integrated_score" %in% colnames(v4)) {
  
  v3_v4_result <- spearman_results %>%
    dplyr::filter(comparison == "V3 all-chrX integrated score vs V4 XCI-linked score")
  
  v3_v4_subtitle <- if (nrow(v3_v4_result) == 1) {
    paste0(
      "Spearman rho = ",
      round(v3_v4_result$spearman_rho, 3),
      ", p = ",
      nice_p(v3_v4_result$p_value),
      ", n = ",
      v3_v4_result$n
    )
  } else {
    ""
  }
  
  p_v3_v4 <- ggplot(
    v4,
    aes(x = v3_integrated_score, y = .data[[score_col]])
  ) +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_bw() +
    labs(
      title = "Comparison of V3 all-chrX score and V4 XCI-linked gene-only score",
      subtitle = v3_v4_subtitle,
      x = "V3 expression + methylation + CNV score",
      y = "V4 XCI-linked gene-only score"
    )
  
  ggsave(
    file.path(figure_dir, "09_v3_vs_v4_integrated_score.png"),
    p_v3_v4,
    width = 7,
    height = 5,
    dpi = 300
  )
}

############################################################
# 20. Plot 9: possible technical artifact checks
############################################################

technical_cols <- intersect(
  c(
    "n_xci_expression_genes_used",
    "n_xci_methylation_genes_used",
    "n_xci_cnv_genes_used",
    "xci_cnv_total_overlap_bp"
  ),
  colnames(v4)
)

if (length(technical_cols) > 0) {
  
  technical_long <- v4 %>%
    dplyr::select(
      patient_barcode,
      dplyr::all_of(score_col),
      dplyr::all_of(technical_cols)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(technical_cols),
      names_to = "technical_metric",
      values_to = "technical_value"
    )
  
  p_technical <- ggplot(
    technical_long,
    aes(x = technical_value, y = .data[[score_col]])
  ) +
    geom_point(alpha = 0.35, size = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    facet_wrap(~ technical_metric, scales = "free_x") +
    theme_bw() +
    labs(
      title = "Technical sanity checks for V4 score",
      subtitle = "Strong trends may suggest technical bias or missingness effects",
      x = "Technical metric",
      y = "V4 XCI-linked gene-only erosion score"
    )
  
  ggsave(
    file.path(figure_dir, "10_v4_technical_artifact_checks.png"),
    p_technical,
    width = 10,
    height = 7,
    dpi = 300
  )
}

############################################################
# 21. Plot 10: XIST expression by V4 high/low group
############################################################

if ("xist_log_cpm" %in% colnames(v4) &&
    "xci_erosion_v4_group" %in% colnames(v4)) {
  
  p_xist_group <- ggplot(
    v4 %>%
      dplyr::filter(
        !is.na(xist_log_cpm),
        !is.na(xci_erosion_v4_group)
      ),
    aes(x = xci_erosion_v4_group, y = xist_log_cpm)
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
    theme_bw() +
    labs(
      title = "XIST expression across V4 high and low erosion groups",
      x = "V4 erosion group",
      y = "XIST expression (logCPM)"
    )
  
  ggsave(
    file.path(figure_dir, "11_xist_expression_by_v4_high_low_group.png"),
    p_xist_group,
    width = 6,
    height = 5,
    dpi = 300
  )
}

############################################################
# 22. save figure list
############################################################

figure_list <- tibble::tibble(
  figure_file = list.files(figure_dir, full.names = FALSE)
)

readr::write_csv(
  figure_list,
  file.path(output_dir, "v4_figure_list.csv")
)

cat("\nSaved outputs to:\n")
cat(output_dir, "\n")

cat("\nSaved figures to:\n")
cat(figure_dir, "\n")

cat("\nMain output files:\n")
cat("1. rq1_v4_plotting_analysis_dataset.csv\n")
cat("2. v4_missingness_summary.csv\n")
cat("3. v4_component_summary.csv\n")
cat("4. v4_spearman_correlation_results.csv\n")
cat("5. v4_kruskal_wallis_results.csv\n")
cat("6. figures/*.png\n")

cat("\nDone.\n")

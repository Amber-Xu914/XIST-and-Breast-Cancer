library(tidyverse)
library(broom)
library(ggplot2)

# Create output directory if it does not exist
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# =========================
# 1. Load XCI erosion score dataset
# =========================

xci_df <- read.csv(
  "outputs/RQ1_V4_XCI_linked/rq1_v4_plotting_analysis_dataset.csv",
  stringsAsFactors = FALSE
)

# Check columns
colnames(xci_df)

# Select the V4 XCI erosion score
xci_df <- xci_df %>%
  mutate(
    patient_barcode = substr(patient_barcode, 1, 12),
    xci_score = xci_erosion_score_v4_xci_linked_gene_only
  ) %>%
  select(
    patient_barcode,
    xci_score
  ) %>%
  distinct(patient_barcode, .keep_all = TRUE)

# =========================
# 2. Load mutation / MAF data
# =========================

maf_df <- read.delim(
  "data/tcga_brca_mutect2.maf.tsv.gz",
  stringsAsFactors = FALSE,
  comment.char = "#"
)

# Check columns
colnames(maf_df)

# Define non-silent mutation types
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

# Identify TP53-mutated patients
tp53_df <- maf_df %>%
  mutate(
    patient_barcode = substr(Tumor_Sample_Barcode, 1, 12)
  ) %>%
  filter(
    Hugo_Symbol == "TP53",
    Variant_Classification %in% non_silent_classes
  ) %>%
  distinct(patient_barcode) %>%
  mutate(
    TP53_mut = 1
  )

cat("Number of TP53-mutated patients in the MAF file:\n")
print(nrow(tp53_df))

head(tp53_df)

# =========================
# 3. Merge XCI score with TP53 mutation status
# =========================

analysis_df <- xci_df %>%
  left_join(
    tp53_df,
    by = "patient_barcode"
  ) %>%
  mutate(
    TP53_mut = ifelse(is.na(TP53_mut), 0, TP53_mut)
  ) %>%
  filter(
    !is.na(xci_score),
    is.finite(xci_score)
  )

cat("\nNumber of patients included in the analysis:\n")
print(nrow(analysis_df))

cat("\nNumber of TP53-mutated patients in the analysis dataset:\n")
print(sum(analysis_df$TP53_mut == 1))

cat("\nNumber of TP53 wild-type patients in the analysis dataset:\n")
print(sum(analysis_df$TP53_mut == 0))

# =========================
# 4. Standardise XCI score and create sextiles
# =========================

analysis_df <- analysis_df %>%
  mutate(
    xci_score_z = as.numeric(scale(xci_score)),
    xci_sextile_number = ntile(xci_score, 6),
    xci_sextile = factor(
      xci_sextile_number,
      levels = c(1, 2, 3, 4, 5, 6),
      labels = c(
        "S1 lowest",
        "S2",
        "S3",
        "S4",
        "S5",
        "S6 highest"
      )
    ),
    high_xci = ifelse(
      xci_sextile == "S6 highest",
      1,
      0
    )
  )

# Check sample sizes across sextiles
print(table(analysis_df$xci_sextile))

# =========================
# 5. Descriptive table:
# TP53 mutation frequency by XCI sextile
# =========================

rq1_table <- analysis_df %>%
  group_by(xci_sextile) %>%
  summarise(
    n = n(),
    TP53_mut_n = sum(TP53_mut == 1),
    TP53_wildtype_n = sum(TP53_mut == 0),
    TP53_mut_percent = round(
      mean(TP53_mut == 1) * 100,
      1
    ),
    .groups = "drop"
  )

print(rq1_table)

write.csv(
  rq1_table,
  "outputs/aim3_rq1_tp53_frequency_by_xci_sextile.csv",
  row.names = FALSE
)

# =========================
# 6. Chi-square test
# =========================

sextile_table <- table(
  analysis_df$xci_sextile,
  analysis_df$TP53_mut
)

print(sextile_table)

chisq_result <- chisq.test(sextile_table)

print(chisq_result)

cat("\nExpected cell counts:\n")
print(chisq_result$expected)

# Save chi-square result
chisq_output <- tibble(
  statistic = unname(chisq_result$statistic),
  degrees_of_freedom = unname(chisq_result$parameter),
  p_value = chisq_result$p.value,
  method = chisq_result$method
)

write.csv(
  chisq_output,
  "outputs/aim3_rq1_chisquare_xci_sextile_tp53.csv",
  row.names = FALSE
)

# Calculate Cramer's V
cramers_v <- sqrt(
  unname(chisq_result$statistic) /
    sum(sextile_table)
)

cat("\nCramer's V:\n")
print(cramers_v)

# =========================
# 7. Univariable logistic regression
# Continuous XCI score model
# =========================

model_rq1 <- glm(
  TP53_mut ~ xci_score_z,
  data = analysis_df,
  family = binomial
)

summary(model_rq1)

model_rq1_result <- tidy(
  model_rq1,
  exponentiate = TRUE,
  conf.int = TRUE
)

print(model_rq1_result)

write.csv(
  model_rq1_result,
  "outputs/aim3_rq1_univariable_logistic_tp53_xci.csv",
  row.names = FALSE
)

# Extract XCI score result for plot annotation
xci_model_result <- model_rq1_result %>%
  filter(term == "xci_score_z")

xci_or <- xci_model_result$estimate
xci_ci_low <- xci_model_result$conf.low
xci_ci_high <- xci_model_result$conf.high
xci_p_value <- xci_model_result$p.value

model_annotation <- paste0(
  "OR per 1-SD increase = ",
  sprintf("%.2f", xci_or),
  " (95% CI ",
  sprintf("%.2f", xci_ci_low),
  "-",
  sprintf("%.2f", xci_ci_high),
  "), p = ",
  format.pval(
    xci_p_value,
    digits = 3,
    eps = 0.001
  )
)

cat("\nContinuous logistic regression result:\n")
cat(model_annotation, "\n")

# =========================
# 8. Trend test across XCI sextiles
# =========================

analysis_df <- analysis_df %>%
  mutate(
    xci_sextile_numeric = as.numeric(xci_sextile)
  )

model_trend <- glm(
  TP53_mut ~ xci_sextile_numeric,
  data = analysis_df,
  family = binomial
)

summary(model_trend)

model_trend_result <- tidy(
  model_trend,
  exponentiate = TRUE,
  conf.int = TRUE
)

print(model_trend_result)

write.csv(
  model_trend_result,
  "outputs/aim3_rq1_logistic_trend_across_xci_sextiles.csv",
  row.names = FALSE
)

# =========================
# 9. Bar plot:
# TP53 mutation frequency across XCI sextiles
# =========================

plot_df <- analysis_df %>%
  group_by(xci_sextile) %>%
  summarise(
    n = n(),
    TP53_mut_n = sum(TP53_mut == 1),
    TP53_mut_percent = mean(TP53_mut == 1) * 100,
    se = sqrt(
      (TP53_mut_percent / 100) *
        (1 - TP53_mut_percent / 100) /
        n
    ) * 100,
    lower = pmax(
      TP53_mut_percent - 1.96 * se,
      0
    ),
    upper = pmin(
      TP53_mut_percent + 1.96 * se,
      100
    ),
    .groups = "drop"
  )

p1 <- ggplot(
  plot_df,
  aes(
    x = xci_sextile,
    y = TP53_mut_percent
  )
) +
  geom_col(
    width = 0.7
  ) +
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.15
  ) +
  geom_text(
    aes(
      label = paste0(
        round(TP53_mut_percent, 1),
        "%"
      )
    ),
    vjust = -0.5,
    size = 3.8
  ) +
  labs(
    title = "TP53 mutation frequency across XCI erosion sextiles",
    subtitle = paste0(
      "Pearson's chi-square test: ",
      "\u03C7\u00B2 = ",
      round(unname(chisq_result$statistic), 2),
      ", df = ",
      unname(chisq_result$parameter),
      ", p = ",
      format.pval(
        chisq_result$p.value,
        digits = 3,
        eps = 0.001
      )
    ),
    x = "XCI erosion sextile",
    y = "TP53 mutation frequency (%)"
  ) +
  expand_limits(
    y = max(plot_df$upper) + 5
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

print(p1)

ggsave(
  "outputs/aim3_rq1_tp53_frequency_by_xci_sextile.png",
  plot = p1,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

# =========================
# 10. Continuous logistic regression curve
# XCI erosion score vs predicted TP53 mutation probability
# =========================

prediction_df <- tibble(
  xci_score_z = seq(
    from = min(
      analysis_df$xci_score_z,
      na.rm = TRUE
    ),
    to = max(
      analysis_df$xci_score_z,
      na.rm = TRUE
    ),
    length.out = 300
  )
)

# Generate predictions on the log-odds scale
continuous_prediction <- predict(
  model_rq1,
  newdata = prediction_df,
  type = "link",
  se.fit = TRUE
)

# Convert log-odds to probabilities
prediction_df <- prediction_df %>%
  mutate(
    predicted_probability = plogis(
      continuous_prediction$fit
    ),
    lower_ci = plogis(
      continuous_prediction$fit -
        1.96 * continuous_prediction$se.fit
    ),
    upper_ci = plogis(
      continuous_prediction$fit +
        1.96 * continuous_prediction$se.fit
    )
  )

write.csv(
  prediction_df,
  "outputs/aim3_rq1_continuous_logistic_predictions.csv",
  row.names = FALSE
)

# Calculate observed TP53 mutation frequency by XCI score decile
# These groups are only used to display the observed data
observed_decile_df <- analysis_df %>%
  mutate(
    xci_decile = ntile(xci_score_z, 10)
  ) %>%
  group_by(xci_decile) %>%
  summarise(
    mean_xci_score_z = mean(
      xci_score_z,
      na.rm = TRUE
    ),
    median_xci_score_z = median(
      xci_score_z,
      na.rm = TRUE
    ),
    observed_tp53_probability = mean(
      TP53_mut == 1
    ),
    TP53_mut_n = sum(
      TP53_mut == 1
    ),
    n = n(),
    .groups = "drop"
  )

print(observed_decile_df)

write.csv(
  observed_decile_df,
  "outputs/aim3_rq1_observed_tp53_frequency_by_xci_decile.csv",
  row.names = FALSE
)

# Continuous logistic regression graph
p_continuous <- ggplot(
  prediction_df,
  aes(
    x = xci_score_z,
    y = predicted_probability
  )
) +
  geom_ribbon(
    aes(
      ymin = lower_ci,
      ymax = upper_ci
    ),
    alpha = 0.20
  ) +
  geom_line(
    linewidth = 1.2
  ) +
  geom_point(
    data = observed_decile_df,
    aes(
      x = mean_xci_score_z,
      y = observed_tp53_probability,
      size = n
    ),
    inherit.aes = FALSE,
    alpha = 0.85
  ) +
  geom_rug(
    data = analysis_df,
    aes(
      x = xci_score_z
    ),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.12
  ) +
  scale_y_continuous(
    labels = scales::label_percent(
      accuracy = 1
    ),
    limits = c(0, 1)
  ) +
  scale_size_continuous(
    name = "Patients per decile"
  ) +
  labs(
    title = "Continuous association between XCI erosion and TP53 mutation",
    subtitle = model_annotation,
    x = "Standardised XCI erosion score (z-score)",
    y = "Predicted probability of TP53 mutation",
    caption = paste(
      "Line: logistic regression prediction.",
      "Shaded region: 95% confidence interval.",
      "Points: observed mutation frequency by XCI score decile."
    )
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    plot.caption = element_text(
      hjust = 0
    )
  )

print(p_continuous)

ggsave(
  "outputs/aim3_rq1_continuous_xci_tp53_probability.png",
  plot = p_continuous,
  width = 8,
  height = 5.8,
  dpi = 300
)

# =========================
# 11. Boxplot:
# XCI score by TP53 mutation status
# =========================

p2 <- ggplot(
  analysis_df,
  aes(
    x = factor(
      TP53_mut,
      levels = c(0, 1),
      labels = c(
        "TP53 wild-type",
        "TP53-mutant"
      )
    ),
    y = xci_score_z
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.25,
    size = 1
  ) +
  labs(
    title = "XCI erosion score by TP53 mutation status",
    x = "TP53 mutation status",
    y = "Standardised XCI erosion score"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

print(p2)

ggsave(
  "outputs/aim3_rq1_xci_score_by_tp53_status.png",
  plot = p2,
  width = 6,
  height = 5,
  dpi = 300
)

# =========================
# 12. Save complete analysis dataset
# =========================

write.csv(
  analysis_df,
  "outputs/aim3_rq1_tp53_xci_complete_analysis_dataset.csv",
  row.names = FALSE
)

cat("\nAnalysis completed successfully.\n")

# =========================
# 13. Classify TP53 mutation types
# =========================

# Identify all patients represented in the MAF file
# This prevents patients without mutation sequencing data
# from being incorrectly classified as TP53 wild-type
maf_patient_df <- maf_df %>%
  mutate(
    patient_barcode = substr(Tumor_Sample_Barcode, 1, 12)
  ) %>%
  distinct(patient_barcode)

# Classify TP53 mutations into broader mutation categories
tp53_type_long <- maf_df %>%
  mutate(
    patient_barcode = substr(Tumor_Sample_Barcode, 1, 12),
    
    mutation_type = case_when(
      Variant_Classification == "Missense_Mutation" ~
        "Missense mutation",
      
      Variant_Classification == "Nonsense_Mutation" ~
        "Nonsense mutation",
      
      Variant_Classification %in% c(
        "Frame_Shift_Ins",
        "Frame_Shift_Del"
      ) ~ "Frameshift insertion/deletion",
      
      Variant_Classification %in% c(
        "In_Frame_Ins",
        "In_Frame_Del"
      ) ~ "In-frame insertion/deletion",
      
      Variant_Classification == "Splice_Site" ~
        "Splice-site mutation",
      
      Variant_Classification == "Translation_Start_Site" ~
        "Translation start-site mutation",
      
      Variant_Classification == "Nonstop_Mutation" ~
        "Nonstop mutation",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    Hugo_Symbol == "TP53",
    !is.na(mutation_type)
  ) %>%
  distinct(
    patient_barcode,
    mutation_type
  )

cat("\nNumber of TP53 mutation records by type:\n")

print(
  tp53_type_long %>%
    count(
      mutation_type,
      sort = TRUE
    )
)

# =========================
# 14. Create one mutation category per patient
# =========================

# A patient may contain more than one TP53 mutation type.
# These patients are classified as "Multiple TP53 mutation types"
# so each patient appears only once in the main plot.

patient_tp53_type <- tp53_type_long %>%
  group_by(patient_barcode) %>%
  summarise(
    number_of_types = n_distinct(mutation_type),
    
    mutation_type_list = paste(
      sort(unique(mutation_type)),
      collapse = "; "
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    mutation_type = ifelse(
      number_of_types == 1,
      mutation_type_list,
      "Multiple TP53 mutation types"
    )
  ) %>%
  select(
    patient_barcode,
    mutation_type,
    number_of_types,
    mutation_type_list
  )

# =========================
# 15. Merge mutation types with XCI erosion score
# =========================

mutation_analysis_df <- analysis_df %>%
  select(
    patient_barcode,
    xci_score,
    xci_score_z
  ) %>%
  
  # Keep only patients represented in the mutation dataset
  inner_join(
    maf_patient_df,
    by = "patient_barcode"
  ) %>%
  
  left_join(
    patient_tp53_type,
    by = "patient_barcode"
  ) %>%
  
  mutate(
    mutation_type = replace_na(
      mutation_type,
      "TP53 wild-type"
    ),
    
    number_of_types = replace_na(
      number_of_types,
      0
    ),
    
    mutation_type = factor(
      mutation_type,
      levels = c(
        "TP53 wild-type",
        "Missense mutation",
        "Nonsense mutation",
        "Frameshift insertion/deletion",
        "In-frame insertion/deletion",
        "Splice-site mutation",
        "Translation start-site mutation",
        "Nonstop mutation",
        "Multiple TP53 mutation types"
      )
    )
  )

cat("\nNumber of patients in each TP53 mutation category:\n")

print(
  mutation_analysis_df %>%
    count(
      mutation_type,
      .drop = FALSE
    )
)

# =========================
# 16. Descriptive summary table
# =========================

mutation_type_summary <- mutation_analysis_df %>%
  group_by(mutation_type) %>%
  summarise(
    n = n(),
    
    mean_xci_score_z = mean(
      xci_score_z,
      na.rm = TRUE
    ),
    
    sd_xci_score_z = sd(
      xci_score_z,
      na.rm = TRUE
    ),
    
    median_xci_score_z = median(
      xci_score_z,
      na.rm = TRUE
    ),
    
    q1_xci_score_z = quantile(
      xci_score_z,
      0.25,
      na.rm = TRUE
    ),
    
    q3_xci_score_z = quantile(
      xci_score_z,
      0.75,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

print(mutation_type_summary)

write.csv(
  mutation_type_summary,
  "outputs/aim3_tp53_mutation_type_xci_summary.csv",
  row.names = FALSE
)

write.csv(
  mutation_analysis_df,
  "outputs/aim3_tp53_mutation_type_xci_patient_dataset.csv",
  row.names = FALSE
)

# =========================
# 17. Jittered scatterplot with boxplots
# =========================

# Add sample size to each category label
mutation_type_counts <- mutation_analysis_df %>%
  count(
    mutation_type,
    .drop = FALSE
  ) %>%
  filter(n > 0)

mutation_type_labels <- setNames(
  paste0(
    as.character(mutation_type_counts$mutation_type),
    " (n = ",
    mutation_type_counts$n,
    ")"
  ),
  as.character(mutation_type_counts$mutation_type)
)

p_mutation_type <- ggplot(
  mutation_analysis_df,
  aes(
    x = mutation_type,
    y = xci_score_z
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.18,
    height = 0,
    alpha = 0.40,
    size = 1.4
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 3
  ) +
  scale_x_discrete(
    labels = mutation_type_labels,
    drop = TRUE
  ) +
  labs(
    title = "XCI erosion score across TP53 mutation types",
    subtitle = paste(
      "Each point represents one patient;",
      "the boxplot shows the median and interquartile range"
    ),
    x = "TP53 mutation category",
    y = "Standardised XCI erosion score (z-score)"
  ) +
  coord_flip() +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )

print(p_mutation_type)

ggsave(
  "outputs/aim3_xci_score_by_tp53_mutation_type.png",
  plot = p_mutation_type,
  width = 9,
  height = 7,
  dpi = 300
)

# =========================
# 18. Kruskal-Wallis test
# Compare XCI erosion between TP53 mutation types
# =========================

mutation_test_df <- mutation_analysis_df %>%
  filter(
    !mutation_type %in% c(
      "TP53 wild-type",
      "Multiple TP53 mutation types"
    )
  ) %>%
  add_count(
    mutation_type,
    name = "group_n"
  ) %>%
  filter(
    group_n >= 5
  ) %>%
  droplevels()

cat("\nMutation types included in the statistical comparison:\n")

print(
  mutation_test_df %>%
    distinct(
      mutation_type,
      group_n
    )
)

if (n_distinct(mutation_test_df$mutation_type) >= 2) {
  
  kruskal_result <- kruskal.test(
    xci_score_z ~ mutation_type,
    data = mutation_test_df
  )
  
  print(kruskal_result)
  
  kruskal_output <- tidy(
    kruskal_result
  )
  
  write.csv(
    kruskal_output,
    "outputs/aim3_kruskal_wallis_tp53_mutation_type_xci.csv",
    row.names = FALSE
  )
  
} else {
  
  cat(
    "\nNot enough mutation categories with at least five patients ",
    "to perform the Kruskal-Wallis test.\n"
  )
}

# =========================
# 19. Pairwise Wilcoxon tests
# Run only if the global test is significant
# =========================

if (
  exists("kruskal_result") &&
  kruskal_result$p.value < 0.05
) {
  
  pairwise_result <- pairwise.wilcox.test(
    x = mutation_test_df$xci_score_z,
    g = mutation_test_df$mutation_type,
    p.adjust.method = "BH",
    exact = FALSE
  )
  
  print(pairwise_result)
  
  pairwise_output <- as.data.frame(
    as.table(pairwise_result$p.value)
  ) %>%
    filter(
      !is.na(Freq)
    ) %>%
    rename(
      mutation_type_1 = Var1,
      mutation_type_2 = Var2,
      adjusted_p_value = Freq
    )
  
  write.csv(
    pairwise_output,
    "outputs/aim3_pairwise_tp53_mutation_type_xci.csv",
    row.names = FALSE
  )
}

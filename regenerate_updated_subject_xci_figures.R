suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, mustWork = TRUE)
setwd(project_dir)
out_dir <- file.path(project_dir, "outputs", "updated_subject_xci_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_thesis <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = "Times") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0, size = base_size),
      legend.title = element_text(face = "bold")
    )
}

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))
}

fmt_num <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

save_both <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot = plot,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot = plot,
         width = width, height = height, bg = "white")
}

merged <- read_csv("results/rq2_v4/rq2_v4_merged_analysis_table.csv", show_col_types = FALSE)
spearman <- read_csv("results/rq2_v4/rq2_unadjusted_spearman_v4_vs_instability.csv", show_col_types = FALSE)
linear <- read_csv("results/rq2_v4/rq2_adjusted_linear_models_genomic_instability.csv", show_col_types = FALSE)
quartile_stats <- read_csv("results/rq2_v4/rq2_quartile_erosion_continuous_outcomes.csv", show_col_types = FALSE)
driver <- read_csv("results/rq2_v4/rq2_quartile_erosion_driver_frequencies.csv", show_col_types = FALSE)
pred <- read_csv("outputs/aim3_rq1_continuous_logistic_predictions.csv", show_col_types = FALSE)
deciles <- read_csv("outputs/aim3_rq1_observed_tp53_frequency_by_xci_decile.csv", show_col_types = FALSE)
logistic <- read_csv("outputs/aim3_rq1_univariable_logistic_tp53_xci.csv", show_col_types = FALSE)
rq3 <- read_csv("outputs/Aim3_RQ3_genomic_instability/rq3_complete_analysis_dataset.csv", show_col_types = FALSE)

quartile_levels <- c(
  "Q1 low erosion", "Q2 low-intermediate erosion",
  "Q3 high-intermediate erosion", "Q4 high erosion"
)
quartile_short <- c("Q1", "Q2", "Q3", "Q4")
quartile_colours <- c("#4E79A7", "#59A14F", "#F28E2B", "#B07AA1")
merged <- merged %>% mutate(xci_erosion_group = factor(xci_erosion_group, levels = quartile_levels))
driver <- driver %>% mutate(xci_erosion_group = factor(xci_erosion_group, levels = quartile_levels))

# Figure 1: continuous FGA association
sp_fga <- spearman %>% filter(outcome == "fraction_genome_altered") %>% slice(1)
p1 <- merged %>%
  filter(!is.na(xci_score), !is.na(fraction_genome_altered)) %>%
  ggplot(aes(xci_score, fraction_genome_altered)) +
  geom_point(size = 1.7, alpha = 0.58, shape = 21, fill = "grey55", colour = "grey20", stroke = 0.25) +
  geom_smooth(method = "lm", se = TRUE, colour = "#2463D4", fill = "#BFC8D8", linewidth = 1) +
  labs(
    title = "Association between the V4 XCI erosion proxy score and fraction of genome altered",
    subtitle = sprintf("Spearman rho = %.3f, P %s; 20 strict subject-to-XCI genes",
                       sp_fga$spearman_rho, fmt_p(sp_fga$p_value)),
    x = "V4 subject-to-XCI composite proxy score",
    y = "Fraction of genome altered"
  ) +
  theme_thesis(13)
save_both(p1, "figure_1_v4_subject_xci_vs_fraction_genome_altered", 8.2, 5.4)

# Figure 2: stage-adjusted genomic-instability forest plot
outcome_labels <- c(
  cnv_altered_segment_count = "Altered copy-number segment count",
  n_nonsynonymous_mutations = "Non-synonymous mutation count",
  log1p_n_nonsynonymous_mutations = "Log-transformed non-synonymous mutation count",
  fraction_genome_altered = "Fraction genome altered",
  cnv_altered_mb = "Altered copy-number burden"
)
forest <- linear %>%
  filter(outcome %in% names(outcome_labels)) %>%
  mutate(
    outcome_label = outcome_labels[outcome],
    detail = sprintf("Estimate = %.3g (95%% CI %.3g to %.3g); P %s; FDR q %s; n = %d",
                     estimate, conf.low, conf.high, fmt_p(p.value), fmt_p(fdr_q), n)
  )
forest_order <- names(outcome_labels)
forest_xlabels <- c(
  cnv_altered_segment_count = "Change in number of altered segments",
  n_nonsynonymous_mutations = "Change in number of mutations",
  log1p_n_nonsynonymous_mutations = "Change in log1p mutation count",
  fraction_genome_altered = "Change in fraction genome altered",
  cnv_altered_mb = "Change in altered megabases"
)
make_forest_panel <- function(outcome) {
  tmp <- forest %>% filter(.data$outcome == .env$outcome)
  ggplot(tmp, aes(estimate, 1)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.14, linewidth = 0.65) +
    geom_point(shape = 21, size = 3.2, fill = "#2F5D9B", colour = "black") +
    scale_y_continuous(NULL, breaks = NULL, limits = c(0.75, 1.25)) +
    labs(title = tmp$outcome_label, subtitle = tmp$detail, x = forest_xlabels[[outcome]]) +
    theme_thesis(9.5) +
    theme(plot.title = element_text(hjust = 0, size = 10.5),
          plot.subtitle = element_text(hjust = 0, size = 8.3),
          axis.title.x = element_text(size = 8.5),
          plot.margin = margin(3, 8, 3, 8))
}
forest_panels <- lapply(forest_order, make_forest_panel)
title2 <- textGrob("Stage-adjusted associations between the V4 XCI erosion proxy score and genomic instability",
                   gp = gpar(fontfamily = "Times", fontface = "bold", fontsize = 16))
sub2 <- textGrob("Linear-model estimates per 1-SD increase in the score; each outcome uses its own horizontal scale",
                 gp = gpar(fontfamily = "Times", fontsize = 10))
g2 <- arrangeGrob(
  grobs = c(list(title2, sub2), forest_panels), ncol = 1,
  heights = c(0.55, 0.38, rep(1, length(forest_panels)))
)
save_both(g2, "figure_2_stage_adjusted_genomic_instability_forest", 10.5, 8.2)

# Figure 3: four quartile panels
panel_specs <- tibble::tribble(
  ~outcome, ~panel, ~label,
  "fraction_genome_altered", "A", "Fraction of genome altered",
  "cnv_altered_mb", "B", "Altered copy-number burden (Mb)",
  "cnv_altered_segment_count", "C", "Altered copy-number segment count",
  "log1p_n_nonsynonymous_mutations", "D", "log1p non-synonymous mutation count"
)
make_quartile_panel <- function(outcome, panel, label) {
  stat <- quartile_stats %>% filter(.data$outcome == .env$outcome) %>% slice(1)
  tmp <- merged %>% filter(!is.na(.data[[outcome]]), !is.na(xci_erosion_group))
  counts <- tmp %>% count(xci_erosion_group, .drop = FALSE)
  x_labels <- setNames(paste0(quartile_short, "\n(n = ", counts$n, ")"), quartile_levels)
  ggplot(tmp, aes(xci_erosion_group, .data[[outcome]])) +
    geom_jitter(width = 0.18, alpha = 0.34, size = 0.75, colour = "grey35") +
    geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey94", linewidth = 0.5) +
    stat_summary(fun = median, geom = "crossbar", width = 0.55,
                 colour = "#2F5D9B", linewidth = 0.8) +
    scale_x_discrete(labels = x_labels) +
    labs(
      title = paste0(panel, ". ", label),
      subtitle = paste0("Kruskal-Wallis P ", fmt_p(stat$p_value), "; FDR q ", fmt_p(stat$fdr_q)),
      x = NULL, y = label
    ) +
    theme_thesis(10) +
    theme(plot.title = element_text(hjust = 0, face = "bold"),
          plot.subtitle = element_text(hjust = 0), axis.text.x = element_text(size = 8.5))
}
quartile_panels <- lapply(seq_len(nrow(panel_specs)), function(i) {
  make_quartile_panel(panel_specs$outcome[i], panel_specs$panel[i], panel_specs$label[i])
})
title3 <- textGrob("Genomic instability across V4 XCI erosion proxy quartiles",
                   gp = gpar(fontfamily = "Times", fontface = "bold", fontsize = 17))
sub3 <- textGrob("Q1: lowest proxy score; Q4: highest proxy score; 20 strict subject-to-XCI genes",
                 gp = gpar(fontfamily = "Times", fontsize = 11))
g3 <- arrangeGrob(
  grobs = c(list(title3, sub3), quartile_panels), ncol = 2,
  layout_matrix = rbind(c(1, 1), c(2, 2), c(3, 4), c(5, 6)),
  heights = c(0.48, 0.34, 3, 3)
)
save_both(g3, "figure_3_genomic_instability_across_quartiles", 11.2, 8.2)

# Figure 4: recurrent driver mutation frequencies
driver_order <- c("AKT1", "BRCA1", "BRCA2", "CDH1", "ERBB2", "ESR1", "GATA3",
                  "MAP2K4", "MAP3K1", "PIK3CA", "PTEN", "RB1", "TP53")
driver_plot <- driver %>%
  mutate(gene = sub("_mut$", "", alteration)) %>%
  filter(gene %in% driver_order) %>%
  mutate(gene = factor(gene, levels = driver_order))
tp53_stat <- driver_plot %>% filter(gene == "TP53") %>% slice(1)
p4 <- ggplot(driver_plot, aes(gene, alteration_frequency, fill = xci_erosion_group)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.76) +
  scale_fill_manual(values = quartile_colours,
                    labels = paste0(c("Q1 Low", "Q2 Low-intermediate", "Q3 High-intermediate", "Q4 High"),
                                    " (n = 195)")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Recurrent driver mutation frequencies across V4 XCI erosion proxy quartiles",
    subtitle = paste0("TP53: Fisher's exact test, P ", fmt_p(tp53_stat$p_value),
                      "; FDR q ", fmt_p(tp53_stat$fdr_q)),
    x = "Recurrently altered gene", y = "Mutation frequency",
    fill = "V4 proxy-score quartile"
  ) +
  theme_thesis(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        legend.position = "right")
save_both(p4, "figure_4_driver_mutation_frequencies_across_quartiles", 10.2, 5.8)

# Figure 5: continuous TP53 probability
or_row <- logistic %>% filter(term == "xci_score_z") %>% slice(1)
rug_df <- merged %>% filter(!is.na(xci_score_z), !is.na(TP53_mut))
p5 <- ggplot(pred, aes(xci_score_z, predicted_probability)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "#C5CBD4", alpha = 0.8) +
  geom_line(colour = "#2463D4", linewidth = 1.2) +
  geom_point(data = deciles, aes(mean_xci_score_z, observed_tp53_probability),
             inherit.aes = FALSE, shape = 21, size = 3.5, fill = "grey35", colour = "black") +
  geom_rug(data = rug_df, aes(x = xci_score_z), inherit.aes = FALSE,
           sides = "b", alpha = 0.16, colour = "grey25") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Continuous association between the XCI erosion proxy score and TP53 mutation probability",
    subtitle = sprintf("OR per 1-SD increase = %.2f (95%% CI %.2f-%.2f), P %s",
                       or_row$estimate, or_row$conf.low, or_row$conf.high, fmt_p(or_row$p.value)),
    x = "Standardised subject-to-XCI composite proxy score (z-score)",
    y = "Predicted probability of TP53 mutation",
    caption = "Line: logistic-regression prediction; shaded region: 95% CI; points: observed frequency by decile."
  ) +
  theme_thesis(12) +
  theme(plot.caption = element_text(hjust = 0, colour = "grey35"))
save_both(p5, "figure_5_continuous_tp53_mutation_probability", 8.5, 5.7)

# Figure 6: formatted four-group summary table
group_levels <- c(
  "TP53 wild-type / lower XCI", "TP53 wild-type / high XCI",
  "TP53-mutant / lower XCI", "TP53-mutant / high XCI"
)
table_summary <- rq3 %>%
  mutate(genomic_group = factor(genomic_group, levels = group_levels)) %>%
  group_by(genomic_group) %>%
  summarise(
    n = n(),
    fga = sprintf("%.3f\n(%.3f-%.3f)", median(fga_autosomal), quantile(fga_autosomal, .25), quantile(fga_autosomal, .75)),
    abs_cnv = sprintf("%.3f\n(%.3f-%.3f)", median(weighted_abs_segment_mean), quantile(weighted_abs_segment_mean, .25), quantile(weighted_abs_segment_mean, .75)),
    segments = sprintf("%.0f\n(%.0f-%.0f)", median(cnv_segment_count), quantile(cnv_segment_count, .25), quantile(cnv_segment_count, .75)),
    tmb = sprintf("%.3f\n(%.3f-%.3f)", median(tmb_per_mb), quantile(tmb_per_mb, .25), quantile(tmb_per_mb, .75)),
    .groups = "drop"
  )
display_table <- data.frame(
  Outcome = c("Autosomal FGA", "Weighted absolute\nCNV magnitude", "Autosomal CNV\nsegment count", "Tumour mutation\nburden"),
  `TP53 WT / lower XCI` = c(table_summary$fga[1], table_summary$abs_cnv[1], table_summary$segments[1], table_summary$tmb[1]),
  `TP53 WT / high XCI` = c(table_summary$fga[2], table_summary$abs_cnv[2], table_summary$segments[2], table_summary$tmb[2]),
  `TP53 mutant / lower XCI` = c(table_summary$fga[3], table_summary$abs_cnv[3], table_summary$segments[3], table_summary$tmb[3]),
  `TP53 mutant / high XCI` = c(table_summary$fga[4], table_summary$abs_cnv[4], table_summary$segments[4], table_summary$tmb[4]),
  check.names = FALSE
)
colnames(display_table)[2:5] <- paste0(colnames(display_table)[2:5],
                                       "\n(n = ", table_summary$n, ")")
table_theme <- ttheme_minimal(
  base_size = 11,
  base_family = "Times",
  core = list(fg_params = list(hjust = 0.5, x = 0.5),
              bg_params = list(fill = rep(c("#F7F8FA", "#EEF2F5"), each = 5), col = "#B8C0CA")),
  colhead = list(fg_params = list(fontface = "bold"), bg_params = list(fill = "#E8EDF3", col = "#7B8DA5")),
  rowhead = list(fg_params = list(fontface = "bold"))
)
tg <- tableGrob(display_table, rows = NULL, theme = table_theme)
tg$heights <- unit(c(0.68, rep(0.52, 4)), "in")
title6 <- textGrob("Genomic-instability outcomes according to combined TP53 mutation and XCI proxy groups",
                   gp = gpar(fontfamily = "Times", fontface = "bold", fontsize = 15))
sub6 <- textGrob("Values are median (interquartile range); high XCI proxy = highest score quartile.",
                 gp = gpar(fontfamily = "Times", fontsize = 11))
foot6 <- textGrob("FGA, fraction genome altered; CNV, copy-number variation; WT, wild-type.",
                  x = 0, hjust = 0, gp = gpar(fontfamily = "Times", fontsize = 9, col = "grey30"))
g6 <- arrangeGrob(
  title6, sub6, tg, foot6,
  ncol = 1,
  heights = unit(c(0.50, 0.30, 2.78, 0.24), "in")
)
save_both(g6, "figure_6_tp53_xci_group_summary_table", 10.8, 4.05)

# Figure 7: horizontal FGA boxplots by combined groups
kw <- kruskal.test(fga_autosomal ~ genomic_group, data = rq3)
rq3_plot <- rq3 %>% mutate(genomic_group = factor(genomic_group, levels = group_levels))
group_labels <- setNames(paste0(group_levels, " (n = ", as.integer(table(rq3_plot$genomic_group)), ")"), group_levels)
p7 <- ggplot(rq3_plot, aes(fga_autosomal, genomic_group)) +
  geom_jitter(height = 0.14, alpha = 0.42, size = 1.05, colour = "grey40") +
  geom_boxplot(outlier.shape = NA, width = 0.58, fill = "grey96", linewidth = 0.65) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.2, fill = "white", colour = "black") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_discrete(labels = group_labels) +
  labs(
    title = "Autosomal fraction of genome altered according to TP53 status and XCI proxy score",
    subtitle = paste0("High XCI proxy: highest score quartile; Kruskal-Wallis P ", fmt_p(kw$p.value)),
    x = "Autosomal fraction of genome altered", y = "Combined TP53 mutation and XCI proxy group",
    caption = "Diamonds represent group means; centre lines represent medians."
  ) +
  theme_thesis(11) +
  theme(plot.caption = element_text(hjust = 0.5), axis.text.y = element_text(size = 9.5))
save_both(p7, "figure_7_fga_by_tp53_and_xci_group", 10.2, 6.0)

manifest <- tibble(
  figure = 1:7,
  png = c(
    "figure_1_v4_subject_xci_vs_fraction_genome_altered.png",
    "figure_2_stage_adjusted_genomic_instability_forest.png",
    "figure_3_genomic_instability_across_quartiles.png",
    "figure_4_driver_mutation_frequencies_across_quartiles.png",
    "figure_5_continuous_tp53_mutation_probability.png",
    "figure_6_tp53_xci_group_summary_table.png",
    "figure_7_fga_by_tp53_and_xci_group.png"
  ),
  score_panel = "HCCS; RAI2; PDK3; DYNLT3; CASK; ELK1; FAM120C; ZXDB; MSN; YIPF6; PGK1; BRWD3; BEX4; ACSL4; UBE2A; XIAP; SLC25A14; MMGT1; FMR1; FLNA"
)
write_csv(manifest, file.path(out_dir, "figure_manifest.csv"))

key_results <- tibble(
  result = c("FGA Spearman rho", "FGA Spearman P", "TP53 OR per SD", "TP53 OR lower CI", "TP53 OR upper CI", "TP53 OR P", "complete V4 n"),
  value = c(sp_fga$spearman_rho, sp_fga$p_value, or_row$estimate, or_row$conf.low,
            or_row$conf.high, or_row$p.value, sum(!is.na(merged$xci_score)))
)
write_csv(key_results, file.path(out_dir, "key_results.csv"))

message("Seven updated figures saved to: ", out_dir)

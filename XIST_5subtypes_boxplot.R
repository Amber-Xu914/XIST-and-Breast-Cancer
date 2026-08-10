############################################################
## TCGA-BRCA XIST expression:
## five PAM50 tumour subtypes plus normal breast tissue
##
## Important: tumour and normal samples are jointly TMM-normalised
## so their logCPM values are directly comparable.
############################################################

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(TCGAbiolinks)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, mustWork = TRUE)
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "outputs", "xist_normal_boxplot")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tumour_rds <- file.path(data_dir, "tcga_brca_star_counts_se.rds")
normal_rds <- file.path(data_dir, "tcga_brca_normal_se.rds")

stopifnot(file.exists(tumour_rds), file.exists(normal_rds))

se_tumour <- readRDS(tumour_rds)
se_normal <- readRDS(normal_rds)

keep_sample_type <- function(se, wanted) {
  meta <- as.data.frame(colData(se))

  if ("sample_type" %in% names(meta)) {
    keep <- meta$sample_type %in% wanted
  } else if ("definition" %in% names(meta)) {
    keep <- meta$definition %in% wanted
  } else if ("shortLetterCode" %in% names(meta)) {
    wanted_code <- if ("Primary Tumor" %in% wanted) "TP" else "NT"
    keep <- meta$shortLetterCode == wanted_code
  } else {
    stop("No recognised sample-type column in colData().")
  }

  se[, !is.na(keep) & keep]
}

se_tumour <- keep_sample_type(
  se_tumour,
  c("Primary Tumor", "Primary solid Tumor")
)
se_normal <- keep_sample_type(
  se_normal,
  c("Solid Tissue Normal", "Solid Tissue Normal, blood derived normal")
)

get_counts <- function(se) {
  assay_name <- if ("unstranded" %in% assayNames(se)) {
    "unstranded"
  } else {
    assayNames(se)[1]
  }
  assay(se, assay_name)
}

tumour_counts <- get_counts(se_tumour)
normal_counts <- get_counts(se_normal)

## Align genes before combining the two count matrices.
common_genes <- intersect(rownames(tumour_counts), rownames(normal_counts))
if (length(common_genes) == 0) {
  stop("Tumour and normal count matrices have no matching row names.")
}

tumour_counts <- tumour_counts[common_genes, , drop = FALSE]
normal_counts <- normal_counts[common_genes, , drop = FALSE]
all_counts <- cbind(tumour_counts, normal_counts)

## Find XIST using tumour row annotation, matched by row name.
row_annot <- as.data.frame(rowData(se_tumour))
symbol_col <- intersect(
  c("gene_name", "external_gene_name", "gene_symbol"),
  names(row_annot)
)[1]
if (is.na(symbol_col)) stop("Could not find a gene-symbol column in rowData().")

gene_symbols <- as.character(row_annot[[symbol_col]])
names(gene_symbols) <- rownames(row_annot)
xist_rows <- common_genes[gene_symbols[common_genes] == "XIST"]
if (length(xist_rows) == 0) stop("XIST was not found in the count matrices.")
xist_row <- xist_rows[1]

## One joint TMM normalisation across every tumour and normal sample.
dge <- DGEList(counts = all_counts)
dge <- calcNormFactors(dge, method = "TMM")
xist_logcpm <- as.numeric(
  cpm(dge, log = TRUE, prior.count = 1)[xist_row, ]
)
names(xist_logcpm) <- colnames(all_counts)

## Obtain PAM50 labels for tumour patients.
subtype_df <- PanCancerAtlas_subtypes() %>%
  filter(cancer.type == "BRCA") %>%
  transmute(
    patient_barcode = substr(pan.samplesID, 1, 12),
    raw_subtype = as.character(Subtype_Selected),
    Group = case_when(
      str_detect(raw_subtype, regex("Basal", ignore_case = TRUE)) ~ "Basal-like",
      str_detect(raw_subtype, regex("Her2", ignore_case = TRUE)) ~ "HER2-enriched",
      str_detect(raw_subtype, regex("LumA", ignore_case = TRUE)) ~ "Luminal A",
      str_detect(raw_subtype, regex("LumB", ignore_case = TRUE)) ~ "Luminal B",
      str_detect(raw_subtype, regex("Normal", ignore_case = TRUE)) ~ "Normal-like",
      TRUE ~ NA_character_
    )
  ) %>%
  distinct(patient_barcode, .keep_all = TRUE)

tumour_barcodes <- colnames(tumour_counts)
normal_barcodes <- colnames(normal_counts)

tumour_df <- tibble(
  sample_barcode = tumour_barcodes,
  patient_barcode = substr(tumour_barcodes, 1, 12),
  XIST_logCPM = xist_logcpm[tumour_barcodes]
) %>%
  left_join(subtype_df, by = "patient_barcode") %>%
  filter(!is.na(Group)) %>%
  select(sample_barcode, patient_barcode, XIST_logCPM, Group)

normal_df <- tibble(
  sample_barcode = normal_barcodes,
  patient_barcode = substr(normal_barcodes, 1, 12),
  XIST_logCPM = xist_logcpm[normal_barcodes],
  Group = "Normal tissue"
)

group_order <- c(
  "Basal-like", "HER2-enriched", "Luminal A",
  "Luminal B", "Normal-like", "Normal tissue"
)

plot_df <- bind_rows(tumour_df, normal_df) %>%
  mutate(Group = factor(Group, levels = group_order))

group_summary <- plot_df %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    median_logCPM = median(XIST_logCPM),
    IQR_logCPM = IQR(XIST_logCPM),
    .groups = "drop"
  )

label_lookup <- setNames(
  paste0(group_summary$Group, "\n(n = ", group_summary$n, ")"),
  group_summary$Group
)

## Overall and normal-vs-each-subtype non-parametric tests.
kruskal_result <- kruskal.test(XIST_logCPM ~ Group, data = plot_df)
normal_tests <- lapply(setdiff(group_order, "Normal tissue"), function(g) {
  test_data <- filter(plot_df, Group %in% c(g, "Normal tissue"))
  wt <- wilcox.test(XIST_logCPM ~ Group, data = test_data, exact = FALSE)
  tibble(comparison = paste("Normal tissue vs", g), p_value = wt$p.value)
}) %>%
  bind_rows() %>%
  mutate(p_adjust_BH = p.adjust(p_value, method = "BH"))

write.csv(
  plot_df,
  file.path(output_dir, "XIST_joint_TMM_plot_data.csv"),
  row.names = FALSE
)
write.csv(
  group_summary,
  file.path(output_dir, "XIST_joint_TMM_group_summary.csv"),
  row.names = FALSE
)
write.csv(
  normal_tests,
  file.path(output_dir, "XIST_normal_vs_PAM50_wilcoxon.csv"),
  row.names = FALSE
)
writeLines(
  sprintf(
    "Kruskal-Wallis chi-squared = %.4f, df = %d, p-value = %.6g",
    unname(kruskal_result$statistic),
    unname(kruskal_result$parameter),
    kruskal_result$p.value
  ),
  file.path(output_dir, "XIST_joint_TMM_kruskal_wallis.txt")
)

set.seed(20260730)
p <- ggplot(plot_df, aes(x = Group, y = XIST_logCPM)) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    fill = "#F2F2F2",
    colour = "#222222",
    linewidth = 0.55
  ) +
  geom_jitter(
    width = 0.14,
    height = 0,
    alpha = 0.34,
    size = 0.75,
    colour = "#4D4D4D"
  ) +
  scale_x_discrete(labels = label_lookup, drop = FALSE) +
  labs(
    title = expression(
      italic("XIST")~
        "normalised expression across normal tissue and PAM50 breast cancer subtypes"
    ),
    x = "Tissue group / PAM50 subtype",
    y = expression(italic("XIST")~"TMM-normalised expression (logCPM)")
  ) +
  theme_classic(base_size = 11, base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 11),
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    panel.grid.major.y = element_line(colour = "#DDDDDD", linewidth = 0.35),
    axis.line = element_line(colour = "#555555", linewidth = 0.45),
    plot.margin = margin(12, 12, 10, 12)
  )

ggsave(
  file.path(output_dir, "XIST_boxplot_with_normal_joint_TMM.png"),
  p,
  width = 10,
  height = 6.2,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(output_dir, "XIST_boxplot_with_normal_joint_TMM.pdf"),
  p,
  width = 10,
  height = 6.2,
  device = "pdf",
  useDingbats = FALSE
)

message("Finished. Outputs saved to: ", normalizePath(output_dir))

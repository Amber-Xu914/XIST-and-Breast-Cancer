############################################################
# RQ1 Version 4:
# Subject-to-XCI gene expression + promoter methylation + CNV proxy score
#
# Output:
# 1. data/rq1_v4_xci_linked_gene_only_integrated_score.csv
# 2. data/rq1_v4_xci_linked_gene_only_component_long.csv
############################################################

############################################################
# 0. packages
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "GenomicRanges",
  "IRanges",
  "dplyr",
  "tidyr",
  "stringr",
  "readr",
  "janitor",
  "tibble",
  "matrixStats",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  for (pkg in missing_packages) {
    if (pkg %in% c(
      "TCGAbiolinks",
      "SummarizedExperiment",
      "GenomicRanges",
      "IRanges",
      "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    )) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(janitor)
  library(tibble)
  library(matrixStats)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

############################################################
# 1. paths
############################################################

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else getwd()
project_dir <- normalizePath(project_dir, mustWork = TRUE)
data_dir <- file.path(project_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

setwd(project_dir)

xci_gene_candidates <- c(
  file.path(project_dir, "resources", "xci_linked_genes.csv"),
  file.path(data_dir, "xci_linked_genes.csv")
)
xci_gene_file <- xci_gene_candidates[file.exists(xci_gene_candidates)][1]

if (is.na(xci_gene_file)) {
  stop("Cannot find resources/xci_linked_genes.csv or data/xci_linked_genes.csv")
}

expr_rds_file <- file.path(data_dir, "tcga_brca_star_counts.rds")
meth_rds_file <- file.path(data_dir, "tcga_brca_methylation_450k.rds")
cnv_rds_file  <- file.path(data_dir, "tcga_brca_cnv_masked_segments.rds")

v4_output_file <- file.path(
  data_dir,
  "rq1_v4_xci_linked_gene_only_integrated_score.csv"
)

v4_long_file <- file.path(
  data_dir,
  "rq1_v4_xci_linked_gene_only_component_long.csv"
)

############################################################
# 2. helper functions
############################################################

z_score <- function(x) {
  as.numeric(scale(x))
}

standardise_chr <- function(x) {
  chr <- as.character(x)
  chr <- stringr::str_replace(chr, "^chr", "")
  chr <- stringr::str_replace(chr, "^CHR", "")
  chr <- toupper(chr)
  
  chr <- dplyr::if_else(
    chr %in% c("23"),
    "X",
    chr,
    missing = NA_character_
  )
  
  return(chr)
}

make_sample_table <- function(sample_barcodes) {
  tibble(
    sample_barcode = as.character(sample_barcodes),
    patient_barcode = stringr::str_sub(sample_barcode, 1, 12),
    sample_type_code = stringr::str_sub(sample_barcode, 14, 15)
  )
}

detect_gene_symbol_col <- function(df) {
  candidates <- c(
    "gene_name",
    "gene_symbol",
    "external_gene_name",
    "symbol",
    "hgnc_symbol"
  )
  detected <- candidates[candidates %in% colnames(df)][1]
  if (is.na(detected)) {
    stop("Cannot detect gene symbol column.")
  }
  detected
}

detect_sample_column <- function(df) {
  candidates <- c(
    "sample",
    "sample_barcode",
    "aliquot_barcode",
    "submitter_id",
    "barcode",
    "gdc_aliquot"
  )
  candidates <- candidates[candidates %in% colnames(df)]
  
  if (length(candidates) == 0) {
    stop("Cannot detect sample barcode column.")
  }
  
  scores <- sapply(candidates, function(col) {
    sum(stringr::str_detect(as.character(df[[col]]), "^TCGA-"), na.rm = TRUE)
  })
  
  candidates[which.max(scores)]
}

############################################################
# 3. load XCI-linked gene list
############################################################

if (!file.exists(xci_gene_file)) {
  template <- tibble(
    gene_symbol = c("XIST"),
    chromosome = c("X"),
    start = NA_real_,
    end = NA_real_
  )
  
  readr::write_csv(template, xci_gene_file)
  
  stop(
    "Created template file: ", xci_gene_file, "\n",
    "Please replace it with your literature-based XCI-linked gene list, then rerun this script."
  )
}

xci_genes <- readr::read_csv(xci_gene_file, show_col_types = FALSE) %>%
  janitor::clean_names()

if (!"gene_symbol" %in% colnames(xci_genes)) {
  stop("xci_linked_genes.csv must contain a gene_symbol column.")
}

xci_genes <- xci_genes %>%
  mutate(
    gene_symbol = stringr::str_trim(gene_symbol),
    gene_symbol_upper = toupper(gene_symbol)
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_symbol_upper, .keep_all = TRUE)

cat("\nNumber of XCI-linked genes provided:\n")
print(nrow(xci_genes))

############################################################
# 4. download or load TCGA-BRCA expression data
############################################################

if (file.exists(expr_rds_file)) {
  cat("\nLoading cached expression object:\n")
  cat(expr_rds_file, "\n")
  expr_se <- readRDS(expr_rds_file)
} else {
  cat("\nDownloading TCGA-BRCA STAR counts from GDC...\n")
  
  query_expr <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  
  GDCdownload(
    query_expr,
    method = "api",
    files.per.chunk = 50
  )
  
  expr_se <- GDCprepare(query_expr)
  saveRDS(expr_se, expr_rds_file)
}

############################################################
# 5. calculate XCI-linked expression score
############################################################

assay_names <- names(SummarizedExperiment::assays(expr_se))
cat("\nExpression assay names:\n")
print(assay_names)

preferred_assays <- c("unstranded", "stranded_first", "stranded_second")
expr_assay <- preferred_assays[preferred_assays %in% assay_names][1]

if (is.na(expr_assay)) {
  expr_assay <- assay_names[1]
}

cat("\nUsing expression assay:\n")
print(expr_assay)

count_mat <- SummarizedExperiment::assay(expr_se, expr_assay)

sample_barcodes <- as.character(colnames(count_mat))
sample_info <- make_sample_table(sample_barcodes)

primary_samples <- sample_info %>%
  filter(sample_type_code == "01") %>%
  arrange(patient_barcode, sample_barcode) %>%
  distinct(patient_barcode, .keep_all = TRUE)

count_mat <- count_mat[, primary_samples$sample_barcode, drop = FALSE]
colnames(count_mat) <- primary_samples$patient_barcode

# Convert counts to logCPM
lib_size <- colSums(count_mat, na.rm = TRUE)
cpm_mat <- t(t(count_mat) / lib_size) * 1e6
logcpm_mat <- log2(cpm_mat + 1)

expr_anno <- as.data.frame(SummarizedExperiment::rowData(expr_se)) %>%
  janitor::clean_names()

gene_col <- detect_gene_symbol_col(expr_anno)

expr_df <- as.data.frame(logcpm_mat) %>%
  tibble::rownames_to_column("row_id") %>%
  mutate(
    gene_symbol = expr_anno[[gene_col]],
    gene_symbol_upper = toupper(gene_symbol)
  ) %>%
  filter(gene_symbol_upper %in% xci_genes$gene_symbol_upper)

cat("\nNumber of XCI-linked genes found in expression data:\n")
print(n_distinct(expr_df$gene_symbol_upper))

if (nrow(expr_df) == 0) {
  stop("No XCI-linked genes found in expression matrix.")
}

expr_gene_level <- expr_df %>%
  select(-row_id, -gene_symbol) %>%
  group_by(gene_symbol_upper) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

expr_gene_mat <- expr_gene_level %>%
  column_to_rownames("gene_symbol_upper") %>%
  as.matrix()

expression_patient_score <- tibble(
  patient_barcode = colnames(expr_gene_mat),
  xci_linked_expression_mean_logcpm = colMeans(expr_gene_mat, na.rm = TRUE),
  xci_linked_expression_median_logcpm = matrixStats::colMedians(expr_gene_mat, na.rm = TRUE),
  n_xci_expression_genes_used = colSums(!is.na(expr_gene_mat))
) %>%
  mutate(
    expression_erosion_score_raw = xci_linked_expression_mean_logcpm
  )

cat("\nExpression patient score dimension:\n")
print(dim(expression_patient_score))

############################################################
# 6. extract gene coordinates for CNV
############################################################

gene_coord_from_file <- xci_genes %>%
  mutate(
    has_coord = all(c("chromosome", "start", "end") %in% colnames(xci_genes))
  )

if (all(c("chromosome", "start", "end") %in% colnames(xci_genes))) {
  
  gene_coord <- xci_genes %>%
    transmute(
      gene_symbol_upper,
      chromosome = standardise_chr(chromosome),
      start = as.numeric(start),
      end = as.numeric(end)
    ) %>%
    filter(
      !is.na(chromosome),
      !is.na(start),
      !is.na(end),
      chromosome == "X",
      end >= start
    )
  
} else {
  
  cat("\nNo chromosome/start/end columns found in xci_linked_genes.csv.\n")
  cat("Trying to extract coordinates from expression rowRanges...\n")
  
  rr <- SummarizedExperiment::rowRanges(expr_se)
  
  if (is.null(rr) || length(rr) == 0) {
    stop(
      "Cannot extract gene coordinates from expression object. ",
      "Please add chromosome, start, and end columns to xci_linked_genes.csv."
    )
  }
  
  if (inherits(rr, "GRangesList")) {
    rr_simple <- GenomicRanges::range(rr)
    rr_simple <- unlist(rr_simple)
  } else {
    rr_simple <- rr
  }
  
  coord_df <- as.data.frame(rr_simple) %>%
    janitor::clean_names()
  
  coord_df$gene_symbol_upper <- toupper(expr_anno[[gene_col]])
  
  gene_coord <- coord_df %>%
    transmute(
      gene_symbol_upper,
      chromosome = standardise_chr(seqnames),
      start = as.numeric(start),
      end = as.numeric(end)
    ) %>%
    filter(
      gene_symbol_upper %in% xci_genes$gene_symbol_upper,
      chromosome == "X",
      !is.na(start),
      !is.na(end),
      end >= start
    ) %>%
    group_by(gene_symbol_upper) %>%
    summarise(
      chromosome = "X",
      start = min(start, na.rm = TRUE),
      end = max(end, na.rm = TRUE),
      .groups = "drop"
    )
}

cat("\nNumber of XCI-linked genes with usable CNV coordinates:\n")
print(nrow(gene_coord))

if (nrow(gene_coord) == 0) {
  stop(
    "No usable XCI-linked gene coordinates for CNV. ",
    "Please add chromosome/start/end to xci_linked_genes.csv."
  )
}

############################################################
# 7. download or load methylation data
############################################################

if (file.exists(meth_rds_file)) {
  cat("\nLoading cached methylation object:\n")
  cat(meth_rds_file, "\n")
  meth_se <- readRDS(meth_rds_file)
} else {
  cat("\nDownloading TCGA-BRCA DNA methylation 450K data from GDC...\n")
  
  query_meth <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "DNA Methylation",
    data.type = "Methylation Beta Value",
    platform = "Illumina Human Methylation 450",
    sample.type = "Primary Tumor"
  )
  
  GDCdownload(
    query_meth,
    method = "client"
  )
  
  meth_se <- GDCprepare(query_meth)
  saveRDS(meth_se, meth_rds_file)
}

############################################################
# 8. calculate XCI-linked promoter methylation score
############################################################

beta_mat <- SummarizedExperiment::assay(meth_se)

meth_sample_info <- make_sample_table(colnames(beta_mat))

meth_primary <- meth_sample_info %>%
  filter(sample_type_code == "01") %>%
  arrange(patient_barcode, sample_barcode) %>%
  distinct(patient_barcode, .keep_all = TRUE)

beta_mat <- beta_mat[, meth_primary$sample_barcode, drop = FALSE]
colnames(beta_mat) <- meth_primary$patient_barcode

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cpg_id") %>%
  janitor::clean_names()

promoter_groups <- c("TSS1500", "TSS200", "5'UTR", "1stExon")

xci_promoter_anno <- anno %>%
  filter(chr == "chrX") %>%
  filter(!is.na(ucsc_ref_gene_name), ucsc_ref_gene_name != "") %>%
  filter(!is.na(ucsc_ref_gene_group), ucsc_ref_gene_group != "") %>%
  tidyr::separate_rows(ucsc_ref_gene_name, sep = ";") %>%
  tidyr::separate_rows(ucsc_ref_gene_group, sep = ";") %>%
  filter(ucsc_ref_gene_group %in% promoter_groups) %>%
  mutate(
    gene_symbol_upper = toupper(ucsc_ref_gene_name)
  ) %>%
  filter(gene_symbol_upper %in% xci_genes$gene_symbol_upper) %>%
  distinct(cpg_id, gene_symbol_upper) %>%
  filter(cpg_id %in% rownames(beta_mat))

cat("\nNumber of XCI-linked promoter CpG-gene mappings:\n")
print(nrow(xci_promoter_anno))

if (nrow(xci_promoter_anno) == 0) {
  stop("No promoter CpGs found for XCI-linked genes.")
}

beta_xci <- beta_mat[unique(xci_promoter_anno$cpg_id), , drop = FALSE]

max_missing_fraction <- 0.20
keep_cpg <- rowMeans(is.na(beta_xci)) <= max_missing_fraction
beta_xci <- beta_xci[keep_cpg, , drop = FALSE]

# Impute remaining missing values by CpG median
for (i in seq_len(nrow(beta_xci))) {
  missing_idx <- is.na(beta_xci[i, ])
  if (any(missing_idx)) {
    beta_xci[i, missing_idx] <- median(beta_xci[i, ], na.rm = TRUE)
  }
}

anno_filtered <- xci_promoter_anno %>%
  filter(cpg_id %in% rownames(beta_xci))

meth_genes <- sort(unique(anno_filtered$gene_symbol_upper))

gene_beta_mat <- matrix(
  NA_real_,
  nrow = length(meth_genes),
  ncol = ncol(beta_xci),
  dimnames = list(meth_genes, colnames(beta_xci))
)

for (gene in meth_genes) {
  gene_cpgs <- anno_filtered %>%
    filter(gene_symbol_upper == gene) %>%
    pull(cpg_id) %>%
    unique()
  
  gene_cpgs <- intersect(gene_cpgs, rownames(beta_xci))
  
  if (length(gene_cpgs) == 1) {
    gene_beta_mat[gene, ] <- beta_xci[gene_cpgs, ]
  } else if (length(gene_cpgs) > 1) {
    gene_beta_mat[gene, ] <- colMeans(beta_xci[gene_cpgs, , drop = FALSE], na.rm = TRUE)
  }
}

methylation_patient_score <- tibble(
  patient_barcode = colnames(gene_beta_mat),
  xci_linked_promoter_methylation_mean_beta = colMeans(gene_beta_mat, na.rm = TRUE),
  xci_linked_promoter_methylation_median_beta = matrixStats::colMedians(gene_beta_mat, na.rm = TRUE),
  n_xci_methylation_genes_used = colSums(!is.na(gene_beta_mat))
) %>%
  mutate(
    methylation_erosion_score_raw = 1 - xci_linked_promoter_methylation_mean_beta
  )

cat("\nMethylation patient score dimension:\n")
print(dim(methylation_patient_score))

############################################################
# 9. download or load CNV data
############################################################

if (file.exists(cnv_rds_file)) {
  cat("\nLoading cached CNV segment data:\n")
  cat(cnv_rds_file, "\n")
  cnv_seg <- readRDS(cnv_rds_file)
} else {
  cat("\nDownloading TCGA-BRCA masked copy number segment data from GDC...\n")
  
  query_cnv <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Copy Number Variation",
    data.type = "Masked Copy Number Segment",
    sample.type = "Primary Tumor"
  )
  
  GDCdownload(
    query_cnv,
    method = "api",
    files.per.chunk = 50
  )
  
  cnv_seg <- GDCprepare(query_cnv)
  saveRDS(cnv_seg, cnv_rds_file)
}

############################################################
# 10. calculate XCI-linked gene CNV score
############################################################

cnv_seg <- as.data.frame(cnv_seg) %>%
  janitor::clean_names()

sample_col <- detect_sample_column(cnv_seg)

chrom_col <- c("chromosome", "chrom", "chr")
chrom_col <- chrom_col[chrom_col %in% colnames(cnv_seg)][1]

start_col <- c("start", "loc_start", "start_position")
start_col <- start_col[start_col %in% colnames(cnv_seg)][1]

end_col <- c("end", "loc_end", "end_position")
end_col <- end_col[end_col %in% colnames(cnv_seg)][1]

segment_mean_col <- c("segment_mean", "seg_mean", "segmentmean")
segment_mean_col <- segment_mean_col[segment_mean_col %in% colnames(cnv_seg)][1]

if (is.na(chrom_col) | is.na(start_col) | is.na(end_col) | is.na(segment_mean_col)) {
  stop("Cannot detect required CNV columns.")
}

cnv_clean <- cnv_seg %>%
  transmute(
    sample_barcode = as.character(.data[[sample_col]]),
    patient_barcode = if_else(
      stringr::str_detect(sample_barcode, "^TCGA-"),
      stringr::str_sub(sample_barcode, 1, 12),
      sample_barcode
    ),
    sample_type_code = if_else(
      stringr::str_detect(sample_barcode, "^TCGA-"),
      stringr::str_sub(sample_barcode, 14, 15),
      NA_character_
    ),
    chromosome = standardise_chr(.data[[chrom_col]]),
    start = as.numeric(.data[[start_col]]),
    end = as.numeric(.data[[end_col]]),
    segment_mean = as.numeric(.data[[segment_mean_col]])
  ) %>%
  filter(
    !is.na(patient_barcode),
    !is.na(chromosome),
    !is.na(start),
    !is.na(end),
    !is.na(segment_mean),
    end >= start
  )

if (any(cnv_clean$sample_type_code == "01", na.rm = TRUE)) {
  cnv_clean <- cnv_clean %>%
    filter(sample_type_code == "01")
}

cnv_keep <- cnv_clean %>%
  distinct(patient_barcode, sample_barcode) %>%
  arrange(patient_barcode, sample_barcode) %>%
  distinct(patient_barcode, .keep_all = TRUE)

cnv_clean <- cnv_clean %>%
  inner_join(cnv_keep, by = c("patient_barcode", "sample_barcode")) %>%
  filter(chromosome == "X")

gene_gr <- GenomicRanges::GRanges(
  seqnames = gene_coord$chromosome,
  ranges = IRanges::IRanges(
    start = gene_coord$start,
    end = gene_coord$end
  ),
  gene_symbol_upper = gene_coord$gene_symbol_upper
)

cnv_gr <- GenomicRanges::GRanges(
  seqnames = cnv_clean$chromosome,
  ranges = IRanges::IRanges(
    start = cnv_clean$start,
    end = cnv_clean$end
  ),
  patient_barcode = cnv_clean$patient_barcode,
  segment_mean = cnv_clean$segment_mean
)

hits <- GenomicRanges::findOverlaps(gene_gr, cnv_gr)

if (length(hits) == 0) {
  stop("No overlap found between XCI-linked genes and CNV segments.")
}

overlap_width <- width(
  GenomicRanges::pintersect(
    gene_gr[queryHits(hits)],
    cnv_gr[subjectHits(hits)]
  )
)

cnv_overlap <- tibble(
  patient_barcode = cnv_gr$patient_barcode[subjectHits(hits)],
  gene_symbol_upper = gene_gr$gene_symbol_upper[queryHits(hits)],
  segment_mean = cnv_gr$segment_mean[subjectHits(hits)],
  overlap_width = overlap_width
)

cnv_gene_score <- cnv_overlap %>%
  group_by(patient_barcode, gene_symbol_upper) %>%
  summarise(
    gene_cnv_segment_mean = weighted.mean(
      segment_mean,
      w = overlap_width,
      na.rm = TRUE
    ),
    gene_overlap_bp = sum(overlap_width, na.rm = TRUE),
    .groups = "drop"
  )

cnv_patient_score <- cnv_gene_score %>%
  group_by(patient_barcode) %>%
  summarise(
    xci_linked_cnv_weighted_segment_mean = weighted.mean(
      gene_cnv_segment_mean,
      w = gene_overlap_bp,
      na.rm = TRUE
    ),
    xci_linked_cnv_median_gene_segment_mean = median(
      gene_cnv_segment_mean,
      na.rm = TRUE
    ),
    xci_linked_cnv_mean_abs_gene_segment_mean = weighted.mean(
      abs(gene_cnv_segment_mean),
      w = gene_overlap_bp,
      na.rm = TRUE
    ),
    n_xci_cnv_genes_used = n_distinct(gene_symbol_upper),
    xci_cnv_total_overlap_bp = sum(gene_overlap_bp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cnv_erosion_score_raw = xci_linked_cnv_weighted_segment_mean,
    cnv_burden_score_raw = xci_linked_cnv_mean_abs_gene_segment_mean
  )

cat("\nCNV patient score dimension:\n")
print(dim(cnv_patient_score))

############################################################
# 11. integrate expression + methylation + CNV
############################################################

v4_data <- expression_patient_score %>%
  full_join(methylation_patient_score, by = "patient_barcode") %>%
  full_join(cnv_patient_score, by = "patient_barcode") %>%
  mutate(
    expression_erosion_z = z_score(expression_erosion_score_raw),
    methylation_erosion_z = z_score(methylation_erosion_score_raw),
    cnv_erosion_z = z_score(cnv_erosion_score_raw),
    
    xci_erosion_score_v4_xci_linked_gene_only =
      rowMeans(
        cbind(
          expression_erosion_z,
          methylation_erosion_z,
          cnv_erosion_z
        ),
        na.rm = FALSE
      ),
    score_gene_panel = "20 strict Balaton-consensus subject-to-XCI genes",
    score_interpretation = paste(
      "Composite proxy: higher subject-gene expression, lower promoter methylation,",
      "and higher copy-number segment mean; not direct allele-specific Xi reactivation"
    )
  )

v4_median <- median(
  v4_data$xci_erosion_score_v4_xci_linked_gene_only,
  na.rm = TRUE
)

v4_data <- v4_data %>%
  mutate(
    xci_erosion_v4_group = case_when(
      is.na(xci_erosion_score_v4_xci_linked_gene_only) ~ NA_character_,
      xci_erosion_score_v4_xci_linked_gene_only >= v4_median ~ "High erosion",
      xci_erosion_score_v4_xci_linked_gene_only < v4_median ~ "Low erosion"
    )
  )

cat("\nNumber of patients with all three V4 components:\n")
print(sum(
  !is.na(v4_data$expression_erosion_score_raw) &
    !is.na(v4_data$methylation_erosion_score_raw) &
    !is.na(v4_data$cnv_erosion_score_raw)
))

############################################################
# 12. add metadata if available
############################################################

possible_metadata_files <- c(
  "outputs/RQ1/rq1_analysis_dataset.csv",
  "outputs/RQ1/RQ1_analysis_dataset.csv",
  "outputs/RQ1/xist_stage_pam50_metadata.csv",
  "data/xist_stage_pam50_metadata.csv"
)

metadata_file <- possible_metadata_files[file.exists(possible_metadata_files)][1]

if (!is.na(metadata_file)) {
  
  metadata <- readr::read_csv(metadata_file, show_col_types = FALSE) %>%
    janitor::clean_names()
  
  if ("patient_id" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(patient_barcode = stringr::str_sub(patient_id, 1, 12))
  } else if ("patient_barcode" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(patient_barcode = stringr::str_sub(patient_barcode, 1, 12))
  } else if ("sample_barcode" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(patient_barcode = stringr::str_sub(sample_barcode, 1, 12))
  }
  
  metadata_keep_cols <- intersect(
    c(
      "patient_barcode",
      "pam50",
      "brca_subtype_pam50",
      "stage_clean",
      "stage4",
      "grade_clean",
      "grade_numeric",
      "age",
      "age_years",
      "sex"
    ),
    colnames(metadata)
  )
  
  metadata <- metadata %>%
    select(all_of(metadata_keep_cols)) %>%
    distinct(patient_barcode, .keep_all = TRUE)
  
  if (!"pam50" %in% colnames(metadata) &&
      "brca_subtype_pam50" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(pam50 = brca_subtype_pam50)
  }
  
  if (!"stage_clean" %in% colnames(metadata) &&
      "stage4" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(stage_clean = stage4)
  }
  
  if (!"age" %in% colnames(metadata) &&
      "age_years" %in% colnames(metadata)) {
    metadata <- metadata %>%
      mutate(age = age_years)
  }
  
  v4_data <- v4_data %>%
    left_join(metadata, by = "patient_barcode")
}

############################################################
# 13. save output CSVs
############################################################

readr::write_csv(v4_data, v4_output_file)

v4_long <- v4_data %>%
  select(
    patient_barcode,
    expression_erosion_z,
    methylation_erosion_z,
    cnv_erosion_z,
    xci_erosion_score_v4_xci_linked_gene_only,
    any_of(c("pam50", "stage_clean", "grade_clean", "grade_numeric"))
  ) %>%
  pivot_longer(
    cols = c(
      expression_erosion_z,
      methylation_erosion_z,
      cnv_erosion_z,
      xci_erosion_score_v4_xci_linked_gene_only
    ),
    names_to = "score_component",
    values_to = "score_value"
  )

readr::write_csv(v4_long, v4_long_file)

cat("\nSaved V4 main CSV to:\n")
cat(v4_output_file, "\n")

cat("\nSaved V4 long-format CSV to:\n")
cat(v4_long_file, "\n")

cat("\nV4 output columns:\n")
print(colnames(v4_data))

cat("\nDone.\n")

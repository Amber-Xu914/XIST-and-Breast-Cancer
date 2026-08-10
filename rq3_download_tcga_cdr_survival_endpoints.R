############################################################
## Download TCGA-CDR survival endpoint table
## Endpoints: OS, DSS, DFI, PFI
############################################################

library(readr)
library(dplyr)
library(stringr)

out_dir <- "data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cdr_url <- "https://pancanatlas.xenahubs.net/download/Survival_SupplementalTable_S1_20171025_xena_sp"

cdr <- read_tsv(cdr_url, show_col_types = FALSE)

cat("TCGA-CDR loaded:\n")
cat("Dimension:", nrow(cdr), "rows x", ncol(cdr), "columns\n")
print(colnames(cdr))

brca_cdr <- cdr %>%
  filter(`cancer type abbreviation` == "BRCA") %>%
  transmute(
    patient_barcode = substr(`_PATIENT`, 1, 12),
    
    age_years = as.numeric(age_at_initial_pathologic_diagnosis),
    
    OS_time = as.numeric(`OS.time`),
    OS_event = as.integer(OS),
    
    DSS_time = as.numeric(`DSS.time`),
    DSS_event = as.integer(DSS),
    
    DFI_time = as.numeric(`DFI.time`),
    DFI_event = as.integer(DFI),
    
    PFI_time = as.numeric(`PFI.time`),
    PFI_event = as.integer(PFI)
  ) %>%
  distinct(patient_barcode, .keep_all = TRUE)

write_csv(
  brca_cdr,
  "data/tcga_brca_cdr_survival_endpoints.csv"
)

cat("\nBRCA TCGA-CDR survival endpoint file saved:\n")
cat("data/tcga_brca_cdr_survival_endpoints.csv\n\n")

cat("Endpoint non-missing summary:\n")
endpoint_summary <- brca_cdr %>%
  summarise(
    n = n(),
    
    OS_n = sum(!is.na(OS_time) & !is.na(OS_event)),
    OS_events = sum(OS_event == 1, na.rm = TRUE),
    
    DSS_n = sum(!is.na(DSS_time) & !is.na(DSS_event)),
    DSS_events = sum(DSS_event == 1, na.rm = TRUE),
    
    DFI_n = sum(!is.na(DFI_time) & !is.na(DFI_event)),
    DFI_events = sum(DFI_event == 1, na.rm = TRUE),
    
    PFI_n = sum(!is.na(PFI_time) & !is.na(PFI_event)),
    PFI_events = sum(PFI_event == 1, na.rm = TRUE)
  )

print(endpoint_summary)

write_csv(
  endpoint_summary,
  "data/tcga_brca_cdr_survival_endpoint_summary.csv"
)
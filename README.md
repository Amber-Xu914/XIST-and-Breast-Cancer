# XIST and Breast Cancer — final analysis scripts

This folder contains the latest script set used for the final thesis analyses as of 10 August 2026. Earlier V1–V3 score scripts, backup files, binary high/low and tertile RQ2 variants, and superseded Aim 3/PFI scripts have been excluded.

## Cohort assumption

The scripts assume that the clinical and molecular inputs have already been restricted to patients recorded as female in the TCGA clinical metadata. Male patients are excluded upstream because the project investigates X-chromosome inactivation in the female biological context.

## Folder structure

```text
scripts/
  00_data_acquisition/       TCGA-CDR survival endpoint download
  01_xist_and_os/            XIST expression and overall-survival analyses
  02_v4_score/               Final 20-gene V4 score construction and validation
  03_genomic_instability/    Final continuous and quartile RQ2 analysis
  04_tp53_and_pfi/           TP53, interaction and PFI analyses
  05_final_figures/          Final thesis figure regeneration
resources/
  xci_linked_genes.csv
  xci_subject_gene_panel_audit.csv
```

## Recommended run order

1. Place the required TCGA-derived inputs in `data/`.
2. Run `scripts/00_data_acquisition/rq3_download_tcga_cdr_survival_endpoints.R` if the TCGA-CDR endpoint table is not already available.
3. Run `scripts/02_v4_score/rq1_version4_xci_linked.R`.
4. Run `scripts/02_v4_score/rq1_plot_v4_xci_linked_score.R`.
5. Run the XIST/OS scripts in `scripts/01_xist_and_os/` as required.
6. Run `scripts/03_genomic_instability/rq2_v4_quartile_analysis_script.R`.
7. Run the scripts in `scripts/04_tp53_and_pfi/`. The `rerun_pfi_restricted_final_analyses.R` script is the final combined workflow for the PFI-restricted TP53 interaction and PFI analyses.
8. Run `scripts/05_final_figures/regenerate_updated_subject_xci_figures.R` after the upstream result tables are available.

Scripts that accept a project-directory argument default to the current working directory. For example:

```sh
Rscript scripts/02_v4_score/rq1_version4_xci_linked.R /path/to/project
```

## Data policy

Raw and large derived TCGA files are intentionally not included. The `.gitignore` file prevents local data, results, R session files and generated figures from being uploaded accidentally. The two small files in `resources/` are included because they define and audit the final 20-gene subject-to-XCI panel.

## Version scope

- Final score: V4 subject-to-XCI composite proxy score.
- Final RQ2 grouping: quartiles, with the continuous score as the primary exposure.
- Final survival endpoint for the V4 score: progression-free interval.
- Final Aim 3 interaction workflow: `rerun_pfi_restricted_final_analyses.R`.

See `MANIFEST.csv` for the retained files and their roles.

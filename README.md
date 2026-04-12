# Hepatocellular-carcinoma-HCC-bulk-rna-seq-DGE

Bulk RNA-seq differential gene expression analysis workflow for hepatocellular carcinoma (HCC).

## Included in this repository
- `r_script/` : R analysis scripts
- `count_matrix_metadata/` : metadata files
- `results_HCC/02_tables/` : DEG, enrichment, GSEA, and WGCNA result tables
- `results_HCC/03_plots/` : PNG and SVG output plots

## Excluded from Git
- raw count TSV files
- `.rds` object files
- TIFF plots
- oversized normalized expression tables

## Main script
- `r_script/liver_lihc_final.R`

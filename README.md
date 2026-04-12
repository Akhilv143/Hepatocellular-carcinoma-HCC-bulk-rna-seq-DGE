# Hepatocellular Carcinoma (HCC) — Bulk RNA-seq Differential Gene Expression Analysis

## Project Overview

This repository contains a comprehensive, multi-cohort bulk RNA-seq bioinformatics pipeline for Hepatocellular Carcinoma (HCC). By integrating four independent GEO datasets aligned to GRCh38.p13, the pipeline models inter-study batch effects within the DESeq2 design formula to identify robust transcriptomic signatures. The workflow covers raw count matrix integration, DESeq2-based differential expression, protein-coding gene filtering, functional enrichment (GO and KEGG), Gene Set Enrichment Analysis (GSEA), and Weighted Gene Co-expression Network Analysis (WGCNA).

**Main Analysis Script:** [`r_script/liver_lihc_final.R`](r_script/liver_lihc_final.R)

---

## Datasets Analyzed

All four datasets were aligned to the **GRCh38.p13 NCBI** reference genome and provided as Entrez gene ID-indexed raw count matrices. Merged by intersecting common genes across all cohorts.

| GEO Accession | Condition | Platform |
|---|---|---|
| [GSE77314](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77314) | Normal liver vs HCC | RNA-seq GRCh38.p13 |
| [GSE124535](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124535) | Normal liver vs HCC | RNA-seq GRCh38.p13 |
| [GSE138485](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138485) | Normal liver vs HCC | RNA-seq GRCh38.p13 |
| [GSE144269](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144269) | Normal liver vs HCC | RNA-seq GRCh38.p13 |

- **Condition labels:** `Control` = normal liver | `Diseased` = HCC tumor
- **Batch variable:** GEO Study ID (GSE accession) — modelled as covariate in DESeq2 design
- **Sample metadata:** [`count_matrix_metadata/HCC_metadata.csv`](count_matrix_metadata/HCC_metadata.csv)

---

## Analytical Pipeline

### 1. Count Matrix Integration
Raw count matrices were loaded with `data.table::fread()`. Common genes across all four datasets were identified using `Reduce(intersect, ...)`, matrices subset and column-bound into a single merged matrix. Rows with zero total counts were removed. Saved as `HCC_merged_raw_counts.tsv`.

### 2. DESeq2 Differential Expression
A `DESeqDataSet` was constructed with design formula `~ batch + condition`. Low-count genes were removed (≥10 counts in ≥10 samples). DESeq2 was run with the negative binomial GLM. Results extracted for `Diseased vs Control` at `alpha = 0.05`. Entrez IDs mapped to HGNC symbols via `org.Hs.eg.db`.

- **Up-regulated:** `padj < 0.05` AND `log2FC > 1`
- **Down-regulated:** `padj < 0.05` AND `log2FC < -1`

### 3. Protein-coding Gene Filter
Significant DEGs were filtered to retain only protein-coding genes by excluding: uncharacterized loci (`LOC*`), miRNAs (`MIR*`), lincRNAs (`LINC*`), mitochondrial genes (`MT-`), snoRNAs (`SNORD*`, `SNORA*`), snRNAs (`RNU*`, `RN7SL`), snoRNA host genes (`SNHG*`, `SCARNA*`), and known lncRNAs (`MALAT1`, `NEAT1`, `H19`, `XIST`, `MEG*`, `PEG*`).

### 4. VST Normalization
Variance Stabilizing Transformation (`vst()`) applied to the full filtered DESeq2 dataset. VST matrix used for PCA, heatmap, and WGCNA. Entrez IDs mapped to gene symbols and saved as `VST_normalized_expression_symbols.csv`.

### 5. Heatmap of Top DEGs
Top 25 upregulated and top 25 downregulated protein-coding DEGs (by adjusted p-value) visualized using `ComplexHeatmap`. VST values row-scaled to Z-scores. Samples ordered Control → Diseased. Column clustering within each condition group. Annotations: condition, GEO batch, regulation direction. Color scale: green → black → red (Z-score −2 to +2).

### 6. Functional Enrichment (ORA)
Over-representation analysis using `clusterProfiler::enrichGO()` (BP, CC, MF) and `enrichKEGG()`. BH correction, `pvalueCutoff = 0.05`, `qvalueCutoff = 0.05`. Separate analyses for Up and Down gene sets saved as `*_with_direction.csv`. Top 10 GO terms per sub-ontology shown as faceted bar plot; top 20 KEGG pathways as dot plot.

### 7. Gene Set Enrichment Analysis (GSEA)
`gseGO()` and `gseKEGG()` run on a pre-ranked gene list sorted by `log2FoldChange`. BH correction at `pvalueCutoff = 0.05`. Positive NES = enriched in HCC tumors; negative NES = enriched in normal liver. Results saved per sub-ontology.

### 8. WGCNA
Performed on the top 5,000 most variable genes from the symbol-indexed VST matrix after `goodSamplesGenes()` QC. Signed co-expression network built with `blockwiseModules()`.

| Parameter | Value |
|---|---|
| Network type | Signed |
| Soft threshold power | 8 (R² ≥ 0.80) |
| Merge cut height | 0.25 |
| Max block size | 10,000 |
| Random seed | 1234 |

Module eigengenes correlated against binary traits (Diseased / Control). Gene Significance (GS) vs Module Membership (MM) scatter plots generated for the top 6 HCC-correlated modules with marginal histograms via `ggExtra::ggMarginal()`.

---

## Visualizations

### 1. Differential Gene Expression

The MA plot shows log2 fold change against mean expression. The Volcano plot labels the top 10 upregulated and top 10 downregulated protein-coding genes using `ggrepel`. Thresholds: `padj < 0.05`, `|log2FC| > 1`.

<p align="center">
  <img src="results_HCC/03_plots/DEG/MA_plot.png" width="48%" alt="MA Plot">
  <img src="results_HCC/03_plots/DEG/Volcano_plot.png" width="48%" alt="Volcano Plot">
</p>

---

### 2. Heatmap of Top 50 DEGs

Hierarchical clustering heatmap of the top 25 upregulated and top 25 downregulated protein-coding DEGs. Rows are Z-score scaled, columns split by condition (Control | Diseased) and annotated by GEO dataset batch.

<p align="center">
  <img src="results_HCC/03_plots/Heatmap/Heatmap_top25up_top25dn.png" width="800" alt="DEG Heatmap">
</p>

---

### 3. Functional Enrichment — GO and KEGG

Over-representation analysis for GO sub-ontologies (top 10 terms each, faceted bar plot) and KEGG pathways (top 20 terms, dot plot sized by gene count, colored by adjusted p-value).

<p align="center">
  <img src="results_HCC/03_plots/enrichment/png/GO_BP.png" width="48%" alt="GO Biological Process">
  <img src="results_HCC/03_plots/enrichment/png/GO_CC.png" width="48%" alt="GO Cellular Component">
</p>
<p align="center">
  <img src="results_HCC/03_plots/enrichment/png/GO_MF.png" width="48%" alt="GO Molecular Function">
  <img src="results_HCC/03_plots/enrichment/png/KEGG.png" width="48%" alt="KEGG Pathways">
</p>

---

### 4. WGCNA — Weighted Gene Co-expression Network Analysis

#### Scale-Free Topology and Mean Connectivity
Soft threshold power selection plots. The red dashed line marks R² = 0.80. Power = 8 was selected as it achieves scale-free topology while maintaining adequate mean connectivity.

> WGCNA plots are saved as high-resolution TIFF files (`results_HCC/03_plots/WGCNA/`). Descriptions below correspond to each output file.

| Plot File | Description |
|---|---|
| `1_Scale_Free_Topology.tiff` | Scale-free topology model fit (R²) vs soft threshold power. Red line at R² = 0.80 marks the selection threshold. |
| `2_Mean_Connectivity.tiff` | Mean network connectivity vs soft threshold power. Used alongside Plot 1 to select optimal power = 8. |
| `3_Cluster_Dendrogram.tiff` | Hierarchical clustering dendrogram of the top 5,000 variable genes with assigned module colors shown below the dendrogram. |
| `4_Module_Correlation_Heatmap.tiff` | Inter-module correlation heatmap. Color scale: red (distance = 0, highly correlated) → white → blue (distance = 1.5, anti-correlated). Hierarchically clustered. |
| `5_Module_Trait_Heatmap.tiff` | Module eigengene–trait correlation heatmap. Each cell shows Pearson r and formatted p-value. Rows = modules, Columns = Diseased / Control traits. |
| `6_GS_MM_Scatter_blue.tiff` | Gene Significance vs Module Membership scatter for the **blue** module. Pearson r and p-value annotated. Marginal histograms shown. |
| `6_GS_MM_Scatter_brown.tiff` | Gene Significance vs Module Membership scatter for the **brown** module. |
| `6_GS_MM_Scatter_green.tiff` | Gene Significance vs Module Membership scatter for the **green** module. |
| `6_GS_MM_Scatter_magenta.tiff` | Gene Significance vs Module Membership scatter for the **magenta** module. |
| `6_GS_MM_Scatter_pink.tiff` | Gene Significance vs Module Membership scatter for the **pink** module. |
| `6_GS_MM_Scatter_yellow.tiff` | Gene Significance vs Module Membership scatter for the **yellow** module. |

---

## Repository Structure

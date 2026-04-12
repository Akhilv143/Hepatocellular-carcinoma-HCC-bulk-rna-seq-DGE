

# step1_Libraries 
library(data.table)
library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(BiocParallel)
library(grid)

register(SerialParam())
options(mc.cores = 1)

# step2_Working directory & folder structure 
setwd("C:/Users/akhil/comet_downlaod/liver_cancer_m/data")

dir.create("results_HCC/01_objects",           recursive = TRUE, showWarnings = FALSE)
dir.create("results_HCC/02_tables/DEG",        recursive = TRUE, showWarnings = FALSE)
dir.create("results_HCC/02_tables/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("results_HCC/03_plots/QC_PCA",      recursive = TRUE, showWarnings = FALSE)
dir.create("results_HCC/03_plots/DEG",         recursive = TRUE, showWarnings = FALSE)
dir.create("results_HCC/03_plots/Heatmap",     recursive = TRUE, showWarnings = FALSE)
dir.create("results_HCC/03_plots/enrichment",  recursive = TRUE, showWarnings = FALSE)

# step3_Load metadata ---------------------------------------------------------
# Condition column: "Control" = normal liver | "Diseased" = HCC tumor
meta           <- read.csv("HCC_metadata.csv", stringsAsFactors = FALSE)
meta$condition <- factor(meta$Condition, levels = c("Control", "Diseased"))
meta$batch     <- factor(meta$GSE)
rownames(meta) <- meta$GSM_ID

message("Metadata: ", nrow(meta), " samples | ",
        sum(meta$condition == "Control"),  " Control | ",
        sum(meta$condition == "Diseased"), " Diseased")
print(table(meta$batch, meta$condition))

# step4_Load count matrices ---------------------------------------------------
load_counts <- function(file) {
  df           <- fread(file, header = TRUE) %>% as.data.frame()
  rownames(df) <- df[[1]]
  df[[1]]      <- NULL
  as.matrix(df)
}

mat_77314  <- load_counts("GSE77314_raw_counts_GRCh38.p13_NCBI.tsv")
mat_124535 <- load_counts("GSE124535_raw_counts_GRCh38.p13_NCBI.tsv")
mat_138485 <- load_counts("GSE138485_raw_counts_GRCh38.p13_NCBI.tsv")
mat_144269 <- load_counts("GSE144269_raw_counts_GRCh38.p13_NCBI.tsv")

# step5_Intersect genes, verify order, merge ----------------------------------
common_genes <- Reduce(intersect, list(
  rownames(mat_77314), rownames(mat_124535),
  rownames(mat_138485), rownames(mat_144269)
))
message("Common genes across all 4 datasets: ", length(common_genes))

mat_77314  <- mat_77314[common_genes, ]
mat_124535 <- mat_124535[common_genes, ]
mat_138485 <- mat_138485[common_genes, ]
mat_144269 <- mat_144269[common_genes, ]

stopifnot(
  all(rownames(mat_77314) == rownames(mat_124535)),
  all(rownames(mat_77314) == rownames(mat_138485)),
  all(rownames(mat_77314) == rownames(mat_144269))
)

raw_counts <- cbind(mat_77314, mat_124535, mat_138485, mat_144269)
raw_counts <- raw_counts[rowSums(raw_counts) > 0, ]
message("Merged matrix: ", nrow(raw_counts), " genes x ", ncol(raw_counts), " samples")

# step6_Save merged raw count matrix 
fwrite(
  as.data.frame(raw_counts) %>% rownames_to_column("GeneID"),
  "HCC_merged_raw_counts.tsv",
  sep = "\t"
)

# step7_Sanity checks 
stopifnot(all(colnames(raw_counts) %in% rownames(meta)))
stopifnot(all(rownames(meta) %in% colnames(raw_counts)))
meta <- meta[colnames(raw_counts), , drop = FALSE]
stopifnot(!any(is.na(meta$condition)))

# step8_DESeq2 object & filter 
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData   = meta,
  design    = ~ batch + condition
)

keep <- rowSums(counts(dds) >= 10) >= 10
dds  <- dds[keep, ]
message("Genes after pre-filtering: ", nrow(dds))

# step9_Run DESeq2 
dds <- DESeq(dds, parallel = FALSE)
summary(results(dds, contrast = c("condition", "Diseased", "Control"), alpha = 0.05))

# step10_Extract results, map gene symbols 
res <- results(dds,
               contrast = c("condition", "Diseased", "Control"),
               alpha    = 0.05)

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  mutate(
    symbol = unname(mapIds(
      org.Hs.eg.db,
      keys      = gene,
      keytype   = "ENTREZID",
      column    = "SYMBOL",
      multiVals = "first"
    )),
    regulation = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >  1 ~ "Up",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(padj)

message("Up: ",   sum(res_df$regulation == "Up",   na.rm = TRUE),
        " | Down: ", sum(res_df$regulation == "Down", na.rm = TRUE))

write.csv(res_df,
          "results_HCC/02_tables/DEG/DESeq2_all_genes_with_symbols.csv",
          row.names = FALSE)

# step11_Non-coding filter 
noncoding_patterns <- c(
  "^LOC[0-9]", "^MIR[0-9]", "^LINC[0-9]", "^MT-",
  "^SNORD[0-9]", "^SNORA[0-9]", "^RNU[0-9]", "^RN7SL",
  "^SNHG[0-9]", "^SCARNA[0-9]", "^MALAT", "^NEAT",
  "^H19$", "^XIST$", "^MEG[0-9]", "^PEG[0-9]"
)

is_noncoding <- function(symbols) {
  pattern       <- paste(noncoding_patterns, collapse = "|")
  result        <- rep(FALSE, length(symbols))
  valid         <- !is.na(symbols)
  result[valid] <- grepl(pattern, symbols[valid], ignore.case = FALSE)
  result
}

res_df$is_noncoding <- is_noncoding(res_df$symbol)

sig_deg <- res_df %>%
  filter(regulation != "NS") %>%
  arrange(padj)

sig_deg_protein <- sig_deg %>%
  filter(!is_noncoding, !is.na(symbol))

message("Protein-coding DEGs: ", nrow(sig_deg_protein),
        " | Up: ",   sum(sig_deg_protein$regulation == "Up"),
        " | Down: ", sum(sig_deg_protein$regulation == "Down"))

write.csv(sig_deg,
          "results_HCC/02_tables/DEG/DESeq2_significant_DEGs_all.csv",
          row.names = FALSE)
write.csv(sig_deg_protein,
          "results_HCC/02_tables/DEG/DESeq2_significant_DEGs_proteinCoding.csv",
          row.names = FALSE)

 
#  CREATE VST WITH GENE SYMBOLS
library(org.Hs.eg.db)
library(tidyverse)

setwd("C:/Users/akhil/comet_downlaod/liver_cancer_m/data")

#  Load saved VST matrix 
vst_data <- readRDS("results_HCC/01_objects/HCC_vst_data.rds")
message("VST loaded: ", nrow(vst_data), " genes x ", ncol(vst_data), " samples")

# Map Entrez to Gene Symbol 
entrez_to_symbol <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(vst_data),
  keytype   = "ENTREZID",
  column    = "SYMBOL",
  multiVals = "first"
)

message("Mapped: ",   sum(!is.na(entrez_to_symbol)), " genes")
message("Unmapped: ", sum(is.na(entrez_to_symbol)),  " genes (will be dropped)")

#  Build symbol VST dataframe 
vst_symbol <- as.data.frame(vst_data) %>%
  rownames_to_column("entrez_id") %>%
  mutate(gene_id = unname(entrez_to_symbol[entrez_id])) %>%
  filter(!is.na(gene_id), gene_id != "") %>%   # drop unmapped
  distinct(gene_id, .keep_all = TRUE) %>%      # drop duplicate symbols (keep first)
  select(gene_id, everything(), -entrez_id)    # gene_id as column 1

# Save 
out_path <- "results_HCC/02_tables/DEG/VST_normalized_expression_geneSymbol.csv"

write.csv(vst_symbol, out_path, row.names = FALSE)

message("✔ Saved: ", out_path)
message("  Genes  : ", nrow(vst_symbol))
message("  Samples: ", ncol(vst_symbol) - 1)
message("\nPreview (first 5 genes, first 4 columns):")
print(vst_symbol[1:5, 1:4])


# step13_Save R objects 
saveRDS(dds,             "results_HCC/01_objects/HCC_dds.rds")
saveRDS(res_df,          "results_HCC/01_objects/HCC_res_df.rds")
saveRDS(sig_deg_protein, "results_HCC/01_objects/HCC_sig_deg_protein.rds")
saveRDS(vst_mat,         "results_HCC/01_objects/HCC_vst_mat.rds")
saveRDS(vst_data,        "results_HCC/01_objects/HCC_vst_data.rds")

# step14_Plot helper: PNG + TIFF 600 DPI, no compression 
save_plot <- function(plot_obj, base_path, w, h) {
  ggsave(paste0(base_path, ".png"),
         plot_obj, width = w, height = h, dpi = 600)
  ggsave(paste0(base_path, ".tiff"),
         plot_obj, width = w, height = h, dpi = 600,
         compression = "none")
}

pal <- c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey75")

# step15_MA plot 
ma_plot <- ggplot(res_df,
                  aes(x = log10(baseMean),
                      y = log2FoldChange,
                      colour = regulation)) +
  geom_point(alpha = 0.4, size = 1.2, stroke = 0) +
  geom_hline(yintercept = c(-1, 1),
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0,
             linetype = "solid", colour = "black", linewidth = 0.3) +
  scale_colour_manual(name = "Regulation", values = pal) +
  labs(x = "Mean expression (log\u2081\u2080)",
       y = "Log\u2082 fold change  (Diseased vs Control)") +
  theme_bw(base_size = 11) +
  theme(
    axis.title            = element_text(face = "bold"),
    axis.text             = element_text(face = "bold"),
    legend.title          = element_text(face = "bold"),
    legend.text           = element_text(face = "bold"),
    legend.position       = "right",
    legend.justification  = "top",
    legend.background     = element_rect(colour = "black", linewidth = 0.4),
    panel.grid.minor      = element_blank()
  )

save_plot(ma_plot, "results_HCC/03_plots/DEG/MA_plot", 8, 6)

# step16_Volcano plot 
top10_up <- res_df %>%
  filter(regulation == "Up", !is.na(symbol), !is_noncoding) %>%
  arrange(padj) %>%
  slice_head(n = 10)

top10_dn <- res_df %>%
  filter(regulation == "Down", !is.na(symbol), !is_noncoding) %>%
  arrange(padj) %>%
  slice_head(n = 10)

top20_labeled <- bind_rows(top10_up, top10_dn)

message("Volcano labels — Up: ",   paste(top10_up$symbol, collapse = ", "))
message("Volcano labels — Down: ", paste(top10_dn$symbol, collapse = ", "))

volcano_plot <- ggplot(res_df,
                       aes(x      = log2FoldChange,
                           y      = -log10(padj),
                           colour = regulation)) +
  geom_point(alpha = 0.5, size = 1.5, stroke = 0) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "black", linewidth = 0.5) +
  ggrepel::geom_text_repel(
    data          = top20_labeled,
    aes(label     = symbol),
    size          = 3.5,
    fontface      = "bold.italic",
    max.overlaps  = 20,
    box.padding   = 0.5,
    point.padding = 0.3,
    show.legend   = FALSE
  ) +
  scale_colour_manual(
    name   = "Regulation",
    values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey70"),
    labels = c("Up" = "Up-regulated", "Down" = "Down-regulated", "NS" = "NS")
  ) +
  labs(
    x = expression(bold(Log[2]~Fold~Change~(Diseased~vs~Control))),
    y = expression(bold(-log[10](italic(p)[adj])))
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text             = element_text(face = "bold"),
    axis.title            = element_text(face = "bold"),
    legend.title          = element_text(face = "bold"),
    legend.text           = element_text(face = "bold"),
    legend.position       = "right",
    legend.justification  = "top",
    legend.key.size       = unit(0.5, "cm"),
    legend.background     = element_rect(colour = "black",
                                         linewidth = 0.5, fill = "white"),
    legend.box.background = element_rect(colour = "black", linewidth = 0.5),
    panel.grid.minor      = element_blank(),
    plot.title            = element_blank()
  )

ggsave("results_HCC/03_plots/DEG/Volcano_plot.png",
       volcano_plot, width = 8, height = 7, dpi = 600)

ggsave("results_HCC/03_plots/DEG/Volcano_plot.tiff",
       volcano_plot, width = 8, height = 7, dpi = 600,
       compression = "none")

message("Volcano plot saved.")



# step17_PCA plots 
pca_df     <- plotPCA(vst_mat, intgroup = c("condition", "batch"),
                      returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

pca_cond <- ggplot(pca_df,
                   aes(PC1, PC2, colour = condition, fill = condition)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.12, colour = NA) +
  stat_ellipse(level = 0.95, linewidth = 1.1) +
  scale_colour_manual(name   = "Condition",
                      values = c("Control" = "#2166AC", "Diseased" = "#D6604D")) +
  scale_fill_manual(name     = "Condition",
                    values   = c("Control" = "#2166AC", "Diseased" = "#D6604D")) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 13) +
  theme(
    axis.title            = element_text(face = "bold"),
    axis.text             = element_text(face = "bold"),
    legend.title          = element_text(face = "bold"),
    legend.text           = element_text(face = "bold"),
    legend.position       = "right",
    legend.background     = element_rect(colour = "black", linewidth = 0.5,
                                         fill = "white"),
    legend.box.background = element_rect(colour = "black", linewidth = 0.5),
    panel.grid.minor      = element_blank()
  )

save_plot(pca_cond, "results_HCC/03_plots/QC_PCA/PCA_by_Condition", 10, 7)

pca_batch <- ggplot(pca_df,
                    aes(PC1, PC2, colour = batch)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(x      = paste0("PC1: ", percentVar[1], "% variance"),
       y      = paste0("PC2: ", percentVar[2], "% variance"),
       colour = "GEO Study") +
  theme_bw(base_size = 13) +
  theme(
    axis.title            = element_text(face = "bold"),
    axis.text             = element_text(face = "bold"),
    legend.title          = element_text(face = "bold"),
    legend.text           = element_text(face = "bold"),
    legend.position       = "right",
    legend.background     = element_rect(colour = "black", linewidth = 0.5,
                                         fill = "white"),
    legend.box.background = element_rect(colour = "black", linewidth = 0.5),
    panel.grid.minor      = element_blank()
  )

save_plot(pca_batch, "results_HCC/03_plots/QC_PCA/PCA_by_Batch", 10, 7)



library(DESeq2)
library(tidyverse)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

setwd("C:/Users/akhil/comet_downlaod/liver_cancer_m/data")

# Load all saved RDS objects
dds      <- readRDS("results_HCC/01_objects/HCC_dds.rds")
vst_data <- readRDS("results_HCC/01_objects/HCC_vst_data.rds")
vst_mat  <- readRDS("results_HCC/01_objects/HCC_vst_mat.rds")

# Load metadata
meta           <- read.csv("HCC_metadata.csv", stringsAsFactors = FALSE)
meta$condition <- factor(meta$Condition, levels = c("Control", "Diseased"))
meta$batch     <- factor(meta$GSE)
rownames(meta) <- meta$GSM_ID

# Verify all loaded correctly
message("dds loaded     : ", class(dds))
message("vst_data dims  : ", nrow(vst_data), " x ", ncol(vst_data))
message("meta samples   : ", nrow(meta))
message("Conditions     : ", paste(levels(meta$condition), collapse = ", "))
message("Batches        : ", paste(levels(meta$batch),     collapse = ", "))
message("✅ All objects loaded — ready to run heatmap script")





# Prerequisites assumed in environment: dds, vst_mat, vst_data, meta 

library(tidyverse)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# Step A: Rebuild res_df with all required columns 
res_raw <- results(dds,
                   contrast = c("condition", "Diseased", "Control"),
                   alpha    = 0.05)

res_df <- as.data.frame(res_raw) %>%
  rownames_to_column("gene") %>%
  mutate(
    symbol = unname(mapIds(
      org.Hs.eg.db,
      keys      = gene,
      keytype   = "ENTREZID",
      column    = "SYMBOL",
      multiVals = "first"
    )),
    regulation = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange >  1 ~ "Up",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(padj)

# Step B: Add is_noncoding flag 
noncoding_patterns <- c(
  "^LOC[0-9]", "^MIR[0-9]", "^LINC[0-9]", "^MT-",
  "^SNORD[0-9]", "^SNORA[0-9]", "^RNU[0-9]", "^RN7SL",
  "^SNHG[0-9]", "^SCARNA[0-9]", "^MALAT", "^NEAT",
  "^H19$", "^XIST$", "^MEG[0-9]", "^PEG[0-9]"
)
is_noncoding <- function(symbols) {
  pattern       <- paste(noncoding_patterns, collapse = "|")
  result        <- rep(FALSE, length(symbols))
  valid         <- !is.na(symbols)
  result[valid] <- grepl(pattern, symbols[valid])
  result
}
res_df$is_noncoding <- is_noncoding(res_df$symbol)

message("res_df columns: ", paste(colnames(res_df), collapse = ", "))
message("Up: ",   sum(res_df$regulation == "Up",   na.rm = TRUE),
        " | Down: ", sum(res_df$regulation == "Down", na.rm = TRUE))

# Step C: Select top 25 Up + top 25 Down 
top25_up <- res_df %>%
  filter(!is.na(padj), !is.na(symbol), !is_noncoding, regulation == "Up") %>%
  arrange(padj) %>% slice_head(n = 25)

top25_dn <- res_df %>%
  filter(!is.na(padj), !is.na(symbol), !is_noncoding, regulation == "Down") %>%
  arrange(padj) %>% slice_head(n = 25)

top_degs <- bind_rows(top25_dn, top25_up)

message("Heatmap genes — Up: ", nrow(top25_up),
        " | Down: ", nrow(top25_dn),
        " | Total: ", nrow(top_degs))

# Step D: Build heatmap matrix 
in_vst   <- top_degs$gene %in% rownames(vst_data)
top_filt <- top_degs[in_vst, ]

stopifnot("No genes found in vst_data" = nrow(top_filt) > 0)

col_order <- c(
  rownames(meta)[meta$condition == "Control"],
  rownames(meta)[meta$condition == "Diseased"]
)

hmat           <- vst_data[top_filt$gene, col_order, drop = FALSE]
rownames(hmat) <- top_filt$symbol
hmat_z         <- t(scale(t(hmat)))
valid          <- apply(hmat_z, 1,
                        function(x) !all(is.na(x)) & !any(is.infinite(x)))
hmat_z         <- hmat_z[valid, , drop = FALSE]

reg_vec <- top_filt$regulation[match(rownames(hmat_z), top_filt$symbol)]
row_ord <- c(which(reg_vec == "Down"), which(reg_vec == "Up"))
hmat_z  <- hmat_z[row_ord, ]
reg_vec <- reg_vec[row_ord]

col_split <- factor(meta[col_order, "condition"],
                    levels = c("Control", "Diseased"))

# Step E: Colours 
condition_colors <- c("Control"  = "#4575B4",
                      "Diseased" = "#FFD700")

n_batch      <- length(levels(meta$batch))
batch_colors <- setNames(
  brewer.pal(max(3, n_batch), "Set2")[1:n_batch],
  levels(meta$batch)
)

col_fun <- colorRamp2(c(-2, 0, 2), c("green3", "black", "red"))

# Step F: Annotations 
ann_df <- data.frame(
  Condition = meta[col_order, "condition"],
  Dataset   = meta[col_order, "batch"],
  row.names = col_order
)

ha_top <- HeatmapAnnotation(
  Condition = ann_df$Condition,
  Dataset   = ann_df$Dataset,
  col = list(
    Condition = condition_colors,
    Dataset   = batch_colors
  ),
  annotation_name_gp   = gpar(fontface = "bold", fontsize = 10),
  show_annotation_name = TRUE,
  show_legend          = FALSE
)

# ── CHANGED: regulation colors only 
ha_left <- rowAnnotation(
  Regulation = reg_vec,
  col = list(
    Regulation = c(
      "Up"   = "#8BC34A",   # lime green
      "Down" = "#CE93D8"    # lavender
    )
  ),
  width                = unit(0.4, "cm"),
  annotation_name_gp   = gpar(fontface = "bold", fontsize = 10),
  annotation_name_side = "top",
  show_annotation_name = TRUE,
  show_legend          = FALSE
)

# Step G: Heatmap object 
ht <- Heatmap(
  hmat_z,
  name                  = "Z-score",
  col                   = col_fun,
  top_annotation        = ha_top,
  left_annotation       = ha_left,
  cluster_rows          = FALSE,
  cluster_columns       = TRUE,
  cluster_column_slices = FALSE,
  column_split          = col_split,
  column_gap            = unit(0, "mm"),
  rect_gp               = gpar(col = NA),
  show_row_names        = TRUE,
  show_column_names     = FALSE,
  row_names_side        = "right",
  row_names_gp          = gpar(fontface = "bold.italic", fontsize = 9),
  show_heatmap_legend   = FALSE,
  row_split             = factor(reg_vec, levels = c("Down", "Up")),
  row_gap               = unit(0, "mm"),
  row_title             = NULL,
  column_title          = NULL,
  border                = TRUE
)

# Step H: Manual legends -
lgd_condition <- Legend(
  title     = "Condition",
  at        = c("Control", "Diseased"),
  legend_gp = gpar(fill = c("#4575B4", "#FFD700")),
  title_gp  = gpar(fontface = "bold", fontsize = 10),
  labels_gp = gpar(fontface = "bold", fontsize = 9)
)


lgd_regulation <- Legend(
  title     = "Regulation",
  at        = c("Down", "Up"),
  labels    = c("Down-regulated", "Up-regulated"),
  legend_gp = gpar(fill = c("#CE93D8", "#8BC34A")),  # Down=lavender, Up=green
  title_gp  = gpar(fontface = "bold", fontsize = 10),
  labels_gp = gpar(fontface = "bold", fontsize = 9)
)

lgd_dataset <- Legend(
  title     = "Dataset",
  at        = levels(meta$batch),
  legend_gp = gpar(fill = batch_colors),
  title_gp  = gpar(fontface = "bold", fontsize = 10),
  labels_gp = gpar(fontface = "bold", fontsize = 9)
)

lgd_zscore <- Legend(
  title         = "Z-score",
  col_fun       = col_fun,
  at            = c(-2, 0, 2),
  labels        = c("-2", "0", "2"),
  legend_height = unit(3.5, "cm"),
  title_gp      = gpar(fontface = "bold", fontsize = 10),
  labels_gp     = gpar(fontface = "bold", fontsize = 9),
  direction     = "vertical"
)

all_legends <- packLegend(
  lgd_condition,
  lgd_regulation,
  lgd_dataset,
  lgd_zscore,
  direction = "vertical",
  gap       = unit(4, "mm")
)


draw_ht <- function() {
  draw(
    ht,
    annotation_legend_list  = all_legends,
    heatmap_legend_side     = "right",
    annotation_legend_side  = "right",
    align_annotation_legend = "heatmap_top",
    merge_legend            = FALSE,
    padding                 = unit(c(2, 2, 2, 2), "mm")
  )
}

dir.create("results_HCC/03_plots/Heatmap", recursive = TRUE, showWarnings = FALSE)

png("results_HCC/03_plots/Heatmap/Heatmap_top25up_top25dn.png",
    width = 14, height = 11, units = "in", res = 600)
draw_ht(); dev.off()

tiff("results_HCC/03_plots/Heatmap/Heatmap_top25up_top25dn.tiff",
     width = 14, height = 11, units = "in", res = 600,
     compression = "none")
draw_ht(); dev.off()

message("Heatmap saved to results_HCC/03_plots/Heatmap/")



# step19_GO and KEGG enrichment 
entrez_input <- unique(as.character(
  sig_deg_protein$gene[!is.na(sig_deg_protein$gene)]
))

wrap_label <- function(x, w = 45) {
  sapply(x, function(s) paste(strwrap(s, w), collapse = "\n"),
         USE.NAMES = FALSE)
}

ego_bp <- tryCatch(
  enrichGO(gene = entrez_input, OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID", ont = "BP",
           pAdjustMethod = "BH", pvalueCutoff = 0.05,
           qvalueCutoff  = 0.05, readable = TRUE),
  error = function(e) NULL)

ego_cc <- tryCatch(
  enrichGO(gene = entrez_input, OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID", ont = "CC",
           pAdjustMethod = "BH", pvalueCutoff = 0.05,
           qvalueCutoff  = 0.05, readable = TRUE),
  error = function(e) NULL)

ego_mf <- tryCatch(
  enrichGO(gene = entrez_input, OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID", ont = "MF",
           pAdjustMethod = "BH", pvalueCutoff = 0.05,
           qvalueCutoff  = 0.05, readable = TRUE),
  error = function(e) NULL)

ekegg <- tryCatch(
  enrichKEGG(gene = entrez_input, organism = "hsa",
             pvalueCutoff = 0.05, pAdjustMethod = "BH"),
  error = function(e) NULL)

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0)
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

if (!is.null(ego_bp)) write.csv(as.data.frame(ego_bp),
                                "results_HCC/02_tables/enrichment/GO_BP_results.csv", row.names = FALSE)
if (!is.null(ego_cc)) write.csv(as.data.frame(ego_cc),
                                "results_HCC/02_tables/enrichment/GO_CC_results.csv", row.names = FALSE)
if (!is.null(ego_mf)) write.csv(as.data.frame(ego_mf),
                                "results_HCC/02_tables/enrichment/GO_MF_results.csv", row.names = FALSE)
if (!is.null(ekegg))  write.csv(as.data.frame(ekegg),
                                "results_HCC/02_tables/enrichment/KEGG_results.csv",  row.names = FALSE)

saveRDS(list(BP = ego_bp, CC = ego_cc, MF = ego_mf, KEGG = ekegg),
        "results_HCC/01_objects/GO_KEGG_objects.rds")

# step20_GO combined barplot 
go_rows <- bind_rows(
  if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0)
    as.data.frame(ego_bp) %>% arrange(p.adjust) %>%
    slice_head(n = 10) %>% mutate(Category = "Biological Process"),
  if (!is.null(ego_cc) && nrow(as.data.frame(ego_cc)) > 0)
    as.data.frame(ego_cc) %>% arrange(p.adjust) %>%
    slice_head(n = 10) %>% mutate(Category = "Cellular Component"),
  if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0)
    as.data.frame(ego_mf) %>% arrange(p.adjust) %>%
    slice_head(n = 10) %>% mutate(Category = "Molecular Function")
)

if (nrow(go_rows) > 0) {
  go_plot_df <- go_rows %>%
    mutate(
      label    = wrap_label(Description, 45),
      Category = factor(Category,
                        levels = c("Biological Process",
                                   "Cellular Component",
                                   "Molecular Function"))
    ) %>%
    arrange(Category, Count) %>%
    mutate(label = factor(label, levels = label))
  
  p_go <- ggplot(go_plot_df,
                 aes(x = Count, y = label, fill = p.adjust)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradient(low = "#D62728", high = "#AEC7E8",
                        name = "Adjusted\np-value") +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Gene count", y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      axis.text          = element_text(face = "bold", size = 9),
      axis.title.x       = element_text(face = "bold", size = 11),
      strip.text         = element_text(face = "bold", size = 10,
                                        colour = "white"),
      strip.background   = element_rect(fill = "#404040"),
      legend.title       = element_text(face = "bold", size = 9),
      legend.text        = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "right",
      legend.background  = element_rect(colour = "black", linewidth = 0.4)
    )
  
  save_plot(p_go, "results_HCC/03_plots/enrichment/GO_combined_barplot", 11, 13)
}

# step21_KEGG dotplot 
if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  kegg_df <- as.data.frame(ekegg) %>%
    arrange(p.adjust) %>% slice_head(n = 20) %>%
    mutate(label = wrap_label(Description, 45)) %>%
    arrange(Count) %>%
    mutate(label = factor(label, levels = label))
  
  p_kegg <- ggplot(kegg_df,
                   aes(x = Count, y = label)) +
    geom_point(aes(size = Count, colour = p.adjust), alpha = 0.85) +
    scale_colour_gradient(low = "#D62728", high = "#AEC7E8",
                          name = "Adjusted\np-value") +
    scale_size_continuous(name = "Gene\ncount", range = c(3, 11)) +
    labs(x = "Gene count", y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      axis.text          = element_text(face = "bold", size = 9),
      axis.title.x       = element_text(face = "bold", size = 11),
      legend.title       = element_text(face = "bold", size = 9),
      legend.text        = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "right",
      legend.background  = element_rect(colour = "black", linewidth = 0.4)
    )
  
  save_plot(p_kegg, "results_HCC/03_plots/enrichment/KEGG_dotplot", 11, 9)
}




# HCC WGCNA ANALYSIS 

library(WGCNA)
library(ggplot2)
library(ggExtra)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

setwd("C:/Users/akhil/comet_downlaod/liver_cancer_m/data")

# Create necessary directories
dir.create("results_HCC/01_objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results_HCC/02_tables/WGCNA", recursive = TRUE, showWarnings = FALSE)
dir.create("results_HCC/03_plots/WGCNA", recursive = TRUE, showWarnings = FALSE)


#  Load metadata 
meta           <- read.csv("HCC_metadata.csv", stringsAsFactors = FALSE)
meta$condition <- factor(meta$Condition, levels = c("Control", "Diseased"))
meta$batch     <- factor(meta$GSE)
rownames(meta) <- meta$GSM_ID

message("Metadata loaded: ", nrow(meta), " samples")


# STEP 1: Prepare expression matrix from symbol VST

vst_symbol <- read.csv(
  "results_HCC/02_tables/DEG/VST_normalized_expression_symbols.csv",
  row.names   = 1,
  check.names = FALSE
)

datExpr_raw <- as.matrix(vst_symbol)
datExpr     <- t(datExpr_raw)

common_samples <- intersect(rownames(datExpr), rownames(meta))
datExpr        <- datExpr[common_samples, ]
meta_wgcna     <- meta[common_samples, ]

message("Matched samples : ", nrow(datExpr))
message("Genes for WGCNA : ", ncol(datExpr))





# STEP 2: Filter low-variance genes (top 5000 most variable)


gene_vars  <- apply(datExpr, 2, var)
top_genes  <- names(sort(gene_vars, decreasing = TRUE))[1:5000]
datExpr    <- datExpr[, top_genes]
message("Genes after variance filter: ", ncol(datExpr))

# Goodness of sample check
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  message("Removed bad samples/genes. Remaining: ",
          nrow(datExpr), " x ", ncol(datExpr))
}


# PLOT 1: Scale-Free Topology — choose soft threshold

powers <- c(1:10, seq(12, 50, by = 2))
sft     <- pickSoftThreshold(datExpr,
                             powerVector  = powers,
                             networkType  = "signed",
                             verbose      = 5)
sft_df  <- sft$fitIndices

p1 <- ggplot(sft_df, aes(Power, SFT.R.sq, label = Power)) +
  geom_point(size = 2.5, color = "black") +
  geom_text(vjust = -0.8, size = 3.5, fontface = "bold") +
  geom_hline(yintercept = 0.8, color = "red", linetype = "dashed",
             linewidth = 0.7) +
  labs(x = "Soft Threshold (Power)",
       y = "Scale-Free Topology Model Fit (signed R²)") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text  = element_text(face = "bold", size = 10))

tiff("results_HCC/03_plots/WGCNA/1_Scale_Free_Topology.tiff",
     width = 6, height = 5, units = "in", res = 600, compression = "none")
print(p1); dev.off()
message("Plot 1 saved: Scale-Free Topology")


# PLOT 2: Mean Connectivity


p2 <- ggplot(sft_df, aes(Power, mean.k., label = Power)) +
  geom_point(size = 2.5, color = "black") +
  geom_text(vjust = -0.8, size = 3.5, fontface = "bold") +
  labs(x = "Soft Threshold (Power)", y = "Mean Connectivity") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text  = element_text(face = "bold", size = 10))

tiff("results_HCC/03_plots/WGCNA/2_Mean_Connectivity.tiff",
     width = 6, height = 5, units = "in", res = 600, compression = "none")
print(p2); dev.off()
message("Plot 2 saved: Mean Connectivity")


# Choose soft power — set based on previous topology analysis
soft_power <- 8


# STEP 3: Build network — blockwise modules

temp_cor <- cor
cor      <- WGCNA::cor

bwnet <- blockwiseModules(
  datExpr,
  maxBlockSize   = 10000,
  TOMType        = "signed",
  power          = soft_power,
  mergeCutHeight = 0.25,
  numericLabels  = FALSE,
  randomSeed     = 1234,
  verbose        = 3
)

cor <- temp_cor                 # restore base cor

moduleColors <- bwnet$colors
MEs          <- bwnet$MEs

message("Modules found: ", length(unique(moduleColors)))
message("Module sizes:")
print(table(moduleColors))

# Save module assignments
module_df <- data.frame(gene = colnames(datExpr),
                        module = moduleColors)
write.csv(module_df,
          "results_HCC/02_tables/WGCNA/Module_Gene_Assignments.csv",
          row.names = FALSE)


# PLOT 3: Hierarchical Clustering Dendrogram

tiff("results_HCC/03_plots/WGCNA/3_Cluster_Dendrogram.tiff",
     width = 12, height = 6, units = "in", res = 600, compression = "none")
par(font.main = 2, font.lab = 2, font.axis = 2,
    cex.main = 1.2, cex.lab = 1.1)
plotDendroAndColors(
  bwnet$dendrograms[[1]],
  moduleColors[bwnet$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  
)
dev.off()
message("Plot 3 saved: Cluster Dendrogram")


# PLOT 4: Module–Module Correlation Heatmap


module_cor  <- cor(MEs, use = "p")
module_dist <- 1 - module_cor
rownames(module_dist) <- gsub("ME", "", rownames(module_dist))
colnames(module_dist) <- gsub("ME", "", colnames(module_dist))

col_fun_mod <- colorRamp2(c(0, 0.75, 1.5),
                          c("red", "white", "blue"))

tiff("results_HCC/03_plots/WGCNA/4_Module_Correlation_Heatmap.tiff",
     width = 8, height = 8, units = "in", res = 600, compression = "none")
ht_mod <- Heatmap(
  module_dist,
  name                        = "Distance",
  col                         = col_fun_mod,
  cluster_rows                = TRUE,
  cluster_columns             = TRUE,
  clustering_distance_rows    = as.dist(module_dist),
  clustering_distance_columns = as.dist(module_dist),
  show_row_names              = TRUE,
  show_column_names           = TRUE,
  row_names_gp                = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp             = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    title      = "Distance",
    title_gp   = gpar(fontsize = 11, fontface = "bold"),
    labels_gp  = gpar(fontface = "bold"),
    direction  = "vertical"
  ),
  border  = TRUE,
  rect_gp = gpar(col = "grey90", lwd = 0.5)
)
draw(ht_mod); dev.off()
message("Plot 4 saved: Module Correlation Heatmap")



# PLOT 5: Module–Trait Heatmap — HCC


# Step 1: Traits matrix
traits <- data.frame(
  Diseased = as.numeric(meta_wgcna$condition == "Diseased"),
  Control  = as.numeric(meta_wgcna$condition == "Control"),
  row.names = rownames(meta_wgcna)
)

# Step 2: Correlations and p-values
moduleTraitCor    <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Step 3: Clean rownames
rownames(moduleTraitCor)    <- gsub("ME", "", rownames(moduleTraitCor))
rownames(moduleTraitPvalue) <- gsub("ME", "", rownames(moduleTraitPvalue))

# Step 4: Order by Diseased correlation (high → low)
ord               <- order(moduleTraitCor[, "Diseased"], decreasing = TRUE)
moduleTraitCor    <- moduleTraitCor[ord, ]
moduleTraitPvalue <- moduleTraitPvalue[ord, ]

# Step 5: Row annotation — actual module colors as colored bar
mod_colors <- rownames(moduleTraitCor)
row_anno <- rowAnnotation(
  Module = mod_colors,
  col    = list(Module = setNames(mod_colors, mod_colors)),
  show_annotation_name = FALSE,
  show_legend          = FALSE,
  simple_anno_size     = unit(0.5, "cm")
)

# Step 6: Cell text — "0.86\n(5e-51)" compact format
cell_text <- matrix("", nrow(moduleTraitCor), ncol(moduleTraitCor))
for (i in seq_len(nrow(moduleTraitCor))) {
  for (j in seq_len(ncol(moduleTraitCor))) {
    r_val <- round(moduleTraitCor[i, j], 2)
    p_raw <- moduleTraitPvalue[i, j]
    p_fmt <- formatC(p_raw, format = "e", digits = 0)
    p_fmt <- sub("e\\+0*(\\d+)", "e+\\1", sub("e-0*(\\d+)", "e-\\1", p_fmt))
    cell_text[i, j] <- paste0(r_val, "\n(", p_fmt, ")")
  }
}

# Step 7: Color — green → white → red (matching reference)
col_fun_trait <- colorRamp2(c(-1, 0, 1),
                            c("#33A02C", "white", "#E31A1C"))

# Step 8: Draw and save
dir.create("results_HCC/03_plots/WGCNA", recursive = TRUE, showWarnings = FALSE)

tiff("results_HCC/03_plots/WGCNA/5_Module_Trait_Heatmap.tiff",
     width = 5, height = 10, units = "in", res = 600, compression = "none")

ht <- Heatmap(
  moduleTraitCor,
  name              = "Correlation",
  col               = col_fun_trait,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  row_names_side    = "left",
  left_annotation   = row_anno,
  column_names_side = "bottom",
  column_names_rot  = 0,
  row_names_gp      = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp   = gpar(fontsize = 12, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(cell_text[i, j], x, y,
              gp = gpar(fontsize = 8, fontface = "bold"))
  },
  heatmap_legend_param = list(
    title      = "Correlation",
    title_gp   = gpar(fontsize = 11, fontface = "bold"),
    labels_gp  = gpar(fontface = "bold"),
    direction  = "vertical",
    at         = c(-1, -0.5, 0, 0.5, 1),
    labels     = c("-1", "-0.5", "0", "0.5", "1")
  ),
  border  = TRUE,
  rect_gp = gpar(col = "white", lwd = 1.5)
)

draw(ht,
     heatmap_legend_side = "right",
     column_title        = "A",
     column_title_gp     = gpar(fontsize = 16, fontface = "bold"))

dev.off()
message("✅ Plot 5 saved: Module-Trait Heatmap")



# PLOT 6: GS vs MM Scatter Plots — top 6 Tumor-correlated modules


geneTraitSignificance <- as.data.frame(
  cor(datExpr, traits$Tumor, use = "p")
)
names(geneTraitSignificance)    <- "GS.Trait"
rownames(geneTraitSignificance) <- colnames(datExpr)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
rownames(geneModuleMembership)  <- colnames(datExpr)

# Pick top 6 modules by absolute Tumor correlation
top6_modules <- rownames(moduleTraitCor)[
  order(abs(moduleTraitCor[, "Tumor"]), decreasing = TRUE)
][1:6]

message("\nGenerating GS vs MM scatter plots for: ",
        paste(top6_modules, collapse = ", "))

for (mod in top6_modules) {
  
  me_col <- paste0("ME", mod)
  if (!me_col %in% colnames(MEs)) {
    message("Module not found: ", mod, " — skipping")
    next
  }
  
  inModule   <- moduleColors == mod
  col_idx    <- match(me_col, colnames(MEs))
  MM_values  <- geneModuleMembership[inModule, col_idx]
  GS_values  <- geneTraitSignificance$GS.Trait[inModule]
  
  plot_data  <- data.frame(MM = MM_values, GS = GS_values) %>%
    filter(complete.cases(.))
  
  if (nrow(plot_data) < 3) next
  
  r_val <- cor(plot_data$MM, plot_data$GS)
  p_val <- cor.test(plot_data$MM, plot_data$GS)$p.value
  
  p_base <- ggplot(plot_data, aes(x = MM, y = GS)) +
    geom_point(color = "navy", alpha = 0.6, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE,
                color = "red", linewidth = 0.7) +
    annotate("text",
             x     = min(plot_data$MM),
             y     = max(plot_data$GS),
             label = paste0("r = ", round(r_val, 2),
                            "\np = ", format(p_val,
                                             scientific = TRUE,
                                             digits     = 2)),
             hjust = 0, vjust = 1,
             size  = 3.5, fontface = "bold") +
    labs(x = paste0("Module Membership (", mod, " module)"),
         y = "Gene Significance for Tumor") +
    theme_classic() +
    theme(
      axis.title   = element_text(size = 10, face = "bold"),
      axis.text    = element_text(size = 9,  face = "bold"),
      panel.border = element_rect(color = "black", fill = NA,
                                  linewidth = 0.8)
    )
  
  p_final <- ggMarginal(
    p_base, type = "histogram", bins = 25, size = 8,
    xparams = list(fill = "orange", color = NA, alpha = 0.7),
    yparams = list(fill = "red",    color = NA, alpha = 0.7)
  )
  
  tiff(paste0("results_HCC/03_plots/WGCNA/6_GS_MM_Scatter_", mod, ".tiff"),
       width = 6, height = 6, units = "in", res = 600,
       compression = "none")
  print(p_final); dev.off()
  
  message("Plot 6 saved: GS-MM Scatter — ", mod, " module")
}

# ============================================================
# Save WGCNA objects for downstream use
# ============================================================

saveRDS(list(bwnet              = bwnet,
             moduleColors       = moduleColors,
             MEs                = MEs,
             moduleTraitCor     = moduleTraitCor,
             moduleTraitPvalue  = moduleTraitPvalue,
             geneModuleMembership   = geneModuleMembership,
             geneTraitSignificance  = geneTraitSignificance),
        "results_HCC/01_objects/WGCNA_objects.rds")







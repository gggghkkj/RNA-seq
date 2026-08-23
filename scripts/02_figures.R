#!/usr/bin/env Rscript
# =============================================================================
# 02_figures.R
# Generate figures for fPSCs transcriptome analysis
# Reference: "Formative pluripotent stem cells show features of epiblast cells
#             poised for gastrulation" (Cell Research, 2021)
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# --- 1. Load data -------------------------------------------------------------
cat("=== Loading data ===\n")

dds <- readRDS("data/processed/dds_object.rds")
vsd <- readRDS("data/processed/vsd_object.rds")

sample_info <- data.frame(
  cell_type = factor(c("EpiSCs", "EpiSCs", "fPSCs", "fPSCs", "fPSCs", "fPSCs", "mESCs", "mESCs"),
                     levels = c("mESCs", "EpiSCs", "fPSCs")),
  row.names = colnames(dds)
)

cell_colors <- c("mESCs" = "#E64B35", "EpiSCs" = "#4DBBD5", "fPSCs" = "#00A087")

# --- 2. Figure 4a: PCA plot ---------------------------------------------------
cat("\n=== Generating Figure 4a: PCA ===\n")

pca_data <- plotPCA(vsd, intgroup = "cell_type", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_type, shape = cell_type)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = cell_colors) +
  scale_shape_manual(values = c("mESCs" = 16, "EpiSCs" = 17, "fPSCs" = 15)) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA - Naive mESCs, EpiSCs, and fPSCs") +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/figures/Figure4a_PCA.pdf", p_pca, width = 7, height = 5)
ggsave("results/figures/Figure4a_PCA.png", p_pca, width = 7, height = 5, dpi = 300)

# --- 3. Figure 4b: Heatmap of top variable genes -----------------------------
cat("\n=== Generating Figure 4b: Heatmap ===\n")

mat <- assay(vsd)
gene_vars <- rowVars(mat)
top_genes <- head(order(gene_vars, decreasing = TRUE), 50)
mat_top <- mat[top_genes, ]
mat_scaled <- t(scale(t(mat_top)))

annotation_col <- data.frame(
  Cell_Type = sample_info$cell_type,
  row.names = rownames(sample_info)
)
ann_colors <- list(
  Cell_Type = c("mESCs" = "#E64B35", "EpiSCs" = "#4DBBD5", "fPSCs" = "#00A087")
)

pdf("results/figures/Figure4b_Heatmap.pdf", width = 8, height = 10)
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = FALSE, show_colnames = TRUE,
  main = "Top 50 Variable Genes",
  fontsize_col = 10
)
dev.off()

png("results/figures/Figure4b_Heatmap.png", width = 800, height = 1000, res = 150)
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = FALSE, show_colnames = TRUE,
  main = "Top 50 Variable Genes",
  fontsize_col = 10
)
dev.off()

# --- 4. Figure 4c: Marker gene expression ------------------------------------
cat("\n=== Generating Figure 4c: Marker genes ===\n")

# Gene name mapping
gtf_lines <- readLines("data/raw/Mus_musculus.GRCm39.113.gtf")
gtf_gene <- gtf_lines[grepl("^.*gene_id.*gene_name.*$", gtf_lines)]
gene_id <- gsub('.*gene_id "([^"]+)".*', "\\1", gtf_gene)
gene_name <- gsub('.*gene_name "([^"]+)".*', "\\1", gtf_gene)
gene_name_map <- data.frame(gene_id = gene_id, gene_name = gene_name,
                            stringsAsFactors = FALSE)
gene_name_map <- gene_name_map[!duplicated(gene_name_map$gene_id), ]

# Marker genes
marker_genes <- list(
  Naive = c("Nanog", "Klf4", "Esrrb", "Tbx3", "Dppa3", "Tfcp2l1"),
  Formative = c("Otx2", "Dnmt3b", "Sox4", "Fgf5", "Zic5", "Utf1"),
  Primed = c("Nodal", "Wnt3", "Eomes", "T", "Gata6", "Foxa2")
)

# Calculate FPKM
counts_norm <- counts(dds, normalized = TRUE)
gene_lengths <- counts_raw[rownames(counts_norm), "Length"]
fpkm <- (counts_norm * 10^9) / (rowSums(counts_norm) * gene_lengths)

# Collect marker data
all_data <- data.frame()
for (category in names(marker_genes)) {
  genes <- marker_genes[[category]]
  for (gene in genes) {
    ensembl_ids <- gene_name_map$gene_id[gene_name_map$gene_name == gene]
    if (length(ensembl_ids) > 0) {
      for (eid in ensembl_ids) {
        if (eid %in% rownames(fpkm)) {
          vals <- as.numeric(fpkm[eid, ])
          temp <- data.frame(
            sample = colnames(fpkm), gene = gene, category = category,
            fpkm = vals, cell_type = sample_info$cell_type
          )
          all_data <- rbind(all_data, temp)
          break
        }
      }
    }
  }
}

# Bar plot (Figure 4c)
if (nrow(all_data) > 0) {
  all_data$gene <- factor(all_data$gene, levels = unlist(marker_genes[levels(all_data$category)]))
  all_data$sample <- factor(all_data$sample,
    levels = c("mESCs_rep1", "mESCs_rep2", "EpiSCs_rep1", "EpiSCs_rep2",
               "fPSCs_P1", "fPSCs_P10", "fPSCs_P20", "fPSCs_P30"))

  p_bar <- ggplot(all_data, aes(x = sample, y = fpkm, fill = cell_type)) +
    geom_bar(stat = "identity", width = 0.7) +
    facet_grid(category ~ gene, scales = "free", space = "free_x") +
    scale_fill_manual(values = cell_colors) +
    labs(x = "", y = "FPKM", title = "Marker Gene Expression") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 7),
          strip.text.x = element_text(face = "italic", size = 9),
          strip.text.y = element_text(face = "bold", size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_blank(),
          panel.spacing = unit(1.5, "lines"))

  ggsave("results/figures/Figure4c_Markers.pdf", p_bar, width = 20, height = 8)
  ggsave("results/figures/Figure4c_Markers.png", p_bar, width = 20, height = 8, dpi = 300)

  # Dot plot (paper style)
  dot_data <- data.frame()
  for (g in unique(all_data$gene)) {
    for (ct in unique(all_data$cell_type)) {
      sub_df <- all_data[all_data$gene == g & all_data$cell_type == ct, ]
      if (nrow(sub_df) > 0) {
        dot_data <- rbind(dot_data, data.frame(
          gene = g, category = sub_df$category[1], cell_type = ct,
          fpkm = mean(sub_df$fpkm), stringsAsFactors = FALSE
        ))
      }
    }
  }

  p_dot <- ggplot(dot_data, aes(x = cell_type, y = gene, size = fpkm, color = cell_type)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(1, 8), name = "FPKM") +
    scale_color_manual(values = cell_colors, name = "Cell Type") +
    labs(x = "", y = "Gene", title = "Gene Expression of Selected Markers") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold"),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major = element_line(color = "grey90"))

  ggsave("results/figures/Figure4c_Markers_dotplot.pdf", p_dot, width = 9, height = 8)
  ggsave("results/figures/Figure4c_Markers_dotplot.png", p_dot, width = 9, height = 8, dpi = 300)
}

cat("\nAll figures saved to results/figures/\n")

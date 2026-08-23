#!/usr/bin/env Rscript
# =============================================================================
# 01_deseq2_analysis.R
# RNA-seq Differential Expression Analysis using DESeq2
# Reference: "Formative pluripotent stem cells show features of epiblast cells
#             poised for gastrulation" (Cell Research, 2021)
# Data: mESCs, EpiSCs, fPSCs (P1, P10, P20, P30)
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

# --- 1. Load count data ------------------------------------------------------
cat("=== Loading count data ===\n")

counts_raw <- read.delim("data/processed/gene_counts.txt",
                         comment.char = "#", row.names = 1)

# Remove annotation columns (Chr, Start, End, Strand, Length)
counts <- counts_raw[, 6:ncol(counts_raw)]
colnames(counts) <- gsub("aligned/|\\.sorted\\.bam", "", colnames(counts))

cat("Count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

# --- 2. Sample metadata ------------------------------------------------------
sample_info <- data.frame(
  sample = colnames(counts),
  cell_type = factor(c("EpiSCs", "EpiSCs", "fPSCs", "fPSCs", "fPSCs", "fPSCs", "mESCs", "mESCs"),
                     levels = c("mESCs", "EpiSCs", "fPSCs")),
  row.names = colnames(counts)
)

cat("\nSample metadata:\n")
print(sample_info)

# --- 3. DESeq2 analysis ------------------------------------------------------
cat("\n=== Running DESeq2 ===\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ cell_type
)

# Pre-filtering: remove genes with very low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# Run DESeq2
dds <- DESeq(dds)

# --- 4. Variance stabilizing transformation ----------------------------------
cat("\n=== VST transformation ===\n")
vsd <- vst(dds, blind = TRUE)

# --- 5. Differential expression - pairwise comparisons -----------------------
cat("\n=== Differential expression results ===\n")

res_epifp <- results(dds, contrast = c("cell_type", "fPSCs", "EpiSCs"))
res_escfp <- results(dds, contrast = c("cell_type", "fPSCs", "mESCs"))
res_escepi <- results(dds, contrast = c("cell_type", "EpiSCs", "mESCs"))

# Save DEG lists
write.csv(as.data.frame(res_epifp[order(res_epifp$padj), ]),
          "results/tables/DEGs_fPSCs_vs_EpiSCs.csv", row.names = TRUE)
write.csv(as.data.frame(res_escfp[order(res_escfp$padj), ]),
          "results/tables/DEGs_fPSCs_vs_mESCs.csv", row.names = TRUE)
write.csv(as.data.frame(res_escepi[order(res_escepi$padj), ]),
          "results/tables/DEGs_EpiSCs_vs_mESCs.csv", row.names = TRUE)

# --- 6. Calculate FPKM -------------------------------------------------------
cat("\n=== Calculating FPKM ===\n")

counts_norm <- counts(dds, normalized = TRUE)
gene_lengths <- counts_raw[rownames(counts_norm), "Length"]
fpkm <- (counts_norm * 10^9) / (rowSums(counts_norm) * gene_lengths)

# Gene name mapping from GTF
gtf_lines <- readLines("data/raw/Mus_musculus.GRCm39.113.gtf")
gtf_gene <- gtf_lines[grepl("^.*gene_id.*gene_name.*$", gtf_lines)]
gene_id <- gsub('.*gene_id "([^"]+)".*', "\\1", gtf_gene)
gene_name <- gsub('.*gene_name "([^"]+)".*', "\\1", gtf_gene)
gene_name_map <- data.frame(gene_id = gene_id, gene_name = gene_name,
                            stringsAsFactors = FALSE)
gene_name_map <- gene_name_map[!duplicated(gene_name_map$gene_id), ]

fpkm_export <- data.frame(
  gene_id = rownames(fpkm),
  gene_name = gene_name_map$gene_name[match(rownames(fpkm), gene_name_map$gene_id)],
  fpkm,
  stringsAsFactors = FALSE
)
write.csv(fpkm_export, "data/processed/FPKM_matrix.csv", row.names = FALSE)

# --- 7. Summary ---------------------------------------------------------------
cat("\n=== Summary ===\n")
cat("Significant DEGs (padj < 0.05, |log2FC| > 1):\n")
cat("  fPSCs vs EpiSCs:", sum(res_epifp$padj < 0.05 & abs(res_epifp$log2FoldChange) > 1, na.rm = TRUE), "\n")
cat("  fPSCs vs mESCs:", sum(res_escfp$padj < 0.05 & abs(res_escfp$log2FoldChange) > 1, na.rm = TRUE), "\n")
cat("  EpiSCs vs mESCs:", sum(res_escepi$padj < 0.05 & abs(res_escepi$log2FoldChange) > 1, na.rm = TRUE), "\n")

# Save RDS for downstream analysis
saveRDS(dds, "data/processed/dds_object.rds")
saveRDS(vsd, "data/processed/vsd_object.rds")

cat("\nAnalysis complete! Results saved to results/tables/\n")

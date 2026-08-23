# Reproducing Figure 4: Formative Pluripotent Stem Cells Transcriptome Analysis

## Project Overview

This project reproduces the core transcriptome analysis (Figure 4a-c) from the paper *"Formative pluripotent stem cells show features of epiblast cells poised for gastrulation"* (Cell Research, 2021). Using RNA-seq data from mouse embryonic stem cells (mESCs), epiblast stem cells (EpiSCs), and formative pluripotent stem cells (fPSCs), we perform differential expression analysis, PCA, and visualization to characterize the transcriptomic signatures of these three pluripotent states.

## Analysis Pipeline

```
Raw Data (FASTQ) 
    ↓
Quality Control (FastQC + MultiQC)
    ↓
Adapter Trimming (fastp)
    ↓
Reference Index Building (HISAT2)
    ↓
Read Alignment (HISAT2 → samtools)
    ↓
Gene Quantification (featureCounts)
    ↓
Differential Expression Analysis (DESeq2)
    ↓
Visualization (PCA, Heatmap, Dot plot)
```

Step-by-step:

1. Quality Control — FastQC assesses raw read quality; MultiQC generates summary reports
2. Adapter Trimming — fastp removes adapters and filters low-quality reads (Q20, length ≥ 50bp)
3. Genome Indexing — HISAT2 builds genome index for GRCm39 (mm10)
4. Read Alignment — HISAT2 aligns paired-end reads to mouse genome (≥95% alignment rate)
5. Gene Quantification — featureCounts counts reads per gene using Ensembl GTF annotation
6. Differential Expression — DESeq2 performs normalization (VST), pairwise comparisons, and identifies DEGs
7. Visualization — PCA plot (Fig 4a), heatmap of top variable genes (Fig 4b), marker gene expression (Fig 4c)

## Environment & Dependencies

| Software | Version | Purpose |
|----------|---------|---------|
| HISAT2 | 2.2.3 | Read alignment |
| samtools | 1.24 | BAM file processing |
| featureCounts | 2.1.1 | Gene quantification |
| fastp | 1.3.6 | Adapter trimming |
| FastQC | 0.12.1 | Quality control |
| MultiQC | 1.35 | QC report aggregation |
| R | 4.5.3 | Statistical analysis |
| DESeq2 | 1.50.2 | Differential expression |
| ggplot2 | 4.0.3 | Visualization |
| pheatmap | 1.0.13 | Heatmap generation |

Conda environment setup:

```bash
# Create conda environment
mamba create -n rnaseq -c bioconda -c conda-forge \
    hisat2 samtools subread fastp fastqc multiqc r-base

# Install R packages
mamba install -n rnaseq -c bioconda -c conda-forge \
    bioconductor-deseq2 bioconductor-annotationdbi \
    r-ggplot2 r-pheatmap r-rcolorbrewer r-circlize

# Activate environment
conda activate rnaseq
```

## How to Run

### Prerequisites

1. Download raw data (fastq.gz) to `data/raw/`
2. Download reference genome and annotation to `data/raw/`:
   - `Mus_musculus.GRCm39.dna.primary_assembly.fa.gz` from [Ensembl](https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/)
   - `Mus_musculus.GRCm39.113.gtf.gz` from [Ensembl](https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/)
3. Decompress reference files:
   ```bash
   gunzip data/raw/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
   gunzip data/raw/Mus_musculus.GRCm39.113.gtf.gz
   ```

### Run Full Pipeline

```bash
# Execute the complete analysis pipeline
bash run_pipeline.sh
```

### Run Individual Steps

```bash
# Step 1-5: Linux commands (QC, trimming, alignment, quantification)
# These are embedded in run_pipeline.sh

# Step 6: DESeq2 analysis (R)
Rscript scripts/01_deseq2_analysis.R

# Step 7: Generate figures (R)
Rscript scripts/02_figures.R
```

### Output Files

| File | Description |
|------|-------------|
| `results/figures/Figure4a_PCA.pdf` | PCA plot showing sample clustering |
| `results/figures/Figure4b_Heatmap.pdf` | Heatmap of top 50 variable genes |
| `results/figures/Figure4c_Markers_dotplot.pdf` | Marker gene expression dot plot |
| `results/tables/DEGs_*.csv` | Differential expression gene lists |
| `data/processed/FPKM_matrix.csv` | Normalized FPKM expression matrix |

## Data Source

RNA-seq data were downloaded from GEO Series GSE154290:

[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154290)

Samples used in this analysis:

| Sample | GEO Accession | SRA Accession | Cell Type |
|--------|---------------|---------------|-----------|
| mESCs rep1 | GSM4668511 | SRR12207279 | Naive mESCs |
| mESCs rep2 | GSM4668512 | SRR12207280 | Naive mESCs |
| EpiSCs rep1 | GSM4668515 | SRR12207283 | Primed EpiSCs |
| EpiSCs rep2 | GSM4668516 | SRR12207284 | Primed EpiSCs |
| fPSCs P1 | GSM4668520 | SRR12207288 | Formative PSCs (Passage 1) |
| fPSCs P10 | GSM4668521 | SRR12207289 | Formative PSCs (Passage 10) |
| fPSCs P20 | GSM4668522 | SRR12207290 | Formative PSCs (Passage 20) |
| fPSCs P30 | GSM4668523 | SRR12207291 | Formative PSCs (Passage 30) |

Download example:

```bash
# Using SRA Toolkit
prefetch SRR12207279
fastq-dump --split-files SRR12207279

# Or using wget (paired-end)
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR12207279/SRR12207279_1.fastq.gz
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR12207279/SRR12207279_2.fastq.gz
```

## Reference Genome

- Assembly: GRCm39 (mm10) — Mus musculus
- Ensembl Release: 113
- Download: [Ensembl FTP](https://ftp.ensembl.org/pub/release-113/)

## Citation

If you use this code or analysis pipeline, please cite:

> Wang X, Xiang Y, Yu Y, et al. Formative pluripotent stem cells show features of epiblast cells poised for gastrulation. *Cell Research*. 2021;31:526-541. doi:[10.1038/s41422-021-00477-x](https://doi.org/10.1038/s41422-021-00477-x)

## License

This project is for educational and research purposes.

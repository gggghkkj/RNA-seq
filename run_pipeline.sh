#!/bin/bash
# =============================================================================
# run_pipeline.sh
# RNA-seq Analysis Pipeline for fPSCs Transcriptome
# Reference: "Formative pluripotent stem cells show features of epiblast cells
#             poised for gastrulation" (Cell Research, 2021)
#
# Usage: bash run_pipeline.sh
#
# Requirements:
#   - conda/mamba with rnaseq environment
#   - Tools: fastp, hisat2, samtools, featureCounts, R
# =============================================================================

set -euo pipefail

# --- Configuration -----------------------------------------------------------
export PATH="/home/gan/rnaseq/miniforge3/envs/rnaseq/bin:$PATH"

RAWDIR="data/raw"
CLEANDIR="data/processed/clean_fastq"
ALIGNDIR="data/processed/aligned"
COUNTDIR="data/processed"
LOGDIR="logs"
THREADS=4

# Reference genome
GENOME_FA="data/raw/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF="data/raw/Mus_musculus.GRCm39.113.gtf"
HISAT2_INDEX="data/processed/hisat2_index/hisat2_index"

# Sample names
SAMPLES=("mESCs_rep1" "mESCs_rep2" "EpiSCs_rep1" "EpiSCs_rep2"
         "fPSCs_P1" "fPSCs_P10" "fPSCs_P20" "fPSCs_P30")

# SRR IDs
declare -A SRR_MAP
SRR_MAP=(
    ["mESCs_rep1"]="SRR12207279"  ["mESCs_rep2"]="SRR12207280"
    ["EpiSCs_rep1"]="SRR12207283" ["EpiSCs_rep2"]="SRR12207284"
    ["fPSCs_P1"]="SRR12207288"    ["fPSCs_P10"]="SRR12207289"
    ["fPSCs_P20"]="SRR12207290"   ["fPSCs_P30"]="SRR12207291"
)

mkdir -p ${CLEANDIR} ${ALIGNDIR} ${LOGDIR} ${HISAT2_INDEX%/*}

# =============================================================================
# Step 1: Quality Control (FastQC)
# =============================================================================
echo "=========================================="
echo "Step 1: Quality Control with FastQC"
echo "=========================================="

for SAMPLE in "${SAMPLES[@]}"; do
    SRR=${SRR_MAP[$SAMPLE]}
    echo "  Running FastQC on ${SAMPLE} (${SRR})..."
    fastqc -o ${LOGDIR} -t ${THREADS} \
        ${RAWDIR}/${SRR}_1.fastq.gz \
        ${RAWDIR}/${SRR}_2.fastq.gz \
        2>&1 | tail -1
done

# Generate MultiQC report
multiqc ${LOGDIR} -o ${LOGDIR} --force 2>&1 | tail -1
echo "  FastQC + MultiQC complete!"

# =============================================================================
# Step 2: Adapter Trimming (fastp)
# =============================================================================
echo ""
echo "=========================================="
echo "Step 2: Adapter Trimming with fastp"
echo "=========================================="

for SAMPLE in "${SAMPLES[@]}"; do
    SRR=${SRR_MAP[$SAMPLE]}
    echo "  Trimming ${SAMPLE} (${SRR})..."

    fastp \
        -i ${RAWDIR}/${SRR}_1.fastq.gz \
        -I ${RAWDIR}/${SRR}_2.fastq.gz \
        -o ${CLEANDIR}/${SAMPLE}_1.fastq.gz \
        -O ${CLEANDIR}/${SAMPLE}_2.fastq.gz \
        --html ${LOGDIR}/${SAMPLE}_fastp.html \
        --json ${LOGDIR}/${SAMPLE}_fastp.json \
        --thread ${THREADS} \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        2>&1 | grep -E "Duplication|Filtering" | head -2
done

echo "  Trimming complete!"

# =============================================================================
# Step 3: Build HISAT2 Index (if not exists)
# =============================================================================
echo ""
echo "=========================================="
echo "Step 3: Build HISAT2 Genome Index"
echo "=========================================="

if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
    echo "  Building HISAT2 index..."

    # Extract splice sites and exons
    hisat2_extract_splice_sites.py ${GTF} > ${HISAT2_INDEX%/*}/splicesites.txt
    hisat2_extract_exons.py ${GTF} > ${HISAT2_INDEX%/*}/exons.txt

    # Build index (basic, without splice sites for lower memory)
    hisat2-build -p ${THREADS} \
        ${GENOME_FA} \
        ${HISAT2_INDEX} \
        2>&1 | tail -3

    echo "  HISAT2 index built!"
else
    echo "  HISAT2 index already exists, skipping..."
fi

# =============================================================================
# Step 4: Alignment (HISAT2)
# =============================================================================
echo ""
echo "=========================================="
echo "Step 4: Alignment with HISAT2"
echo "=========================================="

for SAMPLE in "${SAMPLES[@]}"; do
    echo "  Aligning ${SAMPLE}..."

    hisat2 -x ${HISAT2_INDEX} \
        -1 ${CLEANDIR}/${SAMPLE}_1.fastq.gz \
        -2 ${CLEANDIR}/${SAMPLE}_2.fastq.gz \
        --threads ${THREADS} \
        --dta \
        2> ${LOGDIR}/${SAMPLE}_hisat2.log \
    | samtools sort -@ ${THREADS} -o ${ALIGNDIR}/${SAMPLE}.sorted.bam

    samtools index ${ALIGNDIR}/${SAMPLE}.sorted.bam

    # Print alignment rate
    RATE=$(grep "overall alignment rate" ${LOGDIR}/${SAMPLE}_hisat2.log | head -1)
    echo "    ${SAMPLE}: ${RATE}"
done

echo "  Alignment complete!"

# =============================================================================
# Step 5: Gene Quantification (featureCounts)
# =============================================================================
echo ""
echo "=========================================="
echo "Step 5: Gene Quantification with featureCounts"
echo "=========================================="

featureCounts -T ${THREADS} -p -t exon -g gene_id \
    -a ${GTF} \
    -o ${COUNTDIR}/gene_counts.txt \
    ${ALIGNDIR}/*.sorted.bam \
    2>&1 | grep -E "Successfully|Total"

echo "  Quantification complete!"

# =============================================================================
# Step 6: Differential Expression Analysis (R/DESeq2)
# =============================================================================
echo ""
echo "=========================================="
echo "Step 6: DESeq2 Analysis (R)"
echo "=========================================="

Rscript scripts/01_deseq2_analysis.R
echo "  DESeq2 analysis complete!"

# =============================================================================
# Step 7: Generate Figures (R)
# =============================================================================
echo ""
echo "=========================================="
echo "Step 7: Generate Figures (R)"
echo "=========================================="

Rscript scripts/02_figures.R
echo "  Figure generation complete!"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Figures: results/figures/"
echo "  Tables:  results/tables/"
echo "  Data:    data/processed/"
echo ""
echo "Key outputs:"
echo "  - Figure4a_PCA.pdf: PCA plot"
echo "  - Figure4b_Heatmap.pdf: Top variable genes heatmap"
echo "  - Figure4c_Markers_dotplot.pdf: Marker gene expression"
echo "  - DEGs_*.csv: Differential expression results"
echo "  - FPKM_matrix.csv: Normalized expression matrix"
